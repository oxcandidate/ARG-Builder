/***************************************************************************
 *
 *    arg_builder.cpp
 *
 *    Implementation of logic for building ARGs by joining one sequence at
 *    a time.
 *
 ****************************************************************************/

#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <algorithm>
#include <iostream>
#include <random>
#include <fstream>

#include "arg_builder_logic.h"
#include "vector_set_operations.h"
#include "step_function.h"
#include "clean.h"

bool _multi_roots_given = false;
int _how_verbose = 0;
int _max_number_parents = 2;
float _cost_recurrent_mutation = 1.0;
float _cost_back_mutation = 1.0;
float _cost_recombination = 1.0;
int _recomb_max = INT32_MAX;
int _rm_max = INT32_MAX;
int _bm_max = INT32_MAX;

std::vector<int> _site_multiplicitity;
RunRecord _run_record;
int _run_seed = 0;
int _location_selection_method = 0;
std::string _run_record_file = "";

std::string vector_to_string(const std::vector<int> &v, bool make_negative = false)
{
    std::string s = "";
    bool first = true;
    for (const int i : v)
    {
        if (!first)
            s += ", ";
        if (make_negative)
            s += std::to_string(-i);
        else
            s += std::to_string(i);
        first = false;
    }

    return s;
}

std::string vector_to_string_with_highlights(const std::vector<int> &v, bool make_negative, const std::set<int> &highlights, const std::string highlight_symbol, bool only_show_highlighted)
{
    std::string s = "";
    bool first = true;
    for (const int i : v)
    {
        bool highlight = highlights.find(i) != highlights.end();

        if (highlight || !only_show_highlighted)
        {
            if (!first)
                s += ", ";
            if (make_negative)
                s += std::to_string(-i);
            else
                s += std::to_string(i);

            first = false;

            if (highlight)
                s += highlight_symbol;
        }
    }

    return s;
}

void print_gene(const Gene &g)
{
    if (g.label != "")
    {
        std::cout << g.label << ": ";
    }

    std::cout << "mutations: ";
    for (int i : g.mutations)
    {
        std::cout << i << ", ";
    }
    std::cout << "\n";
}

void print_arg(const ARG &arg)
{
    std::cout << "\n\n--------------printing arg----------------\n\n";
    std::cout << "-----------printing nodes:\n";
    for (auto &node_ptr : arg.nodes)
    {
        if (node_ptr->label != "")
        {
            std::cout << node_ptr->label << ": ";
        }

        std::cout << "mutations: ";
        for (int i : node_ptr->mutations)
        {
            std::cout << i << ", ";
        }
        std::cout << "\n";
    }

    std::cout << "-----------printing edges:\n";
    for (auto &edge_ptr : arg.edges)
    {
        std::cout << edge_ptr->from->label << " -> " << edge_ptr->to->label << " muts: ";
        for (int i : edge_ptr->mutations)
        {
            std::cout << i << ", ";
        }
        std::cout << " back_muts:";
        for (int i : edge_ptr->back_mutations)
        {
            std::cout << i << ", ";
        }
        std::cout << "\n";
    }

    std::cout << "-----------printing edges from nodes:\n";
    for (auto &node_ptr : arg.nodes)
    {
        if (node_ptr->label != "")
        {
            std::cout << node_ptr->label << ": ";
        }

        /* Output edges out of this node */
        if (node_ptr->type == SAMPLE || node_ptr->type == COALESCENCE)
        {
            std::cout << node_ptr->predecessor.one->from->id << " -> " << node_ptr->id;
        }
        else if (node_ptr->type == RECOMBINATION)
        {
            std::cout << node_ptr->predecessor.two.prefix->from->id << " -> " << node_ptr->id << " <- " << node_ptr->predecessor.two.suffix->from->id;
        }
        std::cout << "\n";
    }
}

float get_cost(const int rms, const int bms, const int recombs)
{
    if ((_cost_recurrent_mutation < 0 && rms > 0) || (_cost_back_mutation < 0 && bms > 0) || (_cost_recombination < 0 && recombs > 0))
        return -1.0;

    return static_cast<float>(rms) * _cost_recurrent_mutation +
           static_cast<float>(bms) * _cost_back_mutation +
           static_cast<float>(recombs) * _cost_recombination;
}

float get_cost(const int rms, const int bms)
{
    return get_cost(rms, bms, 0);
}

float get_cost(const ARG &arg)
{
    return get_cost(arg.number_of_recurrent_mutations, arg.number_of_back_mutations, arg.number_of_recombinations);
}

float get_cost_of_addition(const ARG &arg, const int rms, const int bms, const int recombs)
{
    if ((arg.number_of_recombinations + recombs > _recomb_max) || (arg.number_of_recurrent_mutations + rms > _rm_max) || (arg.number_of_back_mutations + bms > _bm_max))
        return -1.0;

    return get_cost(rms, bms, recombs);
}

float get_cost_of_addition(const ARG &arg, const int rms, const int bms)
{
    return get_cost_of_addition(arg, rms, bms, 0);
}

int get_mutations_represented(const std::vector<int> mutations)
{
    int total = 0;
    for (int mut : mutations)
    {
        total += _site_multiplicitity[mut];
    }

    return total;
}

int get_mutations_represented(int mutation)
{
    return _site_multiplicitity[mutation];
}

void arg_output(const ARG &arg, const Genes &genes, FILE *fp,
                ARGOutputFormat format, int how_to_label_edges, ARGOutputLabels node_labels)
{
    int intervals;

    if (fp == NULL)
        fp = stdout;

    /* Output prelude */
    switch (format)
    {
    case ARGDOT:
        fprintf(fp, "digraph ARG {\n");
        break;
    case ARGGDL:
        fprintf(fp, "graph: {\n");
        break;
    case ARGGML:
        fprintf(fp, "graph [\n  directed 1\n");
        break;
    }

    /* Output nodes and their edges */
    for (int i = 0; i < arg.nodes.size(); i++)
    {
        int node_id = arg.nodes[i]->id;

        /* Output node i */
        bool label_printed = false;
        switch (format)
        {
        case ARGDOT:
            fprintf(fp, "  %d [label=\"", node_id);
            break;
        case ARGGDL:
            fprintf(fp, "  node: { title: \"%d\" label: \"", node_id);
            break;
        case ARGGML:
            fprintf(fp, "  node [\n    id %d\n    label \"", node_id);
            break;
        }

        /* Generate node label */
        if (arg.nodes[i]->label != "" && (node_labels == LABEL_BOTH || node_labels == LABEL_LABEL))
        {
            fprintf(fp, "%s", arg.nodes[i]->label.c_str());
            label_printed = true;
        }
        if (node_labels == LABEL_BOTH || node_labels == LABEL_SEQUENCE)
        {
            if (label_printed)
                fprintf(fp, "; %s", vector_to_string(arg.nodes[i]->mutations).c_str());
            else
                fprintf(fp, "%s", vector_to_string(arg.nodes[i]->mutations).c_str());
            label_printed = true;
        }
        if (!label_printed)
        {
            /* No node label printed */
            if ((arg.nodes[i]->type == RECOMBINATION) && (arg.nodes[i]->label != ""))
                /* If a recombination node is still unlabeled, label it (label
                 * should be recombination point) if possible.
                 */
                fprintf(fp, "%s\"", arg.nodes[i]->label.c_str());
            else
            {
                switch (format)
                {
                case ARGDOT:
                    fprintf(fp, "\",shape=point");
                    break;
                case ARGGDL:
                    fprintf(fp, "\" scaling: 0.0");
                    break;
                case ARGGML:
                    fprintf(fp, " \"");
                    break;
                }
            }
        }
        else
            fprintf(fp, "\"");

        switch (format)
        {
        case ARGDOT:
            fprintf(fp, ",color=");
            break;
        case ARGGDL:
            fprintf(fp, " shape: circle bordercolor: ");
            break;
        case ARGGML:
            fprintf(fp, "\n    graphics [\n      outline \"");
            break;
        }
        switch (arg.nodes[i]->type)
        {
        case SAMPLE:
            fprintf(fp, "red");
            break;
        case COALESCENCE:
            fprintf(fp, "green");
            break;
        case RECOMBINATION:
            fprintf(fp, "blue");
            break;
        case ROOT:
            fprintf(fp, "black");
            break;
        case UNSET:
            fprintf(fp, "purple");
            break;
        }
        /* End node description */
        switch (format)
        {
        case ARGDOT:
            fprintf(fp, "];\n");
            break;
        case ARGGDL:
            if (i < genes.genes.size())
                fprintf(fp, " vertical_order: maxlevel");
            fprintf(fp, " }\n");
            break;
        case ARGGML:
            fprintf(fp,
                    "\"\n    ]\n    vgj [\n      type \"Oval\"\n      labelPosition \"in\"\n    ]\n  ]\n");
            break;
        }

        /* Output edges out of this node */
        switch (arg.nodes[i]->type)
        {
        case SAMPLE:
        case COALESCENCE:
        {
            /* One outgoing edge */
            int from_id = arg.nodes[i]->predecessor.one->from->id;
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "  %d -> %d", from_id, node_id);
                break;
            case ARGGDL:
                fprintf(fp, "  edge: { sourcename: \"%d\" targetname: \"%d\"",
                        from_id, node_id);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d",
                        from_id, node_id);
                break;
            }
            /* Output edge label */
            if (how_to_label_edges > 0)
            {
                switch (format)
                {
                case ARGDOT:
                    fprintf(fp, " [label=\"");
                    break;
                case ARGGDL:
                    fprintf(fp, " label: \"");
                    break;
                case ARGGML:
                    fprintf(fp, "\n    label \"");
                    break;
                }
                auto muts = vector_to_string_with_highlights(arg.nodes[i]->predecessor.one->mutations, false, arg.recurrent_sites, "*", how_to_label_edges == 1);
                auto back_muts = vector_to_string(arg.nodes[i]->predecessor.one->back_mutations, true);
                if (muts.length() + back_muts.length() > 0)
                    fprintf(fp, "%s|%s", muts.c_str(), back_muts.c_str());
                fprintf(fp, "\"");
                if (format == ARGDOT)
                    fprintf(fp, "]");
            }
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, ";\n");
                break;
            case ARGGDL:
                fprintf(fp, " }\n");
                break;
            case ARGGML:
                fprintf(fp, "\n  ]\n");
                break;
            }
        }
        break;
        case RECOMBINATION:
        {
            /* Two outgoing edges */
            /* Prefix edge */
            int prefix_id = arg.nodes[i]->predecessor.two.prefix->from->id;
            int suffix_id = arg.nodes[i]->predecessor.two.suffix->from->id;
            int crossover_position = arg.nodes[i]->predecessor.two.position;
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "  %d -> %d [label=\"P %d",
                        prefix_id, node_id, crossover_position);
                break;
            case ARGGDL:
                fprintf(fp,
                        "  edge: { sourcename: \"%d\" targetname: \"%d\" label: \"P %d",
                        prefix_id, node_id, crossover_position);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d\n    label \"P %d",
                        prefix_id, node_id, crossover_position);
                break;
            }
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "\"]\n");
                break;
            case ARGGDL:
                fprintf(fp, "\" }\n");
                break;
            case ARGGML:
                fprintf(fp, "\"\n  ]\n");
                break;
            }
            /* Suffix edge */
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "  %d -> %d [label=\"S %d",
                        suffix_id, node_id, crossover_position);
                break;
            case ARGGDL:
                fprintf(fp,
                        "  edge: { sourcename: \"%d\" targetname: \"%d\" label: \"S %d",
                        suffix_id, node_id, crossover_position);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d\n    label \"S %d",
                        suffix_id, node_id, crossover_position);
                break;
            }
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "\"]\n");
                break;
            case ARGGDL:
                fprintf(fp, "\" }\n");
                break;
            case ARGGML:
                fprintf(fp, "\"\n  ]\n");
                break;
            }
        }
        break;
        case ROOT:
            /* No outgoing edges */
            break;
        }
    }

    /* Output postlude */
    switch (format)
    {
    case ARGDOT:
    case ARGGDL:
        fprintf(fp, "}\n");
        break;
    case ARGGML:
        fprintf(fp, "]\n");
        break;
    }
}

void remove_edge_from_multimap(std::multimap<int, Edge *> &map, int mut, Edge *edge)
{
    auto iterpair = map.equal_range(mut);

    for (auto it = iterpair.first; it != iterpair.second; ++it)
    {
        if (it->second == edge)
        {
            map.erase(it);
            break;
        }
    }
}

std::tuple<std::set<Edge *>, std::vector<int>> find_relevant_edges(ARG &arg, const Gene &g)
{
    std::set<Edge *> relevant_edges;

    std::vector<int> existing_mutations;
    for (int mut : g.mutations)
    {
        bool contained = false;

        if (_multi_roots_given && vector_contains(arg.root_mutations, mut))
            contained = true; // This is only used in the case of multiroots, otherwise everything is flipped so root is all zero.

        // Any edge containing a mutation which is in g may be relevant
        auto range = arg.mutation_to_edges.equal_range(mut);
        for (auto i = range.first; i != range.second; ++i)
        {
            relevant_edges.insert(i->second);
            contained = true;
        }
        if (contained)
            existing_mutations.push_back(mut); // These should be accounted for or else RM

        // Any recombination node which shares mutations with g may be relevant
        range = arg.mutation_to_recombinations.equal_range(mut);
        for (auto i = range.first; i != range.second; ++i)
        {
            relevant_edges.insert(i->second);
        }
    }

    // Any edge with a back mutation for a site not in g may be relevant.
    std::vector<int> relevant_bms;
    std::set_difference(arg.back_mutation_sites.begin(), arg.back_mutation_sites.end(),
                        g.mutations.begin(), g.mutations.end(),
                        std::inserter(relevant_bms, relevant_bms.begin()));
    for (int back_mut : relevant_bms)
    {
        auto range = arg.back_mutation_to_edges.equal_range(back_mut);
        for (auto i = range.first; i != range.second; ++i)
        {
            relevant_edges.insert(i->second);
        }
    }

    // Add all self loops for roots
    for (const auto &self_edge : arg.self_edges)
    {
        if (self_edge->from->type == ROOT)
        {
            relevant_edges.insert(self_edge.get());
        }
    }

    return std::make_tuple(relevant_edges, existing_mutations);
}

/* Splits the edge into two and returns the new node. removable mutations are put below, rest are above */
Node *split_edge(ARG &arg, Edge *edge, std::vector<int> removable_muts, std::vector<int> removable_back_muts)
{
    if (edge->from == edge->to)
    {
        std::cerr << "Trying to split a self loop on " << edge->from->label << "\n";
        return edge->from;
    }

    auto split_node = std::make_unique<Node>(); // Between the top and bottom of best_edge
    auto split_edge = std::make_unique<Edge>(); // From top to split_node

    split_edge->from = edge->from;
    split_edge->to = edge->from = split_node.get();

    split_edge->mutations = vector_difference(edge->mutations, removable_back_muts);
    split_edge->back_mutations = vector_difference(edge->back_mutations, removable_muts);
    edge->mutations = removable_back_muts;
    edge->back_mutations = removable_muts;

    split_node->id = arg.nodes.size();
    split_node->label = "ancestral" + std::to_string(arg.number_of_ancestral_nodes);
    split_node->type = COALESCENCE;
    split_node->predecessor.one = split_edge.get();
    split_node->mutations = vector_difference(vector_union(split_edge->from->mutations, split_edge->mutations), split_edge->back_mutations);

    for (int m : split_edge->mutations)
    {
        arg.mutation_to_edges.insert({m, split_edge.get()});
        remove_edge_from_multimap(arg.mutation_to_edges, m, edge);
    }
    for (int m : split_edge->back_mutations)
    {
        arg.back_mutation_to_edges.insert({m, split_edge.get()});
        remove_edge_from_multimap(arg.back_mutation_to_edges, m, edge);
    }

    Node *new_node = split_node.get();

    arg.edges.push_back(std::move(split_edge));
    arg.nodes.push_back(std::move(split_node));
    arg.number_of_ancestral_nodes += 1;

    return new_node;
}

Node *insert_seq_as_direct_child(ARG &arg, const Gene &g, Node *parent, const std::vector<int> &existing_mutations)
{
    auto new_node = std::make_unique<Node>();
    auto new_edge = std::make_unique<Edge>();

    new_node->id = arg.nodes.size();
    new_node->label = g.label;
    new_node->mutations = g.mutations;
    new_node->type = SAMPLE;
    new_node->predecessor.one = new_edge.get();

    new_edge->from = parent;
    new_edge->to = new_node.get();

    new_edge->mutations = vector_difference(g.mutations, parent->mutations);
    new_edge->back_mutations = vector_difference(parent->mutations, g.mutations);
    auto recurrent_muts = vector_intersect(new_edge->mutations, existing_mutations);

    for (int m : new_edge->mutations)
        arg.mutation_to_edges.insert({m, new_edge.get()});
    for (int m : new_edge->back_mutations)
        arg.back_mutation_to_edges.insert({m, new_edge.get()});

    // Update arg counts
    arg.number_of_recurrent_mutations += get_mutations_represented(recurrent_muts);
    arg.number_of_back_mutations += get_mutations_represented(new_edge->back_mutations);
    for (auto m : recurrent_muts)
    {
        arg.recurrent_sites.insert(m);
    }
    for (auto m : new_edge->back_mutations)
    {
        arg.back_mutation_sites.insert(m);
    }

    Node *return_value = new_node.get();

    arg.edges.push_back(std::move(new_edge));
    arg.nodes.push_back(std::move(new_node));

    return return_value;
}

void pass_sample_to_leaf(ARG &arg, Node *parent)
{
    // Want all samples to be leaves of the graph
    auto new_node = std::make_unique<Node>();
    auto new_edge = std::make_unique<Edge>();

    new_node->id = arg.nodes.size();
    new_node->label = parent->label;
    new_node->mutations = parent->mutations;
    new_node->type = SAMPLE;
    new_node->predecessor.one = new_edge.get();

    new_edge->from = parent;
    new_edge->to = new_node.get();

    new_edge->mutations.clear();
    new_edge->back_mutations.clear();

    arg.edges.push_back(std::move(new_edge));
    arg.nodes.push_back(std::move(new_node));

    parent->label = "ancestral" + std::to_string(arg.number_of_ancestral_nodes);
    parent->type = COALESCENCE;

    arg.number_of_ancestral_nodes += 1;
}

Node *recombine_nodes(ARG &arg, const int pos, Node *prefix, Node *suffix)
{
    if (prefix->type == SAMPLE)
        pass_sample_to_leaf(arg, prefix);
    if (suffix->type == SAMPLE)
        pass_sample_to_leaf(arg, suffix);

    auto recomb_node = std::make_unique<Node>();
    auto prefix_edge = std::make_unique<Edge>();
    auto suffix_edge = std::make_unique<Edge>();

    prefix_edge->from = prefix;
    prefix_edge->to = recomb_node.get();
    suffix_edge->from = suffix;
    suffix_edge->to = recomb_node.get();

    recomb_node->id = arg.nodes.size();
    recomb_node->label = "recomb" + std::to_string(arg.number_of_ancestral_nodes);
    recomb_node->mutations = vector_union(vector_values_below(prefix->mutations, pos),
                                          vector_values_above(suffix->mutations, pos));
    recomb_node->type = RECOMBINATION;
    recomb_node->predecessor.two.position = pos;
    recomb_node->predecessor.two.prefix = prefix_edge.get();
    recomb_node->predecessor.two.suffix = suffix_edge.get();

    auto self_loop = arg.create_self_loop(recomb_node.get());
    for (int m : recomb_node->mutations)
    {
        arg.mutation_to_recombinations.insert({m, self_loop});
    }

    Node *return_value = recomb_node.get();
    arg.recombination_nodes.insert(return_value);

    arg.nodes.push_back(std::move(recomb_node));
    arg.edges.push_back(std::move(prefix_edge));
    arg.edges.push_back(std::move(suffix_edge));
    arg.number_of_ancestral_nodes += 1;
    arg.number_of_recombinations += 1;

    return return_value;
}

std::tuple<Edge *, int, int> find_best_single_parent_location(ARG &arg, const Gene &g, const std::set<Edge *> &related_edges, const std::vector<int> &existing_mutations)
{
    int best_rms = 0;
    int best_bms = 0;
    float best_score = -1.0;
    Edge *best_edge = nullptr;

    for (Edge *edge : related_edges)
    {
        // First find needed mutations relative to target node of edge
        std::vector<int> recurrent_mutations = vector_difference(existing_mutations, edge->to->mutations); // TODO: a symmetric difference is going on
        std::vector<int> back_mutations = vector_difference(edge->to->mutations, existing_mutations);

        // Now consider if some can be removed by splitting edge
        recurrent_mutations = vector_difference(recurrent_mutations, edge->back_mutations);
        back_mutations = vector_difference(back_mutations, edge->mutations);

        int rms = get_mutations_represented(recurrent_mutations);
        int bms = get_mutations_represented(back_mutations);

        float score = get_cost_of_addition(arg, rms, bms);

        if (score >= 0 && (score < best_score || best_score < 0))
        {
            best_score = score;
            best_rms = rms;
            best_bms = bms;
            best_edge = edge;
        }
    }

    return std::make_tuple(best_edge, best_rms, best_bms);
}

/* returns a tuple of the list of single parent locations (edge, score, rms ,bms), and the best score */
std::tuple<std::vector<std::tuple<Edge *, int, int, int>>, int> score_single_parent_locations(ARG &arg, const Gene &g, const std::set<Edge *> &related_edges, const std::vector<int> &existing_mutations)
{
    std::vector<std::tuple<Edge *, int, int, int>> single_parent_locations; // List of (edge, score, rms, bms)

    float best_score = -1.0;

    for (Edge *edge : related_edges)
    {
        // First find needed mutations relative to target node of edge
        std::vector<int> recurrent_mutations = vector_difference(existing_mutations, edge->to->mutations); // TODO: a symmetric difference is going on
        std::vector<int> back_mutations = vector_difference(edge->to->mutations, existing_mutations);

        // Now consider if some can be removed by splitting edge
        recurrent_mutations = vector_difference(recurrent_mutations, edge->back_mutations);
        back_mutations = vector_difference(back_mutations, edge->mutations);

        int rms = get_mutations_represented(recurrent_mutations);
        int bms = get_mutations_represented(back_mutations);

        float score = get_cost_of_addition(arg, rms, bms);

        single_parent_locations.push_back(std::make_tuple(edge, score, rms, bms));

        if (score >= 0 && (score < best_score || best_score < 0))
        {
            best_score = score;
        }
    }

    return std::make_tuple(std::move(single_parent_locations), best_score);
}

void insert_seq_at_edge(ARG &arg, const Gene &g, Edge *best_edge, const std::vector<int> &existing_mutations)
{
    if (best_edge == nullptr)
    {
        std::cerr << "should not have nullptr for best_edge!";
        // Insert at root
        // insert_seq_as_direct_child(arg, g, arg.root_ptr, existing_mutations);
        return;
    }

    // Need to check if best required splitting edge
    std::vector<int> recurrent_mutations = vector_difference(existing_mutations, best_edge->to->mutations);
    std::vector<int> back_mutations = vector_difference(best_edge->to->mutations, existing_mutations);

    auto [removable_rms, required_rms] = vector_split(recurrent_mutations, best_edge->back_mutations);
    auto [removable_back_muts, required_back_muts] = vector_split(back_mutations, best_edge->mutations);
    auto new_mutations = vector_difference(g.mutations, existing_mutations);

    if (removable_rms.size() > 0 || removable_back_muts.size() > 0)
    {
        Node *split_node = split_edge(arg, best_edge, removable_rms, removable_back_muts);

        // Now check if the split node is the node we wish to insert
        // if (split_node->mutations == g.mutations)
        // {
        //     // Then relabel split
        //     split_node->label = g.label;
        //     split_node->type = SAMPLE;
        //     arg.number_of_ancestral_nodes -= 1;
        // }
        // else
        // {
        //     insert_seq_as_direct_child(arg, g, split_node, existing_mutations);
        // }

        insert_seq_as_direct_child(arg, g, split_node, existing_mutations);
    }
    else
    {
        if (best_edge->to->type == SAMPLE)
        {
            pass_sample_to_leaf(arg, best_edge->to);

            insert_seq_as_direct_child(arg, g, best_edge->to, existing_mutations);
        }
        else
        {
            insert_seq_as_direct_child(arg, g, best_edge->to, existing_mutations);
        }
    }
}

void replace_cost_if_better(const ARG &arg, std::map<int, std::tuple<Edge *, int, int>> &costs, int pos, Edge *edge, int rms, int bms)
{
    // Note that we assume that any values current in costs must have non-negative cost
    // add entry to prefix_costs
    auto search = costs.find(pos);
    float new_cost = get_cost_of_addition(arg, rms, bms);

    if (new_cost < 0)
        return;

    if (search != costs.end())
    {
        auto [current_edge, current_rms, current_bms] = search->second;
        if (new_cost < get_cost(current_rms, current_bms))
            costs.insert_or_assign(pos, std::make_tuple(edge, rms, bms));
    }
    else
    {
        costs.insert({pos, std::make_tuple(edge, rms, bms)});
    }
}

std::tuple<int, Edge *, Edge *, int, int> find_best_recomb_location(ARG &arg, const Gene &g, const std::set<Edge *> &related_edges, const std::vector<int> &existing_mutations)
{
    /**
     * This function will return a list of pairs of edges and associated metrics.
     * Each item in list represents a potential recombination join location.
     * Only considers single cross-over recombination.
     */

    /**
     * For each edge we want to calculate the cost of using it as a prefix or suffix.
     * position 0 refers to just before site 0, 1 just before 1 etc.
     * We don't bother considering the position after the end.
     *
     * prefix_costs is a map from cross-over position to a (edge, rms, bms) pair.
     * suffix_costs will be similar and even have the same keys.
     */
    std::map<int, std::tuple<Edge *, int, int>> prefix_costs;
    std::map<int, std::tuple<Edge *, int, int>> suffix_costs;

    for (Edge *edge : related_edges)
    {
        // First find needed mutations relative to target node of edge
        std::vector<int> recurrent_mutations = vector_difference(existing_mutations, edge->to->mutations); // TODO: a symmetric difference is going on
        std::vector<int> back_mutations = vector_difference(edge->to->mutations, existing_mutations);

        // Now consider if some can be removed by splitting edge
        auto req_recurrent_mutations = vector_difference(recurrent_mutations, edge->back_mutations);
        auto req_back_mutations = vector_difference(back_mutations, edge->mutations);

        int i = 0;
        int j = 0;
        int pos = 0;
        int rms = req_recurrent_mutations.size();
        int bms = req_back_mutations.size();

        // set 0 suffix cost and prefix cost
        replace_cost_if_better(arg, suffix_costs, 0, edge, rms, bms);

        while (i < rms && j < bms)
        {
            bool inc_i;
            if (req_recurrent_mutations[i] < req_back_mutations[j])
            {
                pos = req_recurrent_mutations[i];
                inc_i = true;
            }
            else if (req_recurrent_mutations[i] == req_back_mutations[j])
            {
                std::cerr << "Should not be possible to have a recurrent mutation and back mutation at same location";
            }
            else
            {
                pos = req_back_mutations[j];
                inc_i = false;
            }

            // update prefix and suffix costs
            replace_cost_if_better(arg, prefix_costs, pos, edge, i, j);
            replace_cost_if_better(arg, suffix_costs, pos, edge, rms - i, bms - j);

            if (inc_i)
                i += 1;
            else
                j += 1;
        }

        while (i < rms)
        {
            // So j must equal bms
            pos = req_recurrent_mutations[i];

            // update prefix and suffix costs
            replace_cost_if_better(arg, prefix_costs, pos, edge, i, bms);
            replace_cost_if_better(arg, suffix_costs, pos, edge, rms - i, 0);

            i += 1;
        }
        while (j < bms)
        {
            // So j must equal bms
            pos = req_back_mutations[j];

            // update prefix and suffix costs
            replace_cost_if_better(arg, prefix_costs, pos, edge, rms, j);
            replace_cost_if_better(arg, suffix_costs, pos, edge, 0, bms - j);

            j += 1;
        }

        replace_cost_if_better(arg, suffix_costs, pos + 1, edge, 0, 0);
    }

    // Now we can go through the costs (which should be sorted by key)
    // First turn maps into vectors so that we can iterate through both nicely
    std::vector<std::tuple<int, Edge *, int, int>> prefix_costs_vec;
    std::vector<std::tuple<int, Edge *, int, int>> suffix_costs_vec;
    for (const auto &[site, value] : prefix_costs)
    {
        auto [edge, rms, bms] = value;
        prefix_costs_vec.push_back(std::make_tuple(site, edge, rms, bms));
    }
    for (const auto &[site, value] : suffix_costs)
    {
        auto [edge, rms, bms] = value;
        suffix_costs_vec.push_back(std::make_tuple(site, edge, rms, bms));
    }

    // Now go through backwards to create optimal prefix costs, and forwards for optimal suffix costs
    std::vector<std::tuple<int, Edge *, int, int>> optimal_prefix_costs;
    std::vector<std::tuple<int, Edge *, int, int>> optimal_suffix_costs;
    float best_cost = -1.0;
    for (int i = prefix_costs_vec.size() - 1; i >= 0; i--)
    {
        auto [site, edge, rms, bms] = prefix_costs_vec[i];
        float cost = get_cost_of_addition(arg, rms, bms);
        if (cost >= 0 && (cost < best_cost || best_cost < 0))
        {
            optimal_prefix_costs.push_back(prefix_costs_vec[i]);
            best_cost = cost;
        }
    }
    std::reverse(optimal_prefix_costs.begin(), optimal_prefix_costs.end());
    best_cost = -1.0;
    for (int i = 0; i < suffix_costs_vec.size(); i++)
    {
        auto [site, edge, rms, bms] = suffix_costs_vec[i];
        float cost = get_cost_of_addition(arg, rms, bms);
        if (cost >= 0 && (cost < best_cost || best_cost < 0))
        {
            optimal_suffix_costs.push_back(suffix_costs_vec[i]);
            best_cost = cost;
        }
    }

    if (optimal_prefix_costs.size() == 0 or optimal_suffix_costs.size() == 0)
    {
        return std::make_tuple(-1, nullptr, nullptr, 0, 0);
    }

    int p_index = 0, s_index = 0;
    auto [p_site, p_edge, p_rms, p_bms] = optimal_prefix_costs[p_index];
    auto [s_site, s_edge, s_rms, s_bms] = optimal_suffix_costs[s_index];

    int best_pos = -1;
    best_cost = -1.0;
    int best_rms = -1;
    int best_bms = -1;
    Edge *best_prefix = nullptr;
    Edge *best_suffix = nullptr;

    while (p_index < optimal_prefix_costs.size())
    {
        // Now need to find the highest s_index such that
        while (s_index < optimal_suffix_costs.size() - 1)
        {
            int next_site = std::get<0>(optimal_suffix_costs[s_index + 1]);
            if (next_site <= p_site)
            {
                s_index += 1;
                std::tie(s_site, s_edge, s_rms, s_bms) = optimal_suffix_costs[s_index];
            }
            else
            {
                break;
            }
        }

        float current_cost = get_cost_of_addition(arg, p_rms, p_bms) + get_cost_of_addition(arg, s_rms, s_bms);
        if (current_cost >= 0 && (current_cost < best_cost || best_cost < 0))
        {
            best_pos = p_site;
            best_cost = current_cost;
            best_rms = p_rms + s_rms;
            best_bms = p_bms + s_bms;
            best_prefix = p_edge;
            best_suffix = s_edge;
        }

        p_index++;
        std::tie(p_site, p_edge, p_rms, p_bms) = optimal_prefix_costs[p_index];
    }

    if (best_prefix == nullptr)
    {
        if (_how_verbose >= 2)
            std::cout << "Could not find recombination location for g: " << g.label << "\n";
        return std::make_tuple(-1, nullptr, nullptr, -1, -1);
    }
    else if (best_suffix == best_prefix)
    {
        if (_how_verbose >= 2)
            std::cout << "Recombination just used one sequence for g: " << g.label << "\n";
        return std::make_tuple(-1, nullptr, nullptr, -1, -1);
    }
    else
    {
        if (_how_verbose >= 2)
            std::cout << "Best recomb is at pos: " << best_pos << ", cost: " << best_cost << ", prefix: " << best_prefix->to->label << ", suffix: " << best_suffix->to->label << "\n";
        return std::make_tuple(best_pos, best_prefix, best_suffix, best_rms, best_bms);
    }
}

std::tuple<int, std::vector<int>, std::vector<Edge *>, int, int> find_best_multi_recomb_location(ARG &arg, const Gene &g, const std::set<Edge *> &related_edges, const std::vector<int> &existing_mutations)
{
    /**
     * This function will return a tuple: (number of parents, list of cross-over points, list of edges which are parent to each segment, total rms, total bms)
     */

    /**
     * For each edge we want to calculate the cost of using it as a prefix or suffix.
     * position 0 refers to just before site 0, 1 just before 1 etc.
     * We don't bother considering the position after the end.
     *
     * prefix_costs is a map from cross-over position to a (edge, rms, bms) pair.
     * suffix_costs will be similar and even have the same keys.
     */
    std::map<int, std::tuple<Edge *, int, int>> prefix_costs;
    std::map<int, std::tuple<Edge *, int, int>> suffix_costs;

    std::map<Edge *, std::unique_ptr<StepFunction>> edge_cost_functions;

    for (Edge *edge : related_edges)
    {
        // First find needed mutations relative to target node of edge
        std::vector<int> recurrent_mutations = vector_difference(existing_mutations, edge->to->mutations); // TODO: a symmetric difference is going on
        std::vector<int> back_mutations = vector_difference(edge->to->mutations, existing_mutations);

        // Now consider if some can be removed by splitting edge
        auto req_recurrent_mutations = vector_difference(recurrent_mutations, edge->back_mutations);
        auto req_back_mutations = vector_difference(back_mutations, edge->mutations);

        int i = 0;
        int j = 0;
        int pos = 0;
        int total_rms = get_mutations_represented(req_recurrent_mutations);
        int total_bms = get_mutations_represented(req_back_mutations);

        int current_rms = 0;
        int current_bms = 0;

        auto step_function = std::make_unique<StepFunction>(total_rms, total_bms);

        // set 0 suffix cost and prefix cost
        // replace_cost_if_better(arg, suffix_costs, 0, edge, rms, bms);

        while (i < req_recurrent_mutations.size() && j < req_back_mutations.size())
        {
            bool inc_i;
            if (req_recurrent_mutations[i] < req_back_mutations[j])
            {
                pos = req_recurrent_mutations[i];
                inc_i = true;
            }
            else if (req_recurrent_mutations[i] == req_back_mutations[j])
            {
                std::cerr << "Should not be possible to have a recurrent mutation and back mutation at same location";
            }
            else
            {
                pos = req_back_mutations[j];
                inc_i = false;
            }

            // update prefix and suffix costs
            replace_cost_if_better(arg, prefix_costs, pos, edge, current_rms, current_bms);
            replace_cost_if_better(arg, suffix_costs, pos, edge, total_rms - current_rms, total_bms - current_bms);

            float cost = get_cost_of_addition(arg, current_rms, current_bms);
            if (cost >= 0)
                step_function->add_point(pos, current_rms, current_bms, cost);

            if (inc_i)
            {
                current_rms += get_mutations_represented(pos);
                i += 1;
            }
            else
            {
                current_bms += get_mutations_represented(pos);
                j += 1;
            }
        }

        while (i < req_recurrent_mutations.size())
        {
            // So j must equal bms
            pos = req_recurrent_mutations[i];

            // update prefix and suffix costs
            replace_cost_if_better(arg, prefix_costs, pos, edge, current_rms, total_bms);
            replace_cost_if_better(arg, suffix_costs, pos, edge, total_rms - current_rms, 0);

            float cost = get_cost_of_addition(arg, current_rms, total_bms);
            if (cost >= 0)
                step_function->add_point(pos, current_rms, total_bms, cost);

            current_rms += get_mutations_represented(pos);
            i += 1;
        }
        while (j < req_back_mutations.size())
        {
            // So j must equal bms
            pos = req_back_mutations[j];

            // update prefix and suffix costs
            replace_cost_if_better(arg, prefix_costs, pos, edge, total_rms, current_bms);
            replace_cost_if_better(arg, suffix_costs, pos, edge, 0, total_bms - current_bms);

            float cost = get_cost_of_addition(arg, total_rms, current_bms);
            if (cost >= 0)
                step_function->add_point(pos, total_rms, current_bms, cost);

            current_bms += get_mutations_represented(pos);
            j += 1;
        }

        replace_cost_if_better(arg, suffix_costs, pos + 1, edge, 0, 0);
        // replace_cost_if_better(arg, prefix_costs, pos + 1, edge, rms, bms);

        edge_cost_functions.insert({edge, std::move(step_function)});
    }

    // Now we can go through the costs (which should be sorted by key)
    // First turn maps into vectors so that we can iterate through both nicely
    std::vector<std::tuple<int, Edge *, int, int>> prefix_costs_vec;
    std::vector<std::tuple<int, Edge *, int, int>> suffix_costs_vec;
    for (const auto &[site, value] : prefix_costs)
    {
        auto [edge, rms, bms] = value;
        prefix_costs_vec.push_back(std::make_tuple(site, edge, rms, bms));
    }
    for (const auto &[site, value] : suffix_costs)
    {
        auto [edge, rms, bms] = value;
        suffix_costs_vec.push_back(std::make_tuple(site, edge, rms, bms));
    }

    // Now go through backwards to create optimal prefix costs, and forwards for optimal suffix costs
    std::vector<std::tuple<int, Edge *, int, int>> optimal_prefix_costs;
    std::vector<std::tuple<int, Edge *, int, int>> optimal_suffix_costs;
    float best_cost = -1.0;
    for (int i = prefix_costs_vec.size() - 1; i >= 0; i--)
    {
        auto [site, edge, rms, bms] = prefix_costs_vec[i];
        float cost = get_cost_of_addition(arg, rms, bms);
        if (cost >= 0 && (cost < best_cost || best_cost < 0))
        {
            optimal_prefix_costs.push_back(prefix_costs_vec[i]);
            best_cost = cost;
        }
    }
    std::reverse(optimal_prefix_costs.begin(), optimal_prefix_costs.end());
    best_cost = -1.0;
    for (int i = 0; i < suffix_costs_vec.size(); i++)
    {
        auto [site, edge, rms, bms] = suffix_costs_vec[i];
        float cost = get_cost_of_addition(arg, rms, bms);
        if (cost >= 0 && (cost < best_cost || best_cost < 0))
        {
            optimal_suffix_costs.push_back(suffix_costs_vec[i]);
            best_cost = cost;
        }
    }

    if (optimal_prefix_costs.size() == 0 or optimal_suffix_costs.size() == 0)
    {
        return std::make_tuple(0, std::vector<int>(), std::vector<Edge *>(), 0, 0);
    }

    std::vector<int> change_points; // This only needs to be calculated if we have three of more crossover positions
    if (_max_number_parents >= 4)
    {
        int p_index = 0, s_index = 0;
        while (p_index < optimal_prefix_costs.size() && s_index < optimal_suffix_costs.size())
        {
            auto [p_site, x, y, z] = optimal_prefix_costs[p_index];
            auto [s_site, a, b, c] = optimal_suffix_costs[s_index];
            if (p_site < s_site)
            {
                p_index += 1;
                change_points.push_back(p_site);
            }
            else if (p_site == s_site)
            {
                p_index += 1;
                s_index += 1;
                change_points.push_back(p_site);
            }
            else
            {
                s_index += 1;
                change_points.push_back(s_site);
            }
        }

        while (p_index < optimal_prefix_costs.size())
        {
            auto [p_site, x, y, z] = optimal_prefix_costs[p_index];
            change_points.push_back(p_site);
            p_index += 1;
        }
        while (s_index < optimal_suffix_costs.size())
        {
            auto [s_site, x, y, z] = optimal_suffix_costs[s_index];
            change_points.push_back(s_site);
            s_index += 1;
        }
    }

    int best_number_parents = 0;
    std::vector<int> best_positions;
    std::vector<Edge *> best_parents;
    best_cost = -1.0;
    int best_rms = -1;
    int best_bms = -1;

    if (_max_number_parents >= 2)
    {
        int p_index = 0, s_index = 0;
        auto [p_site, p_edge, p_rms, p_bms] = optimal_prefix_costs[p_index];
        auto [s_site, s_edge, s_rms, s_bms] = optimal_suffix_costs[s_index];

        while (p_index < optimal_prefix_costs.size())
        {
            // Now need to find the highest s_index such that s_site <= p_site
            while (s_index < optimal_suffix_costs.size() - 1)
            {
                int next_site = std::get<0>(optimal_suffix_costs[s_index + 1]);
                if (next_site <= p_site)
                {
                    s_index += 1;
                    std::tie(s_site, s_edge, s_rms, s_bms) = optimal_suffix_costs[s_index];
                }
                else
                {
                    break;
                }
            }

            float current_cost = get_cost_of_addition(arg, p_rms + s_rms, p_bms + s_bms, 1);
            if (current_cost >= 0 && (current_cost < best_cost || best_cost < 0))
            {
                best_number_parents = 2;
                best_positions = {p_site};
                best_parents = {p_edge, s_edge};
                best_cost = current_cost;
                best_rms = p_rms + s_rms;
                best_bms = p_bms + s_bms;
            }

            p_index++;
            std::tie(p_site, p_edge, p_rms, p_bms) = optimal_prefix_costs[p_index];
        }
    }

    if (_max_number_parents >= 3 && get_cost_of_addition(arg, 0, 0, 2) >= 0 && (best_cost < 0 || get_cost_of_addition(arg, 0, 0, 2) < best_cost))
    {
        int p_index = 0, s_index = 0;
        auto [p_site, p_edge, p_rms, p_bms] = optimal_prefix_costs[p_index];
        auto [s_site, s_edge, s_rms, s_bms] = optimal_suffix_costs[s_index];

        while (p_index < optimal_prefix_costs.size())
        {
            std::tie(p_site, p_edge, p_rms, p_bms) = optimal_prefix_costs[p_index];

            s_index = optimal_suffix_costs.size() - 1;
            while (s_index >= 0)
            {
                std::tie(s_site, s_edge, s_rms, s_bms) = optimal_suffix_costs[s_index];

                if (s_site <= p_site)
                {
                    break; // This would just be a 2 parent recombination
                }
                else
                {
                    for (Edge *m_edge : related_edges)
                    {
                        auto [m_rms, m_bms] = edge_cost_functions[m_edge]->get_range(p_site, s_site);
                        float current_cost = get_cost_of_addition(arg, p_rms + m_rms + s_rms, p_bms + m_bms + s_bms, 2);
                        if (current_cost >= 0 && (current_cost < best_cost || best_cost < 0))
                        {
                            best_number_parents = 3;
                            best_positions = {p_site, s_site};
                            best_parents = {p_edge, m_edge, s_edge};
                            best_cost = current_cost;
                            best_rms = p_rms + m_rms + s_rms;
                            best_bms = p_bms + m_bms + s_bms;
                        }
                    }
                }

                s_index--;
            }

            p_index++;
        }
    }

    if (_max_number_parents >= 4 && get_cost_of_addition(arg, 0, 0, 3) >= 0 && (best_cost < 0 || get_cost_of_addition(arg, 0, 0, 3) < best_cost))
    {
        int p_index = 0, s_index = 0;
        auto [p_site, p_edge, p_rms, p_bms] = optimal_prefix_costs[p_index];
        auto [s_site, s_edge, s_rms, s_bms] = optimal_suffix_costs[s_index];

        while (p_index < optimal_prefix_costs.size())
        {
            std::tie(p_site, p_edge, p_rms, p_bms) = optimal_prefix_costs[p_index];

            s_index = optimal_suffix_costs.size() - 1;
            while (s_index >= 0)
            {
                std::tie(s_site, s_edge, s_rms, s_bms) = optimal_suffix_costs[s_index];

                if (s_site <= p_site)
                {
                    break; // This would just be a 2 parent recombination
                }
                else
                {
                    for (int midpoint : change_points)
                    {
                        if (midpoint <= p_site)
                            continue;
                        else if (midpoint >= s_site)
                            break;

                        Edge *best_middle1;
                        float best_cost1 = -1.0;
                        Edge *best_middle2;
                        float best_cost2 = -1.0;

                        for (Edge *m_edge : related_edges)
                        {
                            auto [m1_rms, m1_bms] = edge_cost_functions[m_edge]->get_range(p_site, midpoint);
                            float cost1 = get_cost_of_addition(arg, m1_rms, m1_bms);
                            if (cost1 >= 0 && (cost1 < best_cost1 || best_cost1 < 0))
                            {
                                best_cost1 = cost1;
                                best_middle1 = m_edge;
                            }

                            auto [m2_rms, m2_bms] = edge_cost_functions[m_edge]->get_range(midpoint, s_site);
                            float cost2 = get_cost_of_addition(arg, m2_rms, m2_bms);
                            if (cost2 >= 0 && (cost2 < best_cost2 || best_cost2 < 0))
                            {
                                best_cost2 = cost2;
                                best_middle2 = m_edge;
                            }
                        }

                        auto [m1_rms, m1_bms] = edge_cost_functions[best_middle1]->get_range(p_site, midpoint);
                        auto [m2_rms, m2_bms] = edge_cost_functions[best_middle2]->get_range(midpoint, s_site);
                        float current_cost = get_cost_of_addition(arg, p_rms + m1_rms + m2_rms + s_rms, p_bms + m1_bms + m2_bms + s_bms, 3);
                        if (current_cost >= 0 && (current_cost < best_cost || best_cost < 0))
                        {
                            best_number_parents = 4;
                            best_positions = {p_site, midpoint, s_site};
                            best_parents = {p_edge, best_middle1, best_middle2, s_edge};
                            best_cost = current_cost;
                            best_rms = p_rms + m1_rms + m2_rms + s_rms;
                            best_bms = p_bms + m1_bms + m2_bms + s_bms;
                        }
                    }
                }

                s_index--;
            }

            p_index++;
        }
    }

    if (best_cost < 0)
    {
        if (_how_verbose >= 2)
            std::cout << "Could not find recombination location for g: " << g.label << "\n";
        return std::make_tuple(0, std::vector<int>(), std::vector<Edge *>(), 0, 0);
    }
    else
    {
        if (_how_verbose >= 2)
        {
            std::cout << "Best recomb has " << best_number_parents << " parents: ";
            for (auto p : best_parents)
                std::cout << p->to->label << ", ";
            std::cout << ". Cross-over positions: " << vector_to_string(best_positions) << ".";
            std::cout << " Cost was rms: " << best_rms << " bms: " << best_bms << "\n";
        }
        return std::make_tuple(best_number_parents, best_positions, best_parents, best_rms, best_bms);
    }
}

std::tuple<std::vector<std::tuple<int, std::vector<int>, std::vector<Edge *>, int, int, int>>, int> score_multi_recomb_locations(ARG &arg, const Gene &g, const std::set<Edge *> &related_edges, const std::vector<int> &existing_mutations)
{
    /**
     * This function will return a tuple: (number of parents, list of cross-over points, list of edges which are parent to each segment, total rms, total bms)
     */

    std::vector<std::tuple<int, std::vector<int>, std::vector<Edge *>, int, int, int>> locations; // vectors of tuples(#parents, crossovers, parent edges, score, total rms, total bms)
    float best_score = -1.0;

    /**
     * For each edge we want to calculate the cost of using it as a prefix or suffix.
     * position 0 refers to just before site 0, 1 just before 1 etc.
     * We don't bother considering the position after the end.
     *
     * prefix_costs is a map from cross-over position to a (edge, rms, bms) pair.
     * suffix_costs will be similar and even have the same keys.
     */
    std::map<int, std::tuple<Edge *, int, int>> prefix_costs;
    std::map<int, std::tuple<Edge *, int, int>> suffix_costs;

    std::map<Edge *, std::unique_ptr<StepFunction>> edge_cost_functions;

    for (Edge *edge : related_edges)
    {
        // First find needed mutations relative to target node of edge
        std::vector<int> recurrent_mutations = vector_difference(existing_mutations, edge->to->mutations); // TODO: a symmetric difference is going on
        std::vector<int> back_mutations = vector_difference(edge->to->mutations, existing_mutations);

        // Now consider if some can be removed by splitting edge
        auto req_recurrent_mutations = vector_difference(recurrent_mutations, edge->back_mutations);
        auto req_back_mutations = vector_difference(back_mutations, edge->mutations);

        int i = 0;
        int j = 0;
        int pos = 0;
        int total_rms = get_mutations_represented(req_recurrent_mutations);
        int total_bms = get_mutations_represented(req_back_mutations);

        int current_rms = 0;
        int current_bms = 0;

        auto step_function = std::make_unique<StepFunction>(total_rms, total_bms);

        // set 0 suffix cost and prefix cost
        // replace_cost_if_better(arg, suffix_costs, 0, edge, rms, bms);

        while (i < req_recurrent_mutations.size() && j < req_back_mutations.size())
        {
            bool inc_i;
            if (req_recurrent_mutations[i] < req_back_mutations[j])
            {
                pos = req_recurrent_mutations[i];
                inc_i = true;
            }
            else if (req_recurrent_mutations[i] == req_back_mutations[j])
            {
                std::cerr << "Should not be possible to have a recurrent mutation and back mutation at same location";
            }
            else
            {
                pos = req_back_mutations[j];
                inc_i = false;
            }

            // update prefix and suffix costs
            replace_cost_if_better(arg, prefix_costs, pos, edge, current_rms, current_bms);
            replace_cost_if_better(arg, suffix_costs, pos, edge, total_rms - current_rms, total_bms - current_bms);

            float cost = get_cost_of_addition(arg, current_rms, current_bms);
            if (cost >= 0)
                step_function->add_point(pos, current_rms, current_bms, cost);

            if (inc_i)
            {
                current_rms += get_mutations_represented(pos);
                i += 1;
            }
            else
            {
                current_bms += get_mutations_represented(pos);
                j += 1;
            }
        }

        while (i < req_recurrent_mutations.size())
        {
            // So j must equal bms
            pos = req_recurrent_mutations[i];

            // update prefix and suffix costs
            replace_cost_if_better(arg, prefix_costs, pos, edge, current_rms, total_bms);
            replace_cost_if_better(arg, suffix_costs, pos, edge, total_rms - current_rms, 0);

            float cost = get_cost_of_addition(arg, current_rms, total_bms);
            if (cost >= 0)
                step_function->add_point(pos, current_rms, total_bms, cost);

            current_rms += get_mutations_represented(pos);
            i += 1;
        }
        while (j < req_back_mutations.size())
        {
            // So j must equal bms
            pos = req_back_mutations[j];

            // update prefix and suffix costs
            replace_cost_if_better(arg, prefix_costs, pos, edge, total_rms, current_bms);
            replace_cost_if_better(arg, suffix_costs, pos, edge, 0, total_bms - current_bms);

            float cost = get_cost_of_addition(arg, total_rms, current_bms);
            if (cost >= 0)
                step_function->add_point(pos, total_rms, current_bms, cost);

            current_bms += get_mutations_represented(pos);
            j += 1;
        }

        replace_cost_if_better(arg, suffix_costs, pos + 1, edge, 0, 0);
        // replace_cost_if_better(arg, prefix_costs, pos + 1, edge, rms, bms);

        edge_cost_functions.insert({edge, std::move(step_function)});
    }

    // Now we can go through the costs (which should be sorted by key)
    // First turn maps into vectors so that we can iterate through both nicely
    std::vector<std::tuple<int, Edge *, int, int>> prefix_costs_vec;
    std::vector<std::tuple<int, Edge *, int, int>> suffix_costs_vec;
    for (const auto &[site, value] : prefix_costs)
    {
        auto [edge, rms, bms] = value;
        prefix_costs_vec.push_back(std::make_tuple(site, edge, rms, bms));
    }
    for (const auto &[site, value] : suffix_costs)
    {
        auto [edge, rms, bms] = value;
        suffix_costs_vec.push_back(std::make_tuple(site, edge, rms, bms));
    }

    // Now go through backwards to create optimal prefix costs, and forwards for optimal suffix costs
    std::vector<std::tuple<int, Edge *, int, int>> optimal_prefix_costs;
    std::vector<std::tuple<int, Edge *, int, int>> optimal_suffix_costs;
    float best_cost = -1.0;
    for (int i = prefix_costs_vec.size() - 1; i >= 0; i--)
    {
        auto [site, edge, rms, bms] = prefix_costs_vec[i];
        float cost = get_cost_of_addition(arg, rms, bms);
        if (cost >= 0 && (cost < best_cost || best_cost < 0))
        {
            optimal_prefix_costs.push_back(prefix_costs_vec[i]);
            best_cost = cost;
        }
    }
    std::reverse(optimal_prefix_costs.begin(), optimal_prefix_costs.end());
    best_cost = -1.0;
    for (int i = 0; i < suffix_costs_vec.size(); i++)
    {
        auto [site, edge, rms, bms] = suffix_costs_vec[i];
        float cost = get_cost_of_addition(arg, rms, bms);
        if (cost >= 0 && (cost < best_cost || best_cost < 0))
        {
            optimal_suffix_costs.push_back(suffix_costs_vec[i]);
            best_cost = cost;
        }
    }

    if (optimal_prefix_costs.size() == 0 or optimal_suffix_costs.size() == 0)
    {
        return std::make_tuple(std::move(locations), best_score);
    }

    std::vector<int> change_points; // This only needs to be calculated if we have three of more crossover positions
    if (_max_number_parents >= 4)
    {
        int p_index = 0, s_index = 0;
        while (p_index < optimal_prefix_costs.size() && s_index < optimal_suffix_costs.size())
        {
            auto [p_site, x, y, z] = optimal_prefix_costs[p_index];
            auto [s_site, a, b, c] = optimal_suffix_costs[s_index];
            if (p_site < s_site)
            {
                p_index += 1;
                change_points.push_back(p_site);
            }
            else if (p_site == s_site)
            {
                p_index += 1;
                s_index += 1;
                change_points.push_back(p_site);
            }
            else
            {
                s_index += 1;
                change_points.push_back(s_site);
            }
        }

        while (p_index < optimal_prefix_costs.size())
        {
            auto [p_site, x, y, z] = optimal_prefix_costs[p_index];
            change_points.push_back(p_site);
            p_index += 1;
        }
        while (s_index < optimal_suffix_costs.size())
        {
            auto [s_site, x, y, z] = optimal_suffix_costs[s_index];
            change_points.push_back(s_site);
            s_index += 1;
        }
    }

    int best_number_parents = 0;
    std::vector<int> best_positions;
    std::vector<Edge *> best_parents;
    best_cost = -1.0;
    int best_rms = -1;
    int best_bms = -1;

    if (_max_number_parents >= 2)
    {
        int p_index = 0, s_index = 0;
        auto [p_site, p_edge, p_rms, p_bms] = optimal_prefix_costs[p_index];
        auto [s_site, s_edge, s_rms, s_bms] = optimal_suffix_costs[s_index];

        while (p_index < optimal_prefix_costs.size())
        {
            // Now need to find the highest s_index such that s_site <= p_site
            while (s_index < optimal_suffix_costs.size() - 1)
            {
                int next_site = std::get<0>(optimal_suffix_costs[s_index + 1]);
                if (next_site <= p_site)
                {
                    s_index += 1;
                    std::tie(s_site, s_edge, s_rms, s_bms) = optimal_suffix_costs[s_index];
                }
                else
                {
                    break;
                }
            }

            float current_cost = get_cost_of_addition(arg, p_rms + s_rms, p_bms + s_bms, 1);
            locations.push_back(std::make_tuple(2, std::vector<int>({p_site}), std::vector<Edge *>({p_edge, s_edge}), current_cost, p_rms + s_rms, p_bms + s_bms));
            if (current_cost >= 0 && (current_cost < best_cost || best_cost < 0))
            {
                best_number_parents = 2;
                best_positions = {p_site};
                best_parents = {p_edge, s_edge};
                best_cost = current_cost;
                best_rms = p_rms + s_rms;
                best_bms = p_bms + s_bms;
            }

            p_index++;
            std::tie(p_site, p_edge, p_rms, p_bms) = optimal_prefix_costs[p_index];
        }
    }

    if (_max_number_parents >= 3 && get_cost_of_addition(arg, 0, 0, 2) >= 0 && (best_cost < 0 || get_cost_of_addition(arg, 0, 0, 2) < best_cost))
    {
        int p_index = 0, s_index = 0;
        auto [p_site, p_edge, p_rms, p_bms] = optimal_prefix_costs[p_index];
        auto [s_site, s_edge, s_rms, s_bms] = optimal_suffix_costs[s_index];

        while (p_index < optimal_prefix_costs.size())
        {
            std::tie(p_site, p_edge, p_rms, p_bms) = optimal_prefix_costs[p_index];

            s_index = optimal_suffix_costs.size() - 1;
            while (s_index >= 0)
            {
                std::tie(s_site, s_edge, s_rms, s_bms) = optimal_suffix_costs[s_index];

                if (s_site <= p_site)
                {
                    break; // This would just be a 2 parent recombination
                }
                else
                {
                    for (Edge *m_edge : related_edges)
                    {
                        auto [m_rms, m_bms] = edge_cost_functions[m_edge]->get_range(p_site, s_site);
                        int current_rms = p_rms + m_rms + s_rms;
                        int current_bms = p_bms + m_bms + s_bms;
                        float current_cost = get_cost_of_addition(arg, current_rms, current_bms, 2);
                        locations.push_back(std::make_tuple(3, std::vector<int>({p_site, s_site}), std::vector<Edge *>({p_edge, m_edge, s_edge}), current_cost, current_rms, current_bms));
                        if (current_cost >= 0 && (current_cost < best_cost || best_cost < 0))
                        {
                            best_number_parents = 3;
                            best_positions = {p_site, s_site};
                            best_parents = {p_edge, m_edge, s_edge};
                            best_cost = current_cost;
                            best_rms = current_rms;
                            best_bms = current_bms;
                        }
                    }
                }

                s_index--;
            }

            p_index++;
        }
    }

    if (_max_number_parents >= 4 && get_cost_of_addition(arg, 0, 0, 3) >= 0 && (best_cost < 0 || get_cost_of_addition(arg, 0, 0, 3) < best_cost))
    {
        int p_index = 0, s_index = 0;
        auto [p_site, p_edge, p_rms, p_bms] = optimal_prefix_costs[p_index];
        auto [s_site, s_edge, s_rms, s_bms] = optimal_suffix_costs[s_index];

        while (p_index < optimal_prefix_costs.size())
        {
            std::tie(p_site, p_edge, p_rms, p_bms) = optimal_prefix_costs[p_index];

            s_index = optimal_suffix_costs.size() - 1;
            while (s_index >= 0)
            {
                std::tie(s_site, s_edge, s_rms, s_bms) = optimal_suffix_costs[s_index];

                if (s_site <= p_site)
                {
                    break; // This would just be a 2 parent recombination
                }
                else
                {
                    for (int midpoint : change_points)
                    {
                        if (midpoint <= p_site)
                            continue;
                        else if (midpoint >= s_site)
                            break;

                        Edge *best_middle1;
                        float best_cost1 = -1.0;
                        Edge *best_middle2;
                        float best_cost2 = -1.0;

                        for (Edge *m_edge : related_edges)
                        {
                            auto [m1_rms, m1_bms] = edge_cost_functions[m_edge]->get_range(p_site, midpoint);
                            float cost1 = get_cost_of_addition(arg, m1_rms, m1_bms);
                            if (cost1 >= 0 && (cost1 < best_cost1 || best_cost1 < 0))
                            {
                                best_cost1 = cost1;
                                best_middle1 = m_edge;
                            }

                            auto [m2_rms, m2_bms] = edge_cost_functions[m_edge]->get_range(midpoint, s_site);
                            float cost2 = get_cost_of_addition(arg, m2_rms, m2_bms);
                            if (cost2 >= 0 && (cost2 < best_cost2 || best_cost2 < 0))
                            {
                                best_cost2 = cost2;
                                best_middle2 = m_edge;
                            }
                        }

                        auto [m1_rms, m1_bms] = edge_cost_functions[best_middle1]->get_range(p_site, midpoint);
                        auto [m2_rms, m2_bms] = edge_cost_functions[best_middle2]->get_range(midpoint, s_site);
                        int current_rms = p_rms + m1_rms + m2_rms + s_rms;
                        int current_bms = p_bms + m1_bms + m2_bms + s_bms;
                        float current_cost = get_cost_of_addition(arg, current_rms, current_bms, 3);
                        locations.push_back(std::make_tuple(4, std::vector<int>({p_site, midpoint, s_site}), std::vector<Edge *>({p_edge, best_middle1, best_middle2, s_edge}), current_cost, current_rms, current_bms));
                        if (current_cost >= 0 && (current_cost < best_cost || best_cost < 0))
                        {
                            best_number_parents = 4;
                            best_positions = {p_site, midpoint, s_site};
                            best_parents = {p_edge, best_middle1, best_middle2, s_edge};
                            best_cost = current_cost;
                            best_rms = current_rms;
                            best_bms = current_bms;
                        }
                    }
                }

                s_index--;
            }

            p_index++;
        }
    }

    if (best_cost < 0)
    {
        if (_how_verbose >= 2)
            std::cout << "Could not find recombination location for g: " << g.label << "\n";
        return std::make_tuple(std::move(locations), best_score);
    }
    else
    {
        if (_how_verbose >= 2)
        {
            std::cout << "Best recomb has " << best_number_parents << " parents: ";
            for (auto p : best_parents)
                std::cout << p->to->label << ", ";
            std::cout << ". Cross-over positions: " << vector_to_string(best_positions) << ".";
            std::cout << " Cost was rms: " << best_rms << " bms: " << best_bms << "\n";
        }
        return std::make_tuple(std::move(locations), best_score);
    }
}

void insert_seq_at_recomb_location(ARG &arg, const Gene &g, const int pos, Edge *prefix, Edge *suffix, const std::vector<int> &existing_mutations)
{
    auto p_existing_mutations = vector_values_below(existing_mutations, pos);
    auto s_existing_mutations = vector_values_above(existing_mutations, pos);

    // Check if we should split prefix
    std::vector<int> prefix_rms = vector_difference(p_existing_mutations, prefix->to->mutations);
    std::vector<int> prefix_bms = vector_difference(vector_values_below(prefix->to->mutations, pos),
                                                    p_existing_mutations);

    auto [prefix_removable_rms, prefix_required_rms] = vector_split(prefix_rms, prefix->back_mutations);
    auto [prefix_removable_bms, prefix_required_bms] = vector_split(prefix_bms, prefix->mutations);

    Node *prefix_node; // This will either be the target of the prefix edge, or a split node
    if (prefix_removable_rms.size() > 0 || prefix_removable_bms.size() > 0)
    {
        prefix_node = split_edge(arg, prefix, prefix_removable_rms, prefix_removable_bms);
    }
    else
    {
        prefix_node = prefix->to;
    }

    // Now the same for suffix
    std::vector<int> suffix_rms = vector_difference(s_existing_mutations, suffix->to->mutations);
    std::vector<int> suffix_bms = vector_difference(vector_values_above(suffix->to->mutations, pos),
                                                    s_existing_mutations);

    auto [suffix_removable_rms, suffix_required_rms] = vector_split(suffix_rms, suffix->back_mutations);
    auto [suffix_removable_bms, suffix_required_bms] = vector_split(suffix_bms, suffix->mutations);

    Node *suffix_node; // This will either be the target of the suffix edge, or a split node
    if (suffix_removable_rms.size() > 0 || suffix_removable_bms.size() > 0)
    {
        suffix_node = split_edge(arg, suffix, suffix_removable_rms, suffix_removable_bms);
    }
    else
    {
        suffix_node = suffix->to;
    }

    Node *recomb_node = recombine_nodes(arg, pos, prefix_node, suffix_node);
    insert_seq_as_direct_child(arg, g, recomb_node, existing_mutations);
}

void insert_seq_at_multi_recomb_location(ARG &arg, const Gene &g, const int number_parents, std::vector<int> crossovers, std::vector<Edge *> parents, const std::vector<int> &existing_mutations)
{
    std::vector<Node *> parent_nodes;
    for (int i = 0; i < number_parents; i++)
    {
        Edge *parent = parents[i];
        // For each parent decide whether to split
        std::vector<int> existing_mutations_in_range = {};
        std::vector<int> parent_mutations_in_range = {};

        if (i == 0)
        {
            existing_mutations_in_range = vector_values_below(existing_mutations, crossovers[0]);
            parent_mutations_in_range = vector_values_below(parent->to->mutations, crossovers[0]);
        }
        else if (i == number_parents - 1)
        {
            existing_mutations_in_range = vector_values_above(existing_mutations, crossovers.back());
            parent_mutations_in_range = vector_values_above(parent->to->mutations, crossovers.back());
        }
        else
        {
            existing_mutations_in_range = vector_values_between(existing_mutations, crossovers[i - 1], crossovers[i]);
            parent_mutations_in_range = vector_values_between(parent->to->mutations, crossovers[i - 1], crossovers[i]);
        }

        std::vector<int> rms = vector_difference(existing_mutations_in_range, parent->to->mutations);
        std::vector<int> bms = vector_difference(parent_mutations_in_range, existing_mutations_in_range);

        auto [removable_rms, required_rms] = vector_split(rms, parent->back_mutations);
        auto [removable_bms, required_bms] = vector_split(bms, parent->mutations);

        Node *parent_node; // This will either be the target of the prefix edge, or a split node
        if (removable_rms.size() > 0 || removable_bms.size() > 0)
        {
            parent_node = split_edge(arg, parent, removable_rms, removable_bms);
        }
        else
        {
            parent_node = parent->to;
        }

        parent_nodes.push_back(std::move(parent_node));
    }

    Node *recomb_node = parent_nodes[0];
    for (int i = 0; i < number_parents - 1; i++)
    {
        recomb_node = recombine_nodes(arg, crossovers[i], recomb_node, parent_nodes[i + 1]);
    }

    insert_seq_as_direct_child(arg, g, recomb_node, existing_mutations);
}

void add_seq_to_arg(ARG &arg, const Gene &g)
{
    auto [relevant_edges, existing_mutations] = find_relevant_edges(arg, g);

    auto [edge, rms, bms] = find_best_single_parent_location(arg, g, relevant_edges, existing_mutations);
    // nullptr for edge means that no single parent works. A self loop would be returned for a root option.
    float sp_cost = get_cost_of_addition(arg, rms, bms);

    if (_cost_recombination >= 0)
    {
        auto [crossover_pos, prefix_edge, suffix_edge, r_rms, r_bms] = find_best_recomb_location(arg, g, relevant_edges, existing_mutations);
        // crossover_pos = -1 means failure, 0 means all one sequence

        float r_cost = get_cost_of_addition(arg, r_rms, r_bms, 1);

        if (crossover_pos > 0 && r_cost >= 0 && (sp_cost < 0 || r_cost < sp_cost))
        {
            // recombination is valid and better than having single parent
            insert_seq_at_recomb_location(arg, g, crossover_pos, prefix_edge, suffix_edge, existing_mutations);
            return;
        }
    }

    if (edge != nullptr && sp_cost >= 0)
    {
        insert_seq_at_edge(arg, g, edge, existing_mutations);
    }
    else
    {
        std::cerr << "The gene: " << g.label << " could not be added to the ARG so will be ignored.\n";
    }
}

void add_seq_to_arg_with_multi_recombs(ARG &arg, const Gene &g)
{
    auto [relevant_edges, existing_mutations] = find_relevant_edges(arg, g);

    auto [edge, rms, bms] = find_best_single_parent_location(arg, g, relevant_edges, existing_mutations);
    // nullptr for edge means that no single parent works. A self loop would be returned for a root option.
    float sp_cost = get_cost_of_addition(arg, rms, bms);

    if (_cost_recombination >= 0)
    {
        auto [best_number_parents, crossovers, parents, r_rms, r_bms] = find_best_multi_recomb_location(arg, g, relevant_edges, existing_mutations);
        // best_number_parents = 0 means failure

        float r_cost = get_cost_of_addition(arg, r_rms, r_bms, best_number_parents - 1);

        if (best_number_parents >= 2 && r_cost >= 0 && (sp_cost < 0 || r_cost < sp_cost))
        {
            // recombination is valid and better than having single parent
            insert_seq_at_multi_recomb_location(arg, g, best_number_parents, crossovers, parents, existing_mutations);
            return;
        }
    }

    if (edge != nullptr && sp_cost >= 0)
    {
        insert_seq_at_edge(arg, g, edge, existing_mutations);
    }
    else
    {
        std::cerr << "The gene: " << g.label << " could not be added to the ARG so will be ignored.\n";
    }
}

void add_seq_to_arg_with_randomness(ARG &arg, const Gene &g, float delta, int random_number)
{
    auto [relevant_edges, existing_mutations] = find_relevant_edges(arg, g);

    std::vector<int> all_scores;

    auto [scored_single_parent_locations, best_sp_cost] = score_single_parent_locations(arg, g, relevant_edges, existing_mutations);
    int number_of_single_parent_locations = scored_single_parent_locations.size();

    for (const auto &location : scored_single_parent_locations)
    {
        all_scores.push_back(std::get<1>(location));
    }

    std::vector<std::tuple<int, std::vector<int>, std::vector<Edge *>, int, int, int>> scored_multi_parent_locations;
    float best_recomb_cost = -1.0;
    if (_cost_recombination >= 0)
    {
        std::tie(scored_multi_parent_locations, best_recomb_cost) = score_multi_recomb_locations(arg, g, relevant_edges, existing_mutations);
    }

    int number_of_multi_parent_locations = scored_multi_parent_locations.size();
    for (const auto &location : scored_multi_parent_locations)
    {
        all_scores.push_back(std::get<3>(location));
    }

    if (_how_verbose >= 2)
    {
        std::cout << "Neighbourhood created with " << all_scores.size() << " locations.\n";
    }

    if (all_scores.size() == 0)
    {
        std::cerr << "The gene: " << g.label << " could not be added to the ARG so will be ignored.\n";
        return;
    }

    float best_cost = -1.0;
    if (best_sp_cost >= 0)
        best_cost = best_sp_cost;
    if (best_recomb_cost >= 0 && (best_recomb_cost < best_cost || best_cost < 0))
        best_cost = best_recomb_cost;

    std::vector<int> decent_options; // list of indexes for locations within delta of best

    for (int i = 0; i < all_scores.size(); i++)
    {
        if (all_scores[i] <= best_cost + delta)
            decent_options.push_back(i);
    }

    if (_how_verbose >= 2)
    {
        std::cout << "Neighbourhood reduced to " << decent_options.size() << " locations within delta.\n";
    }

    // Now randomly select a location
    int location_index = decent_options[random_number % decent_options.size()];

    if (location_index < number_of_single_parent_locations)
    {
        auto [edge, cost, rms, bms] = scored_single_parent_locations[location_index];
        insert_seq_at_edge(arg, g, edge, existing_mutations);
    }
    else
    {
        location_index -= number_of_single_parent_locations;
        auto [num_parents, crossovers, parents, cost, rms, bms] = scored_multi_parent_locations[location_index];
        insert_seq_at_multi_recomb_location(arg, g, num_parents, crossovers, parents, existing_mutations);
    }
}

/* Functions used to set up before building arg */

/* flips all sites which are mutated in root */
void set_genes_relative_to_root(Genes &genes, const Gene &root)
{
    for (Gene &g : genes.genes)
    {
        g.mutations = vector_symmetric_difference(g.mutations, root.mutations);
    }
}

void set_arg_relative_to_root(ARG &arg, const Gene &root)
{
    for (auto &node : arg.nodes)
    {
        node->mutations = vector_symmetric_difference(node->mutations, root.mutations);
    }
}

Gene create_random_root(int seq_length, int seed)
{
    Gene root;
    root.label = "guessed root";
    auto rng = std::default_random_engine(seed);

    for (int i = 0; i < seq_length; i++)
    {
        if (rng() % 2 == 0)
        {
            root.mutations.push_back(i);
        }
    }

    return std::move(root);
}

ARG _build_arg_from_order(const Genes genes, int roots_given, const std::vector<int> &order)
{
    clock_t tic = clock();

    ARG arg = roots_given > 0 ? ARG(false) : ARG(true);

    for (int root_num = 0; root_num < roots_given; root_num++)
    {
        arg.add_root(genes.genes[root_num].label, genes.genes[root_num].mutations);
    }
    // Note that roots should not be included in order

    auto rng = std::default_random_engine(_run_seed);

    int step = 0;
    double delta = 0;
    if (_location_selection_method == 1)
        delta = 0.0;
    else if (_location_selection_method == 2)
        delta = 0.5;
    else if (_location_selection_method == 3)
        delta = 1.0;

    for (int index : order)
    {
        if (_location_selection_method == 0)
            add_seq_to_arg_with_multi_recombs(arg, genes.genes[index]);
        else
            add_seq_to_arg_with_randomness(arg, genes.genes[index], delta, rng());

        step += 1;
        if (_how_verbose >= 3)
        {
            std::string filename = "network_progress";
            filename += std::to_string(step) + ".dot";
            auto fp = fopen(filename.c_str(), "w");
            arg_output(arg, genes, fp, ARGDOT, true, LABEL_BOTH);
            fclose(fp);
        }
    }

    clock_t toc = clock();
    double timer = (double)(toc - tic) / CLOCKS_PER_SEC;

    // Add arg to record
    _run_record.add_record(arg.number_of_recombinations, arg.number_of_recurrent_mutations, arg.number_of_back_mutations,
                           _run_seed, timer, _cost_recombination, _cost_recurrent_mutation, _cost_back_mutation);

    if (_run_record_file != "")
    {
        Run run = _run_record.runs.back();

        std::fstream fout;
        fout.open(_run_record_file, std::ios::out | std::ios::app);

        fout << run.recombinations << "," << run.back_mutations << "," << run.recurrent_mutations << ","
             << std::to_string(run.seed) << "," << std::to_string(run.build_time) << "," << run.recomb_cost << ","
             << run.rm_cost << "," << run.bm_cost << "\n";
    }

    return std::move(arg);
}

ARG _run_multi_runs(const Genes genes, int number_of_runs, int run_seed, int multi_run_strategy, int roots_given)
{
    // Now proceeding with performing runs
    bool first_run = true;

    std::vector<int> prev_order;
    std::set<int> prev_tricky_sites;

    float best_cost = -1.0;
    ARG best_arg;
    std::vector<int> best_order;
    for (int run_count = 0; run_count < number_of_runs; run_count++)
    {
        auto rng = std::default_random_engine(run_seed);
        _run_seed = run_seed; // Used for record keeping

        std::vector<int> order;
        if (multi_run_strategy == 3)
        {
            // Order from fewest mutation sequences to highest.
            int max_muts = 0;
            for (const auto &g : genes.genes)
            {
                if (g.mutations.size() > max_muts)
                    max_muts = g.mutations.size();
            }

            std::vector<std::vector<int>> num_mutations_to_genes;
            for (int i = 0; i <= max_muts; i++)
            {
                num_mutations_to_genes.push_back(std::vector<int>());
            }
            for (int i = roots_given; i < genes.genes.size(); i++)
            {
                int muts = genes.genes[i].mutations.size();
                num_mutations_to_genes[muts].push_back(i);
            }

            for (int i = 0; i <= max_muts; i++)
            {
                // shuffle all sequences with same number of mutations
                std::shuffle(num_mutations_to_genes[i].begin(), num_mutations_to_genes[i].end(), rng);
                order.insert(order.end(), num_mutations_to_genes[i].begin(), num_mutations_to_genes[i].end());
            }
        }
        else if (first_run || multi_run_strategy == 0)
        {
            for (int i = roots_given; i < genes.genes.size(); i++)
                order.push_back(i);
            if (run_seed != 0) // run seed being zero means don't shuffle
                std::shuffle(order.begin(), order.end(), rng);
        }
        else if (multi_run_strategy == 1 || multi_run_strategy == 2)
        {
            // Deal with tricky sites later
            std::vector<int> tricky_order;
            std::vector<int> easy_order;

            for (int i = roots_given; i < genes.genes.size(); i++)
            {
                auto intersection = vector_set_intersection(genes.genes[i].mutations, prev_tricky_sites);
                if (intersection.size() > 0)
                    tricky_order.push_back(i);
                else
                    easy_order.push_back(i);
            }

            std::shuffle(tricky_order.begin(), tricky_order.end(), rng);
            std::shuffle(easy_order.begin(), easy_order.end(), rng);

            order.clear();
            if (multi_run_strategy == 1)
            {
                order.insert(order.end(), easy_order.begin(), easy_order.end());
                order.insert(order.end(), tricky_order.begin(), tricky_order.end());
            }
            else
            {
                order.insert(order.end(), tricky_order.begin(), tricky_order.end());
                order.insert(order.end(), easy_order.begin(), easy_order.end());
            }
        }

        ARG run_arg = _build_arg_from_order(genes, roots_given, order);

        float run_cost = get_cost(run_arg);
        if (_how_verbose >= 1)
        {
            std::cout << "ARG " << run_count << " (seed: " << run_seed << ") made using " << run_arg.number_of_recurrent_mutations << " recurrent mutations, " << run_arg.number_of_back_mutations;
            std::cout << " back mutations, and " << run_arg.number_of_recombinations << " recombinations\n";
            std::cout << "Made at cost: " << run_cost << ".\n";
        }

        if (run_cost < best_cost || best_cost < 0)
        {
            best_arg = std::move(run_arg);
            best_cost = run_cost;
            best_order = order;
        }

        run_seed = rng();
        first_run = false;
        prev_order = order;
        prev_tricky_sites = run_arg.recurrent_sites;
        for (int m : run_arg.back_mutation_sites)
            prev_tricky_sites.insert(m);
    }

    return std::move(best_arg);
}

ARG _build_arg_multi_runs(const Genes input_genes, bool clean_sequences, int number_of_runs, int run_seed, int roots_given, int multi_run_strategy)
{
    if (roots_given == 1)
    {
        std::cerr << "This should not be called with a single root. This case should be dealth with in main function.\n";
    }

    Genes genes; // copy input genes
    genes.sequence_length = input_genes.sequence_length;
    genes.genes.insert(genes.genes.end(), input_genes.genes.begin() + 1, input_genes.genes.end());

    // Check if cleaning should be done first
    std::vector<HistoryStep> history;
    std::vector<int> site_extent;
    if (roots_given > 1 && clean_sequences)
        std::cerr << "Can't perform clean algorithm with multiple roots given.\n";
    else if (clean_sequences)
        std::tie(history, _site_multiplicitity, site_extent) = clean_genes(genes, _how_verbose);

    ARG best_arg = _run_multi_runs(genes, number_of_runs, run_seed, multi_run_strategy, roots_given);

    if (_how_verbose >= 1)
    {
        std::cout << "Final ARG made using " << best_arg.number_of_recurrent_mutations << " recurrent mutations, " << best_arg.number_of_back_mutations;
        std::cout << " back mutations, and " << best_arg.number_of_recombinations << " recombinations\n";
        std::cout << "Made at cost: " << get_cost(best_arg) << ".\n";
    }

    if (clean_sequences && roots_given <= 1)
        extend_arg_with_history(best_arg, history, site_extent);

    return std::move(best_arg);
}

std::tuple<ARG, RunRecord> build_arg_main(const Genes genes, bool clean_sequences, int how_verbose, int roots_given, int run_seed, int number_of_runs, int multi_run_strategy, int location_selection_method, int find_root_strategy, int find_root_iterations,
                                          int max_number_parents, float cost_rm, float cost_bm, float cost_recomb, int recomb_max, int rm_max, int bm_max,
                                          std::string run_record_file)
{
    _multi_roots_given = roots_given > 1;
    _how_verbose = how_verbose;
    _max_number_parents = max_number_parents;
    _cost_recurrent_mutation = cost_rm;
    _cost_back_mutation = cost_bm;
    _cost_recombination = cost_recomb;
    _recomb_max = recomb_max < 0 ? INT32_MAX : recomb_max;
    _rm_max = rm_max < 0 ? INT32_MAX : rm_max;
    _bm_max = bm_max < 0 ? INT32_MAX : bm_max;

    std::vector<int> default_site_multiplicity(genes.sequence_length, 1);
    _site_multiplicitity = default_site_multiplicity;

    _run_record.clear();
    _location_selection_method = location_selection_method;
    _run_record_file = run_record_file;

    if (how_verbose >= 1)
        std::cout << "Starting to build arg\n";

    if (roots_given == 1)
    {
        Gene root = genes.genes[0];

        Genes genes_copy;
        genes_copy.sequence_length = genes.sequence_length;
        genes_copy.genes.insert(genes_copy.genes.end(), genes.genes.begin() + 1, genes.genes.end());

        set_genes_relative_to_root(genes_copy, root);

        ARG arg = _build_arg_multi_runs(genes_copy, clean_sequences, number_of_runs, run_seed, 0, multi_run_strategy);

        set_arg_relative_to_root(arg, root);

        return std::make_tuple(std::move(arg), std::move(_run_record));
    }
    else if (roots_given == 0 && find_root_strategy == 0)
    {
        ARG arg = _build_arg_multi_runs(genes, clean_sequences, number_of_runs, run_seed, 0, multi_run_strategy);
        return std::make_tuple(std::move(arg), std::move(_run_record));
    }
    else if (roots_given > 1)
    {
        ARG arg = _build_arg_multi_runs(genes, clean_sequences, number_of_runs, run_seed, roots_given, multi_run_strategy);
        return std::make_tuple(std::move(arg), std::move(_run_record));
    }
    else if (find_root_strategy == 1)
    {
        // Simply guess root at random
        auto rng = std::default_random_engine(run_seed);
        Gene best_root;
        float best_cost = -1.0;
        ARG best_arg;
        for (int i = 0; i < find_root_iterations; i++)
        {
            Gene root = create_random_root(genes.sequence_length, rng());
            if (_how_verbose >= 1)
                std::cout << "trying root: " << vector_to_string(root.mutations) << "\n";

            Genes genes_copy;
            genes_copy.sequence_length = genes.sequence_length;
            genes_copy.genes.insert(genes_copy.genes.end(), genes.genes.begin(), genes.genes.end());

            set_genes_relative_to_root(genes_copy, root);

            ARG arg = _build_arg_multi_runs(genes_copy, clean_sequences, number_of_runs, run_seed, 0, multi_run_strategy);

            set_arg_relative_to_root(arg, root);

            float run_cost = get_cost(arg);

            if (run_cost < best_cost || best_cost < 0)
            {
                best_arg = std::move(arg);
                best_cost = run_cost;
                best_root = root;
            }
        }
        return std::make_tuple(std::move(best_arg), std::move(_run_record));
    }
    else if (find_root_strategy == 2 || find_root_strategy == 3)
    {
        // Compare all zero root to root having a 1 at each location

        Genes genes_copy;
        genes_copy.sequence_length = genes.sequence_length;
        genes_copy.genes.insert(genes_copy.genes.end(), genes.genes.begin(), genes.genes.end());
        ARG empty_root_arg = _build_arg_multi_runs(genes_copy, clean_sequences, number_of_runs, run_seed, 0, multi_run_strategy);
        float best_score = get_cost(empty_root_arg);

        Gene result_root;

        for (int i = 0; i < genes.sequence_length; i++)
        {
            Gene root;
            if (find_root_strategy == 2)
                root.mutations = {i};
            else
                root.mutations = vector_union(result_root.mutations, {i});

            if (_how_verbose >= 1)
                std::cout << "trying root: " << vector_to_string(root.mutations) << "\n";

            Genes genes_copy;
            genes_copy.sequence_length = genes.sequence_length;
            genes_copy.genes.insert(genes_copy.genes.end(), genes.genes.begin(), genes.genes.end());

            set_genes_relative_to_root(genes_copy, root);

            ARG arg = _build_arg_multi_runs(genes_copy, clean_sequences, number_of_runs, run_seed, 0, multi_run_strategy);

            set_arg_relative_to_root(arg, root);

            float run_cost = get_cost(arg);

            if (run_cost >= 0 && run_cost < best_score)
            {
                result_root.mutations.push_back(i);
                if (_how_verbose >= 1)
                    std::cout << "root mutation at " << i << " is advantagous.\n";

                if (find_root_strategy == 3)
                    best_score = run_cost;
            }
        }

        genes_copy.genes.clear();
        genes_copy.genes.insert(genes_copy.genes.end(), genes.genes.begin(), genes.genes.end());

        set_genes_relative_to_root(genes_copy, result_root);

        ARG arg = _build_arg_multi_runs(genes_copy, clean_sequences, number_of_runs, run_seed, 0, multi_run_strategy);

        set_arg_relative_to_root(arg, result_root);

        if (_how_verbose >= 1)
            std::cout << "result root: " << vector_to_string(result_root.mutations) << "\n";

        return std::make_tuple(std::move(arg), std::move(_run_record));
    }
}

std::tuple<ARG, RunRecord> build_arg_search(const Genes genes, bool clean_sequences, int how_verbose, int roots_given, int run_seed, int number_of_runs, int multi_run_strategy, int location_selection_method,
                                            int max_number_parents, std::vector<float> costs_rm, std::vector<float> costs_bm, std::vector<float> costs_recomb,
                                            std::string run_record_file)
{
    // This is used for trialling with many different costs. Works only for no root or 1 root.
    // Not really search for best arg, but for run record
    _multi_roots_given = roots_given > 1;
    _how_verbose = how_verbose;
    _max_number_parents = max_number_parents;

    _recomb_max = INT32_MAX;
    _rm_max = INT32_MAX;
    _bm_max = INT32_MAX;

    std::vector<int> default_site_multiplicity(genes.sequence_length, 1);
    _site_multiplicitity = default_site_multiplicity;

    _run_record.clear();
    _location_selection_method = location_selection_method;
    _run_record_file = run_record_file;

    if (how_verbose >= 1)
        std::cout << "Starting arg search\n";

    if (roots_given > 1)
        std::cout << "Can't run a search with more than 1 root";

    Gene root;
    Genes genes_copy;
    genes_copy.sequence_length = genes.sequence_length;

    if (roots_given == 1)
    {
        Gene root = genes.genes[0];
        genes_copy.genes.insert(genes_copy.genes.end(), genes.genes.begin() + 1, genes.genes.end());

        set_genes_relative_to_root(genes_copy, root);
    }
    else
    {
        genes_copy.genes.insert(genes_copy.genes.end(), genes.genes.begin(), genes.genes.end());
    }

    std::vector<HistoryStep> history;
    std::vector<int> site_extent;
    if (clean_sequences)
        std::tie(history, _site_multiplicitity, site_extent) = clean_genes(genes_copy, _how_verbose);

    ARG arg;
    int steps = costs_rm.size() * costs_bm.size() * costs_recomb.size();
    int step = 0;
    int barWidth = 70;
    for (float rm_cost : costs_rm)
    {
        _cost_recurrent_mutation = rm_cost;
        for (float bm_cost : costs_bm)
        {
            if (bm_cost < 0)
                bm_cost = -1;
            else
                _cost_back_mutation = rm_cost + bm_cost;
            for (float recomb_cost : costs_recomb)
            {
                _cost_recombination = recomb_cost;

                arg = _run_multi_runs(genes_copy, number_of_runs, run_seed, multi_run_strategy, 0);
                step += 1;
                if (_how_verbose >= 0)
                {
                    float progress = static_cast<float>(step) / static_cast<float>(steps);
                    std::cout << "[";
                    int pos = barWidth * progress;
                    for (int i = 0; i < barWidth; ++i)
                    {
                        if (i < pos)
                            std::cout << "=";
                        else if (i == pos)
                            std::cout << ">";
                        else
                            std::cout << " ";
                    }
                    std::cout << "] " << int(progress * 100.0) << " %\r";
                    std::cout.flush();
                }
                run_seed = _run_seed;
            }
        }
    }

    if (_how_verbose >= 0)
        std::cout << std::endl;

    if (roots_given == 1)
    {
        set_arg_relative_to_root(arg, root);
    }

    return std::make_tuple(std::move(arg), std::move(_run_record));
}