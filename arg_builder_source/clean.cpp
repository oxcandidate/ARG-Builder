#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <iostream>
#include <random>
#include <algorithm>

#include "clean.h"
#include "arg_builder_logic.h"
#include "vector_set_operations.h"

/* Cleans the genes and returns tuple with history, site multiplicity, and site_extent */
std::tuple<std::vector<HistoryStep>, std::vector<int>, std::vector<int>> clean_genes(Genes &genes, int how_verbose)
{
    std::vector<Gene> gs;
    gs.insert(gs.end(), genes.genes.begin(), genes.genes.end());

    std::vector<HistoryStep> history;

    std::vector<int> remaining_columns;
    std::vector<int> site_extent; // How far does site i cover
    for (int i = 0; i < genes.sequence_length; i++)
    {
        remaining_columns.push_back(i);
        site_extent.push_back(i);
    }
    std::vector<int> column_counts(genes.sequence_length, 0);     // number of mutations per site
    std::vector<int> site_multiplicity(genes.sequence_length, 1); // when merging columns, this tracks their multiplicity

    for (const auto &g : gs)
    {
        for (int site : g.mutations)
            column_counts[site] += 1;
    }

    bool changed = true;
    while (changed)
    {
        changed = false;

        std::vector<int> still_remaining_columns;
        // Check for uninformative columns
        for (int site : remaining_columns)
        {
            if (column_counts[site] == 0)
            {
                if (how_verbose >= 2)
                    std::cout << "site " << site << " is empty\n";
                HistoryStep step;
                step.type = RemoveEmptyCol;
                step.empty_col.site = site;
                history.push_back(step);
                changed = true;
            }
            else if (column_counts[site] == 1)
            {
                if (how_verbose >= 2)
                    std::cout << "site " << site << " has just one mutations\n";
                // Need to find gene with this mutation
                for (auto &g : gs)
                {
                    auto search = std::find(g.mutations.begin(), g.mutations.end(), site);
                    if (search != g.mutations.end())
                    {
                        g.mutations.erase(search);

                        HistoryStep step;
                        step.type = RemoveSingleMutationCol;
                        step.single_mut_col.site = site;
                        step.single_mut_col.sequence_label = g.label;
                        column_counts[site] = 0;
                        history.push_back(step);
                        changed = true;
                        break;
                    }
                }
            }
            else
            {
                still_remaining_columns.push_back(site);
            }
        }
        remaining_columns = std::move(still_remaining_columns);

        // merge identical columns
        if (remaining_columns.size() == 0)
            break;
        still_remaining_columns = {remaining_columns[0]};
        int prev_site = remaining_columns[0];
        for (int i = 1; i < remaining_columns.size(); i++)
        {
            int site = remaining_columns[i];

            if (column_counts[prev_site] != column_counts[site])
            {
                still_remaining_columns.push_back(site);
                prev_site = site;
                continue;
            }

            bool identical = true;
            for (auto &g : gs)
            {
                bool in1 = std::find(g.mutations.begin(), g.mutations.end(), prev_site) != g.mutations.end();
                bool in2 = std::find(g.mutations.begin(), g.mutations.end(), site) != g.mutations.end();

                if (in1 != in2)
                {
                    identical = false;
                    break;
                }
            }

            if (identical)
            {
                if (how_verbose >= 2)
                    std::cout << "merging " << site << " into " << prev_site << "\n";
                // merge site1 into site2
                site_multiplicity[prev_site] += site_multiplicity[site];
                site_multiplicity[site] = 0;
                site_extent[prev_site] = std::max(site_extent[prev_site], site);
                column_counts[site] = 0;

                for (auto &g : gs)
                {
                    g.mutations.erase(std::remove(g.mutations.begin(), g.mutations.end(), site), g.mutations.end());
                }

                HistoryStep step;
                step.type = MergeIdenticalCols;
                step.merge_cols.site_from = site;
                step.merge_cols.site_to = prev_site;
                history.push_back(step);
                changed = true;
            }
            else
            {
                still_remaining_columns.push_back(site);
                prev_site = site;
            }
        }
        remaining_columns = std::move(still_remaining_columns);

        // Remove identical rows
        std::map<std::vector<int>, int> mutations_to_first_seq;
        for (int i = 0; i < gs.size(); i++)
        {
            auto g = gs[i];
            auto search = mutations_to_first_seq.find(g.mutations);
            if (search == mutations_to_first_seq.end())
            {
                mutations_to_first_seq[g.mutations] = i;
            }
            else
            {
                int target = search->second;

                // merge this into target gene
                for (int mut : g.mutations)
                {
                    column_counts[mut] -= 1;
                }

                if (how_verbose >= 2)
                    std::cout << "merging " << g.label << " into " << gs[target].label << "\n";
                HistoryStep step;
                step.type = MergeIdenticalRows;
                step.merge_rows.sequence_label_from = g.label;
                step.merge_rows.sequence_label_to = gs[target].label;
                history.push_back(step);
                changed = true;

                gs.erase(gs.begin() + i);
            }
        }
    }

    if (how_verbose >= 0)
        std::cout << "Cleaning finished, now have " << gs.size() << " genes and " << remaining_columns.size() << " sites.\n";

    genes.genes = std::move(gs);

    return std::make_tuple(std::move(history), std::move(site_multiplicity), std::move(site_extent));
}

/* Adds the steps in the history into the arg. Note that the maps in the arg will no longer be reliable */
void extend_arg_with_history(ARG &arg, const std::vector<HistoryStep> &_history, const std::vector<int> &site_extent)
{
    auto history = _history;
    std::reverse(history.begin(), history.end());

    std::map<std::string, Node *> label_to_node; // Set up a map from node labels to pointers to the nodes in the arg
    for (const auto &node : arg.nodes)
    {
        label_to_node[node->label] = node.get();
    }

    for (auto &step : history)
    {

        if (step.type == RemoveEmptyCol)
        {
            // Nothing to do
        }
        else if (step.type == RemoveSingleMutationCol)
        {
            // Add a mutation at site on the edge into seq
            int site = step.single_mut_col.site;
            std::string seq_label = step.single_mut_col.sequence_label;

            Node *node = label_to_node[seq_label];
            node->mutations.push_back(site);
            node->predecessor.one->mutations.push_back(site);
        }
        else if (step.type == MergeIdenticalCols)
        {
            int site_from = step.merge_cols.site_from;
            int site_to = step.merge_cols.site_to;

            // Update all nodes
            for (auto &node : arg.nodes)
            {
                if (std::find(node->mutations.begin(), node->mutations.end(), site_to) != node->mutations.end())
                {
                    node->mutations.push_back(site_from);
                }
            }

            // Update all edges
            auto range = arg.mutation_to_edges.equal_range(site_to);
            std::vector<Edge *> found_edges;
            for (auto i = range.first; i != range.second; ++i)
            {
                found_edges.push_back(i->second);
            }
            for (Edge *e : found_edges)
            {
                e->mutations.push_back(site_from);
                arg.mutation_to_edges.insert({site_from, e});
            }

            range = arg.back_mutation_to_edges.equal_range(site_to);
            found_edges.clear();
            for (auto i = range.first; i != range.second; ++i)
            {
                found_edges.push_back(i->second);
            }
            for (Edge *e : found_edges)
            {
                e->back_mutations.push_back(site_from);
                arg.back_mutation_to_edges.insert({site_from, e});
            }

            // Update arg statistics
            if (arg.recurrent_sites.count(site_to) > 0)
            {
                arg.recurrent_sites.insert(site_from);
            }
            if (arg.back_mutation_sites.count(site_to) > 0)
            {
                arg.back_mutation_sites.insert(site_from);
            }
        }
        else if (step.type == MergeIdenticalRows)
        {
            std::string seq_from_label = step.merge_rows.sequence_label_from;
            std::string seq_to_label = step.merge_rows.sequence_label_to;

            Node *node = label_to_node[seq_to_label];

            auto new_node = std::make_unique<Node>();
            auto new_edge = std::make_unique<Edge>();

            new_node->id = arg.nodes.size();
            new_node->label = seq_from_label;
            new_node->mutations = node->mutations;
            new_node->type = SAMPLE;
            new_node->predecessor.one = new_edge.get();

            new_edge->from = node;
            new_edge->to = new_node.get();

            new_edge->mutations.clear();
            new_edge->back_mutations.clear();

            label_to_node[seq_from_label] = new_node.get();

            arg.edges.push_back(std::move(new_edge));
            arg.nodes.push_back(std::move(new_node));
        }
    }

    for (auto *node : arg.recombination_nodes)
    {
        int pos = node->predecessor.two.position;
        if (pos > 0 && site_extent[pos-1] != pos-1)
        {
            node->predecessor.two.position = site_extent[pos-1] + 1;
        }
    }

    for (auto &node : arg.nodes)
    {
        std::sort(node->mutations.begin(), node->mutations.end());
    }
    for (auto &edge : arg.edges)
    {
        std::sort(edge->mutations.begin(), edge->mutations.end());
        std::sort(edge->back_mutations.begin(), edge->back_mutations.end());
    }

    std::cout << "Have added history to arg.\n";
}