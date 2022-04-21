/*******************************************************************

    arg_builder.h

    Description of functions for building and outputting an ancestral
    recombination graph

********************************************************************/

#ifndef ARG_BUILDER_H
#define ARG_BUILDER_H

#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <algorithm>

typedef enum
{
    GENE_ANY,
    GENE_BEAGLE,
    GENE_FASTA
} Gene_Format;

typedef enum
{
    GENE_BINARY,
    GENE_NUCLEIC,
    GENE_AMINO
} Gene_SeqType;

typedef enum
{
    ARGDOT, /* Output ARG as ARG in DOT format */
    ARGGDL, /* Output ARG as ARG in GDL format */
    ARGGML  /* Output ARG as ARG in GML format */
} ARGOutputFormat;

typedef enum
{
    LABEL_NONE,     /* Do not label nodes */
    LABEL_BOTH,     /* Label nodes with both label and sequence */
    LABEL_SEQUENCE, /* Label nodes only with sequence */
    LABEL_LABEL     /* Label nodes only with label */
} ARGOutputLabels;

typedef enum
{
    UNSET,         /* Node type has not been set yet */
    SAMPLE,        /* Node represents a sampled sequence */
    COALESCENCE,   /* Node represents a coalescence */
    RECOMBINATION, /* Node represents a recombination */
    ROOT           /* Node representing the root */
} NodeType;

typedef struct _Node Node;

typedef struct _Edge
{
    Node *from;
    Node *to;
    std::vector<int> mutations; /* Vector of sites which mutated along this edge */
    std::vector<int> back_mutations;
} Edge;

typedef struct _Node
{
    int id; // This is used to identify the node, particularly useful for printing graph
    NodeType type;
    std::string label;          // This is a description of the node. Has no effect on logic
    std::vector<int> mutations; // Most sequences have very few mutations, so easier to store as a vector
    union U
    {
        bool none;
        Edge *one;
        struct
        {
            Edge *prefix;
            Edge *suffix;
            int position;
        } two;

        U()
        {
        }

        ~U()
        {
        }
    } predecessor;
} Node;

typedef struct _ARG
{
    // Unique pointers are used to manage memory
    std::vector<std::unique_ptr<Node>> nodes;
    std::vector<std::unique_ptr<Edge>> edges;
    std::multimap<int, Edge *> mutation_to_edges;          // map to all edges which have the mutation
    std::multimap<int, Edge *> back_mutation_to_edges;     // map to all edges which have the back mutation
    std::multimap<int, Edge *> mutation_to_recombinations; // map to all recombination node out-edges containing the mutation

    std::vector<int> root_mutations;

    int number_of_roots = 1;

    int number_of_ancestral_nodes = 0;
    int number_of_back_mutations = 0;
    int number_of_recurrent_mutations = 0;
    int number_of_recombinations = 0;
    std::set<int> recurrent_sites;
    std::set<int> back_mutation_sites;
    std::set<Node *> recombination_nodes;

    std::vector<std::unique_ptr<Edge>> self_edges; // These are used for certain nodes (roots and recombinations) which are always relevant

    Edge *create_self_loop(Node *node_ptr)
    {
        auto self_edge = std::make_unique<Edge>();
        self_edge->from = node_ptr;
        self_edge->to = node_ptr;

        Edge *edge_ptr = self_edge.get();

        self_edges.push_back(std::move(self_edge));

        return edge_ptr;
    }

    Edge *add_root(std::string root_label, std::vector<int> _root_mutations)
    {
        auto new_root = std::make_unique<Node>();

        new_root->type = ROOT;
        number_of_roots += 1;
        new_root->id = -1 * number_of_roots;
        new_root->mutations = _root_mutations;
        new_root->label = root_label;
        new_root->predecessor.none = true;

        create_self_loop(new_root.get());

        root_mutations.insert(root_mutations.end(), _root_mutations.begin(), _root_mutations.end());

        nodes.push_back(std::move(new_root));
    }

    _ARG()
    {
        auto root = std::make_unique<Node>();

        root->type = ROOT;
        root->id = -1;
        root->mutations.clear();
        root->label = "Root";
        root->predecessor.none = true;

        create_self_loop(root.get());

        nodes.push_back(std::move(root));
    }

    _ARG(bool create_empty_root)
    {
        if (create_empty_root)
        {
            auto root = std::make_unique<Node>();

            root->type = ROOT;
            root->id = -1;
            root->mutations.clear();
            root->label = "Root";
            root->predecessor.none = true;

            create_self_loop(root.get());

            nodes.push_back(std::move(root));
        }
    }

    _ARG(std::string root_label, std::vector<int> _root_mutations)
    {
        auto root = std::make_unique<Node>();

        root->type = ROOT;
        root->id = -1;
        root->mutations = _root_mutations;
        root->label = root_label;
        root->predecessor.none = true;

        create_self_loop(root.get());

        root_mutations = _root_mutations;

        nodes.push_back(std::move(root));
    }
} ARG;

typedef struct _GENE
{
    std::string label;
    std::vector<int> mutations;
} Gene;

typedef struct _GENEs
{
    std::vector<Gene> genes;
    int sequence_length;
} Genes;

typedef struct _Run
{
    int recombinations;
    int recurrent_mutations;
    int back_mutations;

    int seed;
    double build_time;
    float recomb_cost;
    float rm_cost;
    float bm_cost;

    _Run(int recombs, int rms, int bms, int run_seed, double run_build_time, float cost_recomb, float cost_rm, float cost_bm)
    {
        recombinations = recombs;
        recurrent_mutations = rms;
        back_mutations = bms;

        seed = run_seed;
        build_time = run_build_time;
        recomb_cost = cost_recomb;
        rm_cost = cost_rm;
        bm_cost = cost_bm;
    }
} Run;

typedef struct _RunRecord
{
    std::vector<Run> runs;
    std::map<int, int> recombs_to_rare_mutations; // Map from # recombinations to best # bms+rms
    // std::vector<std::vector<int>> recombs_and_bms_to_rms;

    void clear()
    {
        runs.clear();
        recombs_to_rare_mutations.clear();
        // recombs_and_bms_to_rms.clear();
    }

    void add_record(int recombs, int rms, int bms, int run_seed, double build_time, float cost_recomb, float cost_rm, float cost_bm)
    {
        Run run(recombs, rms, bms, run_seed, build_time, cost_recomb, cost_rm, cost_bm);
        runs.push_back(run);

        // update recombs_to_rare_mutations
        auto lower_bound = recombs_to_rare_mutations.lower_bound(recombs);
        if (lower_bound->first == recombs)
        {
            // Already have entry for recombs
            if (rms + bms < lower_bound->second)
            {
                recombs_to_rare_mutations[recombs] = rms + bms;
            }
        }
        else if (lower_bound == recombs_to_rare_mutations.begin())
        {
            // recombs is smaller than any entry in map
            recombs_to_rare_mutations[recombs] = rms + bms;
        }
        else
        {
            // recombs is somewhere in the middle. Check it is better than the one before
            auto prev = std::prev(lower_bound);
            if (rms + bms < prev->second)
            {
                recombs_to_rare_mutations[recombs] = rms + bms;
            }
        }
    }
} RunRecord;

float get_cost(const int rms, const int bms, const int rcs);
float get_cost(const int rms, const int bms);
void arg_output(const ARG &arg, const Genes &genes, FILE *fp,
                ARGOutputFormat format, int how_to_label_edges, ARGOutputLabels node_labels);
std::tuple<ARG, RunRecord> build_arg_main(const Genes genes, bool clean_sequences, int how_verbose, int roots_given, int run_seed, int number_of_runs, int multi_run_strategy, int location_selection_method, int find_root_strategy, int find_root_iterations,
                                          int max_number_parents, float cost_rm, float cost_bm, float cost_recomb, int recomb_max, int rm_max, int bm_max, std::string run_record_file);
std::tuple<ARG, RunRecord> build_arg_search(const Genes genes, bool clean_sequences, int how_verbose, int roots_given, int run_seed, int number_of_runs, int multi_run_strategy, int location_selection_method,
                                            int max_number_parents, std::vector<float> costs_rm, std::vector<float> costs_bm, std::vector<float> costs_recomb, std::string run_record_file);
#endif