#ifndef GENE_H
#define GENE_H

#define STORE_FRAGMENT_FUNCTION_TYPE std::function<void(Genes *, RunData)>

#include <stdio.h>

#include <vector>
#include <functional>
#include <memory>
#include <list>

#include "llist.h"
#include "hashtable.h"

#include "common.h"
#include "backtrack.h"

/* Datatypes */


/* Global variable for specifying whether the common ancestral
 * sequence is known or not.
 */
extern int gene_knownancestor;
/* Prototypes */
void free_genes(Genes *g);

// Moved from backtrack to avoid cyclic references


void free_annotatedgenes(AnnotatedGenes *g);
void free_sites(Sites *s);
void free_gene(Gene *g);
void free_site(Site *s);
Genes *new_genes(int n, int length, char **s, Gene_SeqType t);
AnnotatedGenes *read_genes(char *fname, Gene_Format f, Gene_SeqType t);
void output_genes(Genes *g, FILE *fp, char *comment);
void output_labelled_genes(Genes *g, FILE *fp, LList *labels);
void output_genes_indexed(Genes *s, FILE *fp);
void output_annotatedgenes(AnnotatedGenes *a, FILE *fp, char *comment);
char get_genes_character(Genes *g, int seq, int site);
void set_genes_character(Genes *g, int seq, int site, char c);
void swap_genes(Genes *g, int a, int b);
void add_ancestral_sites(Sites *s);
Genes *copy_genes(Genes *g);
Sites *copy_sites(Sites *s);
Genes *copy_allbutone(Genes *g, int a);
Genes *copy_region(Genes *g, int a, int b);
void remove_gene(Genes *g, int a);
void remove_annotatedgene(AnnotatedGenes *g, int a);
void reallocate_genes(Genes *g);
Sites *genes2sites(Genes *g);
char **genes2string(Genes *g);
void output_sites(Sites *s, FILE *fp, char *comment);
int remove_siamesetwins(Genes *g, RunData &run_data);
int remove_uninformative(Genes *g, RunData &run_data);
int remove_nonsegregating(Genes *g, RunData &run_data);
int coalesce_subsumed(Genes *g, RunData &run_data);
int implode_genes(Genes *g, RunData &run_data);
int no_recombinations_required(Genes *g, RunData &run_data);
void force_safeevents(Genes *g, RunData &run_data);
int force_mutations(Genes *g);
int mutate(Genes *g, int pos, int mutant);
std::vector<Genes *> force_mutation(Genes *g, std::vector<Event> &events);
std::vector<Genes *> force_mutation(Genes *g);
int segregating_site(Genes *g, int i);
int compatible(Genes *g, int a, int b);
std::vector<int> incompatible_sites(Genes *g, int a, int b);
int find_safe_coalescence(Genes *g, int a);
int entangled(Genes *g, int a, int b);
void coalesce(Genes *g, int a, int b, RunData &run_data);
std::vector<Genes *> force_coalesce(Genes *g, std::vector<Event> &events, RunData &run_data);
void split(Genes *g, int a, int i, RunData &run_data);
std::vector<Genes *> force_split(Genes *g, int a, RunData &run_data);
std::vector<Genes *> force_split(Genes *g, int a, std::vector<Event> &events, RunData &run_data);
void split_coalesceprefix(Genes *g, int a, int index, int block, int b, RunData &run_data);
void split_coalescepostfix(Genes *g, int a, int index, int block, int b, RunData &run_data);
void splitafter_coalescepostfix(Genes *g, int a, int index, int block, int b, RunData &run_data);
int split_removeprefix(Genes *g, int a, int index, int block, RunData &run_data);
int split_removepostfix(Genes *g, int a, int index, int block, RunData &run_data);
int ancestral_material(Genes *g);
int individual_all_ancestral(Genes *g, int i);
int ancestral_material_overlap(Genes *g);
int minimum_compatiblechops(Genes *g, int a);
Index *maximumsubsumedprefixs(Genes *g);
Index *maximumsubsumedpostfixs(Genes *g);
Index *maximumsubsumedprefix(Genes *g, int s);
Index *maximumsubsumedpostfix(Genes *g, int s);

std::vector<std::unique_ptr<HistoryFragment>> coalesce_compatibleandentangled(Genes *g, const RunData &main_path_data);
void coalesce_compatibleandentangled_map(Genes *g, const RunData &main_run_data, STORE_FRAGMENT_FUNCTION_TYPE f);

std::vector<std::unique_ptr<HistoryFragment>> maximal_prefix_coalesces(Genes *g, Index *a, Index *b, const RunData &main_path_data);
void maximal_prefix_coalesces_map(Genes *g, Index *a, Index *b, const RunData &main_path_data,
                                  STORE_FRAGMENT_FUNCTION_TYPE f);

std::vector<std::unique_ptr<HistoryFragment>> maximal_postfix_coalesces(Genes *g, Index *a, Index *b, const RunData &main_path_data);
void maximal_postfix_coalesces_map(Genes *g, Index *a, Index *b, const RunData &main_path_data,
                                   STORE_FRAGMENT_FUNCTION_TYPE f);

std::vector<std::unique_ptr<HistoryFragment>> maximal_infix_coalesces(Genes *g, Index *a, Index *b, const RunData &main_path_data);
void maximal_infix_coalesces_map(Genes *g, Index *a, Index *b, const RunData &main_path_data,
                                 STORE_FRAGMENT_FUNCTION_TYPE f);

std::vector<std::unique_ptr<HistoryFragment>> maximal_overlap_coalesces(Genes *g, Index *a, Index *b, const RunData &main_path_data);
void maximal_overlap_coalesces_map(Genes *g, Index *a, Index *b, const RunData &main_path_data,
                                   STORE_FRAGMENT_FUNCTION_TYPE f);

void seqerror_flips(Genes *g, const RunData &main_path_data, STORE_FRAGMENT_FUNCTION_TYPE f, const RunSettings &run_settings);
void recmut_flips(Genes *g, const RunData &main_path_data, STORE_FRAGMENT_FUNCTION_TYPE f, const RunSettings &run_settings);

int compare_sites(Sites *s, int a, int b);
int compare_genes(Genes *g, Genes *h);
PackedGenes *pack_genes(Genes *g);
void free_packedgenes(PackedGenes *p);
Genes *unpack_genes(PackedGenes *p);
int compare_packedgenes(PackedGenes *g, PackedGenes *h);
HashTable *new_packedgeneshashtable(int bits);

#endif
