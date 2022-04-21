#ifndef BEAGLE_LOGIC_H
#define BEAGLE_LOGIC_H

#include <stdio.h>
#include "gene.h"
#include "hashtable.h"

extern bool exact_randomise;
int beagle(Genes *g, FILE *print_progress, RunData &run_data);
int beagle_bounded(Genes *g, FILE *print_progress, int lower, int upper, RunData &run_data);
int beagle_reusable(Genes *g, FILE *print_progress, HashTable *t, RunData &run_data);
int beagle_reusable_bounded(Genes *g, FILE *print_progress, int lower,
			    int upper, HashTable *t, RunData &run_data);
EventList beagle_randomised(Genes *g, FILE *print_progress, int r, HashTable *t, RunData &run_data);
HashTable *beagle_allocate_hashtable(Genes *g, int table_size, RunData &run_data);
void beagle_deallocate_hashtable(HashTable *t);

int noexp_rmin(RunData &run_data);
int hb(Genes *g, RunData &run_data);
void reset_beagle_builtins(Genes *g, RunData &run_data);
double get_am();
#endif
