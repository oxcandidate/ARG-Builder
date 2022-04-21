

#ifndef CLEAN_H
#define CLEAN_H

#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>

typedef enum
{
    RemoveEmptyCol,
    RemoveSingleMutationCol,
    MergeIdenticalRows,
    MergeIdenticalCols
} HistoryStepType;

typedef struct _HistoryStep
{
    HistoryStepType type;

    struct
    {
        int site;
    } empty_col;

    struct
    {
        int site;
        std::string sequence_label;
    } single_mut_col;

    struct
    {
        std::string sequence_label_from;
        std::string sequence_label_to;
    } merge_rows;

    struct
    {
        int site_from;
        int site_to;
    } merge_cols;
} HistoryStep;

#include "arg_builder_logic.h"

std::tuple<std::vector<HistoryStep>, std::vector<int>, std::vector<int>> clean_genes(Genes &genes, int how_verbose);

void extend_arg_with_history(ARG &arg, const std::vector<HistoryStep> &_history, const std::vector<int> &site_extent);

#endif