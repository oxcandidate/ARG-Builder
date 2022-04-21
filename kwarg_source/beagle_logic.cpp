/*******************************************************************
 *
 *    beagle_logic.c
 *
 *    Implementation of functions to compute the exact minimum
 *    number of recombinations required under the infinite sites
 *    assumption for an SNP data set (Beagle).
 *
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <memory>
#include <algorithm>

#include "gene.h"
#include "common.h"
#include "mergesort.h"
#include "bounds.h"
#include "llist.h"
#include "beagle_logic.h"
#include "hashtable.h"
#include "bitfunctions.h"
#include "backtrack.h"
#include "common.h"

/* Select method for computing a lower bound on the number of recombinations */
#ifdef BEAGLE_HAPLOTYPE
#define LOWERBOUND(g) haplotype_bound_genes(g)
#elif defined(BEAGLE_HAPLOTYPEHEURISTIC)
#define LOWERBOUND(g) haplotype_heuristic_genes(g, INT_MAX, INT_MAX, 1)
#else
#define LOWERBOUND(g) hudson_kaplan_genes(g)
#endif

typedef struct _BeagleSplitInformation
{
    std::unique_ptr<HistoryFragment> u_hist_ptr;
    int am;
    char splits;
    char representative;

    _BeagleSplitInformation() = default;

    // Can't simply copy due to unique pointer
    _BeagleSplitInformation(_BeagleSplitInformation const &) = delete;
    _BeagleSplitInformation &operator=(_BeagleSplitInformation const &) = delete;

    // But can move
    _BeagleSplitInformation& operator=(_BeagleSplitInformation&& other)
    {
        u_hist_ptr = std::move(other.u_hist_ptr);
        am = other.am;
        splits = other.splits;
        representative = other.representative;

        return *this;
    }

    // move constructor, takes a rvalue reference &&
    _BeagleSplitInformation (_BeagleSplitInformation&& other)
    {
        u_hist_ptr = std::move(other.u_hist_ptr);
        am = other.am;
        splits = other.splits;
        representative = other.representative;
    }
} BeagleSplitInformation;

/* Variables declared globally to avoid parameter clutter */
bool exact_randomise = false;    /* Random evolutionary histories? */
static bool reusable = false;    /* Is hash table reusable? */
static bool skip_lookup = false; /* Should LOOKUP event be allowed? */

/* The smaller element is the one with the smallest number of splits,
 * or if those are equal the element with the smaller amount of
 * ancestral material left.
 */
static int compareboundandancestralmaterial(BeagleSplitInformation *a,
                                            BeagleSplitInformation *b)
{
    if ((a->representative < b->representative) || ((a->representative == b->representative) && ((a->splits > b->splits) || ((a->splits == b->splits) && (a->am < b->am)))))
        return 1;
    else
        return 0;
}

/* Use the random function to draw an integer between 0 and n - 1 */
static int unbiased_random(int n)
{
    long int l = XRAND_MAX / n;
    long int i;

    do
    {
        i = xrandom() / l;
    } while (i >= n); /* i ought to always be at most n, but just to make sure */

    return i;
}

static void permute(BeagleSplitInformation *splits, int n)
{
    int i, j;
    BeagleSplitInformation tmp;

    for (i = n; i > 1;)
    {
        j = unbiased_random(i);
        i--;
        if (j != i)
        {
            memcpy(&tmp, splits + i, sizeof(BeagleSplitInformation));
            memcpy(splits + i, splits + j, sizeof(BeagleSplitInformation));
            memcpy(splits + j, &tmp, sizeof(BeagleSplitInformation));
        }
    }
}

/* Permute the elements in a vector */
template <typename T>
static void permute_vector(std::vector<T> &vec)
{
    int i, j;

    for (i = vec.size(); i > 1;)
    {
        j = unbiased_random(i);
        i--;
        vector_swap_elements(vec, i, j);
    }
}

/* Transfer the ancestral states in genes to the array splits (that
 * should be sufficiently large to hold all transfered ancestral
 * states). For each ancestral state it is checked whether it has
 * already been investigated or whether a computed lower bound on the
 * number of recombinations required plus the base number of
 * recombinations already used exceeds target. The total number of
 * ancestral states transferred is added to n.
 */
static bool transfer2splitinformation(std::vector<BeagleSplitInformation> &alt_splits,
                                      std::vector<std::unique_ptr<HistoryFragment>> &genes, int base, int *n,
                                      HashTable *table, int target, RunData &run_data)
{
    // genes will be cleared out in this process

    /* Insert all feasible branches from genes into splits */
    for (int i = 0; i < genes.size(); i++)
    {
        /* Check whether we have already encountered this ancestral state
         * with the same, or larger, target.
         */
        BeagleSplitInformation new_split;
        new_split.u_hist_ptr = std::move(genes[i]);

        PackedGenes *p = pack_genes(new_split.u_hist_ptr->g);
        void *lookup;
        if (hashtable_lookup((void *)p, table, &lookup))
        {
            /* This ancestral state is already present in the hash table */
            int prevtarget = (intptr_t)lookup;
            new_split.representative = 0;

            if (reusable && (prevtarget < 0) && (base - target <= prevtarget + 1))
            {
                free_packedgenes(p);
                /* We know this set of sequences needs at most as many
                 * recombinations as we have left.
                 */
                if (!skip_lookup)
                {
                    /* And we don't want to search histories */
                    if (g_use_eventlist && run_data.eventlist.in_use)
                    {
                        run_data.eventlist.append(new_split.u_hist_ptr->events);
                        Event e;
                        e.type = LOOKUP;
                        e.event.lookup = -prevtarget - 1;
                        run_data.eventlist.push_back(e);
                    }

                    /* We found a solution - eliminate all candidates */
                    for (i++; i < genes.size(); i++)
                    {
                        genes[i].reset(nullptr);
                    }

                    alt_splits.clear();
                    return true;
                }
            }
            else if ((reusable && ((2 * (target - base) + (base > 0) < prevtarget) || (base - target > prevtarget + 1))) || (!reusable && (target - base <= prevtarget)))
            {
                /* We know this set of sequences needs more recombinations
                 * than we have left.
                 */
                free_packedgenes(p);
                continue;
            }
            else
            {
                /* By the structure of the method, we know that no history
                 * exists with fewer than target - base recombinations.
                 */
                hashtable_update((void *)p,
                                 (void *)((target - base) << reusable), table, NULL);
                free_packedgenes(p);
            }
        }
        else
        {
            new_split.representative = 1;

            /* Compute lower bound for this sequence set */
            int bound = LOWERBOUND(new_split.u_hist_ptr->g);

            /* Is there any hope for this branch */
            if (bound > target - base)
            {
#ifdef BEAGLE_DONOTSTORELEAVES
                free_packedgenes(p);
#else
                if (reusable)
                    hashtable_update((void *)p, (void *)(2 * bound), t, NULL);
                else
                    hashtable_update((void *)p, (void *)(bound - 1), t, NULL);
#endif
                continue; // new_split will be reset
            }

            /* By the structure of the method, we know that no history
             * exists with fewer than target - base recombinations.
             */
            hashtable_update((void *)p,
                             (void *)((target - base) << reusable), table, NULL);
        }

        new_split.am = ancestral_material(new_split.u_hist_ptr->g);
        new_split.splits = base;

        /* Another sequence set inserted into splits */
        alt_splits.push_back(std::move(new_split));
        (*n)++;
    }

    return false;
}

/* Check whether any of the set of sequences in genes can be explained
 * without using recombinations. If so, free the memory used by the
 * sequence sets stored in genes as well as by genes itself and return
 * 1; otherwise return 0.
 */
static int check_for_bottom(const std::vector<std::unique_ptr<HistoryFragment>> &genes, RunData &run_data)
{
    for (int i = 0; i < genes.size(); i++)
    {
        HistoryFragment *history = genes[i].get();
        Genes *g = history->g;
        if (no_recombinations_required(g, run_data))
        {
            /* We can explain the set of sequences with no further recombinations */
#ifdef ENABLE_VERBOSE
            if (verbose())
            {
                printf("Bottom of recursion:\n");
                output_genes_indexed(g, NULL);
            }
#endif
            if (g_use_eventlist && run_data.eventlist.in_use)
            {
                run_data.eventlist.append(history->events);
                history->events.set_null();
            }
            return 1;
        }
    }

    return 0;
}

/* Recursion actually implementing the branch&bound procedure of
 * beagle. The hash table t is used for the dynamic programming and
 * branches known to lead to more recombinations than target are cut
 * off. If try_coalesces is true, all sensible coalesces are tried;
 * otherwise, only events involving at least one split are tried. The
 * rationale behind this, is that calls from an invocation where we
 * try all sensible coalesces should not themselves try all sensible
 * coalesces. If reusable is true, the values with which states are
 * inserted in the hash table are such that the hash table can be
 * reused. Return value just states whether target can be met (or
 * bettered). The number of ancestral states visited is reduced by
 * only attempting splits that will allow coalesces that are in some
 * sense maximal.
 */
static int beagle_recursion(Genes *g, HashTable *t, int target,
                            int try_coalesces, RunData &run_data)
{
    int j, max_possible_splits = 0;
    Index *start, *end;
    std::vector<std::unique_ptr<HistoryFragment>> coalesced, prefix, postfix, infix, overlap;
    std::vector<BeagleSplitInformation> alt_splits;
#ifdef ENABLE_VERBOSE
    int v = verbose();
    set_verbose(0);
#endif
    PackedGenes *p;

    if (try_coalesces)
    {
        /* We have just imploded genes, but we still need to pursue paths
         * coalescing compatible sequences where neither is subsumed in the
         * other but where the ancestral material is still entangled.
         */
        coalesced = coalesce_compatibleandentangled(g, run_data);
        max_possible_splits = coalesced.size();
        if (check_for_bottom(coalesced, run_data))
        {
#ifdef ENABLE_VERBOSE
            if (v)
            {
                printf("Reached by coalescing some compatible sequences from:\n");
                output_genes_indexed(g, NULL);
            }
            set_verbose(v);
#endif
            return 1;
        }
    }
    else
        max_possible_splits = 0;

    int num_splits = 0;
    if (target > 0)
    {
        /* There is room for at least one split */
        /* Determine interesting recombination ranges */
        start = maximumsubsumedprefixs(g);
        end = maximumsubsumedpostfixs(g);

        /* Try all sensible events with one split */
#ifdef ENABLE_VERBOSE
        set_verbose(v - 1);
#endif
        prefix = maximal_prefix_coalesces(g, start, end, run_data);
        if (!exact_randomise && check_for_bottom(prefix, run_data))
        {
            free(start);
            free(end);
            if (try_coalesces)
                coalesced.clear();
#ifdef ENABLE_VERBOSE
            if (v)
            {
                printf("Reached by coalescing a compatible prefix from:\n");
                output_genes_indexed(g, NULL);
            }
            set_verbose(v);
#endif
            return 1;
        }
        postfix = maximal_postfix_coalesces(g, start, end, run_data);
        /* If a random solution is required, we need to randomise the splits */
        if (exact_randomise)
        {
            vector_append(postfix, prefix);
            permute_vector(postfix);
            max_possible_splits += postfix.size();
        }
        else
            max_possible_splits += prefix.size() + postfix.size();
        if (check_for_bottom(postfix, run_data))
        {
            free(start);
            free(end);
            if (try_coalesces)
                coalesced.clear();
            if (!exact_randomise)
                prefix.clear();
#ifdef ENABLE_VERBOSE
            if (v)
            {
                if (exact_randomise)
                    printf("Reached by coalescing a compatible prefix or postfix from:\n");
                else
                    printf("Reached by coalescing a compatible postfix from:\n");
                output_genes_indexed(g, NULL);
            }
            set_verbose(v);
#endif
            return 1;
        }
        if (target > 1)
        {
            /* Try all sensible events with two splits */
            infix = maximal_infix_coalesces(g, start, end, run_data);
            if (!exact_randomise && check_for_bottom(infix, run_data))
            {
                free(start);
                free(end);
                if (try_coalesces)
                    coalesced.clear();
                prefix.clear();
                postfix.clear();
#ifdef ENABLE_VERBOSE
                if (v)
                {
                    printf("Reached by coalescing a compatible infix from:\n");
                    output_genes_indexed(g, NULL);
                }
                set_verbose(v);
#endif
                return 1;
            }
            overlap = maximal_overlap_coalesces(g, start, end, run_data);
            /* If a random solution is required, we need to randomise the splits */
            if (exact_randomise)
            {
                vector_append(overlap, infix);
                permute_vector(overlap);
                max_possible_splits += overlap.size();
            }
            else
                max_possible_splits += infix.size() + overlap.size();

            if (check_for_bottom(overlap, run_data))
            {
                free(start);
                free(end);
                if (try_coalesces)
                    coalesced.clear();
                postfix.clear();
                if (!exact_randomise)
                {
                    prefix.clear();
                    infix.clear();
                }
#ifdef ENABLE_VERBOSE
                if (v)
                {
                    if (exact_randomise)
                        printf("Reached by coalescing compatible overlaps or a compatible infix from:\n");
                    else
                        printf("Reached by coalescing compatible overlaps from:\n");
                    output_genes_indexed(g, NULL);
                }
                set_verbose(v);
#endif
                return 1;
            }

            if (max_possible_splits > 0)
            {
                if (!exact_randomise && transfer2splitinformation(alt_splits, infix, 2, &num_splits, t, target, run_data))
                {
                    free(start);
                    free(end);
                    if (try_coalesces)
                        coalesced.clear();
                    prefix.clear();
                    postfix.clear();
                    overlap.clear();
#ifdef ENABLE_VERBOSE
                    if (v)
                    {
                        printf("Found resolved ancestral state in hash table, reached from coalescing\ncompatible overlaps in:\n");
                        output_genes_indexed(g, NULL);
                    }
                    set_verbose(v);
#endif
                    return 1;
                }
                if (transfer2splitinformation(alt_splits, overlap, 2, &num_splits, t, target, run_data))
                {
                    free(start);
                    free(end);
                    if (try_coalesces)
                        coalesced.clear();
                    if (!exact_randomise)
                        prefix.clear();
                    postfix.clear();
#ifdef ENABLE_VERBOSE
                    if (v)
                    {
                        if (exact_randomise)
                            printf("Found resolved ancestral state in hash table, reached by coalescing a\ncompatible infix or compatible overlaps in:\n");
                        else
                            printf("Found resolved ancestral state in hash table, reached by coalescing a\ncompatible infix in:\n");
                        output_genes_indexed(g, NULL);
                    }
                    set_verbose(v);
#endif
                    return 1;
                }
            }
        }

        if (max_possible_splits > 0)
        {
            if (!exact_randomise && transfer2splitinformation(alt_splits, prefix, 1, &num_splits, t, target, run_data))
            {
                free(start);
                free(end);
                if (try_coalesces)
                    coalesced.clear();
                postfix.clear();
#ifdef ENABLE_VERBOSE
                if (v)
                {
                    printf("Found resolved ancestral state in hash table, reached by coalescing a\ncompatible prefix in:\n");
                    output_genes_indexed(g, NULL);
                }
                set_verbose(v);
#endif
                return 1;
            }
            if (transfer2splitinformation(alt_splits, postfix, 1, &num_splits, t, target, run_data))
            {
                free(start);
                free(end);
                if (try_coalesces)
                    coalesced.clear();
#ifdef ENABLE_VERBOSE
                if (v)
                {
                    if (exact_randomise)
                        printf("Found resolved ancestral state in hash table, reached by coalescing a\ncompatible prefix or postfix in:\n");
                    else
                        printf("Found resolved ancestral state in hash table, reached by coalescing a\ncompatible postfix in:\n");
                    output_genes_indexed(g, NULL);
                }
                set_verbose(v);
#endif
                return 1;
            }
        }

        /* Clean up */
        free(start);
        free(end);
    }

    if (try_coalesces)
    {
        if (coalesced.size() > 0)
        {
            if (transfer2splitinformation(alt_splits, coalesced, 0, &num_splits, t, target, run_data))
            {
#ifdef ENABLE_VERBOSE
                if (v)
                {
                    printf("Found resolved ancestral state in hash table, reached by coalescing some\ncompatible sequences in:\n");
                    output_genes_indexed(g, NULL);
                }
                set_verbose(v);
#endif
                return 1;
            }
        }
    }

    /* We now have all the ancestral states that could sensibly lead to
     * the current ancestral state g stored in splits. Sort them and
     * pursue each branch in sorted order.
     */
#ifdef ENABLE_VERBOSE
    set_verbose(v);
#endif

    fflush(stdout);
    if (num_splits > 0)
    {
        /* ...and there are actually some */
        if (exact_randomise)
        {
            permute_vector(alt_splits);
        }
        else
        {
            struct
            {
                bool operator()(const BeagleSplitInformation &a, const BeagleSplitInformation &b) const
                {
                    if ((a.representative < b.representative) || 
                        ((a.representative == b.representative) && ((a.splits > b.splits) || ((a.splits == b.splits) && (a.am < b.am)))))
                        return true;
                    else
                        return false;
                }
            } splitCmp;
            sort(alt_splits.begin(), alt_splits.end(), splitCmp);
        }
#ifdef ENABLE_VERBOSE
        if (v > 1)
        {
            int j = 0;
            for (const auto &split : alt_splits)
            {
                output_genes_indexed(split.u_hist_ptr->g, NULL);
                printf("is #%d; splits %d, and ancestral material %d\n",
                       j, split.splits, split.am);
                j++;
            }
        }
#endif


        for (const auto &split : alt_splits)
        {
            g = split.u_hist_ptr->g;
            if (beagle_recursion(g, t, target - split.splits,
                                 (split.splits > 0), run_data))
            {
                /* We found a branch leading to at most target recombinations */
#ifdef ENABLE_VERBOSE
                if (v)
                {
                    printf("Obtained from %d splits in:\n", split.splits);
                    output_genes_indexed(g, NULL);
                }
#endif
#ifdef DEBUG
                if (run_data.using_ancestral_state_trace)
                {
                    /* Insert current ancestral state as one visited in the
                     * minimum history we are in the process of returning back
                     * from.
                     */
                    p = pack_genes(g);
                    run_data.ancestral_state_trace.insert(p);
                }
#endif
                if (reusable)
                {
                    p = pack_genes(g);
                    /* We know this set of sequences is present in the hash
                     * table, as it was at the latest inserted in
                     * transfer2splitinformation.
                     */
                    hashtable_update((void *)p,
                                     (void *)(split.splits - target - 1), t, NULL);
                    free_packedgenes(p);
                }
                if (g_use_eventlist && run_data.eventlist.in_use)
                    run_data.eventlist.prepend(split.u_hist_ptr->events);
                
                alt_splits.clear();
                return 1;
            }
            else if (reusable)
            {
                p = pack_genes(g);
                hashtable_update((void *)p,
                                 (void *)(2 * (target - split.splits) + try_coalesces + 1), t, NULL);
                free_packedgenes(p);
            }

            if (g_use_eventlist && run_data.eventlist.in_use)
            {
                split.u_hist_ptr->events.destroy();
            }
        }






    }

    alt_splits.clear();

    /* We cannot explain g with at most target recombinations. Return
     * this informative fact.
     */
    return 0;
}

#ifdef ENABLE_VERBOSE
static void print_gene(Genes *g, va_list args)
{
    output_genes_indexed(g, NULL);
}
#endif

#ifdef HAPLOTYPE_BLOCKS
/* Return maximum of a and b */
static int _intmax(int a, int b)
{
    return (a > b ? a : b);
}
#endif

/* Use branch&bound plus dynamic programming techniques to reconstruct
 * a history for g requiring a minimum number of recombinations. If
 * run_data.eventlist is not NULL, a list of the events leading to this number
 * of recombinations is compiled in run_data.eventlist. If haploblocks is not
 * NULL, local minimum number of recombinations are stored in
 * haploblocks; haploblocks is assumed to be a table initialised to 0s
 * - entry i, j is set to the minimum number of recombinations
 * required in the region from site i to site i + j + 1. It is assumed
 * that lower is a valid lower bound on the number of recombinations,
 * and the search is terminated if it is established that no history
 * with at most upper recombinations exists (and a number larger than
 * upper is returned). If global variable reusable is true, t should
 * not be NULL and it is used as initial hash table, otherwise a local
 * hash table is allocated and later deallocated for keeping track of
 * ancestral states encountered.
 */
static int beagle_core(Genes *g, FILE *print_progress, int lower, int upper,
                       HashTable *t, RunData &run_data)
{
    Genes *h;
    int table_size, bound = 0, r1 = 0, r2 = 0;
    void *lookup;
#ifdef HAPLOTYPE_BLOCKS
    int i, j, n, localbound, oldreusable = reusable;
    LList *l, *tmprep = run_data.representativeness;
    LListCounter *tmprepcount = run_data.representativeness_counter;
    SuperColumn *c;
#endif
    PackedGenes *p;
#ifdef ENABLE_VERBOSE
    int v = verbose();
    set_verbose(0);
#endif

    /* Check whether we can rule out the necessity for recombinations */
    if (no_recombinations_required(g, run_data))
        return 0;

    /* We cannot handle upper bounds larger than half of INT_MAX */
    if (upper >= INT_MAX / 2)
        upper = INT_MAX / 2 - 1;

    /* Check whether sequence set is already present in hash table */
    if (reusable)
    {
        p = pack_genes(g);
        if (hashtable_lookup((void *)p, t, &lookup))
        {
            bound = (intptr_t)lookup;
            free_packedgenes(p);
            if (bound < 0)
            {
                /* We know the true minimum for this sequence set */
                if (!skip_lookup)
                {
                    /* And we are going to use it */
                    if (run_data.do_use_eventlist() && (-bound <= upper))
                    {
                        Event e;
                        e.type = LOOKUP;
                        e.event.lookup = -bound - 1;
                        run_data.eventlist.push_back(e);
                    }
                    return -bound - 1;
                }
            }
            else
                /* We have a lower bound for this sequence set */
                if (2 * upper + 1 < bound)
                    return upper + 1;
        }
        else
            free_packedgenes(p);
    }

#ifdef HAPLOTYPE_BLOCKS
    if (run_data.haploblocks != NULL)
    {
        n = g->length;
        run_data.representativeness = MakeLList();
        for (i = 0; i < g->length; i++)
        {
            c = (SuperColumn *)xmalloc(sizeof(SuperColumn));
            c->left = c->right = i;
            Enqueue(run_data.representativeness, (void *)c);
        }
        run_data.representativeness_counter = MakeCounter(run_data.representativeness, FIRST);
    }
#endif

    /* Compute a good lower bound on the number of recombinations */
    RunData implode_run_data;
    if (run_data.do_use_eventlist())
    {
        implode_run_data.do_track = true;
        implode_run_data.eventlist.reset();
    }
    else
    {
        implode_run_data.do_track = false;
        implode_run_data.eventlist.set_null();
    }

    /* Create working copy of g */
    g = copy_genes(g);
    implode_genes(g, implode_run_data);

#ifdef DEBUG
    /* Insert initial ancestral state as one that was visited */
    if (run_data.using_ancestral_state_trace && (g->n > 0))
    {
        p = pack_genes(g);
        run_data.ancestral_state_trace.insert(p);
    }
#endif

#ifdef HAPLOTYPE_BLOCKS
    if (run_data.haploblocks != NULL)
    {
        l = run_data.representativeness;
        run_data.representativeness = NULL;
        free(run_data.representativeness_counter);
    }
#endif
    if (!no_recombinations_required(g, run_data))
    {
        bound = LOWERBOUND(g);

        if (!reusable)
        {
            /* Determine size of and allocate hash table */
            table_size = msb((g->n - 3) * g->length) + bound;
            t = new_packedgeneshashtable(table_size);
        }

#ifdef HAPLOTYPE_BLOCKS
        if (run_data.haploblocks != NULL)
            reusable = 1;
#endif
#ifdef ENABLE_VERBOSE
        set_verbose(v);
#endif
        /* Determine minimum number of recombinations required for g */
        if (bound < lower)
            bound = lower;
        if (upper < 0)
            upper = INT_MAX;
        if (print_progress != NULL)
            fprintf(print_progress, "At least %d recombination%s required\n", bound,
                    (bound != 1 ? "s" : ""));
        p = pack_genes(g);
        if (hashtable_update((void *)p, (void *)(bound << reusable),
                             t, NULL) < 0)
        {
            r1 = 1;
        }
        for (; (bound <= upper) && !beagle_recursion(g, t, bound, 1, run_data);)
        {
            bound++;
#ifdef ENABLE_VERBOSE
            if (v)
                printf("%d states examined with bound %d\n", hashtable_size(t),
                       bound - 1);
            else if (print_progress != NULL)
                fprintf(print_progress, "At least %d recombinations required\n",
                        bound);
#else
            if (print_progress != NULL)
                fprintf(print_progress, "At least %d recombinations required\n",
                        bound);
#endif
            hashtable_update((void *)p, (void *)(bound << reusable),
                             t, NULL);
        }

        if (!r1)
            free_packedgenes(p);
        if (run_data.do_use_eventlist() && (bound > upper))
            /* We failed to find a valid history */
            implode_run_data.eventlist.destroy();

#ifdef ENABLE_VERBOSE
        if (v)
        {
            printf("%d states examined with true minimum (%d)\n", hashtable_size(t),
                   bound);
            if (v > 1)
            {
                printf("Hash modulo was %ld and number of collisions %d\n", t->size,
                       hashtable_collisions(t));
                hashtable_printlargestbucket(t, (void (*)(void *, va_list))print_gene);
                printf("One largest bucket was %d\n", hashtable_largestbucket(t));
            }
        }
#endif

#ifdef HAPLOTYPE_BLOCKS
        if (run_data.haploblocks != NULL)
        {
            /* We already computed minimum number of recombinations for full
             * sequence set.
             */
            run_data.haploblocks[0][g->length - 2] = bound;
            /* Do not remember events from local bound computations */
            RunData empty_data(true);
            /* Run through regions in order of increasing length */
            for (i = 2; i < g->length; i++)
                for (j = 0; j <= g->length - i; j++)
                {
                    h = copy_region(g, j, j + i);
                    implode_genes(h, empty_data);
                    if (!no_recombinations_required(h, empty_data))
                    {
                        p = pack_genes(h);
                        r2 = 0;
                        lookup = (void *)0;
                        if (!hashtable_lookup((void *)p, t, &lookup) || ((intptr_t)lookup >= 0))
                        {
                            localbound = (intptr_t)lookup / 2;
                            if (i > 2)
                                localbound = _intmax(localbound,
                                                     _intmax(run_data.haploblocks[j][i - 3],
                                                             run_data.haploblocks[j + 1][i - 3]));
                            if (hashtable_update((void *)p, (void *)(localbound << reusable),
                                                 t, NULL) < 0)
                                r2 = 1;
                            for (; (localbound <= upper) && !beagle_recursion(h, t, localbound, 1, empty_data);)
                                hashtable_update((void *)p, (void *)(++localbound << reusable),
                                                 t, NULL);
                            run_data.haploblocks[j][i - 2] = localbound;
#ifdef ENABLE_VERBOSE
                            if (v)
                                printf("%d recombinations required between site %d and site %d\n",
                                       localbound, j, j + i - 1);
#endif
                        }
                        else
                            run_data.haploblocks[j][i - 2] = -(intptr_t)lookup - 1;
                        if (!r2)
                            free_packedgenes(p);
                    }
                    free_genes(h);
                }

            /* So far we have worked on imploded sequence set - reconstruct
             * local bounds for original sequence set.
             */
            explode_local(run_data.haploblocks, l, n);
            /* Clean up */
            while (Length(l) > 0)
                free(Pop(l));
            DestroyLList(l);
            reusable = oldreusable;
        }
#endif
        /* Clean up */
        if (!reusable)
            hashtable_destroy(t, (void (*)(void *))free_packedgenes,
                              NULL, (void (*)(void *))free);
    }
#ifdef HAPLOTYPE_BLOCKS
    else if (run_data.haploblocks != NULL)
    {
        while (Length(l) > 0)
            DestroyLList((LList *)Pop(l));
        DestroyLList(l);
    }
    run_data.representativeness = tmprep;
    run_data.representativeness_counter = tmprepcount;
#endif

    if (g_use_eventlist && run_data.eventlist.in_use)
    {
        if (bound <= upper)
            run_data.eventlist.prepend(implode_run_data.eventlist);
        else
            run_data.eventlist.destroy();
    }

    free_genes(g);
    return bound;
}

/* Compute minimum number of recombinations required by any history for g */
int beagle(Genes *g, FILE *print_progress, RunData &run_data)
{
    return beagle_core(g, print_progress, 0, INT_MAX, NULL, run_data);
}

/* Compute minimum number of recombinations required by any history
 * for g, assuming lower is a lower bound on this number and
 * terminating the search once it has been established that upper does
 * not suffice.
 */
int beagle_bounded(Genes *g, FILE *print_progress, int lower, int upper, RunData &run_data)
{
    return beagle_core(g, print_progress, lower, upper, NULL, run_data);
}

/* Compute minimum number of recombinations required by any history
 * for g, using existing hash table t for checking and storing
 * ancestral states encountered.
 */
int beagle_reusable(Genes *g, FILE *print_progress, HashTable *t, RunData &run_data)
{
    reusable = 1;
    int n = beagle_core(g, print_progress, 0, INT_MAX, t, run_data);
    reusable = 0;

    return n;
}

/* Compute minimum number of recombinations required by any history
 * for g, assuming lower is a lower bound on this number and
 * terminating the search once it has been established that upper does
 * not suffice; t is used as hash table for checking and storing
 * ancestral states encountered.
 */
int beagle_reusable_bounded(Genes *g, FILE *print_progress, int lower,
                            int upper, HashTable *t, RunData &run_data)
{
    int n;

    reusable = 1;
    n = beagle_core(g, print_progress, lower, upper, t, run_data);
    reusable = 0;

    return n;
}

/* Find a random evolutionary history with at most r recombinations.
 * If r is less than the minimum number of recombinations required
 * for g, NULL is returned. If t is not NULL it is assumed to be a
 * reused and reusable hash table of ancestral configurations for
 * computing minimum number of recombinations required for this set of
 * genes. The evolutionary history inferred is returned as an llist of
 * Events.
 */

EventList beagle_randomised(Genes *g, FILE *print_progress, int r, HashTable *t, RunData &run_data)
{
    int reuse = 1, i, j, *permutation;
    Genes *h;
    EventList tmp = std::move(run_data.eventlist), history;
    std::vector<Event> events = {};
    std::vector<Genes *> configurations;
    run_data.eventlist.set_null();

    /* Allocate hash table for ancestral configurations encountered if
     * none such was provided.
     */
    if (t == NULL)
    {
        t = beagle_allocate_hashtable(g, -1, run_data);
        reuse = 0;
    }

    /* Determine number of recombinations required for data set */
    if (beagle_reusable_bounded(g, print_progress, 0, r, t, run_data) > r)
    {
        run_data.eventlist = std::move(tmp);
        history.set_null();
        return history;
    }

    /* Find a history with at most bound recombinations and store events
     * in history.
     */
    history.reset();
    h = copy_genes(g);
    while (h->n > 1)
    {
        /* Non-segregating sites confuse the backtrack - remove */
        /* Removal of nonsegregating sites will eliminate all sequences in
         * a deterministic manner if all sites are eliminated. So remember
         * the number of sites so that we can coalesce them in a random
         * manner if all sites are removed.
         */
        j = h->n;
        remove_nonsegregating(h, run_data);
        /* If this results in all sites being eliminated, we are done (but,
         * possibly, for a few coalescences).
         */
        if (h->length == 0)
        {
            for (i = j; i > 1; i--)
            {
                Event e;
                e.type = COALESCENCE;
                e.event.c.s1 = unbiased_random(i);
                j = unbiased_random(i - 1);
                if (j >= e.event.c.s1)
                    e.event.c.s2 = j + 1;
                else
                {
                    e.event.c.s2 = e.event.c.s1;
                    e.event.c.s1 = j;
                }
                history.push_back(e);
            }
            break;
        }
        /* Determine configurations reachable by mutation */
        configurations = force_mutation(h, events);
        /* Determine configurations reachable by coalescence */

        auto new_configs = force_coalesce(h, events, run_data);
        vector_append(configurations, new_configs);
        /* Determine configurations reachable by recombinations */
        for (i = 0; i < h->n; i++)
        {
            new_configs = force_split(h, i, events, run_data);
            vector_append(configurations, new_configs);
        }

        /* Go through configurations in randomly permuted order and choose
         * first that does not exceed recombination allowance.
         */
        permutation = (int *)xmalloc(configurations.size() * sizeof(int));
        for (i = 0; i < configurations.size(); i++)
            permutation[i] = i;
        for (i = 0; i < configurations.size(); i++)
        {
            j = unbiased_random(configurations.size() - i);
            free_genes(h);
            h = configurations[permutation[i + j]];
            Event e = events[permutation[i + j]];
            if (beagle_reusable_bounded(h, NULL, 0, r - (e.type == RECOMBINATION), t, run_data) <= r - (e.type == RECOMBINATION))
            {
                /* Found next configuration */
                history.push_back(e);
                if (e.type == RECOMBINATION)
                    /* Used one more recombination from our total allowance */
                    r--;
                /* Free all configurations and events not yet inspected */
                permutation[i + j] = permutation[i];
                for (i++; i < configurations.size(); i++)
                {
                    free_genes(configurations[permutation[i]]);
                }
                /* Clean up */
                configurations.clear();
                free(permutation);
                break;
            }
            else
            {
                /* This direction would lead to excess recombinations */
                permutation[i + j] = permutation[i];
                if (i == configurations.size() - 1)
                {
                    fprintf(stderr,
                            "Error in randomised search. Please email error report\n");
                    exit(1);
                }
            }
        }
    }

    /* Clean up */
    free_genes(h);
    run_data.eventlist = std::move(tmp);

    return history;
}

/* Allocate a hash table for storing ancestral states. The table_size
 * should be the logarithm of the number of buckets in the hash table
 * - if no valid number is provided (i.e. table_size is non-positive)
 * and a set of sequences is provided, a reasonable table size is
 * estimated.
 */
HashTable *beagle_allocate_hashtable(Genes *g, int table_size, RunData &run_data)
{
    HashTable *t;

    /* Determine size of and allocate hash table */
    if ((g != NULL) && (table_size <= 0))
    {
        if (no_recombinations_required(g, run_data))
            table_size = 3;
        else
            table_size = msb((g->n - 3) * g->length) + LOWERBOUND(g);
    }
    else if (table_size <= 0)
        table_size = 10;

    t = new_packedgeneshashtable(table_size);

    return t;
}

/* Deallocate a hash table used by beagle and all the elements stored in it */
void beagle_deallocate_hashtable(HashTable *t)
{
    hashtable_destroy(t, (void (*)(void *))free_packedgenes, NULL,
                      (void (*)(void *))free);
}

/* Functions for interfacing with lower bound computations */
static int _greedy_rmin;
static Genes *_greedy_currentstate;
/* Initialise parameters for hashing function call information into m bins */
static int *_greedy_initparam(unsigned long m)
{
    int i, *p = (int *)xmalloc(5 * sizeof(int));

    p[0] = m;
    for (i = 1; i < 5; i++)
        p[i] = xrandom() % m;

    return p;
}

/* Compute minimum number of recombinations needed for current state */
int noexp_rmin(RunData &run_data)
{
    /* Set up hash table for common use for all beagle invocations */
    if (run_data.greedy_beaglereusable == NULL)
    {
        run_data.greedy_beaglereusable = beagle_allocate_hashtable(_greedy_currentstate, -1, run_data);
    }
    if (_greedy_rmin < 0)
    {
        /* We haven't computed r_min for this configuration yet */
        _greedy_rmin = beagle_reusable(_greedy_currentstate, NULL,
                                       run_data.greedy_beaglereusable, run_data);
    }

    return _greedy_rmin;
}

/* Compute haplotype lower bound with the heuristic parameters
 * specified by p.
 */
int hb(Genes *g, RunData &run_data)
{
    int i;
    void **a = (void **)xmalloc(4 * sizeof(void *));

    /* Create array specifying this function call */
    for (i = 0; i < 3; i++)
    {
        a[i + 1] = (void *)INT_MAX;
    }
    a[0] = (void *)hb;

    /* Check whether this function call has previously been invoked; if
     * not, compute and store value.
     */
    if (!hashtable_lookup(a, run_data.greedy_functioncalls, (void **)&i))
    {
        i = haplotype_bound_genes(g);
        hashtable_insert(a, (void *)i, run_data.greedy_functioncalls);
    }

    /* Store lower bound in expression */

    return i;
}

/* Compute lower bound from local exact minimum number of
 * recombinations combined using the composite method.
 */
static int _eagl(Genes *g, RunData &run_data)
{
    void **a = (void **)xmalloc(2 * sizeof(void *));
    int b, **B;
    Sites *s;

    /* Create array for specifying this function call */
    a[0] = (void *)_eagl;
    a[1] = (void *)(int)(10);

    /* Check whether this function call has previously been invoked; if
     * not, compute and store value.
     */
    if (!hashtable_lookup(a, run_data.greedy_functioncalls, (void **)&b))
    {
        s = genes2sites(g);
        B = hudson_kaplan_local(s);
        free_sites(s);
        b = eagl(g, 10, B, NULL, NULL);
        hashtable_insert(a, (void *)b, run_data.greedy_functioncalls);
    }

    return b;
}

/* Hash the function call information stored in elm using the
 * parameters in p.
 */
static unsigned long _greedy_hash(void **elm, int *p)
{
    int i;
    unsigned long v = 0;

    if (*elm == hb)
        for (i = 1; i <= 3; i++)
            v += p[i + 1] * (intptr_t)elm[i];
    else
        v = p[1] + p[2] * (intptr_t)elm[1];

    return v % p[0];
}

/* Compare the two function calls a and b */
static int _greedy_compare(void **a, void **b)
{
    int i;

    if (*a != *b)
        return 0;

    if (*a == hb)
    {
        for (i = 1; i <= 3; i++)
            if ((intptr_t)a[i] != (intptr_t)b[i])
                return 0;
    }
    else
        return (intptr_t)a[1] == (intptr_t)b[1];

    return 1;
}

/* Set current ancestral state to g, update am, seq and len to reflect
 * g, and remove information from previous ancestral state from cache.
 */
static double _am;
void reset_beagle_builtins(Genes *g, RunData &run_data)
{
    _greedy_currentstate = g;

    _am = ancestral_material(g);

    _greedy_rmin = -1;
    if (run_data.greedy_functioncalls != NULL)
    {
        hashtable_cleanout(run_data.greedy_functioncalls, free, NULL);
    }
    else
        run_data.greedy_functioncalls = hashtable_new(6, (unsigned long (*)(void *, void *))_greedy_hash,
                                                      (int (*)(void *, void *))_greedy_compare,
                                                      (void *(*)(unsigned long))_greedy_initparam);
}

double get_am()
{
    return _am;
}
