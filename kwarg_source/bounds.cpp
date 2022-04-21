/***********************************************************************
 *
 *    bounds.c: Implementation of functions to compute bounds of the minimum
 *    number of recombinations required under the infinite sites
 *    assumption for an SNP data set.
 *
 ************************************************************************/

#include <stdio.h>
#include <limits.h>
#include <string.h>

#include "gene.h"
#include "common.h"
#include "mergesort.h"
#include "bounds.h"
#include "beagle_logic.h"
#include "llist.h"
#include "bitfunctions.h"
#ifdef IPBOUNDS
#include "lpinterface.h"
#endif

/* Type for storing a score and sequence number */
typedef int IntPair[2];
/* Compare two IntPairs by their first value */
static int _compare_pairs(IntPair *a, IntPair *b)
{
    return (*a)[0] < (*b)[0];
}

/* Return the smaller of a and b */
static int _intmin(int a, int b)
{
    return (a < b ? a : b);
}

/* Return larger of a and b */
static int _intmax(int a, int b)
{
    return (a > b ? a : b);
}

#ifdef ENABLE_VERBOSE
/* yun is just a wrapper setting up history reconstruction information */
static int _yun_recursion(Genes *g, int best, int *history, int level);
int yun(Genes *g, int best)
{
    int i, j, n, v = verbose(),
                 *history = (int *)xmalloc(2 * g->n * sizeof(int));
    Genes *h;

    /* Find upper bound, generating no output */
    set_verbose(0);
    n = _yun_recursion(g, best, history, 0);

    if (v)
    {
        /* Trace back history giving upper bound */
        set_verbose(v);
        h = copy_genes(g);
        printf("Yun's bound requires %d splits\n", n);
        j = 0;
        for (i = 0; j < n; i++)
        {
            /* Output and pursue i'th event */
            printf(" Sequence %d of remaining sequences is split in %d places\n",
                   history[2 * i], history[2 * i + 1]);
            j += history[2 * i + 1];
            remove_gene(h, history[2 * i]);
            RunData empty_data(true);
            force_safeevents(h, empty_data);
            output_genes_indexed(h, NULL);
        }
        printf("\n");
        free_genes(h);
    }

    free(history);
    return n;
}

/* The real upper bound computation is carried out by yun_recursion,
 * finding the upper bound based on Yun's technique.
 */
static int _yun_recursion(Genes *g, int best, int *history, int level)
#else
/* Compute an upper bound on the number of recombinations required for the
 * data set in g, if this might be smaller than best; otherwise, just
 * return best.
 */
int yun(Genes *g, int best)
#endif
{
    int i, n;
    IntPair *chops;
    Genes *h;
#ifdef ENABLE_VERBOSE
    int y;
#endif

    if (g->n < (gene_knownancestor ? 3 : 4))
        /* Data sets with less than four sequences cannot contain
         * segregating sites.
         */
        return 0;

    n = hudson_kaplan_genes(g);
    if (best <= n)
        /* This branch cannot surpass current best branch */
        return best;

    /* Determine number of segments each sequence needs to be split into
     * to be totally subsumed in other sequences.
     */
    chops = (IntPair *)xmalloc(g->n * sizeof(IntPair));
    for (i = 0; i < g->n; i++)
    {
        chops[i][0] = minimum_compatiblechops(g, i);
        chops[i][1] = i;
    }

    /* Try all possibilities of subsuming one sequence in segments
     * created by recombination into the other sequences. Start with the
     * sequences that need to be chopped into fewest segments.
     */
    merge_sort(chops, g->n, sizeof(IntPair),
               (int (*)(char *, char *))_compare_pairs);
    for (i = 0; i < g->n; i++)
    {
        if (chops[i][0] >= best)
            break;
#ifdef ENABLE_VERBOSE
        if (level < 6)
            printf("Level %d, chopping sequence %d\n", level, chops[i][1]);
#endif
        h = copy_allbutone(g, chops[i][1]);
        /* Reduce data as much as possible */
        RunData empty_data(true);
        force_safeevents(h, empty_data);
#ifdef ENABLE_VERBOSE
        y = _yun_recursion(h, best - chops[i][0], history, level + 1);
        if (chops[i][0] + y < best)
        {
            best = chops[i][0] + y;
            history[2 * level] = chops[i][1];
            history[2 * level + 1] = chops[i][0];
        }
#else
        best = _intmin(best, chops[i][0] + yun(h, best - chops[i][0]));
#endif

        free_genes(h);
    }

    /* Clean up */
    free(chops);

    return best;
}

/* Compute matrix of local bounds using function bound, retrieve
 * global bound by composite method, and free memory used by matrix of
 * local bounds.
 */
static int _global_from_local(Sites *s, int **(*bound)(Sites *s))
{
    int **B = bound(s);
    int i, n;

    if (B == NULL)
        return 0;

    n = local2global(s->length - 1, B);
    for (i = 0; i < s->length - 1; i++)
        free(B[i]);
    free(B);

    return n;
}

/* Compute global bound, first converting g to Sites representaion and
 * then invoking global_from_local.
 */
static int _global_from_local_genes(Genes *g, int **(*bound)(Sites *s))
{
    Sites *s = genes2sites(g);
    int n = _global_from_local(s, bound);

    free_sites(s);

    return n;
}

/* Compute matrix of lower bounds on number of recombinations required
 * within an interval by Hudson and Kaplan's conflicting site
 * technique.
 */
int **hudson_kaplan_local(Sites *s)
{
    unsigned long type00, type01, type10, type11, filter;
    int i, j, k, **B;

    if (s->length < 2)
        /* Too few sites for recombinations to be detectable */
        return NULL;
    B = (int **)xmalloc((s->length - 1) * sizeof(int *));

    for (i = 0; i < s->length - 1; i++)
    {
        /* Find conflicting sites with site i as first site */
        B[i] = (int *)xcalloc(s->length - i - 1, sizeof(int));
        for (j = i + 1; j < s->length; j++)
        {
            /* Compare sites i and j */
            if (gene_knownancestor)
                type00 = ~0;
            else
                type00 = 0;
            type01 = type10 = type11 = 0;
            for (k = 0; k < divblocksize(s->n - 1) + 1; k++)
            {
                filter = s->data[i].ancestral[k] & s->data[j].ancestral[k];
                type00 |= ~s->data[i].type[k] & ~s->data[j].type[k] & filter;
                type01 |= ~s->data[i].type[k] & s->data[j].type[k] & filter;
                type10 |= s->data[i].type[k] & ~s->data[j].type[k] & filter;
                type11 |= s->data[i].type[k] & s->data[j].type[k] & filter;
            }
            if (type00 && type01 && type10 && type11)
                /* Site is segregating */
                B[i][j - i - 1] = 1;
        }
    }

    return B;
}

/* Compute lower bounds on number of recombinations required for s
 * using Hudson and Kaplan's conflicting site technique. The method
 * applies a heuristic part which may reduce the runtime to less than
 * the quadratic worst time.
 */
int hudson_kaplan(Sites *s)
{
    int i, j, k, blocks = divblocksize(s->n - 1) + 1;
    int *B;
    unsigned long type00, type01, type10, type11, filter;

    if (s->length < 2)
        /* No recombinations required */
        return 0;

    B = (int *)xmalloc(s->length * sizeof(int));
    B[0] = 0;
    for (i = 1; i < s->length; i++)
    {
        /* Derermine lower bound for prefix ending at site i */
        j = i - 1;
        B[i] = B[i - 1];
        while ((j >= 0) && (B[j] == B[i - 1]))
        {
            /* We can still improve on bound obtained simply by copying B[i - 1] */
            if (gene_knownancestor)
                type00 = ~0;
            else
                type00 = 0;
            type01 = type10 = type11 = 0;
            for (k = 0; k < blocks; k++)
            {
                filter = s->data[i].ancestral[k] & s->data[j].ancestral[k];
                type00 |= ~s->data[i].type[k] & ~s->data[j].type[k] & filter;
                type01 |= ~s->data[i].type[k] & s->data[j].type[k] & filter;
                type10 |= s->data[i].type[k] & ~s->data[j].type[k] & filter;
                type11 |= s->data[i].type[k] & s->data[j].type[k] & filter;
                if (type00 && type01 && type10 && type11)
                    break;
            }
            if (k < blocks)
            {
                /* Sites j and i are conflicting */
                B[i] = B[j] + 1;
                break;
            }
            j--;
        }
    }

    /* Record lower bound and clean up */
    i = B[s->length - 1];
    free(B);

    return i;
}

int hudson_kaplan_genes(Genes *g)
{
    Sites *s = genes2sites(g);
    //     int n = hudson_kaplan(s);
    int **B;
    int i, n;

    B = hudson_kaplan_local(s);
    n = local2global(g->length - 1, B);

    /* Clean up */
    free_sites(s);
    for (i = 0; i < g->length - 1; i++)
        free(B[i]);
    free(B);

    return n;
}

/* Compare the two sites with indeces a and b in cs_sites */
/* Function for comparing two indeces based on the site they index in
 * cs_sites
 */
static Sites *_cs_sites;
static int _compare_sites(int *i, int *j)
{
    return (compare_sites(_cs_sites, *i, *j) < 0);
}

/* Determine whether elm is among the n elements in set (which is
 * assumed to be sorted).
 */
static int binary_search(int elm, int n, int *set)
{
    if (n == 1)
        return (elm == *set);
    else if (n == 0)
        return 0;
    else if (set[n / 2] < elm)
        return binary_search(elm, n / 2, set);
    else
        return binary_search(elm, (n + 1) / 2, set + n / 2);
}

/* Structure for storing an interval */
typedef struct _Interval
{
    int left;
    int right;
} Interval;

/* Structure for storing a site type; it is represented by the address
 * of the index of its first occurrence and its total number of
 * occurrences (with indeces of remaining occurrences assumed to occupy
 * the following addresses).
 */
typedef struct _Type
{
    int *first;
    int n;
} Type;

/* Find the index of the largest of the n sorted elements specified by
 * a that is smaller than b. Update a to refer to this element. Just
 * having to do this once, the optimal choice would be binary
 * search. But in this setting, if linear search takes a long time we
 * have to do this search a corresponding fewer number of times. So
 * linear search guarantees linear time while binary search might take
 * O(n * log(n)) time.
 */
static void next_interval(int b, Type *a)
{
    int i = 1;

    while ((i < a->n) && (a->first[i] < b))
        i++;

    a->first += i - 1;
    a->n -= i - 1;
}

/* Compare two types by their number of occurrences */
static int compare_typeprevalence(Type *a, Type *b)
{
    /* We want the types sorted in order of decreasing prevalence */
    return (a->n > b->n);
}

/* Compute how many times we have to double m before it becomes at
 * least as large as n.
 */
static int doublings(int n, int m)
{
    int i = msb((unsigned long)n) - msb((unsigned long)m);

    if (m << i < n)
        return i + 1;
    else
        return i;
}

/* Compute matrix of lower bounds on number of recombinations required
 * within an interval by Myers and Griffith's haplotype technique.
 */
int **haplotype_bound_local(Sites *s)
{
    return haplotype_heuristic_local(s, INT_MAX, INT_MAX, 0);
}

int haplotype_heuristic(Sites *s, int maxsetsize, int maxintervallength,
                        int subsetincreasethreshold)
{
    int **B = haplotype_heuristic_local(s, maxsetsize, maxintervallength,
                                        subsetincreasethreshold);
    int i, n;

    if (B == NULL)
        return 0;

    n = local2global(s->length - 1, B);
    for (i = 0; i < s->length - 1; i++)
        free(B[i]);
    free(B);

    return n;
}

int haplotype_heuristic_genes(Genes *g, int maxsetsize, int maxintervallength,
                              int subsetincreasethreshold)
{
    Sites *s = genes2sites(g);
    int n = haplotype_heuristic(s, maxsetsize, maxintervallength,
                                subsetincreasethreshold);

    free_sites(s);

    return n;
}

/* Compute matrix of lower bounds on number of recombinations required
 * within an interval by Myers and Griffith's haplotype technique. Only
 * consider subsets with at most maxsubsetsize types and intervals no
 * longer than maxintervallength (measured as the number of sites
 * contained in the interval). Discard type additions that do not
 * increase the number of subsets with more than subsetincreasethreshold.
 */
int **haplotype_heuristic_local(Sites *s, int maxsetsize,
                                int maxintervallength,
                                int subsetincreasethreshold)
{
    int c, d, h, i, j = 0, k, l, m, n = 0, o,
                    blocks = divblocksize(s->n - 1) + 1, nleft, *typeset, *nintervals,
                    *nsubsets, *sorted, *issegregating, *index, *npartners, **partners, **B;
    Interval **intervals;
    Type *types, current[2];
    unsigned long type00, type01, type10, type11, *filter, **sequences,
        ***subsets;
#ifdef ENABLE_VERBOSE
    int a = 0, b = 0;
#endif

#ifdef ENABLE_VERBOSE
    if (verbose())
        printf("Haplotype bounds\n");
#endif

    if (s->length < 2)
        /* Too few sites for recombinations to be detectable */
        return NULL;

    if (gene_knownancestor)
    {
        /* Add ancestral sequence to data */
        s = copy_sites(s);
        add_ancestral_sites(s);
    }

    /* Create an array of site indeces sorted according to type */
    _cs_sites = s;
    sorted = (int *)xmalloc(s->length * sizeof(int));
    for (i = 0; i < s->length; i++)
        sorted[i] = i;
    merge_sort(sorted, s->length, sizeof(int),
               (int (*)(char *, char *))_compare_sites);

    /* Find number representative of each different types */
    B = (int **)xmalloc((s->length - 1) * sizeof(int *));
    types = (Type *)xmalloc(s->length * sizeof(Type));
    types[0].first = sorted;
    for (i = 1; i < s->length; i++)
    {
        if (_compare_sites(types[n].first, sorted + i))
        {
            /* Strict inequality, hence a new type */
            types[n].n = i - j;
            n++;
            types[n].first = sorted + i;
            j = i;
        }
        /* Also use loop to initialise B */
        B[i - 1] = (int *)xmalloc((s->length - i) * sizeof(int));
    }
    types[n].n = i - j;
    n++;
    merge_sort(types, n, sizeof(Type),
               (int (*)(char *, char *))compare_typeprevalence);

    /* Determine incompatibilities and weed out types that are not
     * incompatible with any other type. This should also weed out
     * uninformative sites.
     */
    l = 0;
    issegregating = (int *)xcalloc(s->length, sizeof(int));
    index = (int *)xmalloc(s->length * sizeof(int));
    npartners = (int *)xcalloc(s->length, sizeof(int));
    filter = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            /* Compare types */
            type00 = type01 = type10 = type11 = 0;
            for (k = 0; k < blocks; k++)
            {
                *filter = s->data[*types[i].first].ancestral[k] & s->data[*types[j].first].ancestral[k];
                type00 |= ~s->data[*types[i].first].type[k] & ~s->data[*types[j].first].type[k] & *filter;
                type01 |= ~s->data[*types[i].first].type[k] & s->data[*types[j].first].type[k] & *filter;
                type10 |= s->data[*types[i].first].type[k] & ~s->data[*types[j].first].type[k] & *filter;
                type11 |= s->data[*types[i].first].type[k] & s->data[*types[j].first].type[k] & *filter;
            }
            /* Use B as temporary storage for number of subsets the two
             * types split the sequences into.
             */
            B[i][j - i - 1] = (type00 != 0) + (type01 != 0) + (type10 != 0) + (type11 != 0);
            if (B[i][j - i - 1] == 4)
            {
                /* Segregating sites partner found for i and j */
                issegregating[i] = 1;
                issegregating[j] = -1;
            }
            if (B[i][j - i - 1] > 2)
            {
                /* It makes sense to have both type i and j in a set of types */
                npartners[i]++;
            }
        }
        if (issegregating[i] == 0)
        {
#ifdef ENABLE_VERBOSE
            if (verbose())
            {
                if (types[i].n > 1)
                {
                    printf(" Sites %d", *types[i].first);
                    for (j = 1; j < types[i].n - 1; j++)
                        printf(", %d", types[i].first[j]);
                    printf(" and %d are compatible with all other sites\n",
                           types[i].first[j]);
                }
                else
                    printf(" Site %d is compatible with all other sites\n",
                           *types[i].first);
            }
#endif
            l++;
        }
        else
        {
            index[i] = i - l;
#ifdef ENABLE_VERBOSE
            if (verbose())
            {
                if (types[i].n > 1)
                {
                    printf(" Sites %d", types[i].first[0]);
                    for (j = 1; j < types[i].n - 1; j++)
                        printf(", %d", types[i].first[j]);
                    printf(" and %d are type %d\n", types[i].first[j], i - l);
                }
                else
                    printf(" Site %d is type %d\n", *types[i].first, i - l);
            }
#endif
        }
    }

    partners = (int **)xmalloc(s->length * sizeof(int *));
    if (l < n)
    {
        /* Move type and number of subsets information to their new
         * positions after all sites that are not incompatible with any
         * other sites are removed, and for each type compile a list of the
         * types together with which it splits the sequences into at least
         * three subsets.
         */
        j = 0;
        for (i = 0; i < n; i++)
        {
            if (issegregating[i] == 1)
            {
                /* Compile list of types it makes sense to have in the same set
                 * as type i.
                 */
                partners[index[i]] = (int *)xmalloc(npartners[i] * sizeof(int));
                k = 0;
                for (j = i + 1; j < n; j++)
                    /* Either i and j should be incompatible, or they should
                     * split the sequences into three sets and j should be
                     * incompatible with a higher indexed type.
                     */
                    if ((B[i][j - i - 1] == 4) || ((issegregating[j] == 1) && B[i][j - i - 1] == 3))
                    {
                        partners[index[i]][k++] = index[j];
                    }
                npartners[index[i]] = k;
            }
            else if (issegregating[i] == -1)
                npartners[index[i]] = 0;
            if (issegregating[i] != 0)
            {
                /* Move type information */
                types[index[i]].first = types[i].first;
                types[index[i]].n = types[i].n;
            }
        }
        n = n - l;

        /* Initialise arrays for storing information about sets of
         * conflicting sites and reset B.
         */
        typeset = (int *)xmalloc(n * sizeof(int));
        intervals = (Interval **)xmalloc((n - 1) * sizeof(Interval *));
        /* The maximum number of minimal intervals is limited by the case
         * where the two most prevalent types occur alternating.
         */
        j = 2 * types[1].n;
        if (types[0].n == types[1].n)
            j--;
        intervals[0] = (Interval *)xmalloc(j * sizeof(Interval));
        for (i = 1; i < n - 1; i++)
        {
            /* For intervals containing more than two types, the limit is the
             * smaller of the number of ways we can combine a set of existing
             * intervals with a new type and the 2 * (|s| - 1) bound obtained
             * from observing that every site can be the leftmost,
             * respectively rightmost, site in at most one minimal
             * interval.
             */
            /* Between each occurence of the new type we can have an interval
             * that gives birth to two new intervals, one by adding the
             * occurence to the left and one by adding the occurence to the
             * right. All other intervals can, at best, be retained if there
             * is an occurence of the new type within the interval. But each
             * occurence of the new type can at most ensure the survival of at
             * most the previous number of types minus 1 existing intervals,
             * as this is the maximum number of minimal intervals that can
             * span a particular site.
             */
            j = _intmin(2 * (s->length - 1), _intmin(types[i + 1].n - 1 + j,
                                                     (i + 2) * types[i + 1].n));
            intervals[i] = (Interval *)xmalloc(j * sizeof(Interval));
        }
        nintervals = (int *)xmalloc((n - 1) * sizeof(Interval *));
        subsets = (unsigned long ***)xmalloc((n - 1) * sizeof(unsigned long **));
        for (i = 0; i < n - 1; i++)
        {
            subsets[i] = (unsigned long **)xmalloc(s->n * sizeof(unsigned long *));
            for (j = 0; j < s->n; j++)
                subsets[i][j] = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        }
        nsubsets = (int *)xmalloc((n - 1) * sizeof(int));
        for (i = 0; i < s->length - 1; i++)
            for (j = 0; j < s->length - i - 1; j++)
                B[i][j] = 0;
        sequences = (unsigned long **)xmalloc((n - 1) * sizeof(unsigned long *));
        for (i = 0; i < n - 1; i++)
            sequences[i] = (unsigned long *)xmalloc(blocks * sizeof(unsigned long));
        /* Start from each pair of types splitting the set of sequences into
         * at least three subsets and keep adding new types to this as long
         * as this makes sense.
         */
        for (i = 0; i < n - 1; i++)
        {
            typeset[0] = i;
            for (j = 0; j < npartners[i]; j++)
            {
                /* Start from type i and i's jth partner */
                typeset[1] = partners[i][j];
                /* Create subsets these two types split the sequences into */
                l = nsubsets[0] = 0;
                /* Determine subset of sequences that have 1 in site i and 1 in
                 * site j, and determine what sequences that have ancestral
                 * material in both sites.
                 */
                for (k = 0; k < blocks; k++)
                {
                    sequences[0][k] = s->data[*types[i].first].ancestral[k] & s->data[*types[typeset[1]].first].ancestral[k];
                    if ((subsets[0][0][k] = s->data[*types[i].first].type[k] & s->data[*types[typeset[1]].first].type[k] & sequences[0][k]) != 0)
                        l = 1;
                }
                nsubsets[0] = l;
                /* Determine subset of sequences that have 0 in site i and 1
                 * in site j.
                 */
                for (k = 0; k < blocks; k++)
                    if ((subsets[0][nsubsets[0]][k] = ~s->data[*types[i].first].type[k] & s->data[*types[typeset[1]].first].type[k] & sequences[0][k]) != 0)
                        l = nsubsets[0] + 1;
                /* Determine subset of sequences that have 1 in site i and 0
                 * in site j.
                 */
                nsubsets[0] = l;
                for (k = 0; k < blocks; k++)
                    if ((subsets[0][nsubsets[0]][k] = s->data[*types[i].first].type[k] & ~s->data[*types[typeset[1]].first].type[k] & sequences[0][k]) != 0)
                        l = nsubsets[0] + 1;
                /* Determine subset of sequences that have 0 in site i and 0
                 * in site j.
                 */
                nsubsets[0] = l;
                for (k = 0; k < blocks; k++)
                    if ((subsets[0][nsubsets[0]][k] = ~s->data[*types[i].first].type[k] & ~s->data[*types[typeset[1]].first].type[k] & sequences[0][k]) != 0)
                        l = nsubsets[0] + 1;
                nsubsets[0] = l;
#ifdef ENABLE_VERBOSE
                if (verbose())
                {
                    printf(" Types %d and %d split the sequence into %d subsets\n",
                           i, typeset[1], nsubsets[0]);
                    if (verbose() > 1)
                        printf(" Minimal intervals are: ");
                }
#endif
                /* Create the minimal intervals spanned by pairs of these
                 * types of segregating sites.
                 */
                nintervals[0] = 0;
                current[0].first = types[i].first;
                current[0].n = types[i].n;
                current[1].first = types[typeset[1]].first;
                current[1].n = types[typeset[1]].n;

                k = (*types[i].first < *types[typeset[1]].first ? 0 : 1);
                /* If i is of type a and j is of type b, then (i, j) is a
                 * minimal interval if there is no i' of type a with i < i' < j
                 * or j' of type b with i < j' < j. If we know there is no j' of
                 * type b with i < j' < j, all we have to do is find the largest
                 * i of type a smaller than j. The next occurrence of type a (if
                 * any) will be in i' where we know there are no occurrences of
                 * type a between j and i'. So now we can continue with the
                 * symmetric situation where we need to find the largest j
                 * smaller than i' and so on.
                 */
                /* Smaller index is at k, larger at 1 - k */
                while (current[1 - k].n > 0)
                {
                    /* Construct next minimal interval */
                    next_interval(*current[1 - k].first, current + k);
                    intervals[0][nintervals[0]].left = *current[k].first;
                    intervals[0][nintervals[0]].right = *current[1 - k].first;
                    if (intervals[0][nintervals[0]].right - intervals[0][nintervals[0]].left < maxintervallength)
                        nintervals[0]++;
                    /* Update corresponding entry of B */
                    if ((nsubsets[0] == 4) && (B[*current[k].first]
                                                [*current[1 - k].first - *current[k].first - 1] == 0))
                        B[*current[k].first]
                         [*current[1 - k].first - *current[k].first - 1] = 1;
                    /* Prepare for next interval */
                    current[k].first++;
                    current[k].n--;
                    k = 1 - k;
#ifdef ENABLE_VERBOSE
                    if (verbose() > 1)
                    {
                        printf("[%d,%d]", intervals[0][nintervals[0] - 1].left,
                               intervals[0][nintervals[0] - 1].right);
                        if (intervals[0][nintervals[0]].right - intervals[0][nintervals[0]].left >= maxintervallength)
                            printf("*");
                    }
#endif
                }
#ifdef ENABLE_VERBOSE
                if (verbose() > 1)
                    printf("\n That is %d interval%s\n", nintervals[0],
                           (nintervals[0] == 1 ? "" : "s"));
#endif
                if (nintervals[0] == 0)
                    /* No intervals of sufficient small length exists */
                    continue;
                /* Keep adding extra types to the initial set of conflicting
                 * sites as long as it seems feasible.
                 */
                /* Use m to keep track of size of current set of types */
                m = 2;
                /* Use l to keep track of the next type to add to the set of types */
                l = typeset[1] + 1;
                while (m > 1)
                {
                    /* Update set of types */
                    /* Try adding an extra type */
                    for (; l < n; l++)
                    {
                        /* Check whether there are enough sequences left after
                         * adding this type to improve bounds.
                         */
                        nleft = 0;
                        for (k = 0; k < blocks; k++)
                            nleft += weight(sequences[m - 1][k] = sequences[m - 2][k] & s->data[*types[l].first].ancestral[k]);
                        if (nleft <= m + 3)
                            /* Too few sequences left after adding l: we already have
                             * m + 1 sites after adding l, and there can be at most nleft
                             * types, so no bound better than nleft - (m + 1) - 1 can be
                             * obtained; any set that gives a non-trivial bound will
                             * contain segregating types that in themselves give a
                             * bound of 1; so if nleft - (m + 1) - 1 is at most 1 there is
                             * no point in continuing down this branch.
                             */
                            continue;
                        /* Check whether l and one of the types already present in
                         * the set splits the remaining sequences into at least
                         * three subsets.
                         */
                        for (h = 0; h < m; h++)
                        {
                            type00 = type01 = type10 = type11 = 0;
                            for (k = 0; k < blocks; k++)
                            {
                                type00 |= ~s->data[*types[typeset[h]].first].type[k] & ~s->data[*types[l].first].type[k] & sequences[m - 1][k];
                                type01 |= ~s->data[*types[typeset[h]].first].type[k] & s->data[*types[l].first].type[k] & sequences[m - 1][k];
                                type10 |= s->data[*types[typeset[h]].first].type[k] & ~s->data[*types[l].first].type[k] & sequences[m - 1][k];
                                type11 |= s->data[*types[typeset[h]].first].type[k] & s->data[*types[l].first].type[k] & sequences[m - 1][k];
                                if ((type00 != 0) + (type01 != 0) + (type10 != 0) + (type11 != 0) > 2)
                                    /* Stop as soon as we know we have found a match */
                                    break;
                            }
                            if (k < blocks)
                                /* Types l and typeset[h] splits the sequences at least
                                 * threeways.
                                 */
                                break;
                        }
                        if (h < m)
                        {
                            /* Types l and typeset[h] splits the sequences at least
                             * threeways. Add l to set of types and stop the search
                             * for the next type to add.
                             */
                            typeset[m] = l;
                            /* Find the new set of subsets the sequences are split into */
                            c = h = nsubsets[m - 1] = 0;
                            for (o = 0; o < nsubsets[m - 2]; o++)
                            {
                                /* Determine new subset where the new type has a 0 */
                                for (k = 0; k < blocks; k++)
                                    if ((subsets[m - 1][nsubsets[m - 1]][k] = sequences[m - 1][k] & subsets[m - 2][o][k] & ~s->data[*types[typeset[m]].first].type[k]) != 0)
                                        h = nsubsets[m - 1] + 1;
                                d = h - nsubsets[m - 1];
                                if ((nsubsets[m - 1] = h) == nleft)
                                    /* All remaining sequences are now their own subset */
                                    break;
                                /* Determine new subset where the new type has a 1 */
                                for (k = 0; k < blocks; k++)
                                    if ((subsets[m - 1][nsubsets[m - 1]][k] = sequences[m - 1][k] & subsets[m - 2][o][k] & s->data[*types[typeset[m]].first].type[k]) != 0)
                                        h = nsubsets[m - 1] + 1;
                                c |= d + h - nsubsets[m - 1];
                                if ((nsubsets[m - 1] = h) == nleft)
                                    /* All remaining sequences are now their own subset */
                                    break;
                            }
                            if ((nsubsets[m - 1] <= nsubsets[m - 2] + subsetincreasethreshold) || ((c & 2) == 0))
                                /* Insufficient increase in number of subsets */
                                continue;
                            /* Compute the best bound ever attainable for extensions
                             * to this set of types.
                             */
                            nleft = nleft - m - doublings(nleft, nsubsets[m - 1]) - 1;
                            if (nleft <= 1)
                                /* Adding the new type can never lead to a bound
                                 * better than 1.
                                 */
                                continue;
                            /* Break out of the search for the next feasible set of types */
                            break;
                        }
                    }
                    if (l >= n)
                    {
                        /* We did not find a feasible type to add to the current set
                         * of types. Remove last type added and continue the search
                         * from there.
                         */
                        m--;
                        l = typeset[m] + 1;
                        continue;
                    }
#ifdef ENABLE_VERBOSE
                    if (verbose())
                    {
                        printf(" Types %d", typeset[0]);
                        for (k = 1; k < m; k++)
                            printf(", %d", typeset[k]);
                        printf(" and %d gives a bound of %d", typeset[m],
                               nsubsets[m - 1] - m - 2);
                        if (verbose() > 1)
                            printf("\n Intervals covered by this set of types are:");
                    }
#endif
                    /* We have found a new feasible set of types. Construct the
                     * minimal intervals spanning this set of types and update B
                     * accordingly.
                     */
                    nintervals[m - 1] = k = h = 0;
                    while ((k < types[typeset[m]].n) && (h < nintervals[m - 2]))
                    {
                        if (types[typeset[m]].first[k] < intervals[m - 2][h].right)
                        {
                            /* Find largest occurrence of new type not to the right
                             * of interval.
                             */
                            while ((k < types[typeset[m]].n - 1) && (types[typeset[m]].first[k + 1] < intervals[m - 2][h].right))
                                k++;
                            if (types[typeset[m]].first[k] < intervals[m - 2][h].left)
                            {
                                /* Combine this occurence with interval to form new interval */
                                intervals[m - 1][nintervals[m - 1]].left = types[typeset[m]].first[k];
                                intervals[m - 1][nintervals[m - 1]].right = intervals[m - 2][h].right;
                                k++;
                            }
                            else
                            {
                                /* j is inside interval - retain interval */
                                intervals[m - 1][nintervals[m - 1]].left = intervals[m - 2][h].left;
                                intervals[m - 1][nintervals[m - 1]].right = intervals[m - 2][h].right;
                                h++;
                            }
                        }
                        else
                        {
                            /* Find rightmost interval to the right of the current
                             * occurrence of the new type.
                             */
                            while ((h < nintervals[m - 2] - 1) && (intervals[m - 2][h + 1].right < types[typeset[m]].first[k]))
                                h++;
                            /* Combine interval with the occurrence of the new type
                             * to form new interval.
                             */
                            intervals[m - 1][nintervals[m - 1]].left = intervals[m - 2][h].left;
                            intervals[m - 1][nintervals[m - 1]].right = types[typeset[m]].first[k];
                            h++;
                        }
#ifdef ENABLE_VERBOSE
                        if (verbose() > 1)
                        {
                            printf(" [%d-%d]", intervals[m - 1][nintervals[m - 1]].left,
                                   intervals[m - 1][nintervals[m - 1]].right);
                            if (intervals[m - 1][nintervals[m - 1]].right - intervals[m - 1][nintervals[m - 1]].left >= maxintervallength)
                                printf("*");
                        }
#endif
                        /* Update B with this interval */
                        if (nsubsets[m - 1] - m - 2 >
                            B[intervals[m - 1][nintervals[m - 1]].left]
                             [intervals[m - 1][nintervals[m - 1]].right - intervals[m - 1][nintervals[m - 1]].left - 1])
                            B[intervals[m - 1][nintervals[m - 1]].left]
                             [intervals[m - 1][nintervals[m - 1]].right - intervals[m - 1][nintervals[m - 1]].left - 1] = nsubsets[m - 1] - m - 2;
                        if (intervals[m - 1][nintervals[m - 1]].right -
                                intervals[m - 1][nintervals[m - 1]].left <
                            maxintervallength)
                            /* We may actually consider intervals that are too long
                             * for updating B above (it would probably take longer
                             * not to do this), but we do not carry them forward.
                             */
                            /* Increase number of intervals */
                            nintervals[m - 1]++;
                    }
#ifdef ENABLE_VERBOSE
                    b += nintervals[m - 1];
                    a++;
                    if (verbose())
                    {
                        if (verbose() > 1)
                            printf(" = %d interval%s\n", nintervals[m - 1],
                                   (nintervals[m - 1] == 1 ? "" : "s"));
                        else
                            printf(" and covers %d minimal interval%s\n", nintervals[m - 1],
                                   (nintervals[m - 1] == 1 ? "" : "s"));
                        fflush(stdout);
                    }
#endif
                    l++;
                    if ((nintervals[m - 1] > 0) && (m < maxsetsize - 1))
                        m++;
                }
            }
        }
#ifdef ENABLE_VERBOSE
        if (verbose())
            printf("Total number of subsets: %d\nTotal number of intervals: %d\n",
                   a, b);
#endif
        /* Clean up */
        for (i = 0; i < n - 1; i++)
        {
            if (npartners[i])
                free(partners[i]);
            free(intervals[i]);
            for (j = 0; j < s->n; j++)
                free(subsets[i][j]);
            free(subsets[i]);
            free(sequences[i]);
        }
        if (npartners[n - 1])
            free(partners[n - 1]);
        free(typeset);
        free(intervals);
        free(nintervals);
        free(subsets);
        free(nsubsets);
        free(sequences);
    }
    else
    {
        /* No recombinations required anywhere */
        for (i = 1; i < s->length; i++)
            for (j = 0; j < s->length - i; j++)
                B[i - 1][j] = 0;
#ifdef ENABLE_VERBOSE
        if (verbose())
            printf("No incompatible segregating sites in data\n");
#endif
    }

    /* Clean up */
    free(sorted);
    free(issegregating);
    free(index);
    free(npartners);
    free(partners);
    free(types);
    free(filter);
    if (gene_knownancestor)
        free_sites(s);

    return B;
}

int haplotype_bound(Sites *s)
{
    return _global_from_local(s, haplotype_bound_local);
}

int haplotype_bound_genes(Genes *g)
{
    return _global_from_local_genes(g, haplotype_bound_local);
}

/* The core computation of exact number of recombinations required for
 * local regions. Parameters are as for eagl, with t being the hash
 * table used for storing ancestral states already visited.
 */
static void eagl_core(Genes *g, int max_length, int **B, HashTable *t,
                      int (*per_region_action)(void),
                      int (*per_increment_action)(int, int, int **), RunData &run_data)
{
    int i, j;
    Genes *h;

    /* Run through regions in order of increasing length */
    for (i = 2; i <= max_length; i++)
    {
        for (j = 0; j <= g->length - i; j++)
        {
            h = copy_region(g, j, j + i);
            implode_genes(h, run_data);
            if (!no_recombinations_required(h, run_data))
            {
                B[j][i - 2] = beagle_reusable_bounded(h, NULL,
                                                      (i > 2 ? _intmax(B[j][i - 3], B[j + 1][i - 3]) : 1),
                                                      INT_MAX, t, run_data);
            }
            free_genes(h);
            /* Check for termination */
            if (per_region_action != NULL)
                if ((*per_region_action)())
                    return;
        }
        /* Check for termination */
        if (per_increment_action != NULL)
            if ((*per_increment_action)(i, g->length - 1, B))
                return;
    }
}

/* Compute a local bound on the number of recombinations necessary to
 * explain g, based on computing the exact minimum number of
 * recombinations necessary to explain all regions of g of at most
 * max_length sites. These local exact values can be combined using
 * local2(locals|global). If B is not NULL, it is assumed to contain
 * existing lower bounds and it is updated with the exact values
 * computed; otherwise, a new matrix is allocated. If
 * per_region_action is not NULL, it is invoked whenever an entry has
 * been computed and if the return value is True the computation is
 * terminated. If per_increment_action is not NULL, it is invoked
 * whenever region size is incremented and if the return value is True
 * the computation is terminated.
 */
int **eagl_local(Genes *g, int max_length, int **B,
                 int (*per_region_action)(void),
                 int (*per_increment_action)(int, int, int **))
{
    int i;
    HashTable *t;

    if (g->length < 2)
        /* Too few sites for recombination to be detectable */
        return NULL;

    if (max_length > g->length)
        max_length = g->length;

    /* Allocate matrix for storing lower bounds and hash table */
    if (B == NULL)
    {
        /* Allocate matrix for local bounds */
        B = (int **)xmalloc((g->length - 1) * sizeof(int *));
        for (i = 0; i < g->length - 1; i++)
            B[i] = (int *)xcalloc((g->length - i - 1), sizeof(int));
        t = new_packedgeneshashtable(msb((g->n - 3) * g->length) + 5);
    }
    else
        t = new_packedgeneshashtable(msb((g->n - 3) * g->length) + B[0][g->length - 2]);

    /* Compute local exact values */
    RunData empty_data(true);
    eagl_core(g, max_length, B, t, per_region_action, per_increment_action, empty_data);

    /* Clean up */
    hashtable_destroy(t, (void (*)(void *))free_packedgenes, NULL,
                      (void (*)(void *))free);

    return B;
}

int eagl(Genes *g, int max_length, int **B,
         int (*per_region_action)(void),
         int (*per_increment_action)(int, int, int **))
{
    int i, n;

    B = eagl_local(g, max_length, B, per_region_action, per_increment_action);
    if (B == NULL)
        return 0;
    n = local2global(g->length - 1, B);
    for (i = 0; i < g->length - 1; i++)
        free(B[i]);
    free(B);

    return n;
}

/* Compute a global lower bound on number of recombinations from the
 * local lower bounds in the n x n matrix B. The i, j entry should
 * hold the local bound on the number of recombinations between site i
 * and site i + j + 1. The first row of the matrix, i.e. the local
 * bounds starting in site 0, will be modified by this function.
 */
// int local2global(int n, int **B)
// {
//     int i, j;
//
//     for (i = 1; i < n; i++){
//         for (j = 0; j < i; j++){
//             if (B[0][i] < B[0][j] + B[j + 1][i - j - 1]) {
//                 B[0][i] = B[0][j] + B[j + 1][i - j - 1];
//             }
//         }
//     }
//
//     return B[0][n - 1];
// }

int local2global(int n, int **B)
{
    int i, j;

    if (gc_enabled)
    {
        if (B[0][0] > 0)
        {
            for (j = 0; j < n - 1; j++)
            {
                // Reduce local bound for all intervals starting with 1 (1->x) by 1.
                if (B[1][j] > 0)
                {
                    B[1][j] = B[1][j] - 1;
                }
            }
        }
    }

    for (i = 1; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            if (B[0][i] < B[0][j] + B[j + 1][i - j - 1])
            {
                B[0][i] = B[0][j] + B[j + 1][i - j - 1];
            }
        }
        if ((B[0][i] > B[0][i - 1]) && gc_enabled && (i < n - 1))
        {
            // Using this to calculate the global bound if we are allowing one-site-long "gene conversions".
            // We have added at least one recombination point just before site i.
            // So, try turning these into one-site-long "gene conversions" to cover (once) the intervals starting at i.
            // Either this will remove one recombination from the lower bound, or leave it unchanged (as intervals further
            // on will now require one more recombination point).
            for (j = 0; j < n - i - 1; j++)
            {
                // Reduce local bound for all intervals starting with i (i->x) by 1.
                if (B[i + 1][j] > 0)
                {
                    B[i + 1][j] = B[i + 1][j] - 1;
                }
            }
        }
    }

    return B[0][n - 1];
}

/* Compute all local lower bounds on number of recombinations from the
 * local lower bounds in the n x n matrix B. The i, j entry should
 * hold the local bound on the number of recombinations between site i
 * and site i + j + 1. The matrix will be modified by this function to
 * reflect the new lower bounds.
 */
void local2locals(int n, int **B)
{
    int i, j, k, bound;

    for (i = 0; i < n - 1; i++)
        /* Compute lower bounds starting in site i */
        for (j = 1; j < n - i; j++)
            /* Compute lower bound ending in site i + j + 1 */
            for (k = 0; k < j; k++)
            {
                /* Compute lower bound from i + k + 1 as intermediate site */
                bound = B[i][k] + B[i + k + 1][j - k - 1];
                if (B[i][j] < bound)
                    B[i][j] = bound;
            }
}
