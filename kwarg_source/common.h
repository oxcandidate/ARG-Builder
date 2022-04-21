/* common.h: header file for common.c; definition of common data types.
 */
#ifndef COMMON_H
#define COMMON_H

#include "llist.h"
#include "hashtable.h"

#include <vector>
#include <list>
#include <set>

#define ENABLE_VERBOSE

// These definitions are passed when kwarg is created using make file
#define HAPLOTYPE_BLOCKS
#define BEAGLE_HAPLOTYPEHEURISTIC
#define BEAGLE_DONOTSTORELEAVES
#define ENUMERATE_DONOTSTORELEAVES
#define ENUMERATE_HAPLOTYPEHEURISTIC

/* Function prototypes; a brief explanation of each function should be
 * provided prior to its implementation in common.c.
 */
#ifdef ENABLE_VERBOSE
int verbose();
void set_verbose(int v);
#endif

typedef struct _Gene
{
    unsigned long *type;      /* Type bit vector */
    unsigned long *ancestral; /* Ancestral site bit vector */
} Gene;

typedef struct _Genes
{
    int n;      /* Number of sequences */
    int length; /* Length of sequences */
    Gene *data; /* Sequences */
} Genes;

typedef struct _AnnotatedGenes
{
    Genes *g;         /* Data */
    LList *positions; /* Site labels */
    LList *sequences; /* Sequence labels */
} AnnotatedGenes;

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

typedef struct _PackedGenes
{
    int n;              /* Number of sequences */
    int length;         /* Length of sequences */
    unsigned int *data; /* Sequences */
} PackedGenes;

typedef struct _Site
{
    unsigned long *type;      /* Type bit vector */
    unsigned long *ancestral; /* Ancestral sequence bit vector */
} Site;

typedef struct _Sites
{
    int n;      /* Number of sequences */
    int length; /* Length of sequences */
    Site *data; /* Sites */
} Sites;

typedef struct _Index
{
    int index;
    int block;
} Index;

typedef enum
{
    SUBSTITUTION,
    COALESCENCE,
    RECOMBINATION,
    REMOVE,
    COLLAPSE,
    SWAP,
    LOOKUP,
    SEFLIP,
    RMFLIP
} EventType;
typedef enum
{
    COAL,
    SE,
    RM,
    RECOMB1,
    RECOMB2
} Action;
typedef struct _Event
{
    EventType type;
    union
    {
        struct
        {
            int seq;
            int site;
        } s;
        struct
        {
            int s1;
            int s2;
        } c;
        struct
        {
            int seq;
            int pos;
        } r;
        struct
        {
            int s1;
            int s2;
        } swap;
        struct
        {
            int seq;
            int site;
        } flip;
        int collapse;
        int remove;
        int lookup;
    } event;

} Event;

typedef struct _EventList
{
    bool in_use; // set to false when not wanting to record some section of computation
    std::list<Event> events;

    _EventList()
    {
        in_use = true;
        events = {};
    }

    int size() const
    {
        if (!in_use)
        {
            fprintf(stderr, "EventList is not in use!");
        }
        return events.size();
    }

    void push_back(const Event &e)
    {
        if (!in_use)
        {
            fprintf(stderr, "EventList is not in use!");
        }
        events.push_back(e);
    }

    void set_null()
    {
        events = {};
        in_use = false;
    }

    void reset()
    {
        events = {};
        in_use = true;
    }

    // splices elements from other list to front
    void prepend(struct _EventList &other)
    {
        events.splice(events.begin(), other.events);
    }

    // splices elements from other list to front
    void append(struct _EventList &other)
    {
        events.splice(events.end(), other.events);
    }

    void destroy()
    {
        in_use = false;
    }

} EventList;

struct HistoryFragment
{
    Genes *g;         /* End configuration */
    EventList events; /* List of events leading from start
                       * configuration to end configuration.
                       */
    double step_cost; /* Number of recombination events */
    std::vector<int> elements;
    std::vector<int> sites;
    Action action;

    HistoryFragment() = default;

    HistoryFragment(HistoryFragment const &) = delete;
    HistoryFragment &operator=(HistoryFragment const &) = delete;

    ~HistoryFragment();
};

typedef struct _RunSettings
{
    // Fixed
    double se_cost;
    double rm_cost;
    double r_cost;
    double rr_cost;
    double temp;
    double run_seed;
    int run_reference;
    int rec_max, rm_max;
} RunSettings;

typedef struct _RunData
{
    // varies during run. Captures data on one path
    bool do_track = true; // Will be used to toggle when event_list etc should be added to

    double current_step_cost; // Note that sequence and site labels aren't used by beagle
    std::vector<int> sequence_labels = {};
    std::vector<int> site_labels = {};
    int seq_numbering;

    EventList eventlist;

    // Following are only used in scoring section
    HashTable *greedy_functioncalls = NULL, *greedy_beaglereusable = NULL;

#ifdef HAPLOTYPE_BLOCKS
    LList *representativeness = nullptr;
    LListCounter *representativeness_counter = nullptr;
    int **haploblocks = nullptr;
#endif

#ifdef DEBUG
    std::set<PackedGenes *> ancestral_state_trace = {};
    bool using_ancestral_state_trace = false;
#endif

    _RunData()
    {
        do_track = true;
        current_step_cost = 0;
        seq_numbering = 0;
        sequence_labels = {};
        site_labels = {};
    }

    _RunData(bool make_empty)
    {
        current_step_cost = 0;
        do_track = !make_empty;
        eventlist.in_use = !make_empty;
    }

    void clear_all();

    /**
     * Create a copy of this but with eventlist reset.
     * The values for do_track and eventlist.in_use will match this
     */
    struct _RunData copy_for_new_fragment() const;

    ~_RunData(); // destructor is in common.cpp so can call functions from elsewhere

    bool do_use_eventlist();

} RunData;

#ifdef HAPLOTYPE_BLOCKS
typedef struct _SuperColumn
{
    int left;
    int right;
} SuperColumn;

void explode_local(int **local, LList *r, int n);
#endif

#include "gene.h"

// These should be global as they apply to all runs/paths
extern bool g_use_eventlist;
extern std::vector<int> g_lookup;
extern int g_howverbose;
extern int gc_enabled; // Used in local2global. Not sure what is does. Is not changed

void *xmalloc(int n);
void *xcalloc(int m, int n);
void *xrealloc(void *oldadr, int n);
#define XRAND_MAX RAND_MAX
void initialise_xrandom();
double initialise_x2random(double seed);
long int xrandom();
long int x2random();
char *i2a(int n);
void pretty_print(FILE *fp, char *s, int l, int i);
void print_option(FILE *fp, char *option, char *description, int l, int i);
void remove_element(int *array, int index, int array_length);
void delete_i(int *array, int i, int array_length);
void delete_by_value(int *array, int v, int array_length);
int max_value(int *array, int array_length);
void print_int_vector(std::vector<int> vec, char *comment);
void set_array(double *a1, double *a2, int a2_length, int b);

/**
 * swaps the values at i and j in the vector
 */
template <typename T>
inline void vector_swap_elements(std::vector<T> &vec, int i, int j)
{
    T temp = std::move(vec[i]);
    vec[i] = std::move(vec[j]);
    vec[j] = std::move(temp);
}

/**
 * append source to destination by moving elements (source left empty)
 */
template <typename T>
inline void vector_append(std::vector<T> &destination, std::vector<T> &source)
{
    if (destination.empty())
        destination = std::move(source);
    else
    {
        std::move(std::begin(source), std::end(source), std::back_inserter(destination));
        // destination.insert(std::end(destination),
        //                    std::make_move_iterator(std::begin(source)),
        //                    std::make_move_iterator(std::end(source)));
    }
}

#endif
