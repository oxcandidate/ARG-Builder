# Tools to use
CC = g++

# DEBUG = yes

# Flags to the C compiler
# If DEBUG is 'yes', the executable is not optimised and compiled with
# information used by the gnu debugging and profiling tools
ifeq (${DEBUG}, yes)
CPPFLAGS = -Wall -Wextra -Wpedantic --std=c++17 -g -pg -DDEBUG -w
else
CPPFLAGS = --std=c++17 -Wall -Wextra -Wpedantic -O3 -w
endif
# If ENABLE_VERBOSE is defined, code for outputting information about the
# progress of the computation at several levels of detail is included
ifdef ENABLE_VERBOSE
CPPFLAGS := ${CPPFLAGS} -DENABLE_VERBOSE
endif
# If BEAGLE_NOHAPLOTYPEBLOCKS is not defined, code for computing all local
# minimum number of recombinations is included
ifndef BEAGLE_NOHAPLOTYPEBLOCKS
CPPFLAGS := ${CPPFLAGS} -DHAPLOTYPE_BLOCKS
endif
# If BEAGLE_HAPLOTYPE is defined, the full composite haplotype bound will be
# used in the branch & bound computation of the minimum number of
# recombinations required.
ifdef BEAGLE_HAPLOTYPE
CPPFLAGS := ${CPPFLAGS} -DBEAGLE_HAPLOTYPE
else
# If BEAGLE_HUDSONKAPLAN is not defined, the composite haplotype bound
# computed using the third heuristic (see HAPLOTYPE_HEURISTICS) with a
# required increase in the number of subsets of at least two is used to
# obtain lower bounds in the branch an bound computation of the minimum
# number of recombinations required. If BEAGLE_HUDSONKAPLAN is defined,
# Hudson & Kaplan's lower bound is used.
ifndef BEAGLE_HUDSONKAPLAN
CPPFLAGS := ${CPPFLAGS} -DBEAGLE_HAPLOTYPEHEURISTIC
endif
endif
# The branch & bound computation of the minimum number of recombinations
# required uses a hash table to store ancestral states already encountered to
# avoid tracing back from the same ancestral state more than once. I.e. a sort
# of dynamic programming to ensure that we only go through the same
# configuration at most once. If BEAGLE_STORELEAVES is not defined, the
# ancestral states that are eliminated directly from the lower bound are not
# stored in the hash table. This leads to a reduction in the amount of
# memory used to store ancestral states already visited, but increases the
# computational time as it is usually faster to look up an ancestral state in
# the hash table than to recompute the lower bound.
ifndef BEAGLE_STORELEAVES
CPPFLAGS := ${CPPFLAGS} -DBEAGLE_DONOTSTORELEAVES
endif
# A minimum number of recombinations computation similar to the branch and
# bound computation of the minimum number of recombinations is part of the
# code to enumerate all ancestral states visitable in minimum, or near minimum,
# histories. This always uses the reduced state space approach, but whether
# the ancestral states that are eliminated directly from the lower bound are
# stored in the hash table is determined by ENUMERATE_STORELEAVES.
ifndef ENUMERATE_STORELEAVES
CPPFLAGS := ${CPPFLAGS} -DENUMERATE_DONOTSTORELEAVES
endif
# The lower bound method used for computing minimum number of recombinations
# required in the visitable ancestral states enumeration is chosen in the
# same manner as for the branch and bound computation of the minimum number
# of recombinations required
ifdef ENUMERATE_HAPLOTYPE
CPPFLAGS := ${CPPFLAGS} -DENUMERATE_HAPLOTYPE
else
ifndef ENUMERATE_HUDSONKAPLAN
CPPFLAGS := ${CPPFLAGS} -DENUMERATE_HAPLOTYPEHEURISTIC
endif
endif
# Flags to the linker
ifeq (${PROFILE}, yes)
LDFLAGS := ${LDFLAGS} -pg
LDLIBS := -lm
ifneq (${DEBUG}, yes)
CPPFLAGS := ${CPPFLAGS} -pg
endif
else
ifeq (${DEBUG}, yes)
LDFLAGS := ${LDFLAGS}
LDLIBS := -lm
else
LDFLAGS := ${LDFLAGS} -O3
LDLIBS := -lm
endif
endif

SRCFILES = beagle.cpp kwarg.cpp simplify.cpp flip.cpp gene.cpp llist.cpp elist.cpp hashtable.cpp getwordsize.cpp common.cpp bounds.cpp beagle_logic.cpp kwarg_logic.cpp mergesort.cpp enumerate.cpp arg.cpp backtrack.cpp topological_sort.cpp mystring.cpp
HDRFILES = gene.h llist.h elist.h hashtable.h common.h bounds.h beagle_logic.h kwarg_logic.h mergesort.h enumerate.h arg.h backtrack.h topological_sort.h mystring.h
EXECUTABLES = beagle kwarg simplify flip getwordsize
DATA = kreitman_snp.txt
SCRIPTS = generatebitfunctions.py
INSTRUCTIONS = 
AUTOSRCFILES = bitfunctions.c
AUTOHDRFILES = bitfunctions.h

all: beagle kwarg simplify flip

beagle: beagle.o gene.o bounds.o beagle_logic.o kwarg_logic.o llist.o elist.o hashtable.o common.o bitfunctions.o mergesort.o arg.o backtrack.o mystring.o

backtrack.o: backtrack.cpp backtrack.h arg.h common.h

arg.o: arg.cpp arg.h llist.h gene.h common.h bitfunctions.h

mystring.o: mystring.cpp mystring.h llist.h

kwarg: kwarg.o gene.o bounds.o common.o beagle_logic.o kwarg_logic.o llist.o elist.o hashtable.o bitfunctions.o mergesort.o arg.o backtrack.o mystring.o

simplify: simplify.o gene.o bounds.o beagle_logic.o llist.o elist.o hashtable.o common.o bitfunctions.o mergesort.o arg.o backtrack.o mystring.o

flip: flip.o gene.o bounds.o beagle_logic.o llist.o elist.o hashtable.o common.o bitfunctions.o mergesort.o arg.o backtrack.o mystring.o

beagle.o: beagle.cpp gene.h bounds.h common.h beagle_logic.h Makefile

kwarg.o: kwarg.cpp gene.h bounds.h common.h beagle_logic.h kwarg_logic.h backtrack.h Makefile

simplify.o: simplify.cpp gene.h bounds.h beagle_logic.h common.h backtrack.h Makefile

flip.o: flip.cpp gene.h bounds.h beagle_logic.h common.h backtrack.h Makefile

gene.o: gene.cpp gene.h common.h bitfunctions.h llist.h elist.h mergesort.h Makefile

bounds.o: bounds.cpp bounds.h beagle_logic.h gene.h common.h mergesort.h Makefile

beagle_logic.o: beagle_logic.cpp beagle_logic.h gene.h common.h mergesort.h bounds.h llist.h hashtable.h bitfunctions.h

kwarg_logic.o: kwarg_logic.cpp kwarg_logic.h beagle_logic.cpp beagle_logic.h gene.h common.h mergesort.h bounds.h llist.h hashtable.h bitfunctions.h

enumerate.o: enumerate.cpp enumerate.h gene.h elist.h common.h bounds.h beagle_logic.h bitfunctions.h mergesort.h

greedy.o: greedy.cpp gene.h common.h bounds.h

llist.o: llist.cpp llist.h Makefile

elist.o: elist.cpp elist.h Makefile

hashtable.o: hashtable.cpp hashtable.h elist.h common.h Makefile

common.o: common.cpp common.h llist.h Makefile

bitfunctions.o: bitfunctions.cpp bitfunctions.h Makefile

bitfunctions.cpp bitfunctions.h: generatebitfunctions.py getwordsize Makefile
	python3 generatebitfunctions.py; rm -f getwordsize

getwordsize: getwordsize.cpp Makefile
	${CC} ${CPPFLAGS} getwordsize.cpp -o getwordsize

mergesort.o: mergesort.cpp mergesort.h Makefile

topological_sort.o: topological_sort.cpp topological_sort.h common.h Makefile

archive: ${SRCFILES} ${HDRFILES} ${DATA} ${SCRIPTS} ${INSTRUCTIONS} Makefile
	tar cvf kwarg.`date +%d%m%y`.tar ${SRCFILES} ${HDRFILES} ${DATA} ${SCRIPTS} ${INSTRUCTIONS} Makefile

# make clean removes all builds and all intermediates
.PHONY: clean
clean:
	rm -f *.o *~ core ${EXECUTABLES} ${AUTOSRCFILES} ${AUTOHDRFILES}

# make semiclean removes all intermediates, like object files, but leaves
# the executables
.PHONY: semiclean
semiclean:
	rm -f *.o ${AUTOSRCFILES} ${AUTOHDRFILES}

