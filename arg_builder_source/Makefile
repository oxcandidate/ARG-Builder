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

SRCFILES =  arg_builder.cpp arg_builder_logic.cpp vector_set_operations.cpp step_function.cpp clean.cpp
HDRFILES =  arg_builder_logic.h vector_set_operations.h step_function.h clean.h
EXECUTABLES = arg_builder
DATA = kreitman_snp.txt
SCRIPTS =
INSTRUCTIONS = 
AUTOSRCFILES =
AUTOHDRFILES =

all: arg_builder

arg_builder: arg_builder.o arg_builder_logic.o vector_set_operations.o step_function.o clean.o

arg_builder.o: arg_builder.cpp arg_builder_logic.h vector_set_operations.h step_function.h clean.h Makefile

arg_builder_logic.o: arg_builder_logic.cpp arg_builder_logic.h vector_set_operations.h step_function.h

vector_set_operations.o: vector_set_operations.h

clean.o: clean.cpp arg_builder_logic.h

step_function.o: step_function.cpp step_function.h

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

