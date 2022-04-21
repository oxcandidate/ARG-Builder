/* mergesort.h: Prototype for mergesort function */

#ifndef MERGESORT_H
#define MERGESORT_H

#include <stdlib.h>

void merge_sort(char *data, size_t n, size_t size,
		int (*less_than)(char *, char *));

inline void merge_sort(void *data, size_t n, size_t size,
		int (*less_than)(char *, char *)) {
  merge_sort(static_cast<char*>(data), n, size, less_than);
}

#endif
