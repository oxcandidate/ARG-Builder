#ifndef VECTOR_SET_OPERATIONS_H
#define VECTOR_SET_OPERATIONS_H

#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <algorithm>

std::vector<int> vector_intersect(const std::vector<int> &u, const std::vector<int> &v);
int vector_intersect_count(const std::vector<int> &u, const std::vector<int> &v);
std::vector<int> vector_difference(const std::vector<int> &u, const std::vector<int> &v);
std::pair<std::vector<int>, std::vector<int>> vector_split(const std::vector<int> &u, const std::vector<int> &v);
std::vector<int> vector_union(const std::vector<int> &u, const std::vector<int> &v);
std::vector<int> vector_values_below(const std::vector<int> &v, const int threshold);
std::vector<int> vector_values_above(const std::vector<int> &v, const int threshold);
std::vector<int> vector_values_between(const std::vector<int> &v, const int min, const int max);
std::vector<int> set_to_vector(const std::set<int> &set);
bool vector_contains(const std::vector<int> &v, const int e);
std::vector<int> vector_symmetric_difference(const std::vector<int> &u, const std::vector<int> &v);
std::set<int> vector_set_intersection(const std::vector<int> &u, const std::set<int> &v);
std::set<int> set_intersection(const std::set<int> &u, const std::set<int> &v);

#endif