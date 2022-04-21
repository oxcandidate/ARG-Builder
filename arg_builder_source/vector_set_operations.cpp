
#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <algorithm>
#include <utility>

/* Return u\v */
std::vector<int> vector_intersect(const std::vector<int> &u, const std::vector<int> &v)
{
    int N = u.size();
    int M = v.size();
    int i = 0, j = 0;
    std::vector<int> intersection = {};

    while (i < N && j < M)
    {
        if (u[i] == v[j])
        {
            intersection.push_back(u[i]);
            i++;
            j++;
        }
        else if (u[i] > v[j])
        {
            j++;
        }
        else
        {
            i++;
        }
    }

    return std::move(intersection);
}

int vector_intersect_count(const std::vector<int> &u, const std::vector<int> &v)
{
    int N = u.size();
    int M = v.size();
    int i = 0, j = 0;
    int count = 0;

    while (i < N && j < M)
    {
        if (u[i] == v[j])
        {
            count++;
            i++;
            j++;
        }
        else if (u[i] > v[j])
        {
            j++;
        }
        else
        {
            i++;
        }
    }

    return count;
}

/* Return u\v */
std::vector<int> vector_difference(const std::vector<int> &u, const std::vector<int> &v)
{
    int N = u.size();
    int M = v.size();
    int i = 0, j = 0;
    std::vector<int> difference = {};

    while (i < N && j < M)
    {
        if (u[i] == v[j])
        {
            i++;
            j++;
        }
        else if (u[i] > v[j])
        {
            j++;
        }
        else
        {
            difference.push_back(u[i]);
            i++;
        }
    }
    while (i < N)
    {
        difference.push_back(u[i]);
        i++;
    }

    return std::move(difference);
}

/* Returns pair with (u intersect v, u\\v) */
std::pair<std::vector<int>, std::vector<int>> vector_split(const std::vector<int> &u, const std::vector<int> &v)
{
    int N = u.size();
    int M = v.size();
    int i = 0, j = 0;
    std::vector<int> intersection = {};
    std::vector<int> rest = {};

    while (i < N && j < M)
    {
        if (u[i] == v[j])
        {
            intersection.push_back(u[i]);
            i++;
            j++;
        }
        else if (u[i] > v[j])
        {
            j++;
        }
        else
        {
            rest.push_back(u[i]);
            i++;
        }
    }
    while (i < N)
    {
        rest.push_back(u[i]);
        i++;
    }

    return std::make_pair(std::move(intersection), std::move(rest));
}

std::vector<int> vector_union(const std::vector<int> &u, const std::vector<int> &v)
{
    int N = u.size();
    int M = v.size();
    int i = 0, j = 0;
    std::vector<int> v_union = {};

    while (i < N && j < M)
    {
        if (u[i] == v[j])
        {
            v_union.push_back(u[i]);
            i++;
            j++;
        }
        else if (u[i] > v[j])
        {
            v_union.push_back(v[j]);
            j++;
        }
        else
        {
            v_union.push_back(u[i]);
            i++;
        }
    }
    while (i < N)
    {
        v_union.push_back(u[i]);
        i++;
    }
    while (j < M)
    {
        v_union.push_back(v[j]);
        j++;
    }

    return std::move(v_union);
}

/* Will return vector containing all element of v which are < threshold */
std::vector<int> vector_values_below(const std::vector<int> &v, const int threshold)
{
    std::vector<int> u;
    for (int a : v)
    {
        if (a >= threshold)
            break;
        u.push_back(a);
    }
    return std::move(u);
}

/* Will return vector containing all element of v which are >= threshold */
std::vector<int> vector_values_above(const std::vector<int> &v, const int threshold)
{
    std::vector<int> u;
    for (int a : v)
    {
        if (a < threshold)
            continue;
        u.push_back(a);
    }
    return std::move(u);
}

/* Will return vector containing all elements x, such that min <= x < max */
std::vector<int> vector_values_between(const std::vector<int> &v, const int min, const int max)
{
    std::vector<int> u;
    for (int a : v)
    {
        if (a < min)
            continue;
        else if (a >= max)
            break;
        u.push_back(a);
    }
    return std::move(u);
}

std::vector<int> set_to_vector(const std::set<int> &set)
{
    std::vector<int> v;
    for (auto i : set)
    {
        v.push_back(i);
    }

    return std::move(v);
}

bool vector_contains(const std::vector<int> &v, const int e)
{
    for (auto i : v)
    {
        if (i == e)
            return true;
        else if (i > e)
            return false;
    }
    return false;
}

/* All elements which are in u or v but not both */
std::vector<int> vector_symmetric_difference(const std::vector<int> &u, const std::vector<int> &v)
{
    int N = u.size();
    int M = v.size();
    int i = 0, j = 0;
    std::vector<int> sym_diff = {};

    while (i < N && j < M)
    {
        if (u[i] == v[j])
        {
            i++;
            j++;
        }
        else if (u[i] > v[j])
        {
            sym_diff.push_back(v[j]);
            j++;
        }
        else
        {
            sym_diff.push_back(u[i]);
            i++;
        }
    }
    while (i < N)
    {
        sym_diff.push_back(u[i]);
        i++;
    }
    while (j < M)
    {
        sym_diff.push_back(v[j]);
        j++;
    }

    return std::move(sym_diff);
}

std::set<int> vector_set_intersection(const std::vector<int> &u, const std::set<int> &v)
{
    std::set<int> intersection;
    std::set_intersection(u.begin(), u.end(),
                          v.begin(), v.end(),
                          std::inserter(intersection, intersection.begin()));

    return std::move(intersection);
}

std::set<int> set_intersection(const std::set<int> &u, const std::set<int> &v)
{
    std::set<int> intersection;
    std::set_intersection(u.begin(), u.end(),
                          v.begin(), v.end(),
                          std::inserter(intersection, intersection.begin()));

    return std::move(intersection);
}
