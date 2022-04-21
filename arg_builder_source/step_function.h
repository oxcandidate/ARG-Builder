#ifndef STEP_FUNCTION_H
#define STEP_FUNCTION_H

#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <algorithm>
#include <iostream>

typedef struct _StepFunction
{
    std::set<int> change_points; // these are the keys for the following map
    std::map<int, std::tuple<int, int, float>> values; // rms, bms, cost
    int max_rms;
    int max_bms;
    int max_cost;

    _StepFunction()
    {
    }

    _StepFunction(int rms, int bms)
    {
        max_rms = rms;
        max_bms = bms;
    }

    void add_point(int pos, int rms, int bms, float cost)
    {
        auto result = change_points.insert(pos);
        if (!result.second)
        {
            std::cerr << "value already assigned at that location";
        }
        values.insert({pos, std::make_tuple(rms, bms, cost)});
    }

    bool add_point_if_better(int pos, int rms, int bms, float cost)
    {
        auto result = change_points.insert(pos);
        if (!result.second)
        {
            auto [current_rms, current_bms, current_cost] = values[pos];
            if (cost < current_cost)
            {
                values[pos] = std::make_tuple(rms, bms, cost);
                return true;
            }
            else
                return false;
        }
        values[pos] = std::make_tuple(rms, bms, cost);
        return true;
    }

    std::tuple<int, int> get_prefix(int pos)
    {
        auto search = change_points.lower_bound(pos); // Either returns pos or the next highest in the set.
        if (search != change_points.end())
        {
            auto [rms, bms, cost] = values[*search];
            return std::make_tuple(rms, bms);
        }
        else
        {
            return std::make_tuple(max_rms, max_bms);
        }
    }

    std::tuple<int, int> get_suffix(int pos)
    {
        auto search = change_points.lower_bound(pos); // Either returns pos or the next highest in the set.
        if (search != change_points.end())
        {
            auto [rms, bms, cost] = values[*search];
            return std::make_tuple(max_rms - rms, max_bms - bms);
        }
        else
        {
            return std::make_tuple(0, 0);
        }
    }

    std::tuple<int, int> get_range(int start, int end)
    {
        auto [start_rms, start_bms] = get_prefix(start);
        auto [end_rms, end_bms] = get_prefix(end);
        return std::make_tuple(end_rms - start_rms, end_bms - start_bms);
    }

} StepFunction;

#endif