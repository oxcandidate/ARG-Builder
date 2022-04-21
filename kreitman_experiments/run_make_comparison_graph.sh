#!/bin/bash

python3 ../python_scripts/graph_results.py \
    -i kwarg_records.csv,threaded_records.csv,kwarg_rooted_records.csv,threaded_rooted_records.csv \
    -b lower_bounds_kreitman.csv \
    -o kreitman_comparison_lb.png \
    -l kwarg,threading,rooted_kwarg,rooted_threading