#!/bin/bash

for j in 25 50 100 200 500 1000
do
    for i in A B C 
    do
        python3 ../python_scripts/graph_results.py -i thread_records_sample${i}${j}.csv -o thread_sample${i}${j}_graph.png -g thread_sample${i}${j}_grid.png
    done
done

for i in A25 B25 C25 A50 B50 C50 A100 B100 C100 A200 C200 A500 
do
    python3 ../python_scripts/graph_results.py -i kwarg_records_sample${i}.csv -o kwarg_sample${i}_graph.png

    python3 ../python_scripts/graph_results.py -i thread_records_sample${i}.csv,kwarg_records_sample${i}.csv \
                                               -o comparison_sample${i}.png -l threaded,kwarg
done