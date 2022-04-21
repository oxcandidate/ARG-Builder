#!/bin/bash

for j in 25 50 100 200 500 1000 all
do
    for i in A B C 
    do
        python3 ../python_scripts/graph_results.py -i thread_records_sample${i}${j}.csv -b lower_bound${i}${j}.csv -o thread${i}${j}_lb_graph.png -t
    done
done

for i in A25 B25 C25 A50 B50 C50 A100 B100 C100 A200 C200 A500 
do
    python3 ../python_scripts/graph_results.py -i kwarg_records_sample${i}.csv -b lower_bound${i}.csv -o kwarg${i}_lb_graph.png -t

    python3 ../python_scripts/graph_results.py -i thread_records_sample${i}.csv,kwarg_records_sample${i}.csv \
                                               -b lower_bound${i}.csv -o comparison${i}lb.png -l threaded,kwarg -t
done