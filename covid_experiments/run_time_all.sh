#!/bin/bash

for j in 25 50 100 200 500 1000
do
    for i in A B C 
    do
        python3 ../python_scripts/sum_run_times.py -i thread_records_sample${i}${j}.csv
    done
done

for j in A25 A50 A100 A200 A500 B25 B50 C25 C50 C100 C200
do
python3 ../python_scripts/sum_run_times.py -i kwarg_records_sample${j}.csv
done