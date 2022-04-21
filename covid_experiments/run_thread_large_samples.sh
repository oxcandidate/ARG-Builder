#!/bin/bash

for j in all
do
    for i in A C B 
    do
        # first attempt with recomb_cost around 1.   2*4*10*2 = 160
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q10 -R0.9,1.0 -B0.1,0.5,1.0,2.0 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q10 -R0.9,1.0 -B0.1,0.5,1.0,2.0 < sample${i}${j}.txt

        # Now try with rm_cost cheaper   4*6*20*2 = 960
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q20 -R2.0,3.0,4.0,5.0 -B0.1,0.3,0.5,1.0,2.0,100.0 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q20 -R2.0,3.0,4.0,5.0 -B0.1,0.3,0.5,1.0,2.0,100.0 < sample${i}${j}.txt

        # Now look for pure recurrent mutations solution  3*5*2 = 30
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q5 -R-1 -B0.1,0.5,1.0 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q5 -R-1 -B0.1,0.5,1.0 < sample${i}${j}.txt

        # Now look for mostly recomb mutations solution 5
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q5 -M50 -B20 < sample${i}${j}.txt
        echo "completed kwarg ${i}${j}"
    done
done
