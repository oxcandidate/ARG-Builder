#!/bin/bash

for j in 200 500
do
    for i in A B C 
    do
        # first attempt with recomb_cost around 1
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q15 -R0.9,1.0,1.1 -B0.0,0.3,0.5,1.0,2.0 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q15 -R0.9,1.0,1.1 -B0.0,0.3,0.5,1.0,2.0 < sample${i}${j}.txt
        
        # Now try with rm_cost cheaper
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q25 -R1.5,2.0,3.0,4.0 -B0.1,0.3,0.5,1.0,2.0,50.0 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q25 -R1.5,2.0,3.0,4.0 -B0.1,0.3,0.5,1.0,2.0,50.0 < sample${i}${j}.txt
        
        # Now look for pure recurrent mutation solution
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q25 -R-1 -B0.1,0.5,1.0 < sample${i}${j}.txt

        # Now look for mostly recomb solution
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q20 -M50 -B10 < sample${i}${j}.txt
        echo "completed kwarg ${i}${j}"
    done
done
