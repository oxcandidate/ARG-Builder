#!/bin/bash

for j in 25 50 100
do
    for i in A B C 
    do
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q25 -R0.5,0.6,0.9,1.0,1.5,1.9,2.0,3.0,4.0,50.0 -B0.0,0.3,0.5,1.0,2.0,50.0 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q25 -R0.5,0.75,1.0,1.5,2.0,3.0,4.0,50.0 -B0.0,0.3,0.5,1.0,2.0,50.0 < sample${i}${j}.txt
        
        # Purely recurrent-mutation
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q25 -R-1 -B0.1,0.5,1.0 < sample${i}${j}.txt

        # Mostly recomb
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q10 -M50 -B10 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q10 -M50 -B10 < sample${i}${j}.txt
        echo "completed kwarg ${i}${j}"
    done
done
