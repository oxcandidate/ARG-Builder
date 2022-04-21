#!/bin/bash

for i in C
do
    for j in 200 
    do
        ./arg_builder -l -V1 -c -orandomselection_sample${i}${j}_records.csv -L1 -r1 -R0.5,0.6,0.9,1.0,1.1,1.5,1.9,2.0,3.0,50.0 -B1.0,1.3,1.5,2.0,3.0,50.0 -Q10 < sampleref${i}${j}.txt
    done
done
