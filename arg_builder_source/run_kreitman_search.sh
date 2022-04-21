#!/bin/bash

for r in 0.5 1.0 1.5 2.0 2.5
do
    for b in 1.0 1.1 1.5 2.0 
    do
        ./arg_builder -l -V1 -c -othreaded_kreitman_records.csv -L1 -Q5 -F3 -R${r} -B${b} < kreitman_labelled.txt
    done
done