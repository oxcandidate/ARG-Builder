#!/bin/bash

for r in 0.5 1.0 1.1 1.5 2.0
do
    for b in 1.0 1.1 1.5 2.0
    do
        ../arg_builder_source/arg_builder -l -V0 -c -othreaded_records.csv -L0 -Q5 -F3 -R${r} -B${b} < kreitman.txt
        ../arg_builder_source/arg_builder -l -V0 -c -othreaded_records.csv -L1 -Q10 -F3 -R${r} -B${b} < kreitman.txt
    done
done