#!/bin/bash

../arg_builder_source/arg_builder -l -V0 -c -r1 -othreaded_rooted_records.csv -L0 -Q5 -R0.5,1.0,1.1,1.5,2.0 -B1.0,1.1,1.5,2.0 < kreitman_rooted.txt
../arg_builder_source/arg_builder -l -V0 -c -r1 -othreaded_rooted_records.csv -L1 -Q10 -R0.5,1.0,1.1,1.5,2.0 -B1.0,1.1,1.5,2.0 < kreitman_rooted.txt