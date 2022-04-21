#!/bin/bash

../arg_builder_source/arg_builder -l -V2 -c -L0 -e1 -r1 -S1751130457 -dthreaded_network.dot -R1.1 -B1 -M1 < kreitman_rooted.txt

../source/kwarg -T50 -Z2093842884 -e -k -S1 -M1.01 -R1 -C2 -dkwarg_network.dot < kreitman_zero_rooted.txt

../source/kwarg -T50 -Z2464303795 -e -S1 -M1.01 -R1 -C2 -dkwarg_alternative_network.dot < kreitman.txt

./dot_to_png.sh