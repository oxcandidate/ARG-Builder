#!/bin/bash

# for i in A C B
# do
# for j in 25 50 100 200
# do
# python3 ../python_scripts/hudson_bound.py -i clean_sample${i}${j}.txt -o lower_bound${i}${j}.csv -h -r4 -s4 -a5 -b10 -m
# done
# done

# for i in B
# do
# for j in 500 1000
# do
# python3 ../python_scripts/hudson_bound.py -i clean_sample${i}${j}.txt -o lower_bound${i}${j}.csv -h -r4 -s4 -a5 -b0 -m
# done
# done

for i in C B
do
for j in all
do
python3 ../python_scripts/hudson_bound.py -i clean_sample${i}${j}.txt -o lower_bound${i}${j}.csv -h -r4 -s0 -a-1 -b-1 -m
done
done