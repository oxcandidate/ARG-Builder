#!/bin/bash

for i in A B C
do
for j in 25 50 100 200 500 1000 all
do
python3 ../python_scripts/clean_binary_sequences.py -i sample${i}${j}.txt -o clean_sample${i}${j}.txt -r
done
done