#!/bin/bash

for i in A B C
do
    for j in 25 50 100 200 500 1000 all
    do 
        python3 ../python_scripts/remove_reference_sequence.py -i sampleref${i}${j}.txt -o sample${i}${j}.txt
    done
done
