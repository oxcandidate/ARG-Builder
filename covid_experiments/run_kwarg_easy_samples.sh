#!/bin/bash


for i in A25 C25 B25 A50 A100
do
    for t in 30 50 100
    do

        ../source/kwarg -T${t} -Q5 -k -r100 -Y150 -Okwarg_records_sample${i}.csv -S-1,3.0,2.0,1.5,1.2,1.1,1.0,0.9,0.7,0.5,0.3,0.1,0.01,1 \
                                                                                        -M-1,3.01,2.01,1.51,1.21,1.11,1.01,0.91,0.71,0.51,0.31,0.11,0.02,1.1 \
                                                                                        -R1,1,1,1,1,1,1,1,1,1,1,1,1,-1 -C2,2,2,2,2,2,2,2,2,2,2,2,2,-1 < sample${i}.txt
    done
    echo "completed kwarg ${i}"
done