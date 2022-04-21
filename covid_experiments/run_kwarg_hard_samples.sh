#!/bin/bash


for i in B100 B200 C500 B500 A1000 C1000 B1000 Aall Call Ball 
do
    for t in 60 100
    do

        ../source/kwarg -T${t} -Q1 -k -r260 -Okwarg_records_sample${i}.csv -S-1,3.0,2.0,1.5,1.1,1.0,0.9,0.7,0.5,0.3,0.1,0.01,1 \
                                                                            -M-1,3.01,2.01,1.51,1.11,1.01,0.91,0.71,0.51,0.31,0.11,0.02,1.1 \
                                                                            -R1,1,1,1,1,1,1,1,1,1,1,1,-1 -C2,2,2,2,2,2,2,2,2,2,2,2,-1 < sample${i}.txt
    done
    echo "completed kwarg ${i}"
done