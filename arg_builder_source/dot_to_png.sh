#!/bin/bash

files="*.dot"
echo $files

for i in *.dot
do 
    name=$(echo "$i" | cut -f 1 -d '.');
    dot -Tpng -o ${name}.png $i;
done