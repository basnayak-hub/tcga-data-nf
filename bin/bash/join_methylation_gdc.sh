#!/bin/bash

outfile="hello.txt"
tmp="tmp.txt"

head -1 ../r/tmpres/TCGA_methylation_paths.txt > "firstfile.txt"

while read p; do
    cat $p | awk '{print $1}' > $outfile
done < firstfile.txt

while read p; do
    echo "$p"
    join $outfile $p > $tmp
    mv $tmp $outfile
done < ../r/tmpres/TCGA_methylation_paths.txt
