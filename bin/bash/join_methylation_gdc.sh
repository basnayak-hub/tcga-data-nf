#!/bin/bash

outfile="hello.txt"
tmp="tmp.txt"

head -1 $2 > "firstfile.txt"

while read p; do
    cat $p | awk '{print $1}' > $outfile
done < firstfile.txt

while read p; do
    echo "$p"
    join $outfile $p > $tmp
    mv $tmp $outfile
done < $2 #tcga_luad_methylation_paths.txt

mv $outfile $1

