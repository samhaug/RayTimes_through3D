#!/bin/bash

infile=$1

csplit --digits=4  --quiet --prefix=outfile $infile "/>/+1" "{*}"
for ii in outfile*;
do
    head -n -1 $ii > temp.txt ; mv temp.txt $ii;
done



