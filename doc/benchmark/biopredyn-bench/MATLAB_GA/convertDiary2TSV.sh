#!/bin/sh

n_repeat=5

for name in B2 B4 B5 B6
do
    for i in `seq $n_repeat`
    do
	infilename=MATLAB_GA_${name}_diary_${i}.txt
	outfilename=MATLAB_GA_${name}_TSV_${i}.dat
	echo Converting $infilename To $outfilename ...
	cat $infilename | sed -e "s/ \+/\t/g" | sed -e "s/^[\t]*//g" | grep -v Gen | grep -v Best | grep -v Optim | grep -v Solution | grep -v Elapsed | sed "/^$/d" | grep -v idum | grep -v '*' > $outfilename
    done
done
