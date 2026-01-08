#! /bin/bash

x_i=0.843
x_f=0.859
x_step=0.002
mkdir -p data
for ll in 24 27 30
do
	temp=$x_i
	while(( $(echo "$temp <= $x_f" | bc ) ))
	do
	    awk -f cal_bins.awk $ll/temp$temp/out.dat >> data/data_$ll
		temp=$(echo "$temp + $x_step" | bc)
	done
	awk -v ll="$ll" '{print ll, $1, $5, $6}' data/data_"$ll" >> data_using
done




