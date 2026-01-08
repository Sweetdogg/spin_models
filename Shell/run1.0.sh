#! /bin/bash
bins=5
mc_eq=1000000
mcs=1000000

j4=0.2

x_i=.40000
x_f=.70000
x_step=0.01000

ifort atv2.f90
#rm -r  8 
for ll in 16 24 32 40 48 56 64
do
	temp=$x_i
	mkdir -p $ll ; cd $ll
	
	while(( $(echo "$temp <= $x_f" | bc ) ))
	do
		echo $(date +%N) > ../seed.in
		mkdir -p temp$temp ; 
		cp ../a.out  ../job ../seed.in  temp$temp ;  # copy a.out to temp folder
		cd temp$temp
		echo "$ll $temp $j4 $mc_eq $bins $mcs" > read.in
		#./a.out  &
		qsub -N $ll-$temp job
		#echo -e "$temp " " \c"
		temp=$(echo "$temp + $x_step" | bc)
		cd ..
		sleep 0.2
	done
	echo "done"
	cd ..
done

