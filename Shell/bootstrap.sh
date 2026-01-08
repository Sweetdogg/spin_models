#! /bin/bash
x_i=1.00
x_f=1.00
x_step=0.002


mkdir -p data_boot
for ll in 32
do
	temp=$x_i
	while(( $(echo "$temp <= $x_f" | bc ) ))
	do
        total_records=$(awk 'END {print NR}' $ll/temp$temp/out.dat)
        bins=$total_records
        samples=$total_records
        for ((i = 1; i <= $samples; i++)); do
            for ((j = 1; j <= $bins; j++)); do
                random_number=$((RANDOM % $samples + 1))
                sed -n "${random_number}p" $ll/temp$temp/out.dat >> $ll/temp$temp/out_bootstrap.dat
            done
            awk '{ sum += $2^2 } END {printf "%5.3f %18.12f\n",$1, sum/NR }' $ll/temp$temp/out_bootstrap.dat >> $ll/temp$temp/bootstrap.dat
            rm $ll/temp$temp/out_bootstrap.dat
        done
	    awk '{ sum += $2; d += $2^2 } END { printf "%5.3f %18.12f %18.12f\n",$1, sum/NR, sqrt(d/NR-(sum/NR)^2)}' $ll/temp$temp/bootstrap.dat >> data_boot/data_$ll
        rm $ll/temp$temp/bootstrap.dat
		temp=$(echo "$temp + $x_step" | bc)
	done
	awk -v ll="$ll" '{print ll, $1, $2, $3}' data_boot/data_"$ll" >> data_using
done
