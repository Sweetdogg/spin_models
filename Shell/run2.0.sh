#! /bin/bash
bins=1
mc_eq=1000000
mcs=1000000

j4=0

x_i=.76000
x_f=.79000
x_step=0.00100
STOP_FILE="stop"

ifort atv2.f90
rm log
start_time=$(date +%s)
until [[ -f "$STOP_FILE" ]]; do
	cur=$(squeue -u $USER|wc -l)
	if (( $cur < 5 )); then
		for ll in 24 32 40 48 56 64 72
		do
			temp=$x_i
			mkdir -p $ll ; cd $ll
			
			while(( $(echo "$temp <= $x_f" | bc ) ))
			do
				echo $(date +%N) > ../seed.in
				mkdir -p temp$temp ; 
				cp ../a.out  ../job ../seed.in  temp$temp
				cd temp$temp
				echo "$ll $temp $j4 $mc_eq $bins $mcs" > read.in
				#./a.out  &
				sbatch -J $ll-$temp job
				#echo -e "$temp " " \c"
				temp=$(echo "$temp + $x_step" | bc)
				cd ..
				sleep 0.1
			done
			echo "done"
			cd ..
		done
	else
		current_time=$(date +%s)
        elapsed_seconds=$(((current_time - start_time)/ 60))
		echo "当前任务数 $cur ，已运行 $elapsed_seconds 分钟..." > log
		sleep 60
	fi
done

echo "检测到停止文件 '$STOP_FILE'，正在退出循环..." >> log
rm -f "$STOP_FILE"

