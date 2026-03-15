#! /bin/bash
bins=10
mc_eq=500000
mcs=100000

j4=0.05

x_i=.23050
x_f=.23250
x_step=0.00010
STOP_FILE="stop"

ifort atv2.f90
rm -f log
start_time=$(date +%s)
for i in {269..500}; do 
	[[ -f "$STOP_FILE" ]] && break

	if [[ ! -d "cl$i" ]]; then
      echo "WARN: cl$i 不存在，跳过" >> log
      continue
    fi

	
	cd cl$i || exit 1
	cur=$(squeue -u $USER|wc -l)
	cp ../a.out  ../job .	
	while true; do
		if (( $cur < 50 )); then
			break
		else
			cur=$(squeue -u $USER|wc -l)
		fi
		sleep 61
	done
	echo "Starting cl$i at $(date)" >> ../log
	for ll in 64 256 576 1024 2304 3136 4096 9216 12544 16384 36864
	do
		temp=$x_i
		mkdir -p $ll ; cd $ll
		
		while(( $(echo "$temp <= $x_f" | bc ) ))
		do
			echo $(date +%N) > ../seed.in
			mkdir -p temp$temp ; 
			cp ../a.out ../nbor${ll}.txt  ../job ../seed.in  temp$temp
			cd temp$temp
      mv nbor${ll}.txt nbor.txt
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
	cd ..
done

echo "检测到停止文件 '$STOP_FILE'，正在退出循环..." >> log
rm -f "$STOP_FILE"
echo "已完成，请检查。" | mail -s "$j4,$x_i to $x_f" xxx@qq.com

