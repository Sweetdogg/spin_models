[toc]

[typora说明文档](https://support.typoraio.cn/zh/Markdown-Reference/)

[markdown官方文档](https://markdown.com.cn/extended-syntax/footnotes.html)

# 每天Linux命令

## usually

- ls -a	ls -R	ls -l	ls -S #按大小排 ls !(||)
- cd - #返回上一次的目录
- {32..88..8} #间隔8,从32到88 
- history命令 !num
- '这是文本'  "\${变量},$(命令)"
- echo -e "$temp " " \c"
- sort -k2 -n #n是按数字大小
- ln         -s #软链接  快捷方式 备份
- mv -v -b -i -f -n

## linux命令行快捷键

- ctrl + w #前面的单词
- ctrl + u 
- ctrl + k
- ctrl + y 粘贴
- ctrl + r  $\rightarrow$选择命令，enter执行命令 
- ctrl + a or ctrl + e

## awk

> awk 'END {print NR}' file
>
> awk '!/^$/' input.txt > out #删除空行 NF#

## cut

> cut -b 13- out.dat > xxx
>
> find -name out.dat -exec sh -c 'cd $(dirname $0);grep "6.00  1.00" out.dat |cut -b 15- >xxx;cd ../.. ' {} \;

## find

- -name
- -iname
- -type
- -maxdepth

> find -name out.dat | xargs rm -v #相当于rm后列举找到所有文件

## tar

- -cvf 压缩
- -xvf 解压

>find -name temp0.8 -type d | xargs tar -cvf file.tar
>
>~~find -name temp0.8 -type d -exec tar -cvf file.tar {} \;~~  

## grep

> grep -v "NaN" out.dat > new.dat #将out.dat中的nan行剔除
>
> grep -rv 0.000000000000 {32..88..16} >out
>
> find -name out.dat -exec sh -c 'cd $(dirname $0);grep -v 0.000000000000 out.dat >outt.dat ; cd ../..' {} \;

## sed

> sed -n '1~3p' data >using
>
> sed -n '1,9p' data >using
>
> awk 'NR % 2 == 0' data_using > using
>
> sed -i '8,99!d' out.dat

## Examples

### 随机数

```sh
$ for((i=1;i<10;i=i+1));do bins=10;echo "$((RANDOM % $bins +1))"; done
```

### 改文件名字

> find -name "a.out" -exec sh -c 'mv "$0" "\${0%/*}/new_filename"' {} \;
>
> find -name spin -exec sh -c 'mv \$0 \$(dirname $0)/s' {} \;

### 计算误差

```sh
#!/bin/bash
binstep=4
comp=0
total_records=$(awk 'END {print NR}' out.dat)
while(( binstep <= 100 ))
do
    step=1
    while ((step <= total_records))
    do
        comp=$((binstep + step ))
        sed -n "${step},${comp}p" out.dat >> new.dat
        step=$((binstep + step + 1))
    done
    awk '{ sum += $2;d+=$2^2 } END { printf "%5.3f %18.12f %18.12f\n", $1, sum/NR , sqrt(d/NR-(sum/NR)^2) }' new.dat >> newnew.dat
    rm new.dat
    binstep=$(binstep + 5)
done
```

### sweep.sh

```sh
#! /bin/bash

x_i=0.847
x_f=0.847
x_step=0.01
mkdir -p data
rm -i data_using
for ll in {32..80..8}
do
	temp=$x_i
	rm -v data/data_$ll
	while(( $(echo "$temp <= $x_f" | bc ) ))
	do
	    awk -f cal_bins.awk $ll/temp$temp/out.dat >> data/data_$ll
		temp=$(echo "$temp + $x_step" | bc)
	done
	awk -v ll="$ll" '{print ll, $1, $5, $6}' data/data_"$ll" >> data_using
done
```

### bootstrap.sh

```sh
#! /bin/bash
x_i=0.841
x_f=0.859
x_step=0.002

mkdir -p data_boot
for ll in 24 27 30
do
	temp=$x_i
	while(( $(echo "$temp <= $x_f" | bc ) ))
	do
		total_records=$(awk 'END {print NR}' out.dat)
		bins=$total_records
		samples=$total_records
        for ((i = 1; i <= $samples; i++)); do
            for ((j = 1; j <= $bins; j++)); do
                random_number=$((RANDOM % 10 + 1))
                sed -n "${random_number}p" $ll/temp$temp/out.dat >> $ll/temp$temp/out_bootstrap.dat
            done
            awk '{ sum += $5 } END {printf "%5.3f %18.12f\n",$1, sum/NR }' $ll/temp$temp/out_bootstrap.dat >> $ll/temp$temp/bootstrap.dat
            rm $ll/temp$temp/out_bootstrap.dat
        done
	    awk '{ sum += $2; d += $2^2 } END { printf "%5.3f %18.12f %18.12f\n",$1, sum/NR, sqrt(d/NR-(sum/NR)^2)}' $ll/temp$temp/bootstrap.dat >> data_boot/data_$ll
        rm $ll/temp$temp/bootstrap.dat
		temp=$(echo "$temp + $x_step" | bc)
	done
	awk -v ll="$ll" '{print ll, $1, $2, $3}' data_boot/data_"$ll" >> data_using
done
```





