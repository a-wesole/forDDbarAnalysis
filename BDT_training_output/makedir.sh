#!/bin/bash

declare -a pt_list=("pt_1_2" "pt_2_3" "pt_3_4" "pt_4_5" "pt_5_6" "pt_6_8" "pt_8_10" "pt_10_15" "pt_15_20" "pt_20_30" "pt_30_50" "pt_50_100" "pt_100_1000")
declare -a cent_list=("cent_0_10" "cent_10_30" "cent_30_50" "cent_50_90")
for val in ${pt_list[@]};do
	for val1 in ${cent_list[@]};do
		mkdir $val\_$val1
	done
#echo $val
done
