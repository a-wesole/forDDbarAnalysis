#!/bin/bash

header_file="header.h"
> $header_file

declare -a pt_list=("pt_1_2" "pt_2_3" "pt_3_4" "pt_4_5" "pt_5_6" "pt_6_8" "pt_8_10" "pt_10_15" "pt_15_20")
declare -a cent_list=("cent_0_10" "cent_10_30" "cent_30_50" "cent_50_90")
for val in ${pt_list[@]};do
	for val1 in ${cent_list[@]};do

    dir_name=${val}_${val1}

    echo "#include \"/home/awesole/forDDbarAnalysis/BDT_training_output/${dir_name}/dataset/weights/TMVAClassification_BDT_${val}_${val1}.class.C\"" >> $header_file

	done
#echo $val
done
