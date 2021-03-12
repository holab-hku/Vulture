#!/bin/bash
echo Output path is:;
echo $1;
echo intermediate files path is:;
echo $2;
echo Cell ranger result path is:;
echo $3;

# mkdir $1
# mkdir $1/STARsolo_outs
# mkdir $1/STARsolo_outs/Gene
mkdir -p $1/STARsolo_outs/Solo.out/Gene/raw
mkdir -p $1/STARsolo_outs/Solo.out/Gene/filtered_feature_bc_matrix

cp $3/outs/raw_feature_bc_matrix/*gz $1/STARsolo_outs/Solo.out/Gene/raw
cp $3/outs/filtered_feature_bc_matrix/*gz $1/STARsolo_outs/Solo.out/Gene/filtered_feature_bc_matrix

cp -r $2 $1


