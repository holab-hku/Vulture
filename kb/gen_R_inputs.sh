#!/bin/bash
#!/bin/bash
if [ "$#" -ne 3 ]; then
  echo "$#";  
  echo "Copy the output of the kb file to fite the formate of the Droputils";
  echo "Usage: $0 <Output_directory> <Intermediate_directory> <KB_directory>" >&2
  exit 1
fi
echo Output path is:;
echo $1;
echo intermediate files path is:;
echo $2;
echo KB result path is:;
echo $3;

# mkdir $1
# mkdir $1/STARsolo_outs
# mkdir $1/STARsolo_outs/Gene
mkdir -p $1/STARsolo_outs/Solo.out/Gene/raw
mkdir -p $1/STARsolo_outs/Solo.out/Gene/filtered

cp $3/outs/counts_unfiltered/cells_x_genes.barcodes.txt $1/STARsolo_outs/Solo.out/Gene/raw/barcodes.tsv;
cp $3/outs/counts_unfiltered/cells_x_genes.genes.txt $1/STARsolo_outs/Solo.out/Gene/raw/genes.tsv;
cp $3/outs/counts_unfiltered/cells_x_genes.mtx $1/STARsolo_outs/Solo.out/Gene/raw/matrix.mtx;
cp -r $2 $1


