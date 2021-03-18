import scanpy as sc
import scipy.io as io
import pandas as pd
import argparse
from pathlib import Path
import logging

def run_main(args):
    # Read h5ad files
    gene_name_dict={"10XV2":"genes.tsv","10XV3":"features.tsv"}
    try:
        adata = sc.read_h5ad(args.input)
    except:
        logging.error("Cannot read h5ad file from:" + args.input)

    barcodes = adata.obs
    genes = adata.var
    
    # Make output directory
    try:
        Path(args.output).mkdir(parents=True, exist_ok=True)
    except:
        logging.error("Cannot make output path:" + args.output)

    #Write matrix 
    io.mmwrite(args.output+"\\matrix.mtx",adata.X.T.astype("int"))
    barcodes.to_csv(args.output+"\\barcodes.tsv",sep="\t",header=None)
    genes.to_csv(args.output+"\\"+gene_name_dict[args.chemistry],sep="\t",header=None)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convser h5ad AnnData file to 10x mtx files.')
    parser.add_argument('--input','-i', type=str, help='an integer for the accumulator')
    parser.add_argument('--output', '-o',  type=str,  help='An output folder')
    parser.add_argument('--chemistry', '-v',  type=str,default="10XV2",  help='10x chemistry vesion')

    args, unknown = parser.parse_known_args()
    run_main(args)