import scanpy as sc
import scipy.io as io
import pandas as pd
import argparse
import numpy as np
from pathlib import Path
import logging

def run_main(args):
    # Read h5ad files
    results_file = args.output
    adata = sc.read_h5ad(args.input) 
    # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx` 
    if (args.vn_mk_unique==True):  
        adata.var_names_make_unique() 
    # Basic filterings
    if (args.min_genes!=0):
        sc.pp.filter_cells(adata, min_genes=200)
    if (args.min_cells!=0):
        sc.pp.filter_genes(adata, min_cells=3) 

    # Assemble some information about mitochondrial genes, which are important for quality control.
    adata.var['mt'] = adata.var_names.str.startswith(args.mito_prefix)  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Basic filterings
    if (args.pct_ngbc!=0):
        adata = adata[adata.obs.n_genes_by_counts < args.pct_ngbc, :]
    if (args.pct_mt!=0):
        adata = adata[adata.obs.pct_counts_mt < args.pct_mt, :]

    # Normalization
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)

    adata.raw = adata
    # Calculate hvgs
    if (args.hvg==True):  
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var.highly_variable]
    # Regress out
    if (args.reg_out==True):  
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

    # Scaling
    sc.pp.scale(adata, max_value=10)
    # PCA and neighbor graphs
    sc.tl.pca(adata, svd_solver='arpack',n_comps=args.n_comps)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_comps)

    # UMAP, 3d umap and 2d umap. X_umap2d == X_umap
    sc.tl.umap(adata,n_components=3)
    adata.obsm["X_umap3d"] = adata.obsm["X_umap"]
    adata.obsm["X_umap2d"] = adata.obsm["X_umap"]

    # Leiden clustering
    sc.tl.leiden(adata,resolution=args.resolution)

    # Write file
    adata.write(results_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process h5ad AnnData file with neccessary processings.')
    # Basic io
    parser.add_argument('--input','-i', required=True,type=str, help='An input path of a h5ad data.')
    parser.add_argument('--output', '-o', default="result.h5ad", type=str,  help='An output path of a h5ad data. DEAFULT: "result.h5ad"')
    # Basic filtering
    parser.add_argument('--vn_mk_unique',action='store_true',default=False, help='Process var_names_make_unique in scanpy.')
    parser.add_argument('--min_genes',action='store',default=0, help='If not 0, filter minimum number of genes expressed required for a cell to pass filtering, see scanpy.pp.filter_cells. DEAFULT: 0 (do not filter)')
    parser.add_argument('--min_cells',action='store',default=0, help='If not 0, filter minimum number of cells expressed required for a gene to pass filtering, see scanpy.pp.filter_genes. DEAFULT: 0 (do not filter)')
    parser.add_argument('--mito_prefix',action='store',default="MT-", help='Mitochondrial prefix to to access info about mitochondrial genes. DEAFULT: "MT-"')
    parser.add_argument('--pct_mt',action='store',default=0, help='If not 0, filter cells by n_genes_by_counts to pass filtering. DEAFULT: 0 (do not filter)')
    parser.add_argument('--pct_ngbc',action='store',default=0, help='If not 0, filter cells by pct_counts_mt to pass filtering. DEAFULT: 0 (do not filter)')
    # Normalization
    parser.add_argument('--target_sum',action='store',default=10000, help='If None, after normalization, each observation (cell) has a total count equal to the median of total counts for observations (cells) before normalization. DEAFULT: 10000')
    # HVG related 
    parser.add_argument('--hvg',action='store_true',default=False, help='Annotate highly variable genes. DEAFULT: "False (not annotate)"')
    parser.add_argument('--hvg_min_mean',action='store',default=0.0125, help='Work if hvg == True, this and all other cutoffs for the means and the normalized dispersions are ignored.. DEAFULT: 0.0125')
    parser.add_argument('--hvg_max_mean',action='store',default=3, help='Work if hvg == True, this and all other cutoffs for the means and the normalized dispersions are ignored. DEAFULT: 3')
    parser.add_argument('--hvg_min_disp',action='store',default=0.3, help='Work if hvg == True, this and all other cutoffs for the means and the normalized dispersions are ignored. DEAFULT: 0.3')
    # Regress out
    parser.add_argument('--reg_out',action='store_true',default=False, help='Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.')
    # PCA     
    parser.add_argument('--n_comps',action='store',default=50, help='Number of principal components to compute. Defaults to 50, or 1 - minimum dimension size of selected representation. DEAFULT: 50')
    parser.add_argument('--n_neighbors','-nn',action='store',default=10, help='Number of principal components to compute. Defaults to 50, or 1 - minimum dimension size of selected representation. DEAFULT: 50')
    # Clustering 
    parser.add_argument('--resolution',action='store',default=0.4, help='A parameter value controlling the coarseness of the clustering. DEAFULT: 0.4')

    args, unknown = parser.parse_known_args()
    run_main(args)