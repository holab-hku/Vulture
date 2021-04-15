#!/usr/bin/env python3
import argparse
from contextlib import ExitStack
import os
import collections
import sys
from subprocess import call

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Generate count matrix based on splitted BUS files')
parser.add_argument('directory', metavar='D', type=str, nargs='+',help=': list of kallisto bus output directories, each of them have a .ec file and a .txt BUS file')
parser.add_argument("--t2g",required=True, help=": a 3-column file with transcript, gene_ID and gene name information")
parser.add_argument("--output","-o",default="./", help=": the output path of the matrix")
parser.add_argument("--ec_name",default='matrix.ec', help=": name of the equivalent to transcripts map file", type=str)
parser.add_argument("--bus_name",default='output.unfiltered.txt', help=": name of the text formated BUS file", type=str)
parser.add_argument("--gene_min",default=200, help=": minimal number of genes detected", type=int)
parser.add_argument("--gene_max", default=5000,help=": maximal number of genes detected", type=int)
parser.add_argument("--len_BC", default=16,help=": bacrode length", type=int)
parser.add_argument("--len_UMI", default=10,help=": UMI length", type=int)
parser.add_argument("--min_g", default=10,help=": minimum genes per cell accepted", type=int)

args = parser.parse_args()

LEN_BC = args.len_BC
LEN_UMI = args.len_UMI
T2G = args.t2g
DIRECTOIRES = args.directory
OUT = args.output
MIN_G = args.min_g

# Generate Transcript to genes
tr2g = {}
trlist = []
with open(T2G) as f:
    for line in f:
        l = line.split()
        tr2g[l[0]] = l[1]
        trlist.append(l[0])

genes = list(set(tr2g[t] for t in tr2g))

# load equivalence classes
ecs = {}
ecs_diff = {}
ec_filenames = [os.path.join(d,args.ec_name) for d in DIRECTOIRES]

# Can when contray
with ExitStack() as stack:
    ec_files = [stack.enter_context(open(i, "r")) for i in ec_filenames]
    zec_files = zip(*ec_files)
    
    for lines in zec_files:
        list_lines = list(lines)
        ls = [l.split() for l in list_lines]

        # Agreements among different ec
        if(all(ls[0] == x for x in ls)==True):
            l=ls[0]
            ec = int(l[0])
            trs = [int(x) for x in l[1].split(',')]
            ecs[ec] = trs

        # If ec is not agreed
        else:
            for i in range(0,len(ls)):
                l = ls[i]
                ec = int(l[0])
                trs = [int(x) for x in l[1].split(',')]
                ecs_diff[i,ec] = trs

print('Finished ec 2 transript calcu')

# Defined my own ec2g if there are different assignment in the ec file        
def ec2g(ec,i=None):
    # ec is in concencus
    if ec in ecs:
        return list(set(tr2g[trlist[t]] for t in ecs[ec]))
    # ec is different among bus files
    elif (i,ec) in ecs_diff:
        return list(set(tr2g[trlist[t]] for t in ecs_diff[i,ec]))
    else:
        return []

########################################################################## Quant the bc-gene matrix ##########################################################################
cell_gene = collections.defaultdict(lambda: collections.defaultdict(float))
pbar=None
pumi=None

filenames = [os.path.join(d,args.bus_name) for d in DIRECTOIRES]
finfile=0
maxbc = 'Z'
with ExitStack() as stack:
    files = [stack.enter_context(open(i, "r")) for i in filenames]
    zfiles = zip(*files)
    rows = list(next(zfiles))

    gs = set()
    while finfile<len(filenames):
        try:
            # Can use min to find the record because it is composed by BC UMI EC COUNT
            min_id = np.nanargmin(rows)
            line = rows[min_id]
            l = line.split()
            barcode,umi,ec,count = line.split()
            ec = int(ec)
            
            if barcode == pbar:
                # same barcode
                if umi == pumi:
                    # same UMI, let's update with intersection of genelist
                    gl = ec2g(ec,min_id)
                    gs.intersection_update(gl)
                else:
                    # new UMI, process the previous gene set
                    for g in gs: 
                        cell_gene[barcode][g] += 1.0/len(gs)
                    # record new umi, reset gene set
                    pumi = umi
                    gs = set(ec2g(ec,min_id))
            else:
                # work with previous gene list
                for g in gs:
                    cell_gene[pbar][g] += 1.0/len(gs)
                
                if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < 10:
                    del cell_gene[pbar]
                
                pbar = barcode
                pumi = umi
                
                gs = set(ec2g(ec,min_id))

            # Iterate over the selec record
            # If ec is in the concecus list, skip all with same ec
            if ec in ecs:
                post_before_cout = line.find("\t",LEN_BC+LEN_UMI+2)
                equal_BC_UMI_EC = [line[:post_before_cout] == r[:post_before_cout] for r in rows]  
                identical_indices = np.where(np.array(equal_BC_UMI_EC)==True)[0]
                for i in identical_indices:  
                    rows[i] = next(files[i])
            else:
                rows[min_id] = next(files[min_id])

        except StopIteration:
            finfile+=1
            rows[min_id] = maxbc
            min_id = np.nanargmin(rows)

        except:
            print('Error when reading from BUS files')

    for g in gs:
        cell_gene[pbar][g] += 1.0/len(gs)
        
    if sum(cell_gene[pbar][g] for g in cell_gene[pbar]) < MIN_G:
        del cell_gene[pbar]

########################################################################## Quant the bc-gene matrix ##########################################################################

barcode_hist = collections.defaultdict(int)
for barcode in cell_gene:
    cg = cell_gene[barcode]
    s = len([cg[g] for g in cg])
    barcode_hist[barcode] += s

#Output a gene count histogram
bcv = [x for b,x in barcode_hist.items() if x > args.gene_min and x < args.gene_max]
plt.switch_backend('agg')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(bcv,bins=100)
ax.set_title("Histogram")
plt.xlabel("number of genes detected")
plt.ylabel("number of barcodes")
fig.savefig(os.path.join(OUT,'gene_hist.png'))

outfile = os.path.join(OUT,'matrix.mtx')

gene_to_id = dict((g,i+1) for i,g in enumerate(genes))
barcodes_to_use = [b for b,x in barcode_hist.items() if x > args.gene_min and x < args.gene_max]

num_entries = 0
for barcode in barcodes_to_use:
    num_entries += len([x for x in cell_gene[barcode].values() if x>0])

with open(outfile, 'w') as of:
    of.write('%%MatrixMarket matrix coordinate real general\n%\n')
    #number of genes
    of.write("%d %d %d\n"%(len(genes), len(barcodes_to_use), round(num_entries)))
    bcid = 0
    for barcode in barcodes_to_use:
        bcid += 1
        cg = cell_gene[barcode]
        gl = [(gene_to_id[g],cg[g]) for g in cg if cg[g] > 0]
        gl.sort()
        for x in gl:
            of.write("%d %d %f\n"%(x[0],bcid,x[1]))

gene_names = {}
with open(T2G) as f:
    f.readline()
    for line in f:
        line_splitted = line.split()
        t,g,gn = line_splitted[0],line_splitted[1],line_splitted[2]
        gene_names[g] = gn

id_to_genes = dict((i,g) for (g,i) in gene_to_id.items())
gl = []
for i in range(1,len(genes)+1):
    g = id_to_genes[i]
    gid = g
#    gid = g[:g.find('.')]
    if gid in gene_names:
        gn = gene_names[gid]
    else:
        gn = ''
    gl.append((g,gn))

with open(os.path.join(OUT,'genes.tsv'),'w') as of:
    for g,gn in gl:
        of.write("%s\t%s\n"%(g,gn))
        
with open(os.path.join(OUT,'barcodes.tsv'),'w') as of:
    of.write('\n'.join(x + '' for x in barcodes_to_use))
    of.write('\n')

