"""
USAGE:-

python bam_concat.py -o <output_file> [input_files]  

"""
import os
import argparse
import pysam

parser = argparse.ArgumentParser(description = 'Generate BAM with CB/UB tag')
parser.add_argument('bams', metavar='B', type=str, nargs='+',help=': list of BAM files')
parser.add_argument('-o', '--output', required=True, help='Output BAM file')

args = parser.parse_args()

bam_list = args.bams

outfile_path = os.path.abspath(args.output)

for infile in bam_list:
    infile_path = os.path.abspath(infile)
    infile = pysam.AlignmentFile(infile_path, "rb")
    outfile = pysam.AlignmentFile(outfile_path, "wb", template=infile)
    iter = infile.fetch(until_eof=True)
    for read in iter:
        outfile.write(read)
    infile.close()

outfile.close()