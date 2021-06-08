# Copyright (c) 2019 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
from pathlib import Path
import time
from Bio.SeqIO.QualityIO import FastqPhredIterator
import logging
from fastqsplitter import key_split_fastqs


# TEST_FILE = Path(__file__).parent / Path("data") / Path("SRR11537950_1_part0.fastq.gz")
# TEST_FILE_2 = Path(__file__).parent / Path("data") / Path("SRR11537950_2_part0.fastq.gz")

# TEST_FILE_INVALID = Path(__file__).parent / Path("data") / Path(
#     "test_invalid.fq.gz")
# number_of_splits = 3

# in this application.
DEFAULT_COMPRESSION_LEVEL = 1
# There is no noticable difference in CPU time between 100 and 1000 group size.
DEFAULT_GROUP_SIZE = 1
# With one thread per file (pigz -p 1) pigz uses way less virtual memory
# (10 vs 300 MB for 4 threads) and total CPU time is also decreased.
# Example: 2.3 gb input file file. Five output files.
# TT=Total cpu time. RT=realtime
# Compression level 1: threads=1 RT=58s TT=3m33, threads=4 RT=57s TT=3m45
# Compression level 5: threads=1 RT=1m48 TT=7m53, threads=4 RT=1m23, TT=8m39
# But this is on a 8 thread machine. When using less threads RT will go towards
# TT. Default compression is 1. So default threads 1 makes the most sense.
DEFAULT_THREADS_PER_FILE = 1
CYTHON_AVAILABLE = False
# This also makes sure we have a valid test file.

# output_files = [Path(__file__).parent / Path("data") / Path("part."+str(_)+".test.fq.gz")
#                 for _ in range(number_of_splits)]
# output_files += [Path(__file__).parent / Path("data") / Path("part."+str(_)+".test2.fq.gz")
#                 for _ in range(number_of_splits)]
# print(output_files)

# start_time = time.time()
# logging.info("--- Start at %s ---" % (start_time))

# #split_fastqs(TEST_FILE, output_files, use_cython=False)
# key_split_fastqs(TEST_FILE,TEST_FILE_2, output_files, use_cython=False)

# end_time = time.time()
# logging.info("--- End at %s  ---" % (end_time))
# logging.info("--- Takes %s seconds ---" % (end_time - start_time))



def main():
    parser = argument_parser()
    parsed_args = parser.parse_args()


    number_of_splits = parsed_args.splits

    output_files = [Path(parsed_args.output) / Path(parsed_args.prefix+"part."+str(_)+".R1."+parsed_args.supfix)
                for _ in range(number_of_splits)]
    output_files += [Path(parsed_args.output) / Path(parsed_args.prefix+"part."+str(_)+".R2."+parsed_args.supfix)
                for _ in range(number_of_splits)]

    key_split_fastqs(parsed_args.input1,parsed_args.input2,
                 output_files,
                 compression_level=parsed_args.compression_level,
                 threads_per_file=parsed_args.threads_per_file,
                 group_size=parsed_args.group_size,
                 bc_offset=parsed_args.bc_offset,
                 bc_length=parsed_args.bc_len,
                 use_cython=False)

def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-r1", "--input1", type=Path, required=True,
                        help="The R1 fastq file to be scattered.")
    parser.add_argument("-r2", "--input2", type=Path, required=True,
                        help="The R2 fastq file to be scattered.")
    parser.add_argument("-o", "--output", action="append", type=Path,
                        required=True,
                        help="Output folder over these output files."
                        )
    parser.add_argument("-s", "--supfix", action="append", type=str, default=".fastq.gz",
                        help="The extensions determine "
                             "which compression algorithm will be used. '.gz' "
                             "for gzip, '.bz2' for bzip2, '.xz' for xz. Other "
                             "extensions will use no compression."
                        )
    parser.add_argument("-p", "--prefix", type=str,
                        required=True, default="output",
                        help="Output file name."
                              "Default=output."
                        )
    parser.add_argument("-n", "--splits", type=int,
                        required=True,
                        help="Number of splits of the output files."
                        )
    parser.add_argument("-c", "--compression-level", type=int,
                        default=DEFAULT_COMPRESSION_LEVEL,
                        help="Only applicable when output files have a '.gz' "
                             "extension. Default={0}"
                        .format(DEFAULT_COMPRESSION_LEVEL))
    parser.add_argument("--bc_offset", type=int, default=0,
                        help="Offset of the barcode in the R1 file"
                             "Default=0")
    parser.add_argument("--bc_len", type=int, default=16,
                        help="Length of the barcode in the R1 file"
                             "extension. Default=16")
    parser.add_argument("-t", "--threads-per-file", type=int,
                        default=DEFAULT_THREADS_PER_FILE,
                        help="Set the number of compression threads per output"
                             " file. NOTE: more threads are only useful when "
                             "using a compression level > 1. Default={0}"
                             "".format(DEFAULT_THREADS_PER_FILE))

    # BELOW ARGUMENTS ARE FOR BENCHMARKING AND TESTING PURPOSES
    parser.add_argument("-g", "--group-size", type=int,
                        default=DEFAULT_GROUP_SIZE,
                        help=argparse.SUPPRESS
                        # help="Specify the group size. This will set how "
                        #      "fine-grained the fastq distribution will be."
                        #      " Default: {0}".format(DEFAULT_GROUP_SIZE)
                        )

    # cython_not_available_text = (
    #     " WARNING: the cython version of the splitting algorithm was not "
    #     "compiled for your system. This is probably due to your system not "
    #     "having a C compiler installed or no pre-built binaries being "
    #     "available for your system. Using fastqsplitter in cython mode will "
    #     "fail.")
    # cython_help = ("Use the cython version of the file splitting algorithm." +
    #                (" (default)" if CYTHON_AVAILABLE else
    #                 cython_not_available_text))
    # python_help = ("Use the python version of the file splitting algorithm." +
    #                ("" if CYTHON_AVAILABLE else " (default)"))
    # cython_group = parser.add_mutually_exclusive_group()
    # cython_group.set_defaults(use_cython=CYTHON_AVAILABLE)
    # cython_group.add_argument("--cython",
    #                           action="store_true",
    #                           dest="use_cython",
    #                           help=cython_help)
    # cython_group.add_argument("--python",
    #                           action="store_false",
    #                           dest="use_cython",
    #                           help=python_help)
    return parser
if __name__ == "__main__":
    main()
