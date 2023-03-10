#!/usr/bin/env python3
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
import contextlib
import io
from pathlib import Path
from shutil import Error
from typing import List
import logging
import os

# xopen opens files as normal files, gzip files, bzip2 files or xz files
# depending on extension.
import xopen

# Import from cython if available
try:
    from key_split_cy import keyfilesplitter as cy_splitter
    #from .split_cy import filesplitter as cy_splitter

    logging.info("Cython version is avaliable, use cython code") 

    CYTHON_AVAILABLE = True
except ImportError as e:

    logging.warning("Cython version is not avaliable, compile cython code") 
    
    try:
        path = os.getcwd()
        code_path = os.path.dirname(os.path.realpath(__file__))
        os.chdir(code_path)
        os.system('python setup.py build_ext --inplace')
        os.chdir(path)

        from key_split_cy import keyfilesplitter as cy_splitter
        #from .split_cy import filesplitter as cy_splitter
        CYTHON_AVAILABLE = True

    except Error as e:      
        CYTHON_IMPORT_ERROR = e
        CYTHON_AVAILABLE = False
        logging.warning("Cython version is not avaliable, use python code") 

from key_split_py import keysplitter as py_splitter

# Choose 1 as default compression level. Speed is more important than filesize
# in this application.
DEFAULT_COMPRESSION_LEVEL = 1
# There is no noticable difference in CPU time between 100 and 1000 group size.
DEFAULT_GROUP_SIZE = 1000
# With one thread per file (pigz -p 1) pigz uses way less virtual memory
# (10 vs 300 MB for 4 threads) and total CPU time is also decreased.
# Example: 2.3 gb input file file. Five output files.
# TT=Total cpu time. RT=realtime
# Compression level 1: threads=1 RT=58s TT=3m33, threads=4 RT=57s TT=3m45
# Compression level 5: threads=1 RT=1m48 TT=7m53, threads=4 RT=1m23, TT=8m39
# But this is on a 8 thread machine. When using less threads RT will go towards
# TT. Default compression is 1. So default threads 1 makes the most sense.
DEFAULT_THREADS_PER_FILE = 1


def key_split_fastqs(input_file_1: Path,input_file_2: Path, output_files: List[Path],
                 compression_level: int = DEFAULT_COMPRESSION_LEVEL,
                 group_size: int = DEFAULT_GROUP_SIZE,
                 threads_per_file: int = DEFAULT_THREADS_PER_FILE,
                 use_cython: bool = CYTHON_AVAILABLE,
                 bc_offset=0,
                 bc_length=10,
                 imp=1):
    # contextlib.Exitstack allows us to open multiple files at once which
    # are automatically closed on error.
    # https://stackoverflow.com/questions/19412376/open-a-list-of-files-using-with-as-context-manager
    with contextlib.ExitStack() as stack:
        input_fastq_1 = stack.enter_context(
            xopen.xopen(input_file_1, mode='rb'))  # type: io.BufferedReader
        input_fastq_2 = stack.enter_context(
            xopen.xopen(input_file_2, mode='rb'))  # type: io.BufferedReader

        output_handles = [
            stack.enter_context(xopen.xopen(
                filename=output_file,
                mode='wb',
                compresslevel=compression_level,
                threads=threads_per_file
            )) for output_file in output_files
        ]  # type: List[io.BufferedWriter]
        if use_cython:
            if CYTHON_AVAILABLE:

                print("Split files using the cython splitter") 
                cy_splitter(input_handle_R1=input_fastq_1,input_handle_R2=input_fastq_2,
                        output_handles=output_handles,
                        lines_per_block=group_size * 4,
                        # 4 lines per fastq record
                        k_offset=bc_offset,
                        k_length=bc_length
                        )
            else:
                raise CYTHON_IMPORT_ERROR
        else:  # Use python fallback
            print("Split files using the python splitter") 
            py_splitter(input_handle_R1=input_fastq_1,input_handle_R2=input_fastq_2,
                        output_handles=output_handles,
                        lines_per_block=group_size * 4,
                        # 4 lines per fastq record
                        key_offset=bc_offset,
                        key_length=bc_length,
                        imp=imp
                        )

def main():
    parser = argument_parser()
    parsed_args = parser.parse_args()


    number_of_splits = parsed_args.splits

    output_files = [Path(parsed_args.output) / Path(parsed_args.prefix+"part."+str(_)+".R1"+parsed_args.supfix)\
                for _ in range(number_of_splits)]
    output_files += [Path(parsed_args.output) / Path(parsed_args.prefix+"part."+str(_)+".R2"+parsed_args.supfix)\
                for _ in range(number_of_splits)]

    key_split_fastqs(parsed_args.input1,parsed_args.input2,
                 output_files,
                 compression_level=parsed_args.compression_level,
                 threads_per_file=parsed_args.threads_per_file,
                 group_size=parsed_args.group_size,
                 bc_offset=parsed_args.bc_offset,
                 bc_length=parsed_args.bc_len,
                 use_cython=parsed_args.use_cython)

def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-r1", "--input1", type=Path, required=True,
                        help="The R1 fastq file to be scattered.")
    parser.add_argument("-r2", "--input2", type=Path, required=True,
                        help="The R2 fastq file to be scattered.")
    parser.add_argument("-o", "--output", type=str,
                        required=True,
                        help="Output folder over these output files."
                        )
    parser.add_argument("-s", "--supfix", type=str, default=".fastq.gz",
                        help="The extensions determine "
                             "which compression algorithm will be used. '.gz' "
                             "for gzip, '.bz2' for bzip2, '.xz' for xz. Other "
                             "extensions will use no compression."
                        )
    parser.add_argument("-p", "--prefix", type=str, default="output",
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
    cython_not_available_text = (
        " WARNING: the cython version of the splitting algorithm was not "
        "compiled for your system. This is probably due to your system not "
        "having a C compiler installed or no pre-built binaries being "
        "available for your system. Using fastqsplitter in cython mode will "
        "fail.")
    cython_help = ("Use the cython version of the file splitting algorithm." +
                   (" (default)" if CYTHON_AVAILABLE else
                    cython_not_available_text))
    python_help = ("Use the python version of the file splitting algorithm." +
                   ("" if CYTHON_AVAILABLE else " (default)"))
    cython_group = parser.add_mutually_exclusive_group()
    cython_group.set_defaults(use_cython=CYTHON_AVAILABLE)
    cython_group.add_argument("--cython",
                              action="store_true",
                              dest="use_cython",
                              help=cython_help)
    cython_group.add_argument("--python",
                              action="store_false",
                              dest="use_cython",
                              help=python_help)
    return parser
if __name__ == "__main__":
    main()
