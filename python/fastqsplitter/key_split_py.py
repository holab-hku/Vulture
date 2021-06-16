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

import io
from typing import List

def keysplitter(input_handle_R1: io.BufferedReader,input_handle_R2: io.BufferedReader,
                 output_handles: List[io.BufferedWriter],
                 lines_per_block=100,
                 key_offset=0,
                 key_length=16,
                 imp=1):

    # Make sure inputs are sensible.
    if len(output_handles) < 1:
        raise ValueError("The number of output files should be at least 1.")
    if lines_per_block < 1:
        raise ValueError("The number of lines per block should be at least 1.")

    n_files = int(len(output_handles)/2)
    number_of_output_files = len(output_handles)
    temp1 = []
    temp2 = []
    blocks = [list() for i in range(number_of_output_files)]
    i = 0
    j = 0
    f = 0
    for line in input_handle_R1:
        temp1.append(line)
        temp2.append(input_handle_R2.readline())

        i += 1
        j += 1
        # Four line per read in fastq
        if j == 4:
            group_no = hash(temp1[1][key_offset:key_length])%n_files
            blocks[group_no]+=(temp1)
            blocks[group_no+n_files]+=(temp2)
            temp1 = []  # reset block.
            temp2 = []
            j = 0
        if i == lines_per_block:
            while f < number_of_output_files:
                output_handles[f].write(b"".join(blocks[f]))
                f += 1
            blocks = [list() for i in range(number_of_output_files)]
            f = 0
            i = 0
    while f < number_of_output_files:
        output_handles[f].write(b"".join(blocks[f]))
        f += 1