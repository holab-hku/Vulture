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

    # group_number = 0
    n_files = int(len(output_handles)/2)
    # number_of_output_files = len(output_handles)
    # An input handle should be an iterable, but since this package uses xopen
    # we cannot be sure.
    # fastq1_iterator = iter(input_handle_R1)
    # fastq2_iterator = iter(input_handle_R2)

    # # Alias the readline method for extra speed.
    # # next(handle) equals handle.readline() but throws a
    # # StopIteration error if the lines is ''.
    # readline1 = fastq1_iterator.__next__
    # readline2 = fastq2_iterator.__next__

    # not_at_end_of_file = True

    # if(imp==2) :
    #     while not_at_end_of_file:
    #         temp1 = []
    #         temp2 = []
    #         # This is the fastest way to read blocks of lines in pure python.
    #         # The method used in split_cy is slower in pure python because it
    #         # stores the counting integer as a PyObject.
    #         # try:
    #         temp1.append(readline1())
    #         temp1.append(readline1())
    #         group_number = hash(temp1[1][key_offset:key_length])%n_files
    #         temp1.append(readline1())
    #         temp1.append(readline1())

    #         temp2.append(readline2())
    #         temp2.append(readline2())
    #         temp2.append(readline2())
    #         temp2.append(readline2())

    #         # Readline should throw a StopIteration at EOF
    #         # except StopIteration:
    #         #     not_at_end_of_file = False
    #         #     if(temp1 == []):
    #         #         return 
    #         #     break

    #         # Select the output handle and write the block to it.
    #         output_handles[group_number].write(b"".join(temp1))
    #         output_handles[group_number+n_files].write(b"".join(temp2))

    # elif(imp==3):
    temp1 = []
    temp2 = []
    i = 0
    for line in input_handle_R1:
        temp1.append(line)
        temp2.append(input_handle_R2.readline())

        i += 1
        if i == lines_per_block:
            group_no = hash(temp1[1][key_offset:key_length])%n_files
            output_handles[group_no].write(b"".join(temp1))
            output_handles[group_no+n_files].write(b"".join(temp2))

            temp1 = []  # reset block.
            temp2 = []
            i = 0
            # This works, if blocksize == 100. Then i will be [0, 1, 2, .., 98, 99]
            # which is exactly 100 numbers.

        # Write remainder to file. Since stuff is only written at i == blocksize
        # there is a remainder that needs to be written.
        # group_no = hash(temp1[1][key_offset:key_length])%n_files
        # output_handles[0].write(b"".join(temp1))
        # output_handles[0+n_files].write(b"".join(temp2))

    # else:
    #     while not_at_end_of_file:
    #         temp1 = []
    #         temp2 = []


    #         # This is the fastest way to read blocks of lines in pure python.
    #         # The method used in split_cy is slower in pure python because it
    #         # stores the counting integer as a PyObject.
    #         for j in range(lines_per_block):
    #             try:
    #                 temp1.append(readline1())
    #                 temp2.append(readline2())

    #             # Readline should throw a StopIteration at EOF
    #             except StopIteration:
    #                 not_at_end_of_file = False
    #                 if(temp1 == []):
    #                     return 
    #                 break

    #         # Select the output handle and write the block to it.
    #         group_number = hash(temp1[1][key_offset:key_length])%n_files
    #         output_handles[group_number].write(b"".join(temp1))
    #         output_handles[group_number+n_files].write(b"".join(temp2))