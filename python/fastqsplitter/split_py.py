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


def filesplitter(input_handle: io.BufferedReader,
                 output_handles: List[io.BufferedWriter],
                 lines_per_block=100):
    print(lines_per_block)
    # Make sure inputs are sensible.
    if len(output_handles) < 1:
        raise ValueError("The number of output files should be at least 1.")
    if lines_per_block < 1:
        raise ValueError("The number of lines per block should be at least 1.")

    group_number = 0
    number_of_output_files = len(output_handles)
    # An input handle should be an iterable, but since this package uses xopen
    # we cannot be sure.
    fastq_iterator = iter(input_handle)
    # Alias the readline method for extra speed.
    # next(handle) equals handle.readline() but throws a
    # StopIteration error if the lines is ''.
    readline = fastq_iterator.__next__
    not_at_end_of_file = True

    while not_at_end_of_file:
        block = []
        # This is the fastest way to read blocks of lines in pure python.
        # The method used in split_cy is slower in pure python because it
        # stores the counting integer as a PyObject.
        for j in range(lines_per_block):
            try:
                block.append(readline())
            # Readline should throw a StopIteration at EOF
            except StopIteration:
                not_at_end_of_file = False
                break

        # Select the output handle and write the block to it.
        output_handles[group_number].write(b"".join(block))

        # Set the group number for the next group to be written.
        group_number += 1
        if group_number == number_of_output_files:
            group_number = 0

def keysplitter(input_handle_R1: io.BufferedReader,input_handle_R2: io.BufferedReader,
                 output_handles: List[io.BufferedWriter],
                 lines_per_block=100,
                 key_offset=0,
                 key_length=16):

    # Make sure inputs are sensible.
    if len(output_handles) < 1:
        raise ValueError("The number of output files should be at least 1.")
    if lines_per_block < 1:
        raise ValueError("The number of lines per block should be at least 1.")

    group_number = 0
    n_files = int(len(output_handles)/2)
    number_of_output_files = len(output_handles)
    # An input handle should be an iterable, but since this package uses xopen
    # we cannot be sure.
    fastq1_iterator = iter(input_handle_R1)
    fastq2_iterator = iter(input_handle_R2)

    # Alias the readline method for extra speed.
    # next(handle) equals handle.readline() but throws a
    # StopIteration error if the lines is ''.
    readline1 = fastq1_iterator.__next__
    readline2 = fastq2_iterator.__next__

    not_at_end_of_file = True
    while not_at_end_of_file:
        blocks1 = [ [] for i in range(n_files)]
        blocks2 = [ [] for i in range(n_files)]


        # This is the fastest way to read blocks of lines in pure python.
        # The method used in split_cy is slower in pure python because it
        # stores the counting integer as a PyObject.
        for j in range(lines_per_block):
            temp1 = []
            temp2 = []
            for i in range(0,4):
                try:
                    temp1.append(readline1())
                    temp2.append(readline2())

                # Readline should throw a StopIteration at EOF
                except StopIteration:
                    not_at_end_of_file = False
                    if(temp1 == []):
                        return 
                    break

            group_number = hash(temp1[1][key_offset:key_length])%n_files
            blocks1[group_number].append(temp1)
            blocks2[group_number].append(temp2)

        for f in range(n_files):
            # Select the output handle and write the block to it.
            output_handles[f].write(b"".join(blocks1[f][0]))
            output_handles[f+n_files].write(b"".join(blocks2[f][0]))
            
    # while not_at_end_of_file:
    #     block1 = []
    #     block2 = []
    #     #table = [[0]*2 for i in range(n_files)]


    #     # This is the fastest way to read blocks of lines in pure python.
    #     # The method used in split_cy is slower in pure python because it
    #     # stores the counting integer as a PyObject.
    #     for j in range(lines_per_block):
    #         try:
    #             block1.append(readline1())
    #             block2.append(readline2())

    #         # Readline should throw a StopIteration at EOF
    #         except StopIteration:
    #             not_at_end_of_file = False
    #             if(block1 == []):
    #                 return 
    #             break

    #     # Select the output handle and write the block to it.

    #     group_number = hash(block1[1][key_offset:key_length])%n_files
    #     output_handles[group_number].write(b"".join(block1))
    #     output_handles[group_number+n_files].write(b"".join(block2))