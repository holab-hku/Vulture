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

# cython: language_level=3


import io
from typing import List


def keyfilesplitter(input_handle_R1: io.BufferedReader,input_handle_R2: io.BufferedReader,
                 output_handles: List[io.BufferedWriter],
                 lines_per_block=100,
                 k_offset=0,
                 k_length=16):
    # Make sure inputs are sensible.

    if len(output_handles) < 1:
        raise ValueError("The number of output files should be at least 1.")
    if lines_per_block < 1:
        raise ValueError("The number of lines per block should be at least 1.")

    # cdef bytes line  # Faster to not type.
    temp1 = []
    temp2 = []

    cdef unsigned int blocksize = lines_per_block
    cdef unsigned int i = 0
    cdef unsigned int j = 0
    cdef unsigned int f = 0
    cdef unsigned int group_no = 0
    cdef unsigned int n_files = int(len(output_handles)/2)
    cdef unsigned int n_outputs = len(output_handles)

    cdef unsigned int key_offset = k_offset
    cdef unsigned int key_length = k_length

    blocks = [[] for i in range(n_outputs)]

    # for line in handle is the fastest way to read lines in python that I
    # know of. Implementations with next(handle) or handle.readline are
    # slower.
    for line in input_handle_R1:
        temp1.append(line)
        temp2.append(input_handle_R2.readline())

        i += 1
        j += 1

        if j == 4:
            group_no = hash(temp1[1][key_offset:key_length])%n_files
            blocks[group_no]+=(temp1)
            blocks[group_no+n_files]+=(temp2)
            temp1 = []  # reset block.
            temp2 = []
            j = 0

        if i == blocksize:
            while f < n_outputs:
                output_handles[f].write(b"".join(blocks[f]))
                f += 1
            blocks = [[] for i in range(n_outputs)] # reset block.
            f = 0
            i = 0
        # This works, if blocksize == 100. Then i will be [0, 1, 2, .., 98, 99]
        # which is exactly 100 numbers.

    while f < n_outputs:
        output_handles[f].write(b"".join(blocks[f]))
        f += 1
    f = 0
    # Write remainder to file. Since stuff is only written at i == blocksize
    # there is a remainder that needs to be written.

    #output_handles[0].write(b"".join(temp1))
    #output_handles[0+n_files].write(b"".join(temp2))

    # The resetting of i at blocksize accomplishes two things:
    # 1. It will never be larger than blocksize. We do not have to keep
    # counting to infinity, which will lead to an integer overflow.
    # 2. We do not have to use modulo (i % blocksize) == 0 to determine
    # whether the blocksize is reached. Direct comparison is faster than
    # modulo.
