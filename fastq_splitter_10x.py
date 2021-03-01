#! /usr/bin/python3

import io
import gzip
import subprocess
import os
#import boto3
import argparse
import sys
import shutil

global record_separator, pair_separator, parser_result, max_split_size
max_split_size = 1024 * 10 # 10G
record_separator = "\n"
#pair_separator = "\t"

TEMP_OUTPUT_FOLDER = "/home/junyichen/temp/output"


def split_reads_10x(file_list):
    # We assume that 10x has two reads and the size R2>R1
    global parser_result

    #if len(file_list) == 2:
    file_r1, file_r2 = file_list

    # File one name is the last name in the path name and the file
    file_r2_name = file_r2.split("/")[-1]
    file_r1_name = file_r1.split("/")[-1]

    gzipped_file = False
    if file_r2.endswith(".gz"):
        gzipped_file = True
        output_prefix_1, suffix_extension_1 = file_r1_name[:-3].rsplit(".", 1)
        output_prefix_2, suffix_extension_2 = file_r2_name[:-3].rsplit(".", 1)

    else:
        output_prefix_1, suffix_extension_1 = file_r1.rsplit(".", 1)
        output_prefix_2, suffix_extension_2 = file_r2.rsplit(".", 1)

    output_dir_1 = TEMP_OUTPUT_FOLDER + "/processing_" + output_prefix_1
    output_dir_2 = TEMP_OUTPUT_FOLDER + "/processing_" + output_prefix_2


    try:
        os.makedirs(output_dir_1, exist_ok=True)
        os.makedirs(output_dir_2, exist_ok=True)

        #os.mkdir(output_dir)
    except:
        shutil.rmtree(output_dir_1, ignore_errors=True)
        os.makedirs(output_dir_1, exist_ok=True)

        shutil.rmtree(output_dir_2, ignore_errors=True)
        os.makedirs(output_dir_2, exist_ok=True)
        print('Output directory {}/{} exist.'.format(output_dir_1,output_dir_2), file=sys.stderr)

    split_file_names = file_list
    
    file_part = 0
    output_name_1 = "{}/{}_part{}.{}".format(output_dir_1, output_prefix_1, str(file_part), suffix_extension_1)
    output_name_gzip_1 = output_name_1 + ".gz"

    output_name_2 = "{}/{}_part{}.{}".format(output_dir_2, output_prefix_2, str(file_part), suffix_extension_2)
    output_name_gzip_2 = output_name_2 + ".gz"


    if gzipped_file:
        file_r2_reader = io.TextIOWrapper(gzip.open(split_file_names[1]))
        file_r1_reader = io.TextIOWrapper(gzip.open(split_file_names[0]))

    else:
        file_r2_reader = open(split_file_names[1])
        file_r1_reader = open(split_file_names[0])

    output_writer_1 = open(output_name_1, 'w')
    output_writer_2 = open(output_name_2, 'w')


    r2_record = []
    r1_record = []
    index = 0

    for read_line in file_r2_reader:
        # After collecting one FASTQ record, write it to files
        if index % 4 == 0 and index != 0:

            record_output_2 = record_separator.join(r2_record) + "\n"
            record_output_1 = record_separator.join(r1_record) + "\n"

            r2_record = []
            r1_record = []

            output_writer_2.write(record_output_2)
            output_writer_1.write(record_output_1)


            if output_writer_2.tell() > (max_split_size * 1024 * 1024):

                output_writer_2.close()
                #subprocess.call(["gzip", "-f", "--fast", output_name_2])
                subprocess.Popen(["gzip", "-f", "--fast", output_name_2])

                output_writer_1.close()
                #subprocess.call(["gzip", "-f", "--fast", output_name_1])
                subprocess.Popen(["gzip", "-f", "--fast", output_name_2])

                # If there is still more record to process, open a new file to write to
                if file_r2_reader.buffer.peek(10):
                    file_part += 1
                    output_name_2 = "{}/{}_part{}.{}".format(output_dir_2, output_prefix_2, str(file_part), suffix_extension_2)
                    output_name_gzip_2 = output_name_2 + ".gz"
                    output_writer_2 = open(output_name_2, 'w')

                    output_name_1 = "{}/{}_part{}.{}".format(output_dir_1, output_prefix_1, str(file_part), suffix_extension_1)
                    output_name_gzip_1 = output_name_1 + ".gz"
                    output_writer_1 = open(output_name_1, 'w')

        r2_record.append(read_line.strip())
        r1_record.append(file_r1_reader.readline().strip())

        index += 1

    # Close the file if it's not yet closed and write the last line!
    if not output_writer_2.closed:

        record_output_2 = record_separator.join(r2_record) + "\n"
        record_output_1 = record_separator.join(r1_record) + "\n"

        output_writer_2.write(record_output_2)
        output_writer_2.close()
        subprocess.Popen(["gzip", "-f", "--fast", output_name_2])

        output_writer_1.write(record_output_1)
        output_writer_1.close()
        subprocess.Popen(["gzip", "-f", "--fast", output_name_1])

    return "Success"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='File splitter for spark-based RNA-seq Pipeline')
    parser.add_argument('--input', '-i', action="store", dest="input_dir",
                        help="Input location (S3 bucket)", required=True)
    parser.add_argument('--output', '-o', action="store", dest="output_dir", help="output location", required=True)
    parser.add_argument('--region', '-r', action="store", dest="s3_region", help="Region for S3 bucket",
                        nargs="?", default="us-east-1")
    parser.add_argument('--size', '-s', action="store", dest="split_size", help="Split size of each file (GigaByte)", type=int,
                        nargs="?", default=10)   
    parser_result = parser.parse_args()

    TEMP_OUTPUT_FOLDER = parser_result.output_dir

    parser_result.output_dir = parser_result.output_dir.strip().rstrip("/")
    file_list = parser_result.input_dir.strip().split(",")

    max_split_size = parser_result.split_size * 1024


    print(file_list)
    split_reads_10x(file_list)
    # for line in sys.stdin:
    #     file_list = line.strip().split()
    #     split_reads(file_list)
