#!/bin/bash
echo input path $1;
echo output path $2;
kb ref \
-i=$2 \
-f1=$1 \
-g=$2 \