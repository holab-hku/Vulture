#!/bin/bash
echo Arguments input to the container are: "$@"
echo "$@" >> ~/args.txt;
echo Arguments are stored to "~/args.txt";
/usr/bin/tini --;