#!/usr/bin/env bash

# simulate fastq_pass during ONT run
# $1 is source barcode folder with fastq files
# $2 is cp destination
# $3 is number of files to cp
# $4 is sleep time

while [[ true ]] 
do 
    files=$(find $1 -type f -name "*.fastq" | shuf -n $3)
    rsync -R $files $2
    printf "copied to $2:\n$files\nsleeping for $4 sec\n"
    sleep $4
done