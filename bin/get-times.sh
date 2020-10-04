#!/usr/bin/env bash
# return min and max time stamps from fastq headers

# $1 is min or max
# $2 is a ONT fastq file (will not work on fastq.gz)
# 

if [ "$1" == min ]; then
sed -n 's/.*start_time=//p' $2 | cut -f 1 -d " " | sort | head -n 1 
elif [ "$1" == max ]; then
sed -n 's/.*start_time=//p' $2 | cut -f 1 -d " " | sort | tail -n 1
else 
echo "Usage: $0 min fastq_file"
exit
fi
