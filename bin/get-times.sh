#!/usr/bin/env bash
# return min and max time stamps from fastq headers

# $1 is min or max
# input is stdin (fastq string, eg, head or tail)
# so it is fast also on big fastq files!!

if [ "$1" == min ]; then
sed -n 's/.*start_time=//p' | cut -f 1 -d " " | sort | head -n 1 
elif [ "$1" == max ]; then
sed -n 's/.*start_time=//p' | cut -f 1 -d " " | sort | tail -n 1
else 
echo "Usage: $0 {min,max} -"
exit
fi
