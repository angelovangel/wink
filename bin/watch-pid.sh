#! /usr/bin/env sh

while [[ true ]]
do
    ps -p $1 -o pid,ppid,start,etime,%mem,%cpu,state
    sleep 2
done