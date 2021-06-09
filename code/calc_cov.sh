#!/bin/bash

stats=$1
fai=$2

cat $1 | grep "total length" | cut -f 3 | while read total ; do cat $2 | tail -n 1 | awk '{print "'"$total"'"/$3}' ; done
