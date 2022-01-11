#!/bin/bash

grep -v "^>" $1 | wc -m | awk '{print $1*"'"$2"'"}'