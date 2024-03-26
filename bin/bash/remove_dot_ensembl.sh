#!/bin/bash

cat $1 | awk 'BEGIN { OFS=FS="\t" } {  sub(/\..*$/, "", $1); print  }' >> $2
