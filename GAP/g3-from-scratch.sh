#!/bin/bash

if [ $# -lt 2 ]; then
  echo "Usage:  g3-from-scratch.sh  <genome> <tag>"

  exit;
fi


genome=$1
tag=$2

awkpath=/bin
glimmerpath=/usr/global/blp/bin

# add/change glimmer options here
glimmeropts="-o 50 -g 100 -t 30 -l"

$glimmerpath/long-orfs -n -t 1.15 $genome $tag.longorfs
$glimmerpath/extract -t $genome $tag.longorfs > $tag.train
$glimmerpath/build-icm -r $tag.icm < $tag.train
$glimmerpath/glimmer3 $glimmeropts $genome $tag.icm $tag


