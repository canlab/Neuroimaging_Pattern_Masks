#!/bin/bash

[ $# -lt 1 ] && { echo 'Usage :
  $1 = folder path
  [Optional] $2 = roi1 suffix
  [Optional] $3 = roi2 suffix';
  exit 1; }

./mrstats_median.sh $1 median_ad.csv ad $2 $3
./mrstats_median.sh $1 median_rd.csv rd $2 $3
./mrstats_median.sh $1 median_fa.csv fa $2 $3
./mrstats_median.sh $1 median_adc.csv adc $2 $3
