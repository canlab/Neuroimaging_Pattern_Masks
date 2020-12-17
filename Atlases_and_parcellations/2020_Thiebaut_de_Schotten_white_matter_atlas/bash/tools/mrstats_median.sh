#!/bin/bash

[ $# -lt 3 ] && { echo 'Usage :
  $1 = folder path
  $2 = spreadsheet name(it will not overwrite)
  $3 = input suffix (ad/ar/fa etc ...)
  [Optional] $4 = roi1 suffix
  [Optional] $5 = roi2 suffix';
  exit 1; }

folder=$1

name=$(basename $1)

suffix=$3

roi1="6roi"
roi2="7roi"
if [[ $4 != '' ]];
then
  roi1=$4
fi;
if [[ $5 != '' ]];
then
  roi2=$5
fi;

val1=`mrstats $1/${name}_${suffix}.mif -mask $1/${name}_${roi1}.mif -output median`
echo $name,$val1 >> $2
val2=`mrstats $1/${name}_${suffix}.mif -mask $1/${name}_${roi2}.mif -output median`
echo $name,$val2 >> $2
