#! /bin/bash
[ $# -lt 3 ] && { echo 'Usage :
  $1 = MNI152 image
  $2 = mask
  $3 = resultFolder
  $4 = lower threshold for mask
  $5 = upper threshold for mask'; exit 1; }


fileName() {
  name=$(basename $1)
  name=${name%%.*}
  echo -n $name
}

low_thr=$4

if [[ $4 == '' ]];
then
  low_thr=0
fi;

low_thr=$4

if [[ $5 == '' ]];
then
  up_thr=1
fi;

mkdir -p $3

fsleyes render \
  -no \
  -sz 600 750 \
  -of=$3/`fileName $2`.png \
  -hc \
  -xh \
  -yh \
  -p 3 \
  -hl \
  $1 \
  $2 \
  -cm=render3 \
  -dr $low_thr $up_thr
