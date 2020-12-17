#! /bin/bash
[ $# -lt 5 ] && { echo 'Usage :
  $1 = MNI152 image
  $2 = patient prefix
  $3 = lesions folder
  $4 = disconnectomes folder
  $5 = resultFolder
  $6 = lower threshold for disconnectomes'; exit 1; }


echo $@

fileName() {
  name=$(basename $1)
  name=${name%%.*}
  echo -n $name
}

disco=`ls $4/*$2*`
les=`ls $3/*$2*`
thr=$6

if [[ $6 == '' ]];
then
  thr=0
fi;

mkdir -p $5
echo "fsleyes render \
  -no \
  -sz 1800 800 \
  -of=$5/`fileName $2`.png \
  -hc \
  --scene=lightbox \
  -nr=3 \
  -nc=8 \
  -ss=5.6 \
  -zr 16 149 \
  $1 \
  $disco \
  -cm=hot \
  -dr $thr 1 \
  $les \
  -ot=mask \
  -mc 0 0 1"
fsleyes render \
  -no \
  -sz 1800 800 \
  -of=$5/`fileName $2`.png \
  -hc \
  --scene=lightbox \
  -nr=3 \
  -nc=8 \
  -ss=5.6 \
  -zr 16 149 \
  $1 \
  $disco \
  -cm=hot \
  -dr $thr 1 \
  $les \
  -ot=mask \
  -mc 0 0 1
  # -mc 0.25 1 0.25

# -cm=hot
# -cm=blue
