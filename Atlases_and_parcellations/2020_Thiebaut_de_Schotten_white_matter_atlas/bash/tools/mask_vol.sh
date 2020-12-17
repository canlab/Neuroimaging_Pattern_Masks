# /bin/bash

[ $# -lt 2 ] && { echo "Usage : $0 LesionFolder result_file(csv)"; exit 1; }

echo -n "" > $2

for f in $1/*nii*;
do
  vol=`fslstats $f -V | awk '{print $1}'`;
  echo $(basename $f),$vol >> $2
done;
