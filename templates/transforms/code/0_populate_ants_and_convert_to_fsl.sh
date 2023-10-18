#!/bin/bash

mkdir -p ../ants ../spm ../fsl
./fsl_to_fmriprep.sh
./fmriprep_to_fsl.sh
