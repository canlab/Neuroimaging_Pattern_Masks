#!/bin/bash

matlab -nojvm -nodesktop -r "m=importdata('../ants/00_fsl_to_fmriprep_subctx_AffineTransform.mat'); mm=ea_antsmat2mat(m.AffineTransform_double_3_3,m.fixed); csvwrite('../spm/00_fsl_to_fmriprep_subctx_AffineTransform.csv',mm);\
m=importdata('../ants/01_fmriprep_to_fsl_subctx_AffineTransform.mat'); mm=ea_antsmat2mat(m.AffineTransform_double_3_3,m.fixed); csvwrite('../spm/01_fmriprep_to_fsl_subctx_AffineTransform.csv',mm); exit()"

