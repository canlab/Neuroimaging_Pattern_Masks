The ANTs (*.h5) transforms were obtained from leads-dbs via figshare,
https://figshare.com/articles/dataset/MNI_T1_6thGen_NLIN_to_MNI_2009b_NLIN_ANTs_transform/3502238

For more details see here,
https://www.lead-dbs.org/about-the-mni-spaces/

These were converted into FSL format using the following sequence of commands

1) in the fmriprep 20.0.3 singularity container:
for file in MNI_6thgen_2_MNI2009b.h5; do CompositeTransformUtil --disassemble $file ${file%%.*}; done

2) Using wb_command v. 1.5.0 (from HCP's connectome workbench)
wb_command -convert-warpfield -from-itk 01_MNI_6thgen_2_MNI2009b_DisplacementFieldTransform.nii.gz \
    -to-fnirt 01_MNI_6thgen_2_MNI2009b_DisplacementFieldTransform_fsl.nii.gz \
              ../../templates/MNI152NLin6Asym_T1_1mm.nii.gz
wb_command -convert-warpfield -from-itk 00_MNI_6thgen_2_MNI2009b_Inverse_DisplacementFieldTransform.nii.gz \
    -to-fnirt 00_MNI_6thgen_2_MNI2009b_Inverse_DisplacementFieldTransform_fsl.nii.gz \
              ../../templates/MNI152NLin2009cAsym_T1_1mm.nii.gz
