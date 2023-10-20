### MNIColin27v1988

Transformation matrices were computed by running the Colin27 template through fmriprep 20.0.3 with surface reconstruction
enabled with MNI152NLin6Asym and MNI152NLin2009cAsym selected as output spaces. There are two versions of Colin27, one from 
1998 and another that's higher resolution and includes T2 and PD images from 2008. They are not coregistered. This was run 
on the 1998 version. In theory it should be possible to obtain better transforms to MNI152NLin2009cAsym space in particular
by performing a rigid body transfrom from the 1998 to 2008 data and then running multimodal alignment of T1, T2 and PD data
to the MNI152NLin2009cAsym equivalents, but this seemed like more work than it was worth.


### MNI152NLin6Asym to MNI152NLin2009cAsym

Transformation matrices were computed by aligning the templateFlow MNI152NLin6Asym 1mm T1 image to MNI152NLin2009cAsym
1mm T1 image using fmriprep 20.0.3 with surface reconstruction enabled. fsl_to_fmriprep.sh and fmriprep_to_fsl.sh then 
used the resultant h5 files to populate the ants and fsl sister folders here and fsl_to_fmriprep.py and fmriprep_to_fsl.py 
were then run to convert ants to SPM format.

use antsApplyTransform to use the files in the ants directory. The order of application is important. antsApplyTransform
can apply multiple transforms at once, but the arguments need to be specified in the first in last out order, which is
counterintuitive (the last transform you specify is applied first). If in doubt use antsApplyTransform in multiple
separate invocations applying one transform each time. The prefixed numbers indicate the order in which a transform should
be applied. But only use this to double check a multi-transform invocation since each application resamples your data.
Better to resample once than multiple times like with the individual invocations for each separate transform. Refer to the 
antsApplyTransform help for more details

use applywarp with the files in the fsl directory. Ideally you apply both transforms simultaneously, like below, but make
sure you specify premat or postmat based on whether or not your prefixed index for the *.mat file is before that for the 
warp field. If you want to apply these transforms sequentially (bad idea, double resampling is less accurate than single
resampling) you need to think through the original transform. Going from fsl to fmriprep space the affine transform also
results in the change of space, so that the warp occurs entirely in fmriprep space. When inverting this that means that the
reference template for the fmriprep to fsl warp is the fmriprep template while the affine transform that follows is when
the conversion to fsl space happens. Save yourself the headache and apply them in one shot.
applywarp -i MNI152NLin6Asym_T1_1mm.nii \
          -r MNI152NLin2009cAsym_T1_1mm.nii \
          -o MNI152NLin6Asym_to_MNI152NLin2009cAsym_fsl.nii \
          --premat=../../_resources/transforms/fsl/00_fsl_to_fmriprep_AffineTransform_mod.mat \
          -w ../../_resources/transforms/fsl/01_fsl_to_fmriprep_DisplacementFieldTransform_mod.nii.gz

To apply transformations in matlab using SPM utilities try this,
apply_transforms('MNI152NLin2009cAsym_T1_1mm.nii', ...                                              % moving image
                 'MNI152NLin6Asym_T1_1mm.nii',...                                                   % fixed image
                 [],...                                                                             % premat
                 'spm/01_fmriprep_to_fsl_AffineTransform.csv',  ...                                 % postmat
                 'spm/y_00_fmriprep_to_fsl_DisplacementFieldTransform.nii',...                      % warp file
                 'MNI152NLin2009cAsym_to_MNI152NLin6Asym_spm.nii',...                               % output
                 1)                                                                                 % interpolation
See help for details
