## Overview

Although all data provided by neuroimaging pattern masks should be in MNI 
space, different reference templates may have been used for spatial 
normalization for generating specific atlases and maps. This is probably
not consistently documented, but some comment templates are provided here
for those cases in which the standard space template is described.

The Montreal Neurological Institute originally created their standard space
coordinates for use with MNI305, which was the result of linearly coregistering
305 individuals to 241 brains that had been coregistered to the Taliarch atlas.
Subsequent refinements have been made, resulting in a variety of standard
spaces that all use MNI coordinates but differ in their particular geometry due
to differences in registration approaches and subject samples used. Differences 
between templates are generally small, and likely don't matter if you're doing
traditional fMRI analysis with low resolution acquistions and spatial smoothing
but might be important for achieving voxel level precision in various more
exacting applications.

For a helpful discussion on templates you can refer here,
https://www.lead-dbs.org/about-the-mni-spaces/

Most official MNI atlases are found here: https://nist.mni.mcgill.ca/atlases/

This index is also a useful record of the different templates and what 
software uses what templates:
https://bids-specification.readthedocs.io/en/stable/appendices/coordinate-systems.html

Add any you find useful here, with comments in this document regarding when to use them.
The reference template, spatial resolution and imaging modality should be clear from the
filename, similar to those already included. If you need it spelled out,

\<ref_space\>\_\<modality\>\_\<resolution\>\.\<extension\>

## Templates

### HCP
HCP uses the assymetric version of the MNI152 nonlinear 6th generation template, aka MNI152NLin6ASym. This is different from the official atlas. Rather than the MNI, HCP templates can be obtained from 
https://github.com/Washington-University/HCPpipelines/tree/master/global/templates

Source files are called, <br />
MNI152NLin6ASym_T1_1mm.nii.gz <br />
MNI152NLin6ASym_T1_2mm.nii.gz

### FSL

As of the time of this writing FSL uses the same template as HCP above. You can
confirm by verifying that \<FSLROOT\>/data/standard/MNI152_T1_2mm.nii.gz is 
identical to MNI152NLin6ASym_T1_2mm.nii.gz if in doubt with future releases of
FSL.

### SPM 12
SPM uses the MNI152 linear template with some slight modifications. It is 
distributed with spm12 and is available from canonical/avg152T1.nii, but is 
provided again here for convenience for non SPM users that are working with 
SPM derivatives. 

Source file is called, <br />
IXI549Space_T1_2mm.nii.gz

### SPM99 to SPM8
These use MNI152Lin.

Source file is called, <br />
MNI152Lin_T1_1mm.nii.gz

It was obtained here:
https://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152Lin

### SPM96 and BioImage Suite
These packages use MNIColin27. This is a template based on one subject (Colin) as opposed to the more typical N=152 sample.
It benefits from more detail at the expense of being less representative. In particular it features a somewhat squat 
appearance with more laterally displaced ventral structures. Colin was scanned 27 times, and those 27 images were 
coregistered and averaged to produce this template, hence Colin27. This template seems to be used most often when there's
a need for exquisite detail. An example of such a scenario would be coregistration of histological sections, which have
micron resolution.

Source file is called, <br />
MNIColin27_T1_1mm.nii.gz

It comes from https://nist.mni.mcgill.ca/colin-27-average-brain/

There are two versions of this, one from 1998 and another from 2008. 
This is the 1998 version. The 2008 version has slightly different geometry.

### fMRIPrep (and QSIPrep)
fMRIPrep and QSIPrep use MNI152 Nonlinear 2009c Asymmetric, aka MNI152NLin2009CAsym. This is newer than what FSL or SPM use by default as of the time of this writing. Use of this template by fMRIPrep is documented here, 

https://fmriprep.org/en/stable/spaces.html#standard-spaces

And to ensure consistency with the fMRIPrep repo I've pulled a copy from the TemplateFlow repository linked to in the fmriprep docs above, via this page,

https://www.templateflow.org/browse/

QSIprep doesn't document it's defaults in as explicit a form, but there are repeated references throughout to this being its default template, and most importantly for this repo, this is the default template that atlases need to be registered to for use in qsiprep.

The source file is called, <br />
MNI152NLin2009CAsym_T1_1mm.nii.gz

##
Bogdan Petre 10/10/23
