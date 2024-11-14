PainGen Participants:

Participants were recruited by telephone from the Colorado Community Twin Sample, which is
derived from the Colorado Twin Registry, a population-based registry which has been run by the
Institute of Behavioral Genetics (IBG) at the University of Colorado since 1984. The study
includes a larger sample of participants, but the data here are from a subset of 248 unrelated 
participants. Subjects from the same families, or subjects with image quality metrics greater
than 3SD away form the mean were excluded. MRIQC 0.16.1 was was used to compute image quality
metrics. Participants were excluded from participation in the study if they did not pass MRI 
screening (e.g., had metals in their body) or had a history of liver disease/damage, allergies 
to local analgesics, or were breastfeeding. Of those 248 that remained, 241 had T1w images
and cortical reconstructions available for this analysis.

In the full dataset (including siblings and participants with bad image quality metrics) 
Participants’ age ranged from approximately 30 to 43 years (M = 35, SD = 2.6), and the sample 
included approximately 100 men and 140 women. All participants provided their informed consent 
at the beginning of the experiment. The experiment was approved by the institutional review board 
of the University of Colorado Boulder. 


BMRK5 participants:

Participants. Participants are part of a larger sample of 97 individuals 
(47 male, age 19–54 years; M = 28.98, s.d. = 5.56) 33 AA individuals (15 male), 
32 non-Hispanic WA individuals (16 male) and 32 HA individuals (16 males), based 
on self-reported ethnicity. 88 had data available for this analysis.

Participants were recruited from the greater Denver area through Craigslist or
one of three different participant pools from the University of Colorado Boulder
Institute for Behavioral Genetics (IBG) to capitalize on existing genetic data on
these participants in future analyses. Participants reported no current or recent
(past 6 months) neurological or psychiatric diagnosis and reported no current use
of psychoactive or pain medications. Participants also reported no pain-related
medical conditions.

The study was approved by the University of Colorado Boulder Institutional Review Board 
and we complied with all relevant ethical regulations when carrying out the study. Written 
informed consent was obtained from all participants, who were financially compensated 
for their participation.

SpaceTop Participants:

Participants were recruited from the Dartmouth College and surrounding community. 

All participants provided their informed consent at the beginning of the study. The study was
approved by the institutional review board of Dartmouth College. 113 subjects were scanned, 
of which 112 had T1w images available with cortical reconstructions for use in this analysis.


PainGen and SpaceTop Imaging data acquisition:

Imaging data were acquired using a 3T Siemens Prisma MRI scanner with a 32-channels head
coil, at the University of Colorado at Boulder Center for Innovation and Creativity (PainGen) 
or Dartmouth College (SpaceTop). A T1-weighted structural scan was acquired using a magnetization 
prepared rapid gradient echo (MPRAGE) pulse sequence with parallel imaging factor (iPAT) of 3, 
TR=2000 ms, TE=2.11ms, flip angle=8 degrees, FOV=256mm, resolution=0.8 x 0.8 x 0.8mm. Additional 
functional and structural scans were collected but are not relevant to this dataset.


BMRK5 Data acquisition:

Data were collected on a 3 Tesla Siemens Trio MRI scanner at the University of 
Colorado Boulder Center for Innovation and Creativity. A high-resolution T1-weighted 
magnetization-prepared rapid gradient echo (MPRAGE) structural scan (1x1x1 mm3 voxels,
repetition time (TR): 2,530 ms, echo time (TE): 1.64 ms, flip angle: 7*, inversion
time (Tl): 1,200 ms, field of view (FoV) read: 256 mm, echo spacing: 12.2 ms,
bandwidth: 651 Hz Px-1, time: 6:03) was performed on each participant.


Anatomical data preprocessing (identical for PainGen and BMRK5): 

BMRK5 and PainGen data were preprocssed by running fmriprep 20.2.3 while SpaceTop was preprocessed by
fmriprep 21.0.3, in all cases with surface reconstruction enabled. What follows is the fmriprep 
boilerplate lightly edited to indicate differences between fmriprep versions:

The T1-weighted (T1w) image was corrected for intensity nonuniformity (INU) with N4BiasFieldCorrection 
(Tustison et al. 2010), distributed with ANTs 2.3.3 (RRID:SCR_004757), and used as T1w-reference throughout 
the workflow. The T1w-reference was then skullstripped with a Nipype implementation of the 
antsBrainExtraction.sh workflow (from ANTs), using OASIS30ANTs as target template. Brain tissue segmentation 
of cerebrospinal fluid (CSF), whitematter (WM) and gray-matter (GM) was performed on the brain-extracted 
T1w using fast (BMRK5 & PainGen: FSL 5.0.9, RRID:SCR_002823; SpaceTop: FSL 6.0.5.1:57b01774, RRID:SCR_002823, 
Zhang, Brady, and Smith 2001). Brain surfaces were reconstructed using recon-all (FreeSurfer 6.0.1, 
RRID:SCR_001847, Dale, Fischl, and Sereno 1999), and the brain mask estimated previously was refined with a 
custom variation of the method to reconcile ANTs-derived and FreeSurfer-derived segmentations of the cortical 
gray-matter of Mindboggle (RRID:SCR_002438, Klein et al. 2017). Volume-based spatial normalization to two 
standard spaces (MNI152NLin2009cAsym, MNI152NLin6Asym) was performed through nonlinear registration with 
antsRegistration (ANTs 2.3.3), using brain extracted versions of both T1w reference and the T1w template. 
The following templates were selected for spatial normalization: ICBM 152 Nonlinear Asymmetrical template 
version 2009c (Fonov et al. [2009], RRID:SCR_008796; TemplateFlow ID: MNI152NLin2009cAsym], FSL's MNI ICBM 152 
nonlinear 6th Generation Asymmetric Average Brain Stereotaxic Registration Model (Evans et al. [2012], 
RRID:SCR_002823; TemplateFlow ID: MNI152NLin6Asym).


Data analysis:

Registration fusion was performed using the approach detailed by Wu, Ngo, Greve, et al. (2018), except 
that all transformations used were computed by fmriprep as detailed above. 
Registration fusion was implemented on an HPC system with one job instance per participant. An example
script is included with this data release. See single_subject_registration_fusion.sh for comprehensive
methods. In brief, a segmentation was estimated by freesurfer as part of the standard surface reconstruction
pipeline. This was projected into subject specific volumetric space for each hemisphere separately using
mri_surf2vol and then aligned to one of two MNI template spaces (MNI152NLin6Asym and MNI152NLin2009cAsym)
using ANTs and transformation matrices computed as part of the standard fmriprep anatomical pipeline.
The results for the left and right hemisphere concatenated across subjects for each space, study and hemisphere.


References:
BMRK5 has been published here,
Losin, Woo, Medina, Andrews, Eisenbarth, Wager. (2020). Nature Human Behavior 4(5) 517-530.

PainGen has had data partially published here,
Botvinik-Nezer, Petre, Ceko, Lindquist, Friedman, Wager. (2023). BioRxiv. DOI:https://doi.org/10.1101/2023.09.21.558825

Registration fusion is described here:
Wu, Ngo, Greve. (2018). Human Brain Mapping 39(9) 37930-3808. DOI: 10.1002/hbm.24213
