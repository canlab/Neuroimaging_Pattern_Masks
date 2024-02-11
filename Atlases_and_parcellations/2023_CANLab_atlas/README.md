## Overview

This is a full brain atlas mashup. It draws from the following,
* Cortex: HCP Multimodal Parcellation (Glasser et al. + Petre 2023 volumetric projection)
* Thalamus: cytoarchitecture (Morel, Licensing restriction)
* Basal Ganglia (except GP): Functional gradients (Tian)
* Basal Ganglia, Globus Palidus: T1/T2 contrast (CIT168 subcortical parcellation)
* Hippocampal formation: cytoarchitecture (Julitch)
* Amygdala: T1/T2 contrast (CIT168 Amygdalar parcellation)
* Cerebellum: Deidrichsen SUIT
* Brainstem: Multimodal gray matter parcellation (Bianciardi, Licensing restriction)
* Brainstem: resting state gross parcellation (Shen)
* Midbrain PAG: BOLD 7T (Kragel et al. 2019)
* Midbrain SN, RN, STH: T1/T2 contrast (CIT168 amygdala parcellation, Pauli 2018)

There were two goals which motivated atlas construction. In order of priority
* Provide a probablistic spatial reference for functional localization of (mainly) group level results in multiple references spaces
* Provide parcels for automated full brain parcelwise analysis
The first requires a higher level of spatial detail than the second, so a fine scale and coarse scale version of 
the atlas are provided to facilitate both goals. The coarse atlas is designed for applications that need some slack, 
and differs from the fine atlas in the following ways.
* fine thalamic parcels are merged into larger subdivisions in the coarse atlas
* bianciardi atlas is coarse version rather than fine. Bianciardi offers a nested subdivision of the coarse atlas
  which includes subnuclei, while the coarse version does not take things as far. Aadditionally the coarse CANLab2023
  has bianciardi voxel probabilities increased to be given priority over surrounding (mostly white matter) Shen parcels. 
  This results in wider boundaries which are still informed by the original parcellation data.
* Amygdala is subdivided into Central, Centromedian and Basolateral instead of 10 subnuclei

These atlases are available in multiple spaces. The projection from one space to another differs depending on
the nature of the underlying data. See section below for details.
* MNI152NLin2009cAsym (fmriprep)
* MNI152NLin6Asym (fsl)
* LAS MNI152NLin2009cAsym (qsiprep)
* HCP Grayordinate

## Licensing

Use of this atlas is restricted according to the licenses of the constituent atlases. The constituent atlases are 
indicated in the labels_5 field of the atlas object and in the references property of the atlas object. Licenses
(when available) are included in Atlases_and_parcellations/2023_CANLab_atlas/licenses. For the most part the licenses are
unrestrictive except for Bianciardi and Morel. If you intend to use any parts of this atlas for any commercial applications
you would be wise to also contact the original authors of any constituent atlas for which no license is provided. The
absence of a license file here should not be taken as an indication of the absence of a usage license.

Bianciardi has a distribution restriction. CANLab2023 cannot be distributed directly as a result. Instead, we technically 
only distribute a partial atlas which includes the entire atlas less the bianciardi parcels and each user must assemble 
the remainder locally. In theory the abridged atlases could be used directly, but code to automate the construction process 
is provided and invoked automatically when you load the atlas.

Additionally, note that	the morel atlas	is also	restricted, and this usage restriction seems more severe than Bianciardi's.
Not only can we not distribute it, we also can't provide you with a legitimate resource for obtaining it short of contacting
the original host institution (see Atlases_and_parcellations/2023_CANLab_atlas/licenses/morel_copyright). Consequently, this
is the only atlas not directly available in Neuroimaging_Pattern_Masks (it's available to CANLab members in MasksPrivate, a
private repo). However, due to an oversight it was incorporated into canlab2018, which was distributed publically years ago, 
so that cat's out of the bag. Consequently, I (BP) didn't bother to separate the morel parcels out the way I did with 
Bianciardi's parcels. It's integrated into the 'abridged' version of this atlas directly distributed by this repo. Take
note if using this atlas in contexts where open licensing may be important to you.


## SETUP

To use this atlas with canlab tools you don't need to do anything, but taking the steps below will make subsequent usage flow
more smoothly. If you ignore this step then you'll simply have a slower evaluation the first time you run load_atlas on a 
new version of the atlas, nothing more. If you're using this atlas outside matlab you'll need these setup utilities.

This atlas technically has 10 different versions. 2 parcel granularities x 2 spatial resolutions x 2 sampling spaces + a
qsiprep compatible version and a CIFTI version. The CIFTI version cannot be built in matlab, and requires connectome 
workbench and (most likely) a traditional *nix environment to build. I'm not sure if it will work on Mac and will require
modifications to run on windows since the install script is written in the bash shell scripting language and would need to 
be ported to work in the windows powershell. If you're using Windows you'll likely encounter errors when trying to follow the
instructions below. Ignore these. The other 9 atlases should still be correctly built.

If you're not using windows modify the first couple of lines of the setup.sh script to point to your local copy of wb_command 
(provided by connectome workbench, available here: https://www.humanconnectome.org/software/get-connectome-workbench; I'm
currently using v1.5.0).

Run setup.m from matlab to construct all the necessary files. You will need an internet connection. Your PC should also 
ideally have at least 32G of memory and disk space, since the uncompressed probability maps are quite large. Once assembled,
the atlas files only take up a couple hundred MB. CIFTI file assembly is handled at the end by a system call to setup.sh.



## Usage

Most use is anticipated to be for analysis of functional data preprocessed by fmriprep data, in which case you
can load a useful version of the atlas like so,

canlab2023 = load_atlas('canlab2023_coarse_fmriprep20_2mm').threshold(0.2);

See "help load_atlas" in matlab for details (make sure CanlabCore is in your path).


### Bespoke applications

Although versions of this atlas have been provided that should be useful to most users most of the time if you have
exacting applications you may want generate atlases specific to your data and usage intentions. For localization purposes
probability maps are helpful, but if hard borders are desired probablistic thresholding and parcel defragmentation
may be desired. This can be accomplished using the atlas/threshold function. If you 
have specific data you intend to use this atlas with it may also make sense to resample it to your target space before 
thresholding. Presumably your data is already in alignment with the desired MNI template, but may be resampled to a 
different resolution, origin or orientation. You can potentially resample this atlas to a new space like so,

canlab = fmri_data(<path totemplate>)
canlab = load_atlas('canlab2023_coarse_fmriprep20'); % replace with desired parcellation and fsl6 version if appropriate
canlab = canlab.resample_space(template_2mm);
canlab = canlab.threshold(0.2); % adjust parameters as desired

Parcel defragmentation in particular may take some time. You may want to save the result for reuse. Saving and
reloading can be done like so,

save(<path>, 'canlab');
load(<path>, 'canlab');

Small regions will be severely affected by partial volume effects during resampling. In particular the brainstem. Your
best bet is to regenerate brainstem nuclei from source using the scripts in the 2023_Bianciardi* sister folder to this
directory and incorporate you target space into a custom apply_spm_warps.m script, similar to the existing one in the
templates/transforms/code folder. The idea is to accomplish all transformations, including your own resampling, in a
single step to avoid compounding partial volume erors across interpolations. This will probably require some substantial
engineering on your part to achieve, but to put these into perspective, this atlas took 4-5 weeks of dedicated work to 
assemble. It is likely worthwhile to build off of it rather than starting from scratch even if it takes a day or two of 
work.


### Grayordinates

The CIFTI files are not designed for use with CANLab tools, which as of this time does not support CIFTI formated
data. They can be used with connectome workbench or other toolboxes though. Run the setup scripts to generate them 
(see SETUP section above).


### QSIPrep

QSIprep is to DWI analysis what fMRIPrep is to fMRI analysis and some. Whereas fMRIprep only performs preprocessing, QSIPrep 
also performs what might be a the DWI analog of a first level analysis. I mean this in the sense that its qsirecon workflows 
automated common aspects of post-preprocessing analysis. In the case of tractography, the question of interest is what 
regions are connected to what other regions. The regions are defined by atlas files. Several default atlases are provided 
with qsiprep, but none are more detailed or exhaustive than the CANLab atlas. Thankfully qsiprep allows for use of custom 
atlases if appropriately formatted and accompanied by the necessary metadata (specified in a json file). To use the CANLab 
2023 atlas with qsiprep you need to pass a path to a directory containing the atlas and its atlas_config.json file in during 
the QSIprep run. If running qsiprep in a singularity container you do this by mounting this folder at a particular location 
within the virtual environment by using the --bind <full_atlas_directory_path>:/atlas/qsirecon_atlases. For use without a 
singularity container refer to the QSIprep docs,
https://qsiprep.readthedocs.io/en/latest/reconstruction.html

For the most basic usage (e.g. you only want to use the CANLab 2023 atlas in your qsirecon pipelines, and no other atlases) 
you can simply mount this directory in your singularity container after running the setup script (see SETUP above). You 
will only need to manually configure anything if you want to use your own custom atlases or use custom atlases in 
combination with qsiprep stock atlases. This will require merging the atlas_config.json with the stock atlas_config.json 
file to combine multiple atlases into a single reconstruction workflow. You will need the stock qsiprep atlases' files. 
They're linked to in the qsiprep docs, but for your convenience they can be found here,
https://upenn.box.com/shared/static/8k17yt2rfeqm3emzol5sa0j9fh3dhs0i.xz

Note, this has only been minimally tested. The parcellations in canlab2023 are in many cases poorly suited to tractography,
especially subcortical regions. I'm freezing development of CANLab2023 to focus on CANLab2024 though, and am not investing 
more time and effort into resolving this problem here. To understand some of the kinds of problems that you might run into 
though have a look at attempts to map hippocampal tractography in Dalton et al. (2022) eLife. Basically you'd likely need to
further modify the ROI masks to make allow tracts to penetrate into subcortical structures. Cortical structures should be
handled just fine though, and I've taken care to make sure that the nifti file is correctly formatted for QSIprep (correctly
oriented, sampled, right header information, etc).


### Usage Restrictions

There is also a .gitignore file in the relevant folders to prevent reupload of Bianciardi atlas files and their derivatives. 
Do not change this behavior without permission from Bianciardi.


## Mappings between spaces

Source parcellations are generally in idiosyncratic spaces and most needed to be projected into a new space for
any particlar space of interest. The following projections were used:

Registration fusion using fmriprep output to project labels from fsaverage space into MNI space.
* Glasser parcelation. See README in 2016_Glasser* sister folder of this one for details

Nonlinear registrations of MNI152NLin6Asym to MNI152NLin2009cAsym using ANTs SyN. Inverse transform was estimated and 
used for the MNI152NLin2009cAsym to MNI152NLin6Asym alignment. Some of the source atlases exclusively contained
subcortical structures, in which case a subcortical weighting mask was used when transformations were computed. Full
brain alignments were estimated using fmriprep 20.2.3 LTS with cortical reconstruction enabled. Subcortically weighted
alignments were run using a similar ANTs parameterization as fmriprep's, but with a subcortical mask for the cost 
function. The details of this alignment are documented in templates/transforms/code/subctx_alignment.sh.
* Hippocampal formation (full brain alignment), Julich atlas
* Basal Ganglia (subcortically weighted alignment), Tian atlas
* Thalamus (subcortically weighted alignment), Morel atlas
* Bianciardi brainstem parcels (subcortically weighted alignment), Bianciardi Atlas
* HCP grayordinate volumetric mask (subcortically weighted alignment; used as a mask, not directly incorporated)

Nonlinear registration of MNIColin27v1998 to MNI152NLin2009cAsym and MNI152NLin6Asym using ANTs SyN and implemented
by fmriprep with surface reconstruction enabled.
* Shen

CIT was already available in all desired spaces and was not transformed.

Alignments used linear interpolation of probability maps when available, or nearest neighbor interpolation otherwise 
(e.g. morel, shen).

For details regarding any particular regions transformation you should refer to the README file in that atlases 
directory in this repo. All transformations used are in the templates/transforms folder. They are available in
multiple formats for use with ANTs, FSL and SPM but were all computed using ANTs.


## Parcel Probabilities

* Cortex: likelihood a voxel will be circumscribed by a label in surface space, based on registration fusion
* Basal Ganglia: Same as cortex. Obtained by personal correspondence with the authors, not available elsewhere to my 
knowledge.  See Tian atlas directory README for details.
* MTL, Brainstem nuclei, CIT regions and cerebellum: Probablistic likelihoods that a labeled region was found in a participant 
at a particular voxel after alignment to standard space. Sample sizes are small, so don't expect them to be well
calibrated.
* Thalamus: bogus values I assigned since there weren't any natively
* Brainstem background regions (shen parcels): bogus values I asigned to be 0.35 to provide minimal constraints on nuclear
parcel boundaries. Anywhere these probabilities exceed those of Bianciardi atlas regions I impose a value 0.1 lower than 
bianciardi's labels. In tihs way I use a greedy algorithm to asign voxels to Bianciardi's brainstem nuclei when available.


## CANLab2018 Comparison

This atlas was created as a drop in replacement for CANLab2018. Taken individually the differences are relatively 
minor but extensive, and all together represent a substantial change. Differences are as follows
* Internally consistent grayordinate and volume formated versions available in register with multiple standard templates
* Probablistic cortical parcels obtained through registration fusion
* Removal of redundant Glasser cortical hippocampal segmentation. It doesn't comprehensively cover the hippocampal volume and intersects comprehensive cytoarchitectonic atlases.
* Different basal ganglia segmentation which better respects gross anatomical subdivisions and is probablistic
* New brainstem nuclear segmentation that is substantially more accurate and probablistic
* Correction of misalignments of MTL (severe), thalamus (severe), brainstem and cerebellum.
* (mostly) Probablistic (exceptions: thalamus and some brainstem filler regions)
* no RVM or trigeminal analog. The canlab2018 areas were not credible when compared against Duvernoy's Atlas (which is the
authoritative reference), so they were not carried over. Most other brainstem regions have analogs here, although potentially
under a different name (for instance the dorsal motor nucleus of the vagus, or DMNX, is now the viscero-sensory-motor nuclei,
or VSM).
* new PAG label, based on Phil's 7T columnar subdivisions
* New 10 region Amygdala labels instead of Julich 3 label subdivision.
* Perfect hemispheric subdivision (midline is assigned to x+ direction). Canlab2018 had irregular midline boundaries across regions.

## Methods

Coming soon? You could piece something together by referring to the README.md files in the constituent atlases' directories
and to the source scripts in the src subfolder here. Here's a bit re PAG though,

### PAG

Phil Kragel's PAG parcellation was redone to provide probablistic labels. 19/24 participants had good parcellations (1,2,
4-10,12,14,15,17,19,20-24). These participants were reprojected into their target space using transformations obtained
from Phil's dropbox into the same target space as the 2019 paper (IXI549) except linear interpolation was used instead of
cubic splines to avoid gibbs ringing. The results are saved in the source subfolder here as KragelPAG_MNI152NLin6Asym.nii.gz.
The space designation is justified because the IXI sample was registered to MNI152NLin6Asym before generating the IXI549
template used by Dartel to produce the warps used. Although there are differents between these templates the location and
orientation of the cerebral aqueduct is the same, so there's no need for further alignment to MNI152NLin6Asym space. Individual
subject alignments (partial volume effects and all) were averaged to produce a probablistic PAG map.

This procedure did not reproduce the PAG columns. These were not derived on a per subject level though so no probablistic
delineation between columns can be made. Instead we simply diluted the existing kragel2019pag atlas from this repository
to span a mask defined by the probablistic labels derived above and used nearest neighbor interpolation to label the newly
identified voxels within the dilution mask. These were then used to asign voxel probabilities to each of the individual
columns. Because we do not have subject specific probabilities the intercolumn probabilities are nonintersecting, but the
exterior margin of each column adopts the newly derived probablistic values.

## Parcel Discussion

Many of the parcels available here are quite small. In many cases you shouldn't consider their localization authoritative.
Take small brianstem nuclei or thalamic parcels for instance. It would be inappropriate to conclude involvement of a structure
based solely on whether or not your data overlaps said structure in this atlas. Even taking the probability maps into account 
you may want to take localization with a grain of salt and use multiple methods to establish the identity of a region. This 
atlas is best used in a blind and dumb fashion at coarse scale, or alternatively at a fine grained level as a rough guide 
indicating possible structural associations rather than definitive ones. Despite substantial effort devoted to alignment the 
templates used by the source datasets for generating many of these atlases are intrinsically limited. For instance, the 
MNI152NLin6Asym template used by the bianciardi atlas natively is exclusively a T1 template, but many of the structures in 
question have minimal to no T1w contrast. This leaves few landmarks for ANTs to exploit when computing mappings between 
structures. Additionally, the parcellation provided here, particularly for the brainstem, is incomplete. Other regions not 
yet included may be involved, which if involved might have a much more probable boundary for a region than any yet included 
here. For instance, we don't have an RVM or a trigeminal nucleus in this parcellation, but they certainly exist. 

Additionally, some of these regions are derived based on criteria that have only been investigated by single studies. In 
particular, many cortical parcels are based solely on Glasser's parcellation criteria, many basal ganglia parcels are based
solely on Tian's parcellation criteria, etc. Contrast that with the Julich atlas' cytoarchitectonic areas. Cytoarchitecture
has been successfully used as a way of demarcating regions for decades, while multimodal imaging metrics (Glasser) lack that
kind of validation to say nothing of novel criterial like functional connectivity streamline watersheds (Tian). Using these
kinds of novel criteria are helpful from the perspective of ontology: we can label brain regions in some consistent fashion.
We can't presuppose physiological significance of the boundaries however, and future updates to this atlas should revisit
such delineations.

Finally, note that this parcellation and associated probabilities are derived from a particular population of subjects that
may not be representative of yours. Most studies are disproportionately conducted on individuals between 25-35, who are
healthy and don't suffer from any pathological conditions that might affect brain morphology. The exception to this are
histlogical sections (morel and julich) which are obviously conducted on cadavers. Although the doners did not have brain
pathologies, they were nevertheless old enough to have anticipated their deaths and consented to the donation in the first 
place, so use your best judgement here. Finally, these are all studies conducted in western developed nations and are likely
biased towards associated demographics. The more your sample deviates from these characteristics the more you should take
in interpreting probability maps and label assignments.

## Developer Notes
Version tracking is handled by the *latest files which contain hashes of the 
associated objects. Running the create_CANlab2023_atlas.m script updates these
hashes automatically. If you modify the script or any other upstream atlases
that you want incorporated into a new version of the canlab2023 atlas you will
need to rerun create_CANLab2023_atlas.m for each of the 8 variations of the
canlab2023 atlas (two spaces x 2 resolutions x 2 granularities of parcellations)
so that the associated hashes are updated. If you do not do this then users
who pull updated git repos won't know that those atlas spaces/resolutions/parcellations
have been updated and load_atlas() won't know to update their atlas *.mat files.
For convenience I've included a setup.m script that runs the creation script for all 
atlas versions, so you can just run that.

## References

### README references

Dalton M, D'souza A, LV J, Calamante F. New insights into anatomical connectivity along the anterior-posterior axis of the human hippocampus using in vivo quantitative fibre tracking. eLife (2022) 11.

### atlas source references 

* Amunts et al (2023) [Dataset v3.0.3] DOI:10.25493/56EM-75H
* Amunts, K., Kedo, O., Kindler, M., Pieperhoff, P., Mohlberg, H., Shah, N. J., Habel, U., Schneider, F., & Zilles, K. (2005). Cytoarchitectonic mapping of the human amygdala, hippocampal region and entorhinal cortex: intersubject variability and probability maps. Anatomy and Embryology, 210(5–6), 343–352. https://doi.org/10.1007/s00429-005-0025-5 DOI: 10.1007/s00429-005-0025-5                                                      "
* Amunts, K., Mohlberg, H., Bludau, S., & Zilles, K. (2020). Julich-Brain: A 3D probabilistic atlas of the human brain s cytoarchitecture. Science, 369(6506), 988–992. https://doi.org/10.1126/science.abb4588 DOI: 10.1126/science.abb4588
* Diedrichsen, Jörn, Joshua H. Balsters, Jonathan Flavell, Emma Cussans, and Narender Ramnani. 2009. A Probabilistic MR Atlas of the Human Cerebellum. NeuroImage 46 (1): 39?46.
* Glasser, Matthew F., Timothy S. Coalson, Emma C. Robinson, Carl D. Hacker, John Harwell, Essa Yacoub, Kamil Ugurbil, et al. 2016. A Multi-Modal Parcellation of Human Cerebral Cortex. Nature 536 (7615): 171?78.
* Petre B, Ceko M, Friedman NP, Keller MC, Lindquist MA, Losin EAR, Wager T. [Dataset v9] DOI: 10.6084/m9.figshare.24431146
* Jakab A, Blanc R, Berényi EL, Székely G. (2012) Generation of Individualized Thalamus Target Maps by Using Statistical Shape Models and Thalamocortical Tractography. AJNR Am J Neuroradiol. 33: 2110-2116, doi: 10.3174/ajnr.A3140.
* Kedo, O., Zilles, K., Palomero-Gallagher, N., Schleicher, A., Mohlberg, H., Bludau, S., & Amunts, K. (2017). Receptor-driven, multimodal mapping of the human amygdala. Brain Structure and Function. https://doi.org/10.1007/s00429-017-1577-x DOI: 10.1007/s00429-017-1577-x
* Kragel, P. A., Bianciardi, M., Hartley, L., Matthewson, G., Choi, J. K., Quigley, K. S., ... & Satpute, A. B. (2019). Functional involvement of human periaqueductal gray and other midbrain nuclei in cognitive control. Journal of Neuroscience, 2043-18.
* Krauth A, Blanc R, Poveda A, Jeanmonod D, Morel A, Székely G. (2010) A mean three-dimensional atlas of the human thalamus: generation from multiple histological data. Neuroimage. 2010 Feb 1;49(3):2053-62. Jakab A, Blanc R, Berényi EL, Székely G. (2012) Generation of Individualized Thalamus Target Maps by Using Statistical Shape Models and Thalamocortical Tractography. AJNR Am J Neuroradiol. 33: 2110-2116, doi: 10.3174/ajnr.A3140
* Pauli, Wolfgang M., Amanda N. Nili, and J. Michael Tyszka. 2018. ?A High-Resolution Probabilistic in Vivo Atlas of Human Subcortical Brain Nuclei.? Scientific Data 5 (April): 180063.
* Shen X, Tokoglu F, Papademetris X, Constable R. Groupwise whole-brain parcellation from resting-state fMRI data for network node identification. Neuroimage 82, 403-415, 2013.
* Tian Y, Margulies D, Breakspear M, Zalesky A (2020). Nature Neuroscience. 23(11) 1421-1432.
* Tyszka, J. M. & Pauli, W. M. In vivo delineation of subdivisions of the human amygdaloid complex in a high-resolution group template. Hum. Brain Mapp. 37, 3979–3998 (2016).
* Wu J, Ngo GH, Greve D, Li J, He T, Fischl B, Eickhoff SB, Yeo T. Accurate nonlinear mapping between MNI volumetric and FreeSurfer surface coordinate systems. 2018. Human Brain Mapping 39(9) 3793-3808. DOI: 10.1002/hbm.
