## Tian 2020

This atlas was adapted from a different repo,

https://github.com/yetianmed/subcortex/tree/e80ee787732536e9e89534c9a623b10aff7928f4

Group-Parcellation files were copied here for posterity. These may now be obsolete though.
While the first draft of this atlas was derived from these parcels, subsequently the
authors shared their individual participant segmentations with me (BP). I used these
to create probablistic labels which now form the basis of this atlas.

This atlas was generated in a hierarchical fashion with finer and finer parcellations
at each level until a natural stopping condition in the parcellation algorithm
was reached (see Tian et al \[2020\] Nat Neuro for details). We store this information
using labels_2, labels_3, labels_4, fields of the atlas. By default, the voxel map is 
indexed at the finest scale available (54 parcels), and each index label corresponds
to an incrementally more coarse parcellation. To obtain a voxel map of one of these 
coarser parcellations you can invoke atlas/get_coarser_parcellation. E.g. to reduce 
the 54 parcels to their coarsest parent parcellation (labeled by labes_4) you can use,

fine_tian = load_atlas('tian_fmriprep20')
reduced_tian = fine_tian.get_coarser_parcellation('labels_4')

Furthermore, This atlas is provided registered to two different reference templates 
for use with FSL/HCP or fMRIprep data as the case may be, which is likewise specified 
by FSL6 or fmriprep20 (the latest LTS) in the name when invoking load_atlas.

## Atlas space

The authors of this atlas are professionally invested in deep brain stimulation 
treatments and require greater precision from their segmentations than is typical
in fMRI research. Consequently they make a point of specifying their reference spaces.
I've made efforts to retain this specificity here, so multiple versions of this atlas
are provided which correspond to multiple common reference spaces.

As of this writing, most wagerlab users are relying on fMRIprep for spatial normalization.
This aligns all data to the MNI152NLin2009cAsym reference space. To use this atlas please
use the Tian_3T_S4_2mm_fmriprep20 segmentation, and use the appropriate underlay
for visualization. Please see the "templates" folder of this repo for more details.

Note: fmriprep20 refers to version 20 of fmriprep, which as of this writing is the 
most recent "long term release" of this software. fsl6 likewise refers to fsl version 6.

## Methods

I combined the individual subject segmentation provided by Tian with with familial information	
from the HCP restricted	data to	derive probability maps from unrelated individuals. For	
convenience I reused a list of i.i.d. subjects I'd generated for other purposes. Those other 
purposes required first excluding any participants that lacked any GLM task contrasts or 
either resting state sessions. The intersection of this iid list with the participants Tian had
used resulted in a sample of 377 subject's segmentations. The total number of iid participants
in the HCP dataset is in the 450 range, but they don't all have resting state data, which is
needed to estimate this parcellation in the first place, so the maximum number of i.i.d.
participants I might have used falls somewhere below 450, but above 377.

The full parcellation that Tian et al release has very clean delineations of large scale
structures, but the probablistic map I have is not so clean. For instance putamen can
bleed into accumbens, accumbens can bleed into caudate, thalamus can bleed into the fornix, 
etc. To achieve a cleaner parcellation I've modified the probability to prevent parcel
probability maps from overlapping inappropriate structures. Each of the 54 parcels was 
assigned to a CIFTI structure based on which subcortical CIFTI structure had the largest 
intersection. Any probability values that overlapped with other CIFTI structures were then
set to zero.


## Comparisons with Pauli 2016

Pauli 2016 is a basal ganglia specific atlas which predates this one in the literature and in canlab tools. It 
was developed in collaboration with Tor, and is based on clustering of coactivation maps on Neurosynth. 

The Tian atlas is based on gradients in functional connectivity which are evaluated using a novel tractography based method.
Spectral methods are used to define gradients of change in voxel-level functional connectivity fingerprints. Each
spectral basis constitutes a gradient, which together can be used to define tensors based on their direction of
maximal change. Tractography methods can be applied to these tensors. The authors define boundaries between regions
based on 'watersheds' in the tractograms. Consider how water flows in a mountain range. The watersheds are the 
ridgelines, like the continental divide in the American Rockies. The reason for this complicated approach is that there's
debate in the field regarding whether or not the basial ganglia can be divided into discrete parcells or instead exist
along a gradient. The idea is that there are cortico-basal ganglia loops which are topographically organized and define
distinct zones, but at the same time there's also integration within the basal ganglia, which is import for integration
action selection, motor function, and reward learning (Haber & Knutson [2010] Neuropsychopharmacology). Tian's gradient+watershed
based approach formally models all of this and is able to both respect the gradient like organization of the subcortex
while simultaneously providing a map of subdivisions defined by subtle fluctuations in those gradients. This leads 
naturally to a hierarchy of parcellation with increasingly 'soft' watersheds as the parcellations becomes more detailed.

In my opinion the Tian atlas provides a better parcellation than Pauli for two reasons,
1) it's directly derived form data based based on a large and high quality sample (1000 HCP participants). This gives it 
greater precision than is possible with neurosynth meta-analysis maps.
2) Tian delineates accumbens core and shell, which is something pauli doesn't do, but which I think is especially relevant 
for a lab interested in affect and pain like the Wager lab. Additionally, I find this particularly impressive because I've 
seen people try and fail to do this with BOLD in the past, and only succeed using white matter tractography. Tractography 
appears to be the gold standard in general (Haber & Knutson [2010] Neuropsychopharmacology), so this achievement lends credence 
to the superiority of Tian et al's approach.

## References

* Tian Y, Margulies D, Breakspear M Zalesky A. Topographic organization of the human subcortex unveiled with functional
connectivity gradients. Nature Neuroscience 23(11), 1421-1432, 2020.
* Haber S & Knutson B. The reward circuit: Linking primate anatomy and human imaging. Neuropsychopharmacology 35(1), 4-26, 2010.

## 

Bogdan Petre <br />
10/10/23
