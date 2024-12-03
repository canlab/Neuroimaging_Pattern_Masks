The spatial maps here are the principal gradients of transcriptomic variation computed 
from 6 postmortem brains registered to MNI152NLin6Asym space. The gradients were
generated as part of the following paper:

Vogel et al. "Deciphering the functional specialization of whole-brain spatiomolecular gradients in the adult brain" (2024) PNAS

The transcritomic characterization (and pinricpal paper in some sense) is here:

Hawrylycz et al. An anatomically comprehensive atlas of the adult human brain transcriptome. (2012) Nature

The code saved as create_gradients.py is a condensed version of the code provided by Vogel et al. at

https://github.com/PennLINC/Vogel_PLS_Tx-Space/

After generating transcriptomic_gradients.nii.gz the results were smoothed using the following command

fslmaths transcriptomic_gradients.nii.gz -dilall -s 2.548 -mas ~/.matlab/canlab/Neuroimaging_Pattern_Masks/templates/MNI152NLin6Asym_T1_1mm.nii.gz transcriptomic_gradients.nii.gz

It isn't entirely clear what template the coordinates from Hawrylycz et al. are registered to in Vogel et al. It's
not the original space but rather a nonlinear one but there are actually two available. One is available here 

https://github.com/gdevenyi/AllenHumanGeneMNI/

And another is in the poster and gitrepo from 

Gorgolewski KJ, Fox AS, Chang L, Schäfer A, Arélin K, Burmann I, Sacher J, Margulies DS. "Tight fitting genes: finding relations between statistical maps and gene expression patterns" (2014) OHBM

Neurovault also implements a spatial correlation analysis with these genetic maps using the nonlinear coordinates, and
uses the FSL template, which is MNI152Nlin6Asym.

See here how get_standard_mask is called here, and this provides the affine matrix they use:

https://github.com/NeuroVault/NeuroVault/blob/413e8ff64717bcba4e53b59f6f46b2f4a61c70a8/scripts/preparing_AHBA_data.py

get_standard_mask comes from pybraincompare

https://github.com/vsoch/pybraincompare/blob/master/pybraincompare/mr/datasets.py

And they use FSL's MNINlin6Asym

A comparison of the mni_nlin_x/y/z coordinates from 
https://github.com/gdevenyi/AllenHumanGeneMNI/blob/master/transformed-points/H0351.2001_SampleAnnot_RAS_mni_nonlin.csv
with those obtained from data/gxp_correlation_wholebrain_results_NEW.csv from https://github.com/PennLINC/Vogel_PLS_Tx-Space
suggests that these are the same though, and the AllenHumanGeneMNI repo uses MNI152NLin2009cSym, so I think that's the
original template here.

Based on visual inspection the results of the gradients seem to align better with the boundaries of MNI152NLin6Asym 
than the slightly larger MNI152Nlin2009cSym. At the end of the day it's splitting hairs though. With the level of 
uncertainty there is in the source spaces and the level of precision you might expect from slices in 6 participants, 
it doesn't really matter. I've smoothed the data by a 6mm kernel regardless to account for this variation.

If trying to figure this out for yourself, note that subject IDs in these datasets come in two formats, one is the API
format, the other is the web format. The abagen code has the mapping, which you can also see here:
https://abagen.readthedocs.io/en/stable/_modules/abagen/datasets/fetchers.html


The CIFTI version of the gradients was created using registration fusion. The AHBA subject freesurfer surfaces from 
the abagen repo and https://github.com/gdevenyi/AllenHumanGeneMNI/ MNI to native T1 transformations were used.

Bogdan
11/28/2024
