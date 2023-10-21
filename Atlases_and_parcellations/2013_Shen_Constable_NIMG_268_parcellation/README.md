## Source Notes

The original source of the shen atlas source files is unknown. It predates me (BP). I downloaded a copy from NITRC though,
https://www.nitrc.org/frs/?group_id=51
It overlaps perfectly with the legacy shen atlas parcelation, so I'm pretty sure this is the same file, albeit maybe not
the same source. I haven't modified the 2mm version, but I wanted a 1mm version.

## Space Notes

The original parcellation appears to be in Colin27 space. This is based on the fact that the reference paper says it used
BioImage Suite. If you download bioimage suite and look at its source image it looks identical to a skull stripped version
of Colin27. When you then overlay the shen parcels on the Colin27 it shows perfect alignment, but misalignment when used
on MNI152NLin6Asym (FSL6 standard template) or MNI152NLin2009cAsym (fmriprep standard template). This is most obvious in
the ventral parts of the brain where Colin27 is most atypical. You can also find references to Colin27 in the contemporary
bioimage suite software docs (https://medicine.yale.edu/bioimaging/suite/bioimagesuite_manual_95522_2907_5_v1.pdf).

Data was aligned from MNI Colin 27 v 1998 (there's also the 2008 version, but here we used the 1998 one) space to other
MNI spaces using ANTs with MultiLabel interopolation (sigma = 0.5). The sigma defines a smoothing kernel. 0.5 was selected
from among 0, 0.5, 1 and 2 based on comparing the results to the original data. sigma=0 (nearest neighbor) was more jagged
than the original, while 1 was smoother than the original. 0.5 seemed to produce parcels which approximated the smoothness
of the original best.

Bogdan Petre
10/12/2023
