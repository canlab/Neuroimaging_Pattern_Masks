## Overview

This atlas is based on histology, immunochemistry and ex vivo DWI tractography 
in 3 postmortem individuals. The resulting parcels have been extensively 
characterized in terms of their cortical, subcoritcal and commisural 
structural connectivity, and functional connectivity with resting state 
functional connectivity, specifically the DMN (Edlow, 2023).

Although the measurement techniques are extensively documented (Edlow 2012, 
2023), technical details of the projection to MNI space are obscure.  The 
atlas is in MNI 1mm space but it's not clear what template was used to achieve 
this projection. Based on the publication date of the original atlas (v1.0), 
the use of FSL tools for tractography, and under the assumption that the same 
templates were used for the current atlas version (v2.0, a safe assumption: 
differences with version 1.0 are detailed here: 
https://datadryad.org/stash/dataset/doi:10.5061/dryad.zw3r228d2 and they are
minor as far as spatial characteristics are concerned), I've assumed they 
used the FSL standard 1mm template, which is MNI152NLin6Asym, and used the
subcortically weighted alignment from MNI152NLin6Asym to MNI152NLin2009cAsym
space in the templates/transforms folder to project it to MNI152NLin2009cAsym
space. That said, given the obscurities associated with the template selection 
the boundaries of these regions should be considered approximate.

Consider comparing regions with the Bianciardi brainstem nuclei. All of these
are also available there, albeit under a more restrictive distribution license.

## License

This dataset is licensed under creative commons (CC0) 1.0 Universal and is in
the public domain.

## External Links

* https://datadryad.org/stash/dataset/doi:10.5061/dryad.zw3r228d2
* https://www.nmr.mgh.harvard.edu/resources/aan-atlas


## References

Note: Edlow 2023a is the recommended citation

* Edlow, BL; Kinney, HC (2023a), Harvard Ascending Arousal Network Atlas – Version 2.0, Dryad Digital Repository, https://doi.org/10.5061/dryad.zw3r228d2
* Brian L Edlow, Mark Olchanyi, Holly J Freeman, Jian Li, Chiara Maffei, Samuel B Snider, Lilla Zollei, Juan Eugenio Iglesias, Jean Augustinack, Yelena G. Bodien, Robin L Haynes, Douglas N Greve, Bram R Diamond, Allison Stevens, Joseph T Giacino, Christophe Destrieux, Andre van der Kouwe, Emery N Brown, Rebecca D Folkerth, Bruce Fischl, Hannah C Kinney. (2023b), Sustaining wakefulness: Brainstem connectivity in human consciousness bioRxiv 2023.07.13.548265; doi: https://doi.org/10.1101/2023.07.13.548265
* Edlow BL, Takahashi E, Wu O, Benner T, Dai G, Bu L, Grant PE, Greer DM, Greenberg SM, Kinney HC, Folkerth RD. Neuroanatomic connectivity of the human ascending arousal system critical to consciousness and its disorders. (2012) Journal of Neuropathology and Experimental Neurology. 71:531-546.  PMCID: PMC3387430.
