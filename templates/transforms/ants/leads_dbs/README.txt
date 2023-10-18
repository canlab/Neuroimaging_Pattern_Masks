The ANTs (*.h5) transforms were obtained from leads-dbs via figshare,
https://figshare.com/articles/dataset/MNI_T1_6thGen_NLIN_to_MNI_2009b_NLIN_ANTs_transform/3502238

For more details see here,
https://www.lead-dbs.org/about-the-mni-spaces/

It seems that they map from MNI152NLin6Sym to MNI152NLin2009bSym, but what we need is a mapping
from MNI152NLin6Asym to MNI152NLin2009cAsym, the assymetric spaces. The leads_dbs transforms
had the involvement of someone from MNI, and are likely to be better mappings that what I've
managed to generate myself, so which mapping you use depends on what you think matters more:
a quality inter template matching, or for the templates in question to be the asymetric templates
rather than the symmetric templates.

refer to ../../code for code that will help you convert this h5 data to fsl and spm formats
if you wish
