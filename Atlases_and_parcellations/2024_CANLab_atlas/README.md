## Usage

### Parcellation granularity
labels - "Fine" granularity. Designed for use in BOLD signal attribution, 
under assumptions that all probable regions can be considered simultaneously.
labels_2 - "Coarse" granularity. Designed for BOLD signal extraction under 
assumption that a winner-takes-all parcellation will be used. Neighboring regions 
that are unlikely to be distinguishable a-priori based on atlas labels are merged 
relative to labels_1
labels_3 - Even coarser granularity deisnged for prototyping.


## Subregion composition

### Amygdala

CIT168 amygdala atlas. 
labels - Intercalated nuclei are also subdivided into subregions for grouping at coarser scales
labels_2 - intercalated nuclei are assigned to their nearest other nucleus. Central nucleus (CEN) and Crotical median nuclei (CMN) are also combined.
labels_3 - Nuclei are subdivided into lateral, basolateral and centromedian (superficial) nuclear subdivisions.
labels_4 - amygdala


## Comparison with Canlab2023

### Major changes
* Tian NAc shell/core like regions replaced with Cartmell DWI based Core/Shell
parcellation
* Hypothalamus and Mammillary bodies replaced with Billot/Iglesias/Freesurfer 
hypothalamic subnuclei.
* Morel Thalamic parcellation replaced with Iglesias 2018/Freesurfer thalamic
parcellation. This includes LGN/MGN, so we lose parvocellular/magnocellular
LGN subdivisions.
* Bianciardi Dorsal Raphe, Locus Coeruleus, and Vestibulo-sensory-motor 
regions replaced by Levinson-Bari Dorsal Raphe, Locus Coeruleus and nucleus
tractus solitarius, respectively (not necessarily better or worse, just have a
more generous license).

### Minor changes
#### Amygdala
* Subdivision of intercalated nuclei of amygdala. The original authors could not 
identify these with the contrasts available and assigned these labels by default
to a meshwork of voxels surrounding the other nuclei. This leads to an awkwardly
shaped and fragmented parcel, and doesn't play well with downsampling across
granularities. To address this in CANLab2024 I (BP) have subdivided the IC nuclei
using a nearest neighbor parcellation based on the other 9 larger nuclei. These 
are available at the highest level of granularity ('labels') of the 1mm atlas. 
At labels_2 in the 1mm map and both labels and labels_2 of the 2mm map IC nuclei 
are subsumed by their nearest larger nucleus.
* Amygdala labels 2 now includes 9 regions instead of 3 (the 9 regions + neighboring
voxels of intercalated nuclei). The 3 nuclei subdivision (CE, BLA and CM) remains
available at labels_3 which in canlab2023 was redundant with labels_2 anyway.
* Brain mask was slightly dilated in certain areas to provide better coverage of 
BST_SLEA, VeP and anterior hypothalamic areas. This results in a slight inconsistency 
between grayordinate and volumetric versions of the atlas.
#### Hippocampus
* CA2/CA3 and CA1 are kept separate in coarse parcellation as two distinct
regions rather than being merged into one.
#### Brainstem
* Cerebral aqueduct was masked in Shen filler parcels


## openCANLab2024

openCANLab2024 is variation on CANLab2024 that only uses source atlases with open
usage and distribution licenses. This amounts to substitution of regions from the
Harvard Ascending Arousal Network Atlas, the Levinson-Bari Limbic Brainstem Atlas,
and the CIT168 atlas for a number of regions otherwise provided by the Bianciardi
Atlas. The Bianciardi Atlas appears to be superior to these, but has distribution
restrictions. CANLab2024 has some workarounds for these distribution restrictions
but the end result is ineviably less portable or versatile, hence the need for
an open software variant.

The following substitutions were made:

Bianciardi -> open Atlas region (Atlas Name) -> New Name
DR_B7 -> DR (Levinson-Bari) -> DR_B7 <br />
LC -> LC (Levinson-Bari) -> LC <br />
VSM -> NTS (Levinson-Bari) -> NTS <br />
LDTg_CGPn -> LDTg (Harvard AAN) -> LDTg <br />
MPB/LPB -> PBC (Harvard AAN) -> MPB_LPB <br />
PTg -> PTg (Harvard AAN) -> PTg <br />
PnO_PnC_B5 -> PnO (Harvard AAN) -> PnO <br />
isRt -> mRt (Harvard AAN) -> isRt <br />
MnR/PMnR -> MnR (Harvard AAN) -> MnR <br />

The following regions are missing from openCANLab2024, and do not appear to be
available anywhere besides the Bianciardi atlas at the moment:


Of the substitute atlases in openCANLab2024, the Harvard AAN atlas is not 
probablistic. We set synthetic probabilities to 0.8 for all labeled voxels and
dilated the very small regions (MPB_LPB, PTg and PnO) with a 3mm fwhm gaussian 
smoothing kernel to provide some more consistent behavior between atlases. In 
all cases these regions remain more conservatively sized than in Bianciardi.

The substitute regions do not all perfectly overlap with Bianciardi regions. The
region plots below illustrate the degree of overlap.
