# Templates

Reference T1-weighted templates in common MNI / stereotactic spaces.
These are the canonical target spaces that atlases and signatures in
this repo are aligned to. Different software packages and atlas builds
use different MNI templates — collect them here so you can verify
alignments at a glance and use the right one for your pipeline.

> Although all data provided by **Neuroimaging Pattern Masks** should
> be in MNI space, different reference templates may have been used
> for spatial normalisation for generating specific atlases and maps.
> This is probably not consistently documented, but some common
> templates are provided here for those cases in which the standard-
> space template is described.

For a helpful discussion of MNI templates, see
[lead-dbs.org / About the MNI spaces](https://www.lead-dbs.org/about-the-mni-spaces/).
Most official MNI atlases are at
[nist.mni.mcgill.ca/atlases](https://nist.mni.mcgill.ca/atlases/),
and the
[BIDS coordinate-systems table](https://bids-specification.readthedocs.io/en/stable/appendices/coordinate-systems.html)
is a useful cross-reference of which software uses which template.

Naming convention: `<ref_space>_<modality>_<resolution>.<extension>`.

## Templates included

### HCP / FSL

HCP and current FSL use the asymmetric MNI152 nonlinear 6th-generation
template — **MNI152NLin6Asym**. Originally distributed via
[`HCPpipelines/global/templates`](https://github.com/Washington-University/HCPpipelines/tree/master/global/templates).

- `MNI152NLin6Asym_T1_1mm.nii.gz`
- `MNI152NLin6Asym_T1_2mm.nii.gz`

You can confirm FSL still tracks this build by checking that
`$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz` matches the 2 mm file
above.

### fMRIPrep / QSIPrep

**MNI152NLin2009cAsym** (used by [fMRIPrep](https://fmriprep.org/en/stable/spaces.html#standard-spaces)
and QSIPrep). Atlases destined for these preprocessing pipelines
should be in this space.

- `MNI152NLin2009cAsym_T1_1mm.nii.gz`
- `MNI152NLin2009cAsym_T1_2mm.nii.gz`
- `MNI152NLin2009cAsym_T1_1mm_brainmask.nii.gz` — brain mask.
- `MNI152NLin2009cAsym_1mm_t1s_lps.nii.gz` — same volume in LPS orientation.
- `MNI152NLin2009bAsym_T1_0.5mm.README` — note on the 2009b variant.

Pulled from the [TemplateFlow](https://www.templateflow.org/browse/) repo.

### SPM12

SPM12 uses a slightly modified MNI152 linear template, distributed as
`canonical/avg152T1.nii` in SPM and provided here for SPM-derivative
work outside of SPM:

- `IXI549Space_T1_2mm.nii.gz`

### SPM99–SPM8

These use **MNI152Lin**. Source:
[bic.mni.mcgill.ca / ICBM152Lin](https://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152Lin).

- `MNI152Lin_T1_1mm.nii.gz`

### SPM96 / BioImage Suite

**MNIColin27** — a single-subject template (Colin scanned 27 times,
co-registered and averaged). High detail but less representative.
Used when extreme spatial detail is needed (e.g., histology
registration). Source:
[nist.mni.mcgill.ca / Colin 27 average brain](https://nist.mni.mcgill.ca/colin-27-average-brain/).

- `MNIColin27_T1_1mm.nii.gz`

This is the 1998 version; the 2008 version has slightly different
geometry.

## Sub-folders

| Folder | Description |
| --- | --- |
| [cerebellum](./cerebellum) | Cerebellar atlas / template resources. |
| [transforms](./transforms) | Cross-template transforms (e.g., MNI152NLin6Asym ↔ MNI152NLin2009cAsym). |

## Contributing

Add any template you find useful here, with a comment in this document
on when to use it. The reference template, spatial resolution, and
imaging modality should be clear from the filename.

— *Bogdan Petre, 10/10/23 (extended into landing-page form here).*
