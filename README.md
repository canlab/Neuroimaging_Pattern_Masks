# Neuroimaging_Pattern_Masks

This repository contains pre-defined brain "signatures" (multivariate predictive patterns), atlases of local regions and networks, and masks and regions derived from published meta-analyses of neuroimaging data. It includes a fairly comprehensive set of such resources developed by the Cognitive and Affective Neuorscience Lab (Tor Wager, PI) and our collaborators, and also includes some products from other groups shared publically or by permission from the creators.  

Some of these resources are used in other toolboxes, particularly the CAN Lab’s <a href = "https://github.com/canlab/CANlab_help_examples">Help Examples and Batch Scripts</a> repository. They are also very useful when doing interactive analysis with the CAN lab's object-oriented neuroimaging toolbox, <a href = "https://github.com/canlab/CanlabCore">Canlab Core Tools</a>. 

The three types of brain maps included are:
- Pre-defined brain "signatures" (aka multivariate predictive patterns, brain biomarkers, or "neuromarkers") that can be applied to new individual participants to generate predictions and validate predictive models. Most CANlab signatures are publically available and can be downloaded here. A few, the Neurologic Pain Signature (NPS) and fibromyalgia-predictive patterns, are available for research use upon request (contact Prof. Tor Wager). 

- Atlases with pre-defined brain parcels (regions) and networks. This can reduce brain space to a smaller set of (hopefully) meaningful units of analysis. These are saved as Analyze (.img) or NIFTI (.nii) files, and also as "atlas"-type objects, an object type defined in <a href = "https://github.com/canlab/CanlabCore">Canlab Core Tools</a> that facilitates working with brain atlases.

- Brain maps from published meta-analyses of neuroimaging data, which define consensus regions across studies for multiple psychological/task categories -- e.g., emotion, working memory, PTSD, and more. These masks can be used to specify a priori regions of interest or as "patterns of interest" in new studies.

Getting help and additional information:
------------------------------------------------------------
We have several sources of documentation for this repository. A website with references and pictures is here: https://sites.google.com/dartmouth.edu/canlab-brainpatterns/home

1.  You can use any software to load, view, and apply the models. We have an object-oriented Matlab toolbox that makes it easy to load and apply the models. For the philosophy behind the  object-oriented toolbox and code walkthroughs in Matlab, see canlab.github.io. For function-by-function help documents on the Core Tools objects and functions, see the <a href = http://canlabcore.readthedocs.org/en/latest/>help pages on Readthedocs</a>. The code to create the walkthroughs, and a batch script system that uses the CanlabCore object-oriented tools for second-level neuroimaging analysis, is here: <a href='https://github.com/canlab/CANlab_help_examples'>CANlab_help_examples github repository</a>

2. The CANlab website is https://sites.google.com/dartmouth.edu/, and we also maintain a WIKI with more information on some of our toolboxes and fMRI analysis more generally, which is <a href = "https://canlabweb.colorado.edu/wiki/doku.php/help/fmri_tools_documentation">here</a>.  For more information on fMRI analysis generally, see <a href = "https://leanpub.com/principlesoffmri">Martin and Tor's online book</a> and our free Coursera videos and classes <a href = "https://www.coursera.org/learn/functional-mri">Principles of fMRI Part 1</a> and <a href = "https://www.coursera.org/learn/functional-mri-2">Part 2 </a>.

Dependencies: These should be installed to use this toolbox
------------------------------------------------------------
Matlab www.mathworks.com

<recommended> Matlab statistics toolbox
  
<recommended> Matlab signal processing toolbox
  
<recommended> Statistical Parametric Mapping (SPM) software https://www.fil.ion.ucl.ac.uk/spm/

<recommended> the CANlab Core Tools repository https://github.com/canlab/CanlabCore
  
<recommended> the canlab_help_examples repository  https://github.com/canlab/CANlab_help_examples
  
Want to contribute? Please get in touch
------------------------------------------------------------  
Most of the maps/models come from our lab, but some are those that other groups have agreed to share and which we find particularly useful. We'd love to grow the collection, so if you want to contribute an atlas, meta-analysis results, or multivariate model, please get in touch! Maryam Amini is a good contact <maryam.amini@dartmouth.edu>.
