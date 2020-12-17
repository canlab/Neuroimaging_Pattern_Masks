# -*-coding:Utf-8 -*

import os

import nipype.interfaces.ants as ants
import BCBlib.tools.constants as cst

proba = os.path.join(cst.get_ants_priors_folder(), "brainPrior.nii.gz")
temp = os.path.join(cst.get_ants_priors_folder(),
                           "brainWithSkullTemplate.nii.gz")

def skull_strip(dim, anatomical, brain_proba=proba, template=temp):
    be = ants.BrainExtraction()
    be.inputs.dimension = dim
    be.inputs.anatomical_image = anatomical
    be.inputs.brain_template = brain_proba
    be.inputs.brain_probability_mask = template
    be.cmdline

# antsBrainExtraction.sh -d 3 -a /data/Chrisbrain/S10_T1.nii.gz -e
# Tools/extraFiles/Priors/brainWithSkullTemplate.nii.gz -m
# Tools/extraFiles/Priors/brainPrior.nii.gz -o /data/proper_atropos/Chris_
