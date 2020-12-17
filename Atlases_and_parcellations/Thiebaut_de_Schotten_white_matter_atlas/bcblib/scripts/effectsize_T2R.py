"""
Convert T maps into R maps for effect size

Authors: Chris Foulon & Michel Thiebaut de Scotten
"""

#see http://www.bcblab.com/BCB/Coding/Coding.html for a step by step installation of python
#don't forget to pip install nibabel :)

### modules ###
import os
from os import system as oss
import sys
import pandas as pd
import numpy as np
import nibabel as nib
print("python Effectsize_T2R.py numberofparticipants tvaluemap.nii.gz output.nii.gz")


# Defines the variable
n = int(sys.argv[1])
tval = sys.argv[2]
Rpower = sys.argv[3]

#import the nifti
i_tval = nib.load(tval)

# extract the data from the images
data_tval = i_tval.get_fdata()

#create df n - 2
df = np.subtract(n,2)
#Create t^2
step1 = np.square(data_tval)
#Create t^2 +44
step2 = np.add(step1,df)
#create (t^2 + 44)^1/2
step3 = pow(step2,0.5)
#create R = t / (t^2 + 44)^1/2
step4 = np.divide(data_tval,step3)

#save the data
i_Rpower = nib.Nifti1Image(step4, i_tval.affine)
nib.save(i_Rpower, Rpower)
