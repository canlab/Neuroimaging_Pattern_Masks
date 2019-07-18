please find the attached S1.zip file, which includes SVM weight map
(MNI space, 2-mm resolution)and its intercept.

this weight map was calculated using paired-SVM method which highlights
individual difference in the PRE and POST-maneuver data. so if you have
REST before the pain exacerbation (using cuff bladder!) and after, that
would be proper dataset to test. to apply the weight, the data need to be
prepared into two forms: POST-PRE set and PRE-POST set. then apply
classification for 'POST-PRE (e.g., positive class)' vs 'PRE-POST
(negative)' (i.e., paired-SVM).

if you are comparing PRE- and POST-therapy, hopefully POST has less pain
levels compared to the PRE. then, probably i would set 'PRE-POST' as
positive class and 'POST-PRE' as negative :) let me know what you think!