This package provides code and data for users who are interested in investigating functional network reorganization during pain, and also for further replication of our findings.

Reference: 
Zheng, W. et al. (2019). Pain-evoked reorganization in functional brain networks, Cerebral Cortex.

Code:
This package contains code for network construction. For graph theoretical analysis, the Brain Connectivity Toolbox (BCT) from the Sporns Lab or another toolbox with similar functions is needed. The paper uses BCT functions. The pipeline for using this code is below:

Step 1. 
Make sure you have the CanlabCore toolbox (from Github) and canlab_datasets_and_metadata.mat (included here), which provides needed functions and datasets.

Step 2. 
Run preparation_extract_data.m using canlab_datasets_and_metadata.mat as the input. This function will segment the brain using a prior template (the brain connectome atlas we privide in Template folder), and extract the regional activity of each brain region. The outputs are organized according to the datasets we used, and provide the meta information and regional data for each dataset.

Step 3.
 Run mean_data_cross_trails.m, using the outputs of step 2 as the inputs. Because the datasets we used are heterogeneous data acquired on different MRI machines, and include psychological manipulations (e.g., placebo treatment), this function can remove the placebo/high-vif (high Variance Inflation Factor) trials for relevant datasets, and extract the mean activity of each brain region across trials for each individual. Signals related to mean grey matter, white matter and CSF are then regressed out from the cross-trial mean activity of brain regions within each dataset. Median absolute deviation (MAD) is used to rescale the residuals of each dataset to decrease the difference between datasets.

Step 4. 
Run CRsubj_net.m using the output of Step 3 as the input. This function calculate the cross-subject network based on average trials of subjects under similar stimuli. We use 45.3 ¡æ as threshold in this work. In a subject has multiple stimulation level above (or below) the threshold, we randomly selected one stimulation intensity for network construction. 

Module:
The module folder provides the nifty-format files for each module for noxious, innocuous, painful, and non-painful. The module segmentation index in the template is also provided.

Network:
The network folder includes network matrices and regional data that we used in our work.

Template:
Nifit file for the Brainnetome brain connectome atlas including 273 brain regions and their coordinates and names


   


