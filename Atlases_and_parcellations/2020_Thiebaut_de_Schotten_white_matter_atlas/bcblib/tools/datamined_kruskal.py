
import numpy as np
import pandas as pd
import scipy.stats as st



def nii_suffix(string):
    if string.endswith(".nii.gz"):
        return string
    elif string.endswith(".nii"):
        return string + '.gz'
    else:
        return string + '.nii.gz'

def del_suffix(string):
    if string.endswith(".nii.gz"):
        return string[:-7]
    elif string.endswith(".nii"):
        return string[:-4]
    else:
        return string

# controls = pd.read_csv("/data/Chris/res_anacom_disco_mw/copycontrols.csv",
#                        header=None)
#
# patients = pd.read_csv("/data/Chris/res_anacom_disco_mw/copypatients.csv",
#                        header=None)

pat = {'patient1.nii.gz': 32.0,
 'patient2.nii.gz': 32.0,
 'patient3.nii.gz': 27.0,
 'patient4.nii.gz': 26.0,
 'patient5.nii.gz': 28.0,
 'patient6.nii.gz': 40.0,
 'patient7.nii.gz': 22.0,
 'patient8.nii.gz': 21.0,
 'patient9.nii.gz': 37.0,
 'patient10.nii.gz': 20.0,
 'patient11.nii.gz': 28.0,
 'patient12.nii.gz': 20.0,
 'patient13.nii.gz': 21.0,
 'patient14.nii.gz': 16.0,
 'patient15.nii.gz': 31.0,
 'patient16.nii.gz': 38.0,
 'patient17.nii.gz': 29.0,
 'patient18.nii.gz': 33.0,
 'patient19.nii.gz': 33.0,
 'patient20.nii.gz': 33.0,
 'patient21.nii.gz': 31.0,
 'patient22.nii.gz': 29.0,
 'patient23.nii.gz': 23.0,
 'patient24.nii.gz': 19.0,
 'patient25.nii.gz': 38.0,
 'patient26.nii.gz': 10.0,
 'patient27.nii.gz': 27.0,
 'patient29.nii.gz': 42.0,
 'patient30.nii.gz': 33.0,
 'patient31.nii.gz': 26.0,
 'patient32.nii.gz': 43.0,
 'patient33.nii.gz': 35.0,
 'patient34.nii.gz': 38.0,
 'patient35.nii.gz': 20.0,
 'patient36.nii.gz': 45.0,
 'patient37.nii.gz': 29.0,
 'patient38.nii.gz': 34.0}

# pat
# for i in np.arange(len(patients.index)):
#         pat[patients[0][i]] = patients[1][i]

ctr = [37,47,35,21,38,37,62,46,57,36,23,39,34,45,39,41,30,42,34,45,37,35,38,36,
       29,41,33,47,23,47,29,38,48,29,36,40,37,38,33,41,45,29,32,50,44,33,32,32,
       53,33,43,59,26,25]

disco_pat = input("Disconnected patients: ")
clu_name = input("Cluster name: ")
disco = disco_pat.split(',')


disco = [nii_suffix(p) for p in disco]

dis_scores = [pat[p] for p in disco]
spared = [p for p in pat.keys() if p not in disco]
spa_sco = [pat[p] for p in pat.keys() if p not in disco]
print(disco)
print(spared)

h, p = st.kruskal(dis_scores, spa_sco, ctr)


dct_h, dct_p = st.mannwhitneyu(dis_scores, ctr)
cct_h, cct_p = st.mannwhitneyu(spa_sco, ctr)
cd_h, cd_p = st.mannwhitneyu(spa_sco, dis_scores)

bonf_corr = 276
print("Cluster: " + clu_name)
print("Pvalue: " + str(p))
print("Stat (H): " + str(h))
print("Pval bonf corr: " + str(p * bonf_corr))

f = open("/data/Chris/BCBTK_paper/kw_clu.csv", 'a')
s = str(clu_name + ',' + str(p) + ',' + str(h) + ',' + str(p*bonf_corr) + '\n')
f.write(s)
f.close()
