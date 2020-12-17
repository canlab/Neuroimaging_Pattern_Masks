import sys
import os

import numpy as np
from statsmodels.stats.multitest import multipletests
import scipy.stats as st
import pandas as pd

def nii_suffix(string):
    if string.endswith(".nii.gz"):
        return string
    elif string.endswith(".nii"):
        return string + '.gz'
    else:
        return string + '.nii.gz'

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

log_file = sys.argv[1]
res_folder = sys.argv[2]

with open(log_file, "r") as text_file:
    lines = text_file.readlines()
indices = [ind for ind,l in enumerate(lines) if l.startswith("+(335)")]
clu_pat_list = [l[14:-1] for ind, l in enumerate(lines) if ind in indices]

clu_names = []
for i in indices:
    tmp_i = i - 1
    line = lines[tmp_i]
    while not line.startswith('+(322)'):
        tmp_i -= 1
        line = lines[tmp_i]
    clu_names.append(os.path.basename(line.split(' ')[1]))
# clu_names =  [l.split(' ')[1] for ind, l in enumerate(
#     lines) if (ind + 7) in indices]
# print(clu_names)

    # print(filtered)
disco = [[nii_suffix(p) for p in d.split(',')] for d in clu_pat_list]
dis_scores = [[pat[p] for p in clu] for clu in disco]
spared = [[p for p in pat.keys() if p not in clu] for clu in disco]
spa_sco = [[pat[p] for p in pat.keys() if p not in clu] for clu in disco]

kruskal = []
for ind, dis in enumerate(dis_scores):
    h, p = st.kruskal(dis_scores[ind], spa_sco[ind], ctr)
    spa_deco_h, spa_deco_p = st.mannwhitneyu(
        spa_sco[ind], dis_scores[ind], alternative='two-sided')
    deco_ctr_h, deco_ctr_p = st.mannwhitneyu(
        dis_scores[ind], ctr, alternative='two-sided')
    spa_ctr_h, spa_ctr_p = st.mannwhitneyu(
        spa_sco[ind], ctr, alternative='two-sided')
    kruskal.append([h, p, spa_deco_h, spa_deco_p, deco_ctr_h, deco_ctr_p,
                    spa_ctr_h, spa_ctr_p])

kw_pvals = [clu[1] for clu in kruskal]
spa_deco_pvals = [clu[3] for clu in kruskal]
deco_ctr_pvals = [clu[5] for clu in kruskal]
spa_ctr_pvals = [clu[7] for clu in kruskal]
# test = np.array([0.0003, 0.000023,0.004,0.026, 0.03])
reject, kw_corr, alphacSidak, alphacBonf = multipletests(
    pvals=kw_pvals, alpha=0.05, method='holm')
reject, spa_deco_corr, alphacSidak, alphacBonf = multipletests(
    pvals=spa_deco_pvals, alpha=0.05, method='holm')
reject, deco_ctr_corr, alphacSidak, alphacBonf = multipletests(
    pvals=deco_ctr_pvals, alpha=0.05, method='holm')
reject, spa_ctr_corr, alphacSidak, alphacBonf = multipletests(
    pvals=spa_ctr_pvals, alpha=0.05, method='holm')

for ind,v in enumerate(kw_corr):
    kruskal[ind].append(v)
    kruskal[ind].append(spa_deco_corr[ind])
    kruskal[ind].append(deco_ctr_corr[ind])
    kruskal[ind].append(spa_ctr_corr[ind])
# for ind,v in enumerate(kw_corr):
#     kruskal[ind] = np.insert(kruskal[ind], 2, v)
#     kruskal[ind] = np.insert(kruskal[ind], 5, spa_deco_corr[ind])
#     kruskal[ind] = np.insert(kruskal[ind], 8, deco_ctr_pvals[ind])
#     kruskal[ind] = np.insert(kruskal[ind], -1, spa_ctr_corr[ind])
# print(kruskal)
print(deco_ctr_pvals[268])
final_dict = {}

for ind, clu in enumerate(kruskal):
    final_dict[clu_names[ind]] = clu
col = ['kw_st','kw_p','spa_deco_st', 'spa_deco_p',
       'deco_ctr_st', 'deco_ctr_p', 'spa_ctr_st','spa_ctr_p','kw_p_corr',
       'spa_deco_p_corr', 'deco_ctr_p_corr', 'spa_ctr_p_corr']
# out_res = pd.DataFrame.from_dict(final_dict)
out_res = pd.DataFrame(data=kruskal, index=clu_names, columns=col)
# out_res.rename()

writer = pd.ExcelWriter(os.path.join(res_folder, 'kw_ph_old_anacom.xlsx'))
# out_test.to_excel(writer, 'tests')
out_res.to_excel(writer, 'FluA_kw_ph', columns=['kw_st','kw_p','kw_p_corr',
                                                'spa_deco_st', 'spa_deco_p',
                                                'spa_deco_p_corr','deco_ctr_st',
                                                'deco_ctr_p','deco_ctr_p_corr',
                                                'spa_ctr_st','spa_ctr_p',
                                                'spa_ctr_p_corr'])
writer.save()

# print(final_dict)
# print(np.array([clu[2] for clu in kruskal]) < 0.05)
