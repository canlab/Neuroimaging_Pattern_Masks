"""One-time catalog seed from repository documentation; review data/catalog.json after edits."""
import json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'Multivariate_signature_patterns'
# Folder prefix, display title, domains, modalities, targets. Map-level tags refine these.
STUDIES=[
('2011_Wager','Placebo analgesia','Pain','Somatosensory','Placebo response'),
('2015_Chang','PINES · Negative affect','Aversive / negative affect','Visual','Negative affect'),
('2015_Kragel','Emotion categories · BPLS','Aversive / negative affect|Appetitive / reward','Visual|Auditory','Emotion categories'),
('2015_Woo','Romantic rejection & physical pain','Cognitive & social|Pain','Visual|Somatosensory','Social rejection'),
('2016_Eisenbarth','Autonomic responses · GSR & HR','Physiology','Multimodal','Skin conductance|Heart rate|Social threat'),
('2016_Krishnan','VPS · Vicarious pain','Pain','Visual','Vicarious pain'),
('2017_Ashar','Empathic care & distress','Cognitive & social|Aversive / negative affect|Appetitive / reward','Auditory','Empathic care|Empathic distress'),
('2017_Woo','SIIPS1 · Cerebral contributions to pain','Pain','Somatosensory','Pain intensity|Expectancy|Perceived control'),
('2018_Kragel','Medial frontal cortex patterns','Pain|Aversive / negative affect|Cognitive & social','Multimodal','Pain|Negative affect|Cognitive control'),
('2018_Reddan','ImEx · Conditioned threat','Aversive / negative affect','Auditory','Fear|Threat conditioning|Imagined extinction'),
('2019_Kragel','Emotion schemas','Aversive / negative affect|Appetitive / reward','Visual','Emotion categories'),
('2019_Lee','Chronic back pain · Perfusion','Clinical|Pain','Somatosensory','Chronic back pain'),
('2019_Yu','Interpersonal guilt','Cognitive & social|Aversive / negative affect','Multimodal','Guilt'),
('2020_Geuter','Pain mediation · PDM','Pain','Somatosensory','Pain mediation'),
('2020_Silvestrini','Pain & cognitive control','Pain|Cognitive & social','Somatosensory|Visual','Pain|Cognitive control|Stroop interference'),
('2020_Zhou','General vicarious pain','Pain','Visual','Vicarious pain|Facial expressions'),
('2021_Ceko','MPA2 · Multi-aversive patterns','Aversive / negative affect|Pain','Visual|Auditory|Somatosensory','Negative affect|Mechanical pain|Thermal pain'),
('2021_Zhou','VIFS · Subjective fear','Aversive / negative affect','Visual','Fear'),
('2021_vantHoff','BASIC · Sexual-image classifier','Appetitive / reward','Visual','Sexual processing'),
('2022_Koban','NCS · Craving','Appetitive / reward','Visual','Food craving|Drug craving'),
('2022_coll','Decision value · Pain & money','Pain|Appetitive / reward','Somatosensory|Visual','Pain decision value|Monetary reward'),
('2023_Speer','BRS · Brain reward signature','Appetitive / reward','Visual','Monetary reward'),
('2024_FEPS','FEPS · Facial expression of pain','Pain','Somatosensory','Facial expression of pain'),
('2026_Acil','Mentalizing · Self & other','Cognitive & social','Visual','Self-related processing|Other-related processing|Mentalizing'),
('2026_Murillo','PiFoneM · Fear of neck movement','Clinical|Aversive / negative affect|Pain','Visual','Fear of pain|Neck movement'),
]
NAMES={
'Rating_Weights_LOSO_2':'PINES · full predictive weights',
'nonnoc_v11_4_137subjmap_weighted_mean':'SIIPS1 · full predictive weights',
'Mentalization_Boot_Unthr_11-Jun-2024':'MS · Mentalizing',
'Self_Boot_Unthr_11-Jun-2024':'MS-Self · Self-related mentalizing',
'Other_Boot_Unthr_11-Jun-2024':'MS-Other · Other-related mentalizing',
'SvO_Boot_Unthr_11-Jun-2024':'MS-SvO · Self vs. other',
'craving_wmapN99_boot10K_02-May-2022':'NCS · Drug and food craving',
'wmap_onlyFOOD_l2nGM_N99_20220428':'Food craving · weights',
'wmap_onlyDRUGS_l2nGM_N99_20220428':'Drug craving · weights',
'PiFoneM_unthresholded':'PiFoneM · full predictive weights',
'Reward_Signature_bootstrapped_0.5':'BRS · supplied bootstrapped map',
'Geuter_2020_cPDM_combined_pain_map':'cPDM · combined pain mediation',
'feps_mean_xval_weights':'FEPS · mean cross-validated weights',
'feps_z_unthresholded':'FEPS · standardized unthresholded map',
'VIFS':'VIFS · Subjective fear',
}
def plain(s):
 s=re.sub(r'\[([^\]]+)\]\([^)]*\)',r'\1',s)
 return re.sub(r'\s+',' ',s.replace('**','').replace('*','').replace('`','')).strip()
def stem(p):
 return re.sub(r'\.(nii|img)(\.gz)?$','',p.name)
items=[]
for prefix,title,domains,modalities,targets in STUDIES:
 folder=next(BASE.glob(prefix+'*')); text=(folder/'contents_description.md').read_text()
 overview=text.split('## Overview',1)[1].split('**Primary reference',1)[0].strip()
 ref=re.search(r'\*\*Primary reference.*?\*\*\s*(.*?)(?:\n\[doi|\n\n)',text,re.S)
 doi=re.search(r'https://doi.org/([^\s)]+)',text)
 citation=plain(ref.group(1)) if ref else title
 year=re.search(r'\((20\d\d)\)',citation)
 files={}
 for p in sorted(folder.rglob('*')):
  if not p.name.endswith(('.nii','.nii.gz','.img','.img.gz')): continue
  if prefix=='2019_Lee' and 'PCASL' not in str(p): continue
  if stem(p) in ('mask',) or any(s in p.name for s in ['pvalue','subcluster_maps','dACC_rf','dACC_hw']): continue
  if p.name.endswith('.gz') and p.with_suffix('').exists(): continue
  files[str(p.relative_to(folder)).removesuffix('.gz')]=p
 maps=[]
 for p in files.values():
  s=stem(p); low=s.lower()
  support=any(k in low for k in ['fdr','thr','bootz','bsz','importance','pruned','p005']) and not any(k in low for k in ['unthr','nothresh'])
  role='Statistical support' if support else 'Predictive weights'
  if prefix=='2015_Kragel': role='Bootstrap z map'
  if prefix=='2023_Speer': role='Bootstrapped map'
  name=NAMES.get(s,s.replace('PLS_betas_','Emotion schema · ').replace('bPLS_','MFC · ').replace('_',' '))
  mt=targets.split('|')
  if prefix=='2019_Kragel': mt=[s.replace('PLS_betas_','')]
  if prefix=='2015_Kragel': mt=[s.split('_')[2].capitalize()]
  if prefix=='2026_Acil': mt={'Self': ['Self-related processing'],'Other':['Other-related processing'],'SvO':['Self-related processing','Other-related processing']}.get(s.split('_')[0],['Mentalizing'])
  if prefix=='2021_Ceko': mt=[s.split('_')[0]+' aversion']
  if prefix=='2022_Koban' and 'onlyFOOD' in s: mt=['Food craving']
  if prefix=='2022_Koban' and 'onlyDRUGS' in s: mt=['Drug craving']
  maps.append(dict(id=re.sub('[^a-z0-9]+','-',s.lower()).strip('-'),name=name,source=str(p.relative_to(ROOT)),role=role,targets=mt,description=name+'. '+('Supplied statistical/support map; use the full predictive weights for pattern expression.' if support or role!='Predictive weights' else 'Unthresholded map as supplied in the repository. Display thresholds do not alter the downloaded data.')))
 maps.sort(key=lambda m:(m['role']!='Predictive weights', 'full predictive' not in m['name'],m['name']))
 item=dict(id=folder.name.lower().replace('_','-'),folder=str(folder.relative_to(ROOT)),name=title,year=int(year.group(1)) if year else int(prefix[:4]),description=plain(overview.split('\n\n')[0]),citation=citation,paper='https://doi.org/'+doi.group(1) if doi else '',domains=domains.split('|'),modalities=modalities.split('|'),targets=targets.split('|'),space='MNI coordinates; exact template variant not documented',maps=maps)
 if prefix=='2016_Eisenbarth': item['description']='Brain patterns predicting heart rate and skin conductance responses to social threat.'
 if prefix=='2026_Acil': item.update(paper='https://www.nature.com/articles/s41467-026-73945-w',citation='Açıl, D., et al. (2026). Brain neuromarkers predict self- and other-related mentalizing across adult, clinical, and developmental samples. Nature Communications, 17, 7229.')
 if prefix=='2024_FEPS': item.update(paper='https://elifesciences.org/articles/87962',citation='Picard, F., et al. (2023). A distributed brain response predicting the facial expression of acute nociceptive pain. eLife, 12, RP87962.',year=2023)
 items.append(item)
(ROOT/'site/data/catalog.json').write_text(json.dumps({'version':1,'studies':items},indent=2,ensure_ascii=False)+'\n')
print(len(items),'studies;',sum(len(s['maps']) for s in items),'map files')
