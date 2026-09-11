"""Regenerate the descriptive lookup from repository atlas dictionaries (no network).
The generated JSON/CSV are committed; deployment does not require source atlases.
"""
from pathlib import Path
import csv,json,re,hashlib
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'site/data/shared'
BASE='https://github.com/canlab/Neuroimaging_Pattern_Masks/blob/master/'
CAN='Atlases_and_parcellations/2024_CANLab_atlas/src/openCANLab2024_MNI152NLin6Asym_labels.csv'
BIO='Atlases_and_parcellations/2023_Bianciardi_BrainstemNavigatorV0.9/source_files/bianciardi_coarse_labels.csv'
GLASSER='https://www.humanconnectome.org/storage/app/media/documentation/AABC2/AreaNamesAndIndices%20-%20NIHMS68870-supplement-Neuroanatomical_Supplementary_Results.pdf'
rows=list(csv.DictReader((ROOT/CAN).open())); brain={}
for r in csv.DictReader((ROOT/BIO).open()):
 for side in ['left','right']:
  if r[side]: brain[re.sub(r'_([lr])$',lambda m:'_'+m[1].upper(),r[side].strip())]=r['full_label'].replace('prabigeminal','parabigeminal')
def clean(s):
 s=s.replace('_',' ').replace('NucleusAccumbens','Nucleus accumbens').replace('triansition','transition').replace('Anteroir','Anterior').replace('cudal','caudal')
 return re.sub(r'\s*\((left|right)\)|^(Left|Right)\s+', '', s, flags=re.I).strip()
labels=json.loads((OUT/'atlas-labels.json').read_text()); result=[]
for code,label in enumerate(labels,1):
 side={'L':'Left','R':'Right'}.get(label.split('_')[-1]); suffix=', '+side if side else ''
 matching=[r for r in rows if r['labels_2']==label] or [r for r in rows if r['labels']==label]
 components=list(dict.fromkeys(clean(r['label_descriptions']) for r in matching))
 source=BASE+CAN; note=''; name='; '.join(components)
 if label.startswith('Ctx_'):
  name='Cortex: '+name; source=GLASSER+' ; '+source
 elif label.startswith('Cblm_'):
  part=label.removeprefix('Cblm_'); part=re.sub(r'_[LR]$','',part)
  name='Cerebellum: '+part.replace('Vermis_','Vermis, ').replace('CrusII','Crus II').replace('CrusI','Crus I').replace('I_IV','lobules I–IV')
  if 'Crus' not in name and 'lobules' not in name: name=name.replace('Cerebellum: ','Cerebellum: lobule ')
  source='https://www.diedrichsenlab.org/imaging/propatlas.htm ; '+source
 elif re.match(r'BG_(CAU|PUT)_',label):
  _,structure,part,_=label.split('_'); name='Basal ganglia: '+{'CAU':'Caudate','PUT':'Putamen'}[structure]+', '+{'DA':'dorsal anterior','VA':'ventral anterior','DP':'dorsal posterior','VP':'ventral posterior','body':'body','tail':'tail'}[part]
  source='https://github.com/yetianmed/subcortex ; '+source
 elif label.startswith('BStem_'):
  key=label[6:]; plain=re.sub(r'_[LR]$','',key)
  overrides={'LC+':'Locus coeruleus and subcoeruleus','RObPaMg':'Raphe obscurus, pallidus and magnus','OC':'Inferior olivary nucleus and superior olivary complex','PAG':'Periaqueductal gray and merged cuneiform nucleus','STH':'Subthalamic nucleus'}
  if plain in overrides:
   name=overrides[plain]; source=BASE+'Atlases_and_parcellations/2024_CANLab_atlas/create_CANLab2024_atlas.m'; note='Composite atlas grouping; see construction code.' if plain!='STH' else 'Hemisphere follows the displayed label; source fine/coarse naming is inconsistent.'
  elif key in brain: name=brain[key]; source=BASE+BIO
  if key.startswith('Shen_') and side: name=re.sub(r'\b'+side+r'\b','',name,flags=re.I).strip()
  name='Brainstem: '+name
 elif label.startswith('BG_'): name='Basal ganglia: '+name
 elif label.startswith('MTL_'): name='Medial temporal lobe: '+name
 elif label.startswith('Thal_'): name='Thalamus: '+name
 elif label.startswith('hypothalamus_'): name='Hypothalamus: '+re.sub(r'_[LR]$','',label[len('hypothalamus_'):]).replace('_',' ')
 if len(matching)>1: note=(note+' Coarse parcel merges the listed source components.').strip()
 assert name and not name.endswith(': '),label
 result.append(dict(code=code,short_name=label,full_name=name+suffix,components=components,source=source,note=note))
(OUT/'atlas-descriptions.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
with (OUT/'atlas-descriptions.csv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=list(result[0]),lineterminator='\n');writer.writeheader()
 for r in result: writer.writerow({**r,'components':'; '.join(r['components'])})
assert len(result)==518 and result[0]['code']==1
print(f'Wrote {len(result)} descriptive labels')

provenance_path=OUT/'provenance.json'
provenance=json.loads(provenance_path.read_text())
for filename in ['atlas-descriptions.json','atlas-descriptions.csv']:
 provenance['files'][filename]=hashlib.sha256((OUT/filename).read_bytes()).hexdigest()
provenance_path.write_text(json.dumps(provenance,indent=2)+'\n')
