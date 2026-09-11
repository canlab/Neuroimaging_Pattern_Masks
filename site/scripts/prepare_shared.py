"""Vendor shared anatomy/atlas and HCP meshes from a specified CANlab checkout."""
import ast,base64,json,re,sys,shutil,hashlib,subprocess
from pathlib import Path
import nibabel as nib
HERE=Path(__file__).resolve().parents[1]
core=Path(sys.argv[1]); viewer=core/'CanlabCore/Visualization_functions/canlab_niivue'
out=HERE/'data/shared'; out.mkdir(exist_ok=True)
text=(viewer/'sample/emotionreg_ttest.html').read_text()
for key,name in [('UNDERLAY','underlay'),('ATLAS','atlas')]:
 data=base64.b64decode(re.search(r"const "+key+r"_B64 = '([^']+)'",text).group(1))
 (out/(name+'.nii.gz')).write_bytes(data)
labels=ast.literal_eval(re.search(r'atlasLabels: (\[.*?\])',text).group(1))
(out/'atlas-labels.json').write_text(json.dumps(labels))
surfs=core/'CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces'
for side in ['L','R']:
 shutil.copyfile(surfs/f'S1200.{side}.midthickness_MSMAll.32k_fs_LR.surf.gii',out/f'{side}.surf.gii')
info={'canlab_commit':subprocess.check_output(['git','-C',str(core),'rev-parse','HEAD'],text=True).strip(),'atlas':'CANlab2024 (exported CANlab demo label volume)','anatomy':'CANlab viewer demo MNI underlay','surface':'HCP S1200 MSMAll fsLR 32k midthickness; nominal MNI coordinates','files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir()) if p.name!='provenance.json'}}
(out/'provenance.json').write_text(json.dumps(info,indent=2)+'\n')
print(info['canlab_commit'],len(labels),'atlas labels')
