"""Validate scientific assets and deployment integrity against source maps."""
import gzip,json,struct
from pathlib import Path
import nibabel as nib
import numpy as np
from build import ROOT,OUT,load_source,sign_stats
from scipy.ndimage import map_coordinates

def main():
 d=json.loads((OUT/'catalog.json').read_text());ids=set();n=0
 for s in d['studies']:
  assert s['id'] not in ids;ids.add(s['id']);assert s['paper'].startswith('https://'),s['name']
  assert (OUT/'marker'/s['id']/'index.html').exists()
  mids=set()
  for m in s['maps']:
   assert m['id'] not in mids;mids.add(m['id'])
   src=load_source(ROOT/m['source']);a=src.get_fdata();a=a[...,m['frame']] if m['frame'] is not None else a
   exported=nib.load(OUT/m['url']);b=exported.get_fdata()
   np.testing.assert_allclose(a,b,rtol=0,atol=0,equal_nan=True)
   np.testing.assert_allclose(src.affine,exported.affine,rtol=0,atol=1e-5)
   stats=sign_stats(a.astype(np.float32));assert stats==m['stats']
   for side,path in m['surfaces'].items():
    raw=gzip.decompress((OUT/path).read_bytes());magic,attr,faces,verts,skip=struct.unpack('<HHIII',raw[:16]);assert magic==23117 and attr==8 and faces==0 and verts==32492
    values=np.frombuffer(raw,dtype='<f4',offset=16);assert len(values)==verts and np.isfinite(values).all()
   assert (OUT/m['preview']).exists();n+=1
 # Independent interpolation check for one representative map (physical world-coordinate transform).
 m=d['studies'][1]['maps'][0];src=load_source(ROOT/m['source']);g=nib.load(OUT/d['assets']['shared']/'L.surf.gii');pts=g.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
 coords=nib.affines.apply_affine(np.linalg.inv(src.affine),pts)
 expected=map_coordinates(np.nan_to_num(src.get_fdata(dtype=np.float32)),coords.T,order=1,mode='constant',cval=0)
 actual=np.frombuffer(gzip.decompress((OUT/m['surfaces']['L']).read_bytes()),dtype='<f4',offset=16)
 np.testing.assert_allclose(actual,expected)
 print(f'PASS: {n} map exports match source values/affines; thresholds, links, surface files and projection verified.')
if __name__=='__main__':main()
