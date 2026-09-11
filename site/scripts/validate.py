"""Validate scientific assets and deployment integrity against source maps."""
import csv,gzip,hashlib,json,struct,zlib
from pathlib import Path
import nibabel as nib
import numpy as np
from build import ROOT,OUT,load_source,sign_stats
from scipy.ndimage import map_coordinates

def main():
 d=json.loads((OUT/'catalog.json').read_text());ids=set();n=0
 shared=OUT/d['assets']['shared']; labels=json.loads((shared/'atlas-labels.json').read_text()); descriptions=json.loads((shared/'atlas-descriptions.json').read_text())
 assert len(descriptions)==len(labels)==518
 assert [r['code'] for r in descriptions]==list(range(1,519))
 assert [r['short_name'] for r in descriptions]==labels
 assert all(r['full_name'].strip() and r['source'].startswith('https://') for r in descriptions)
 assert next(r['full_name'] for r in descriptions if r['short_name']=='Ctx_p24_L')=='Cortex: Area posterior 24, Left'
 assert len(list(csv.DictReader((shared/'atlas-descriptions.csv').open())))==518
 for name,checksum in json.loads((shared/'provenance.json').read_text())['files'].items(): assert hashlib.sha256((shared/name).read_bytes()).hexdigest()==checksum,name
 assert set(np.unique(nib.load(shared/'atlas.nii.gz').get_fdata())).issubset(set(range(519)))
 brainmask=nib.load(shared/'preview-brainmask.nii.gz');brain=brainmask.get_fdata(dtype=np.float32)
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
   # Independently inspect exported PNG pixels: no colored overlay may escape the
   # nearest-neighbor template brain mask, including padding and raster scaling.
   ras=nib.as_closest_canonical(nib.Nifti1Image(a.astype(np.float32),src.affine));v=ras.get_fdata(dtype=np.float32)
   coords=nib.affines.apply_affine(np.linalg.inv(brainmask.affine)@ras.affine,np.indices(v.shape).reshape(3,-1).T)
   inside=map_coordinates(brain,coords.T,order=0,mode='constant',cval=0).reshape(v.shape)>0
   z=int(np.argmax(np.sum(np.abs(np.nan_to_num(v))*inside,axis=(0,1))))
   mask=np.rot90(inside[:,:,z]);height,width=mask.shape;scale=min(340/width,240/height);h,w=max(1,int(height*scale)),max(1,int(width*scale))
   ys=np.minimum((np.arange(h)/scale).astype(int),height-1);xs=np.minimum((np.arange(w)/scale).astype(int),width-1)
   allowed=np.zeros((260,360),bool);allowed[(260-h)//2:(260-h)//2+h,(360-w)//2:(360-w)//2+w]=mask[ys[:,None],xs[None,:]]
   png=(OUT/m['preview']).read_bytes();offset=8;compressed=b''
   while offset<len(png):
    size=struct.unpack('>I',png[offset:offset+4])[0];kind=png[offset+4:offset+8]
    if kind==b'IDAT':compressed+=png[offset+8:offset+8+size]
    offset+=12+size
   scanlines=np.frombuffer(zlib.decompress(compressed),dtype=np.uint8).reshape(260,1081);assert np.all(scanlines[:,0]==0)
   rgb=scanlines[:,1:].reshape(260,360,3);colored=(rgb[:,:,0]!=rgb[:,:,1])|(rgb[:,:,1]!=rgb[:,:,2])
   assert not np.any(colored&~allowed),m['id']
 # Independent interpolation check for one representative map (physical world-coordinate transform).
 m=d['studies'][1]['maps'][0];src=load_source(ROOT/m['source']);g=nib.load(OUT/d['assets']['shared']/'L.surf.gii');pts=g.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
 coords=nib.affines.apply_affine(np.linalg.inv(src.affine),pts)
 expected=map_coordinates(np.nan_to_num(src.get_fdata(dtype=np.float32)),coords.T,order=1,mode='constant',cval=0)
 actual=np.frombuffer(gzip.decompress((OUT/m['surfaces']['L']).read_bytes()),dtype='<f4',offset=16)
 np.testing.assert_allclose(actual,expected)
 print(f'PASS: {n} map exports match source values/affines; thresholds, links, surface projection, 518 atlas names, and all thumbnail brain masks verified.')
if __name__=='__main__':main()
