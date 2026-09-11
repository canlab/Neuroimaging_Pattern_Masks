"""Reproducible static gallery: no MATLAB, network, or server runtime required."""
from pathlib import Path
import copy,gzip,hashlib,html,json,os,re,shutil,struct,tempfile,zlib
import nibabel as nib
import numpy as np
from scipy.ndimage import map_coordinates
ROOT=Path(__file__).resolve().parents[2]; SITE=ROOT/'site'; OUT=SITE/'dist'
def load_source(source):
 if source.name.endswith('.img.gz'):
  with tempfile.TemporaryDirectory() as temp:
   local=Path(temp)/source.name[:-3]; local.write_bytes(gzip.decompress(source.read_bytes()))
   hdr=source.with_suffix('').with_suffix('.hdr')
   if hdr.exists(): shutil.copyfile(hdr,local.with_suffix('.hdr'))
   else: local.with_suffix('.hdr').write_bytes(gzip.decompress(Path(str(hdr)+'.gz').read_bytes()))
   loaded=nib.load(local)
   # Materialize before temporary source disappears; preserve affine from Analyze origin.
   return nib.Nifti1Image(loaded.get_fdata(),loaded.affine)
 return nib.load(source)
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write_json(p,data): p.write_text(json.dumps(data,ensure_ascii=False,separators=(',',':'),allow_nan=False))
def sign_stats(data):
 result={}
 for name,v in [('positive',data[data>0]),('negative',-data[data<0])]:
  v=v[np.isfinite(v)]
  result[name]={'count':int(v.size),'max':float(v.max()) if v.size else 0,'cutoffs':[float(np.quantile(v,1-p/100)) if v.size else 0 for p in range(101)]}
 return result

def save_preview(rgb,path):
 """Write a deterministic RGB PNG without an image-decoding dependency."""
 height,width=rgb.shape[:2]; scale=min(340/width,240/height)
 h,w=max(1,int(height*scale)),max(1,int(width*scale))
 ys=np.minimum((np.arange(h)/scale).astype(int),height-1)
 xs=np.minimum((np.arange(w)/scale).astype(int),width-1)
 canvas=np.zeros((260,360,3),dtype=np.uint8)
 canvas[(260-h)//2:(260-h)//2+h,(360-w)//2:(360-w)//2+w]=rgb[ys[:,None],xs[None,:]]
 def chunk(kind,data): return struct.pack('>I',len(data))+kind+data+struct.pack('>I',zlib.crc32(kind+data)&0xffffffff)
 raw=b''.join(b'\0'+row.tobytes() for row in canvas)
 path.write_bytes(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',360,260,8,2,0,0,0))+chunk(b'IDAT',zlib.compress(raw,9))+chunk(b'IEND',b''))

def main():
 # Remove only generated output, avoiding stale assets and stale catalog routes.
 if OUT.exists(): shutil.rmtree(OUT)
 OUT.mkdir(exist_ok=True)
 for folder in ['assets','maps','previews','marker']: (OUT/folder).mkdir(exist_ok=True)
 for p in (SITE/'src').iterdir():
  if p.is_file(): shutil.copyfile(p,OUT/p.name)
 version=hashlib.sha256(b''.join(p.read_bytes() for p in sorted((SITE/'vendor').iterdir()))).hexdigest()[:12]
 vendor=OUT/'assets'/version; shutil.copytree(SITE/'vendor',vendor,dirs_exist_ok=True)
 shared=SITE/'data/shared'; shash=digest(shared/'provenance.json')[:12]
 dest=OUT/'assets'/shash; shutil.copytree(shared,dest,dirs_exist_ok=True)
 assets={'viewer':f'assets/{version}/canlab_niivue_viewer.js','niivue':f'assets/{version}/niivue.js','css':f'assets/{version}/canlab_niivue.css','shared':f'assets/{shash}/'}
 # Sample anatomy on each map grid for faithful slice previews.
 anatomy=nib.load(shared/'underlay.nii.gz'); anatomical=anatomy.get_fdata(dtype=np.float32)
 meshes={}
 for side in ['L','R']:
  g=nib.load(shared/f'{side}.surf.gii'); meshes[side]=g.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
 catalog=json.loads((SITE/'data/catalog.json').read_text()); emitted=[]; source_paths=set()
 for study in catalog['studies']:
  s=copy.deepcopy(study); s['maps']=[]
  for entry in study['maps']:
   source=ROOT/entry['source']; source_paths.add(entry['source'])
   img=load_source(source); data=img.get_fdata(dtype=np.float32)
   if data.ndim not in [3,4]: raise ValueError(f'Unsupported dimensions: {source}: {data.shape}')
   if not np.isfinite(img.affine).all() or abs(np.linalg.det(img.affine[:3,:3]))<1e-8: raise ValueError(f'Bad affine: {source}')
   for frame in range(data.shape[3] if data.ndim==4 else 1):
    m=copy.deepcopy(entry); arr=data[...,frame] if data.ndim==4 else data
    if not np.any(np.isfinite(arr)&(arr!=0)):
     if m['role']=='Predictive weights': raise ValueError(f'Empty predictive map: {source} frame {frame}')
     m['description']+=' This supplied support component has no surviving nonzero voxels.'
    if data.ndim==4: m['id']+=f'-frame-{frame+1}'; m['name']+=f' · component {frame+1}'
    m['frame']=frame if data.ndim==4 else None
    key=hashlib.sha256(source.read_bytes()+img.affine.tobytes()+img.header.binaryblock+str(frame).encode()+b'v3').hexdigest()[:20]
    url=f'maps/{key}.nii.gz'; target=OUT/url
    if not target.exists():
     # Only container/frame conversion; no reorientation, resampling, normalization, or display threshold.
     if data.ndim==3 and source.name.endswith('.nii.gz'): shutil.copyfile(source,target)
     elif data.ndim==3 and source.name.endswith('.nii'): target.write_bytes(gzip.compress(source.read_bytes(),mtime=0))
     else:
      full=img.get_fdata(); full=full[...,frame] if data.ndim==4 else full
      outimg=nib.Nifti1Image(full,img.affine); outimg.header.set_xyzt_units('mm'); outimg.set_sform(img.affine,4)
      nib.save(outimg,target)
    m.update(url=url,sha256=digest(target),source_sha256=digest(source),shape=list(arr.shape),voxel_size=[float(v) for v in nib.affines.voxel_sizes(img.affine)],orientation=''.join(nib.aff2axcodes(img.affine)),bytes=target.stat().st_size,stats=sign_stats(arr))
    # Trilinear midthickness sampling in world coordinates. This is a display derivative,
    # not a surface-based predictive model or a template registration.
    m['surfaces']={}
    for side,points in meshes.items():
     surl=f'maps/{key}-{side}.mz3'; starget=OUT/surl
     if not starget.exists():
      coords=nib.affines.apply_affine(np.linalg.inv(img.affine),points)
      values=map_coordinates(np.nan_to_num(arr),coords.T,order=1,mode='constant',cval=0).astype('<f4')
      raw=struct.pack('<HHIII',23117,8,0,len(values),0)+values.tobytes()
      starget.write_bytes(gzip.compress(raw,mtime=0))
     m['surfaces'][side]=surl
    # One orthographic axial preview: actual map, shared anatomy, per-sign 20%.
    preview=f'previews/{key}.png'; pp=OUT/preview
    if not pp.exists():
     ras=nib.as_closest_canonical(nib.Nifti1Image(arr,img.affine)); v=ras.get_fdata(dtype=np.float32)
     z=int(np.argmax(np.sum(np.abs(np.nan_to_num(v)),axis=(0,1))))
     grid=np.indices(v.shape[:2]); ijk=np.stack([grid[0].ravel(),grid[1].ravel(),np.full(grid[0].size,z)],axis=1)
     coords=nib.affines.apply_affine(np.linalg.inv(anatomy.affine)@ras.affine,ijk)
     bg=map_coordinates(anatomical,coords.T,order=1,mode='constant').reshape(v.shape[:2])
     bg=np.clip(bg/(np.percentile(anatomical[anatomical>0],99) or 1),0,1)
     rgb=np.repeat((bg*165)[...,None],3,axis=2); sl=v[:,:,z]
     pos=(sl>0)&(sl>=m['stats']['positive']['cutoffs'][20]); neg=(sl<0)&(-sl>=m['stats']['negative']['cutoffs'][20])
     rgb[pos]=[245,134,50]; rgb[neg]=[62,169,234]
     save_preview(np.rot90(rgb.astype(np.uint8)),pp)
    m['preview']=preview; s['maps'].append(m)
   print('.',end='',flush=True)
  emitted.append(s)
 catalog['studies']=emitted; catalog['assets']=assets
 catalog['notes']=['Connectivity-edge models excluded.','Fibromyalgia redistribution status requires reconciliation with repository README; not packaged.','Matthewson: no local maps. Van Oudenhove: network coefficients only; not packaged.']
 write_json(OUT/'catalog.json',catalog)
 # Real static detail routes supply metadata to social crawlers without JavaScript.
 template=(OUT/'index.html').read_text(); origin=os.environ.get('GALLERY_ORIGIN','https://canlab.github.io/Neuroimaging_Pattern_Masks').rstrip('/')
 for s in emitted:
  p=OUT/'marker'/s['id']; p.mkdir(exist_ok=True)
  markup=template.replace('<head>','<head>\n<base href="../../">').replace('<title>Neuromarker Gallery · CANlab</title>',f'<title>{html.escape(s["name"])} · CANlab</title>')
  meta=f'<meta property="og:title" content="{html.escape(s["name"],quote=True)}"><meta property="og:description" content="{html.escape(s["description"],quote=True)}"><meta property="og:image" content="{origin}/{s["maps"][0]["preview"]}"><meta property="og:url" content="{origin}/marker/{s["id"]}/"><meta name="twitter:card" content="summary_large_image">'
  (p/'index.html').write_text(markup.replace('</head>',meta+'</head>'))
 (OUT/'.nojekyll').touch()
 report={'studies':len(emitted),'maps':sum(len(s['maps']) for s in emitted),'unique_sources':len(source_paths),'total_bytes':sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())}
 write_json(OUT/'build-report.json',report)
 print('\n'+json.dumps(report,indent=2))
 if report['total_bytes']>900_000_000: raise ValueError('Site exceeds 900 MB deployment budget')
if __name__=='__main__': main()
