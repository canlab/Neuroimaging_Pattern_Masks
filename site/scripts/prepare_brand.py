"""Offline scientific rendering of CANlab's left cutaway with NPS weights.
Usage: python prepare_brand.py CANONICAL_CUTAWAY.mat NPS_WEIGHTS.img.gz
Requires matplotlib in addition to the gallery build environment. Generated PNGs
are committed; deployment does not run this preparation step.
"""
import sys,json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.ndimage import map_coordinates
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from build import load_source
OUT=Path(__file__).resolve().parents[1]/'src/brand';OUT.mkdir(exist_ok=True)
anatomy=loadmat(sys.argv[1],simplify_cells=True);img=load_source(Path(sys.argv[2]));values=img.get_fdata(dtype=np.float32)
cut={sign:np.quantile(np.abs(values[values*sign>0]),.65) for sign in [1,-1]}
triangles=[];base=[];paint=[]
for kind in ['isosurf','isocap']:
 for mesh in anatomy[kind]:
  vertices=mesh['vertices'];faces=mesh['faces'].astype(int)-1
  # Spatial vertex clustering reduces only the display mesh, never the model.
  cells=np.round(vertices/1.4).astype(int);_,inv=np.unique(cells,axis=0,return_inverse=True);n=inv.max()+1
  counts=np.bincount(inv);v=np.stack([np.bincount(inv,weights=vertices[:,i])/counts for i in range(3)],axis=1)
  f=inv[faces];f=f[(f[:,0]!=f[:,1])&(f[:,1]!=f[:,2])&(f[:,0]!=f[:,2])];f=np.unique(f,axis=0)
  tri=v[f];centers=tri.mean(axis=1);normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-8)
  light=np.array([.4,.2,1]);light/=np.linalg.norm(light);shade=.56+.40*np.abs(normal@light)
  gray=np.repeat(shade[:,None],3,axis=1)*np.array([.91,.94,.96])
  if kind=='isocap':
   intensity=np.bincount(inv,weights=mesh['facevertexcdata'])/counts
   gray*=np.clip((intensity[f].mean(axis=1)-100)/140,.35,1)[:,None]
  coordinates=np.linalg.inv(img.affine)@np.c_[centers,np.ones(len(centers))].T
  weight=map_coordinates(values,coordinates[:3],order=0,mode='constant',cval=0)
  colored=gray.copy()
  for sign in [1,-1]:
   selected=weight*sign>=cut[sign];mag=np.clip(np.abs(weight[selected])/np.quantile(np.abs(values[values*sign>0]),.99),0,1)
   cmap=plt.get_cmap('inferno' if sign==1 else 'winter');colored[selected]=cmap(.3+.65*mag)[:,:3]*(.75+.25*shade[selected,None])
  triangles.append(tri);base.append(gray);paint.append(colored)
triangles=np.concatenate(triangles);base=np.concatenate(base);paint=np.concatenate(paint)
def brain(ax,colors):
 ax.add_collection3d(Poly3DCollection(triangles,facecolors=colors,edgecolors='none',linewidths=0,antialiased=False,zsort='average'))
 ax.set(xlim=(-83,83),ylim=(-110,78),zlim=(-78,89));ax.set_box_aspect((166,188,167));ax.view_init(elev=15,azim=45);ax.set_proj_type('ortho');ax.set_axis_off();ax.patch.set_alpha(0)
for name,colors in [('brain-anatomy',base),('brain-pattern',paint)]:
 fig=plt.figure(figsize=(6.4,6),dpi=140);ax=fig.add_axes([-.13,-.13,1.26,1.26],projection='3d');brain(ax,colors)
 fig.savefig(OUT/(name+'.png'),transparent=True);plt.close(fig)
fig=plt.figure(figsize=(5.12,5.12),dpi=100);ax=fig.add_axes([-.25,-.31,1.50,1.50],projection='3d');brain(ax,paint)
fig.savefig(OUT/'brain-icon.png',transparent=True);plt.close(fig)
fig=plt.figure(figsize=(12,6.3),dpi=100,facecolor='#faf8f4');ax=fig.add_axes([-.045,-.03,.62,1.1],projection='3d');brain(ax,paint)
fig.text(.55,.63,'Neuromarkers',fontsize=37,fontfamily='DejaVu Serif',color='#1a2332')
fig.text(.555,.48,'Population-level brain patterns',fontsize=18,color='#3f6390')
fig.text(.555,.34,'Browse · Explore · Download',fontsize=17,color='#3a4049')
fig.text(.555,.18,'CANlab  |  Brain models for new studies',fontsize=12,color='#646b75')
fig.savefig(OUT/'neuromarkers-social.png',facecolor=fig.get_facecolor());plt.close(fig)
(OUT/'provenance.json').write_text(json.dumps({'anatomy':'CANlab canonical left-cutaway mesh; display-only vertex clustering at 1.4 mm','anatomy_source':'https://github.com/canlab/CanlabCore/blob/master/CanlabCore/canlab_canonical_brains/Canonical_brains_surfaces/canlab_canonical_brain_surface_left_cutaway.mat','pattern':'Neurologic Pain Signature (NPS), Wager et al. 2013; strongest 35% per sign, nearest-neighbor sampling at displayed triangle centers','pattern_source':'https://github.com/canlab/Paradigms_Public/tree/master/2022_Han_Kragel_RTNF_real_time_neurofeedback/updated_code','paper':'https://doi.org/10.1056/NEJMoa1204471','note':'Illustrative static rendering; not a surface predictor. Only rendered images are distributed here.'},indent=2)+'\n')
print('Rendered',len(triangles),'triangles into branding assets')
