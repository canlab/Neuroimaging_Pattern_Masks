import {initHero} from './hero.js';
import {renderGraph} from './graph.js';
import {cutoff, maskValues, summarize, matches, parseDisplay} from './state.js';
const $=id=>document.getElementById(id), root=new URL('./',document.baseURI), url=p=>new URL(p,root).href;
// Keep static logo/icon URLs stable while SPA navigation changes the document URL.
const documentBase=document.querySelector('base')||document.head.appendChild(document.createElement('base'));documentBase.href=root.href;
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fetchJSON=async p=>{const r=await fetch(url(p));if(!r.ok)throw Error(`${p}: HTTP ${r.status}`);return r.json()};
let catalog, study, map, nv, surface, original, surfaceOriginal=[], display=parseDisplay(new URLSearchParams(location.search)), busy=false, localFile=null, filters={domains:[],modalities:[],targets:[]};
let modules, atlasLabels;
const params=()=>new URLSearchParams(location.search);
function tags(tags){return tags.map(t=>`<span class="tag">${esc(t)}</span>`).join('')}
function setBusy(value){busy=value; for(const id of ['map-select','reset','local-file','png'])$(id).disabled=value;}
function dispose(){for(const n of [nv,surface]){if(n){n.cleanup();n.gl?.getExtension('WEBGL_lose_context')?.loseContext()}}nv=null;surface=null;for(const id of ['gl','surface-gl']){const old=$(id),fresh=old.cloneNode();old.replaceWith(fresh)}$('canlab-niivue-controls').replaceChildren();$('canlab-niivue-readout').replaceChildren();original=null;surfaceOriginal=[];}
function displayURL(embed=false){const u=new URL(`marker/${study.id}/`,root);const p=u.searchParams;p.set('map',map.id);p.set('pos',display.positive);p.set('neg',display.negative);p.set('mode',display.mode);p.set('sync',display.sync?'1':'0');p.set('sign',display.sign);p.set('opacity',display.opacity);p.set('anatomy',display.anatomy);p.set('negativeColor',$('negative-color').value);p.set('hemisphere',$('hemisphere').value);if(nv){p.set('color',nv.volumes.at(-1)?.colormap||'inferno');p.set('layout',nv.opts.sliceType);p.set('xyz',Array.from(nv.scene.crosshairPos).map(x=>x.toFixed(5)).join(','));const sels=$('canlab-niivue-controls').querySelectorAll('select');if(sels.length>1)p.set('atlas',sels[sels.length-1].value)}if(surface){p.set('azimuth',surface.scene.renderAzimuth);p.set('elevation',surface.scene.renderElevation)}if(embed)p.set('embed','1');return u;}
function remember(){if(study&&!localFile)history.replaceState(null,'',displayURL(params().get('embed')==='1'))}
function renderCatalog(){
 const q=$('search').value;let shown=catalog.studies.filter(s=>matches(s,q,filters));
 shown.sort((a,b)=>$('sort').value==='name'?a.name.localeCompare(b.name):($('sort').value==='newest'?b.year-a.year:a.year-b.year));
 $('count').textContent=`${shown.length} studies · ${shown.reduce((n,s)=>n+s.maps.length,0)} maps`;
 $('results').innerHTML=shown.length?shown.map(s=>`<a class="card" href="${url(`marker/${s.id}/`)}"><img src="${url(s.maps[0].preview)}" loading="lazy" alt="Axial preview: ${esc(s.name)}" width="360" height="260"><div class="card-body"><h3>${esc(s.name)}</h3><div class="tags">${tags(s.targets.slice(0,3))}</div><p>${esc(s.description.slice(0,150))}${s.description.length>150?'…':''}</p></div><div class="card-meta">${s.year} · ${s.maps.length} maps · ${esc(s.citation.split(' (')[0].slice(0,50))}</div></a>`).join(''):'<p>No matches. Try fewer filters or a broader search.</p>';
 renderGraph($('signature-graph'),shown,url);
 const u=new URL(root);u.hash=location.hash;if(q)u.searchParams.set('q',q);for(const [axis,selected]of Object.entries(filters))for(const t of selected)u.searchParams.append(axis,t);u.searchParams.set('sort',$('sort').value);history.replaceState(null,'',u);
}
function renderFilters(){
 $('filters').innerHTML=Object.entries({domains:'Domain',modalities:'Sensory modality',targets:'Target'}).map(([axis,label])=>{
 const values=[...new Set(catalog.studies.flatMap(s=>[...s[axis],...(axis==='targets'?s.maps.flatMap(m=>m.targets):[])]))].sort();
 return `<details class="filter-group" ${axis!=='targets'||filters[axis].length?'open':''}><summary>${label}${filters[axis].length?' · '+filters[axis].length+' selected':''}</summary><fieldset aria-label="${label}">${values.map(t=>`<label class="filter-option"><input type="checkbox" data-axis="${axis}" value="${esc(t)}" ${filters[axis].includes(t)?'checked':''}><span>${esc(t)}</span><small>${catalog.studies.filter(s=>s[axis].includes(t)||(axis==='targets'&&s.maps.some(m=>m.targets.includes(t)))).length}</small></label>`).join('')}</fieldset></details>`}).join('');
}
function syncControls(){
 $('threshold-mode').value=display.mode;$('sync').checked=display.sync;$('sign').value=display.sign;$('opacity').value=display.opacity;$('anatomy-opacity').value=display.anatomy;
 for(const side of ['positive','negative']){const max=display.mode==='percent'?100:Math.max(map.stats.positive.max,map.stats.negative.max,1e-9);for(const suffix of ['range','number']){const el=$(side+'-'+suffix);el.max=max;el.step=display.mode==='percent'?1:'any';if(document.activeElement!==el||el.type!=='number')el.value=display[side];el.disabled=!map.stats[side].count;} }
 document.querySelectorAll('.unit').forEach(e=>e.textContent=display.mode==='percent'?'%':'');
}
function applyThreshold(){
 if(!nv||!original)return;
 const pc=cutoff(map.stats.positive,display.positive,display.mode),nc=cutoff(map.stats.negative,display.negative,display.mode),ov=nv.volumes.at(-1);
 const positive=display.sign==='negative'?Infinity:pc,negative=display.sign==='positive'?Infinity:nc;
 // Work in physical values even when a source uses NIfTI scaling.
 const physical=original.physical; const masked=maskValues(physical,positive,negative,display.sign);
 // Converted web volumes are floating point, and the viewer copy may safely use physical units.
 ov.img=masked;ov.hdr.scl_slope=1;ov.hdr.scl_inter=0;ov.hdr.datatypeCode=16;
 ov.cal_min=Number.isFinite(pc)?pc:map.stats.positive.max+1;ov.cal_max=Math.max(map.stats.positive.max,ov.cal_min+1e-9);ov.cal_minNeg=-(Number.isFinite(nc)?nc:map.stats.negative.max+1);ov.cal_maxNeg=-Math.max(map.stats.negative.max,-ov.cal_minNeg+1e-9);ov.colormapType=1;ov.colormapNegative=$('negative-color').value;
 nv.setOpacity(nv.volumes.length-1,display.opacity);nv.setOpacity(0,display.anatomy);nv.updateGLVolume();drawColorbar();
 for(const side of ['positive','negative']){const c=side==='positive'?positive:negative;$(side+'-cutoff').textContent=map.stats[side].count?`${map.stats[side].count.toLocaleString()} voxels · ${Number.isFinite(c)?'≥ '+c.toPrecision(3):'hidden'}`:'No values of this sign';}
 if(surface){surface.meshes.forEach((mesh,i)=>{const layer=mesh.layers[0];if(!layer)return;layer.values=maskValues(surfaceOriginal[i],positive,negative,display.sign);Object.assign(layer,{colormap:ov.colormap,colormapNegative:ov.colormapNegative,cal_min:ov.cal_min,cal_max:ov.cal_max,cal_minNeg:ov.cal_minNeg,cal_maxNeg:ov.cal_maxNeg,useNegativeCmap:true,isTransparentBelowCalMin:true,opacity:display.opacity});mesh.updateMesh(surface.gl)});surface.drawScene();}
}
function updateDetails(){
 $('map-role').textContent=localFile?'Local overlay':map.role;$('map-title').textContent=map.name;$('map-description').textContent=map.description;
 $('study-description').textContent=localFile?'This local file stays in your browser. Its publisher citation and template space are unknown. Restore a catalog map to share or embed a reproducible view.':study.description;
 $('citation').textContent=localFile?'No publication associated with this local overlay.':study.citation;$('citation').href=study.paper||'#';$('citation').hidden=!!localFile;
 $('source-link').href='https://github.com/canlab/Neuroimaging_Pattern_Masks/blob/master/'+map.source?.split('/').map(encodeURIComponent).join('/');$('source-link').hidden=!!localFile;
 $('share').disabled=$('embed').disabled=!!localFile;
 $('technical').innerHTML=[['Map type',map.role],['Space',localFile?'Not verified':study.space],['Dimensions',map.shape?.join(' × ')||'Local volume'],['Voxel size',map.voxel_size?.map(x=>x.toFixed(2)).join(' × ')+' mm'],['Orientation',map.orientation||'From header'],['Frame',map.frame==null?'Single 3D map':`${map.frame+1} (exported as 3D)`],['Source file',map.source?.split('/').at(-1)||localFile?.name]].map(([k,v])=>`<dt>${esc(k)}</dt><dd>${esc(v)}</dd>`).join('');
 $('fallback').src=url(map.preview||study.maps[0].preview);
}
async function loadMap(selected){
 if(busy)return;setBusy(true);dispose();localFile=null;map=selected;$('viewer').hidden=false;$('fallback').hidden=true;$('viewer-status').textContent='Loading map and anatomy…';$('surface-status').textContent='Loading cortical surfaces…';$('action-status').textContent='';updateDetails();syncControls();
 try{
  if(!modules){const css=document.createElement('link');css.rel='stylesheet';css.href=url(catalog.assets.css);document.head.prepend(css);modules=await Promise.all([import(url(catalog.assets.viewer)),import(url(catalog.assets.niivue))]);atlasLabels=(await fetchJSON(catalog.assets.shared+'atlas-descriptions.json')).map(r=>r.short_name+' — '+r.full_name);}
  const p=params();nv=await modules[0].canlabNiivue('gl',{underlay:url(catalog.assets.shared+'underlay.nii.gz'),atlas:url(catalog.assets.shared+'atlas.nii.gz'),atlasLabels,overlay:url(map.url),cal_min:0,cal_max:Math.max(map.stats.positive.max,map.stats.negative.max),showOpacity:false,allowLoadOverlay:false,showRenderInMultiplanar:false,colormap:p.get('color')||'inferno',backColor:[.025,.04,.055,1]});
  nv.volumes[0].colorbarVisible=false;nv.opts.isColorbar=false;
  const ov=nv.volumes.at(-1);original={physical:Float32Array.from(ov.img,x=>x*(ov.hdr.scl_slope||1)+(ov.hdr.scl_inter||0))};
  $('negative-color').value=['winter','cool','blue'].includes(p.get('negativeColor'))?p.get('negativeColor'):'winter';
  const layout=Number(p.get('layout'));if(p.has('layout')&&layout>=0&&layout<=4)nv.setSliceType(layout);
  const xyz=(p.get('xyz')||'').split(',').map(Number);if(xyz.length===3&&xyz.every(x=>Number.isFinite(x)&&x>=0&&x<=1))nv.scene.crosshairPos=xyz;
  const selects=$('canlab-niivue-controls').querySelectorAll('select');if(selects.length>1){selects[selects.length-1].value=['outline','shaded','off'].includes(p.get('atlas'))?p.get('atlas'):'off';selects[selects.length-1].dispatchEvent(new Event('change'));}
  applyThreshold();$('viewer-status').textContent='';
  try{
   surface=new modules[1].Niivue({backColor:[.025,.04,.055,1],isColorbar:false,show3Dcrosshair:false});await surface.attachTo('surface-gl');surface.setSliceType(4);
   await surface.loadMeshes(['L','R'].map(side=>({url:url(catalog.assets.shared+side+'.surf.gii'),rgba255:[175,187,196,255],layers:[{url:url(map.surfaces[side]),colormap:'inferno',colormapNegative:'winter',useNegativeCmap:true,opacity:display.opacity,cal_min:0,cal_max:Math.max(map.stats.positive.max,map.stats.negative.max)}]})));
   surface.meshes.forEach(m=>surface.setMeshShader(m.id,'Matte'));surfaceOriginal=surface.meshes.map(m=>m.layers[0].values.slice());surface.setRenderAzimuthElevation(p.has('azimuth')?Number(p.get('azimuth')):110,p.has('elevation')?Number(p.get('elevation')):15);$('hemisphere').value=['L','R'].includes(p.get('hemisphere'))?p.get('hemisphere'):'both';hemispheres();applyThreshold();$('surface-status').textContent='';
  }catch(e){$('surface-status').textContent='Surface unavailable. Orthographic views remain interactive.';console.error(e);}
  remember();
 }catch(e){$('viewer-status').textContent='The interactive viewer could not load. A static preview and map download are available. '+e.message;$('viewer').hidden=true;$('fallback').hidden=false;console.error(e)}finally{setBusy(false)}
}
function hemispheres(){if(surface){surface.meshes.forEach((m,i)=>m.visible=$('hemisphere').value==='both'||$('hemisphere').value===['L','R'][i]);surface.drawScene()}}
async function openStudy(s){study=s;document.body.classList.toggle('embed',params().get('embed')==='1');$('browse').hidden=true;document.querySelector('.intro').hidden=true;$('detail').hidden=false;$('study-title').textContent=s.name;$('study-kicker').textContent=`${s.year} · ${s.maps.length} maps`;$('study-tags').innerHTML=tags(s.domains);$('map-select').innerHTML=s.maps.map(m=>`<option value="${esc(m.id)}">${esc(m.name)} — ${esc(m.role)}</option>`).join('');const selected=s.maps.find(m=>m.id===params().get('map'))||s.maps[0];$('map-select').value=selected.id;await loadMap(selected);}
function browse(){dispose();study=null;document.body.classList.remove('embed');$('detail').hidden=true;$('browse').hidden=false;document.querySelector('.intro').hidden=false;renderCatalog();if(location.hash==='#browse'){$('browse').scrollIntoView({block:'start'});$('search').focus({preventScroll:true})}}
async function route(){if(busy)return;const s=catalog.studies.find(s=>location.pathname.includes('/marker/'+s.id+'/'));if(s){display=parseDisplay(params());await openStudy(s)}else browse()}
let raf;function change(side,value){if(String(value).trim()==='')return;const max=display.mode==='percent'?100:Math.max(map.stats.positive.max,map.stats.negative.max,1e-9);const v=Math.max(0,Math.min(max,Number(value)));if(!Number.isFinite(v))return;display[side]=v;if(display.sync)display[side==='positive'?'negative':'positive']=v;syncControls();cancelAnimationFrame(raf);raf=requestAnimationFrame(()=>{applyThreshold();remember()})}
async function copy(text){try{await navigator.clipboard.writeText(text);$('action-status').textContent='Copied to clipboard.'}catch{$('action-status').textContent=text}}
function downloadBlob(blob,name){const u=URL.createObjectURL(blob),a=document.createElement('a');a.href=u;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(u),10000)}
async function main(){
 initHero();
 catalog=await fetchJSON('catalog.json');$('inventory').textContent=`${catalog.studies.length} studies · ${catalog.studies.reduce((n,s)=>n+s.maps.length,0)} maps`;
 for(const axis of Object.keys(filters))filters[axis]=params().getAll(axis);$('search').value=params().get('q')||'';$('sort').value=params().get('sort')||'name';renderFilters();
 for(const view of ['graph','tile'])$(view+'-view').onclick=()=>{const graph=view==='graph';document.querySelector('.skip').href=graph?'#signature-graph':'#results';$('graph-panel').hidden=!graph;$('results').hidden=graph;$('graph-view').setAttribute('aria-pressed',String(graph));$('tile-view').setAttribute('aria-pressed',String(!graph))};
 $('search').oninput=renderCatalog;$('sort').onchange=renderCatalog;$('filters').onchange=e=>{const {axis}=e.target.dataset;if(!axis)return;filters[axis]=Array.from($('filters').querySelectorAll(`input[data-axis="${axis}"]:checked`),x=>x.value);renderCatalog()};$('clear').onclick=()=>{filters={domains:[],modalities:[],targets:[]};$('search').value='';renderFilters();renderCatalog()};
 $('browse').onclick=async e=>{const a=e.target.closest('a.card,a.graph-node');if(!a)return;e.preventDefault();history.pushState(null,'',a.href);display=parseDisplay(params());await route();window.scrollTo(0,0)};$('back').onclick=()=>{if(!busy){history.pushState(null,'',root);browse()}};window.onpopstate=()=>busy?location.reload():route();
 $('map-select').onchange=()=>{display=parseDisplay(new URLSearchParams());loadMap(study.maps.find(m=>m.id===$('map-select').value))};
 for(const side of ['positive','negative'])for(const suffix of ['range','number'])$(side+'-'+suffix).oninput=e=>change(side,e.target.value);
 $('sync').onchange=()=>{display.sync=$('sync').checked;if(display.sync)change('positive',display.positive);remember()};
 $('threshold-mode').onchange=()=>{const next=$('threshold-mode').value;if(next==='absolute'){display.positive=cutoff(map.stats.positive,display.positive,display.mode);display.negative=cutoff(map.stats.negative,display.negative,display.mode);for(const side of ['positive','negative'])if(!Number.isFinite(display[side]))display[side]=map.stats[side].max;display.sync=false;}else{display.positive=display.negative=35}display.mode=next;syncControls();applyThreshold();remember()};
 $('sign').onchange=()=>{display.sign=$('sign').value;applyThreshold();remember()};for(const [id,key] of [['opacity','opacity'],['anatomy-opacity','anatomy']])$(id).oninput=()=>{display[key]=Number($(id).value);applyThreshold();remember()};$('negative-color').onchange=()=>{applyThreshold();remember()};
 $('canlab-niivue-controls').addEventListener('change',()=>{applyThreshold();remember()});$('canlab-niivue-controls').addEventListener('click',()=>setTimeout(remember,0));$('canlab-niivue-canvas-container').addEventListener('pointerup',()=>setTimeout(remember,0));$('surface-container').addEventListener('pointerup',()=>setTimeout(remember,0));
 $('hemisphere').onchange=()=>{hemispheres();remember()};$('surface-reset').onclick=()=>{surface?.setRenderAzimuthElevation(110,15);remember()};$('reset').onclick=()=>{history.replaceState(null,'',new URL(`marker/${study.id}/?map=${map.id}`,root));display=parseDisplay(new URLSearchParams());loadMap(study.maps.find(m=>m.id===$('map-select').value))};
 $('share').onclick=()=>copy(displayURL().href);$('embed').onclick=()=>copy(`<iframe src="${displayURL(true).href}" title="${esc(study.name)}" width="100%" height="1000" loading="lazy" style="border:0"></iframe>`);
 $('download').onclick=async()=>{try{let blob=localFile||await fetch(url(map.url)).then(r=>{if(!r.ok)throw Error('Download failed');return r.blob()});const bytes=new Uint8Array(await blob.slice(0,2).arrayBuffer());if(bytes[0]===31&&bytes[1]===139)blob=await new Response(blob.stream().pipeThrough(new DecompressionStream('gzip'))).blob();downloadBlob(blob,map.id+'.nii');$('action-status').textContent='Downloaded the complete selected map, without applying display thresholds.'}catch(e){$('action-status').textContent=e.message}};
 $('png').onclick=()=>{try{if(!nv)throw Error('Interactive viewer unavailable.');nv.drawScene();surface?.drawScene();const a=$('gl'),b=$('surface-gl'),out=document.createElement('canvas');out.width=a.width+(surface?b.width:0);out.height=Math.max(a.height,surface?b.height:0)+152;const c=out.getContext('2d');c.fillStyle='#080d12';c.fillRect(0,0,out.width,out.height);c.drawImage(a,0,0);if(surface)c.drawImage(b,a.width,0);c.drawImage($('weight-colorbar'),0,a.height,a.width,72);c.fillStyle='white';c.font='20px sans-serif';c.fillText(map.name,20,out.height-45);c.font='16px sans-serif';c.fillText(`${display.mode}: positive ${Number(display.positive).toPrecision(4)}, negative ${Number(display.negative).toPrecision(4)} · ${display.sign}`,20,out.height-18);out.toBlob(blob=>{if(blob)downloadBlob(blob,map.id+'.png')})}catch(e){$('action-status').textContent=e.message}};
 $('local-file').onchange=async e=>{const file=e.target.files[0];if(!file||!nv)return;try{setBusy(true);const volume=await modules[1].NVImage.loadFromFile({file});if(volume.hdr.dims[4]>1)throw Error('Please select a single 3D local overlay.');const old=nv.volumes.at(-1);nv.removeVolume(old);nv.addVolume(volume);original={physical:Float32Array.from(volume.img,x=>x*(volume.hdr.scl_slope||1)+(volume.hdr.scl_inter||0))};localFile=file;map={id:'local-overlay',name:file.name,role:'Local overlay',description:'Locally loaded image. Display controls operate on this image; the file is not uploaded.',stats:summarize(original.physical),shape:Array.from(volume.hdr.dims).slice(1,4),voxel_size:Array.from(volume.hdr.pixDims).slice(1,4)};if(surface){surface.meshes.forEach(m=>m.visible=false);surface.drawScene()}$('surface-status').textContent='No precomputed surface for local overlays.';surfaceOriginal=[];surface?.meshes.forEach(m=>m.layers=[]);display=parseDisplay(new URLSearchParams());updateDetails();syncControls();applyThreshold();}catch(e){$('action-status').textContent=e.message}finally{setBusy(false)}};
 await route();
 const context=document.modelContext;if(context?.registerTool){try{await context.registerTool({name:'search_neuromarkers',description:'Search the catalog and show matching map-based neuromarkers.',inputSchema:{type:'object',properties:{query:{type:'string'}},required:['query'],additionalProperties:false},execute:({query})=>{if(typeof query!=='string'||query.length>500)throw Error('Invalid query');if(busy)throw Error('Viewer is loading');browse();$('search').value=query;renderCatalog();return catalog.studies.filter(s=>matches(s,query,filters)).map(s=>({id:s.id,name:s.name,maps:s.maps.length}))}})}catch(e){console.warn('Optional WebMCP unavailable',e)}}
}
main().catch(e=>{$('results').textContent='Could not load the catalog. '+e.message;$('inventory').textContent='Catalog unavailable';console.error(e)});

function drawColorbar(){
 const canvas=$('weight-colorbar'),ctx=canvas.getContext('2d'),ov=nv.volumes.at(-1);
 ctx.clearRect(0,0,canvas.width,canvas.height);ctx.fillStyle='#121b22';ctx.fillRect(0,0,1000,72);
 const sides=['negative','positive'].filter(side=>display.sign==='both'||display.sign===side);
 const fmt=v=>v===0?'0':Number(v.toPrecision(3)).toString();
 sides.forEach((side,j)=>{
  const negative=side==='negative',x=28+j*(944/sides.length),width=944/sides.length-32;
  const lut=nv.colormap(negative?ov.colormapNegative:ov.colormap),max=map.stats[side].max,cut=cutoff(map.stats[side],display[side],display.mode);
  ctx.fillStyle='#c6d5dc';ctx.font='19px system-ui';ctx.textAlign='left';ctx.fillText(negative?'Negative weights':'Positive weights',x,20);
  if(!max||!Number.isFinite(cut)||cut>max){ctx.fillStyle='#90a4b0';ctx.fillText('Hidden / no values',x,55);return;}
  for(let i=0;i<width;i++){const f=i/(width-1),value=negative?max-(max-cut)*f:cut+(max-cut)*f,index=Math.min(255,Math.round(value/max*255))*4;ctx.fillStyle=`rgb(${lut[index]},${lut[index+1]},${lut[index+2]})`;ctx.fillRect(x+i,28,1,13)}
  ctx.fillStyle='#d9e5eb';ctx.font='18px system-ui';ctx.textAlign='left';ctx.fillText(fmt(negative?-max:cut),x,65);ctx.textAlign='right';ctx.fillText(fmt(negative?-cut:max),x+width,65);
 });
}
