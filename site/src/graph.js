// Editorial taxonomy, not a learned embedding or a similarity network.
const groups=[
 ['Pain','Pain & nociception','#a85c2c'],
 ['Aversive / negative affect','Threat & negative affect','#7855a3'],
 ['Cognitive & social','Self & social cognition','#3f6390'],
 ['Appetitive / reward','Reward & desire','#287b69'],
 ['Clinical','Clinical targets','#a6506b'],
 ['Physiology','Body & physiology','#7c7528']
];
const aliases={
 '2011-wager-jneuro-placebo-prediction':'Placebo analgesia',
 '2015-chang-plosbiology-pines':'PINES · Negative affect',
 '2015-kragel-emotionclassificationbpls':'BPLS · Emotion categories',
 '2015-woo-naturecomms-rejection':'Rejection & physical pain',
 '2016-eisenbarth-jneuro-autonomic-patterns':'Autonomic · GSR & HR',
 '2016-krishnan-elife-vps':'VPS · Vicarious pain',
 '2017-ashar-care-distress':'Empathic care & distress',
 '2017-woo-siips1':'SIIPS1 · Pain',
 '2018-kragel-mfc-generalizability':'Medial frontal patterns',
 '2018-reddan-threat-conditioning-imex':'ImEx · Conditioned threat',
 '2019-kragel-emotion-schemas':'Emotion schemas',
 '2019-lee-jpain-backpain':'Back pain · Perfusion',
 '2019-yu-koban-guilt':'Interpersonal guilt',
 '2020-geuter-pain-multivariate-mediation-pdm':'PDM · Pain mediation',
 '2020-silvestrini-rainville-pain-cogcontrol-interaction-amcc':'Pain & cognitive control',
 '2020-zhou-general-vicarious-pain':'General vicarious pain',
 '2021-ceko-mpa2-multiaversive':'MPA2 · Multi-aversive',
 '2021-zhou-subjective-fear':'VIFS · Subjective fear',
 '2021-vanthoff-basic-sexual-image-classifier':'BASIC · Sexual images',
 '2022-koban-ncs-craving':'NCS · Craving',
 '2022-coll-pain-monetary-reward-decision-value':'Decision value · Pain & money',
 '2023-speer-brain-reward-signature-brs':'BRS · Brain reward',
 '2024-feps-facial-expressions-of-pain-signature':'FEPS · Pain expression',
 '2026-acil-mentalizing-self-other':'Mentalizing · Self & other',
 '2026-murillo-pifonem':'PiFoneM · Movement fear'
};
const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
document.addEventListener('keydown',e=>{if(e.key==='Escape')document.querySelector('#graph-popup')?.remove()});
document.addEventListener('scroll',()=>document.querySelector('#graph-popup')?.remove(),{capture:true,passive:true});
export function renderGraph(container,studies,url){
 document.querySelector('#graph-popup')?.remove();
 if(!studies.length){container.innerHTML='<p>No matches. Try fewer filters or a broader search.</p>';return;}
 container.innerHTML='<div class="tree-origin"><span class="origin-dot"></span> NEUROMARKERS <span class="origin-line"></span></div><div class="branch-grid"></div>';
 const grid=container.querySelector('.branch-grid');
 for(const [domain,title,color] of groups){
  const members=studies.filter(s=>s.domains[0]===domain);if(!members.length)continue;
  const height=members.length*48+52, mid=height/2;
  const branch=document.createElement('section');branch.className='graph-branch';branch.style.setProperty('--branch',color);branch.style.height=height+'px';
  branch.innerHTML=`<h3>${esc(title)} <small>${members.length}</small></h3><svg viewBox="0 0 400 ${height}" preserveAspectRatio="none" aria-hidden="true"><path class="branch-stem" d="M 0 0 L 16 0 Q 30 0 30 20 L 30 ${height-24}"/>${members.map((s,i)=>{const y=55+i*48,x=68+(i%3)*12;return `<path d="M 30 ${Math.min(mid,y)} C 30 ${y}, ${x-25} ${y}, ${x} ${y}"/>`}).join('')}</svg>`;
  members.forEach((s,i)=>{
   const a=document.createElement('a');a.className='graph-node';a.href=url(`marker/${s.id}/`);a.style.top=(34+i*48)+'px';a.style.left=(17+(i%3)*3)+'%';a.style.right='2%';
   a.innerHTML=`<span class="node-dot"></span><span class="node-label">${esc(aliases[s.id]||s.name)}<small>${s.year} <span>· ${s.maps.length} maps</span></small></span>`;a.setAttribute('aria-label',`${s.name}, ${s.maps.length} maps. Open study.`);
   const show=()=>{
    document.querySelector('#graph-popup')?.remove();const popup=document.createElement('div');popup.id='graph-popup';popup.className='graph-popup';popup.setAttribute('role','tooltip');
    popup.innerHTML=`<p class="eyebrow">${esc(s.domains.join(' · '))}</p><h3>${esc(s.name)}</h3><p>${esc(s.description)}</p><div class="tags">${s.targets.slice(0,4).map(t=>`<span class="tag">${esc(t)}</span>`).join('')}</div><p class="popup-meta">${s.year} · ${s.maps.length} available maps <span>Open study ↗</span></p>`;
    document.body.append(popup);a.setAttribute('aria-describedby',popup.id);const box=a.getBoundingClientRect();popup.style.left=Math.max(12,Math.min(box.left,innerWidth-popup.offsetWidth-12))+'px';popup.style.top=Math.max(12,box.bottom+8+popup.offsetHeight<innerHeight?box.bottom+8:box.top-popup.offsetHeight-8)+'px';
   };
   const hide=()=>{document.querySelector('#graph-popup')?.remove();a.removeAttribute('aria-describedby')};
   a.onpointerenter=show;a.onfocus=show;a.onpointerleave=()=>{if(document.activeElement!==a)hide()};a.onblur=hide;a.onclick=hide;a.onkeydown=e=>{if(e.key==='Escape')hide()};branch.append(a);
  });grid.append(branch);
 }
}
