const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
try{
 const catalog=await fetch('catalog.json').then(r=>r.json()),base=catalog.assets.shared;
 for(const type of ['csv','json']){const a=document.getElementById(type);a.href=base+'atlas-descriptions.'+type;a.download='atlas-descriptions.'+type}
 const rows=await fetch(base+'atlas-descriptions.json').then(r=>r.json());
 function render(){const q=document.getElementById('atlas-search').value.toLowerCase(),shown=rows.filter(r=>JSON.stringify(r).toLowerCase().includes(q));document.getElementById('atlas-count').textContent=shown.length+' of '+rows.length+' parcels';document.getElementById('atlas-rows').innerHTML=shown.map(r=>`<tr><td><small>${r.code}</small>${esc(r.short_name)}</td><td>${esc(r.full_name)}${r.components.length?`<details><summary>Source components</summary>${esc(r.components.join('; '))}</details>`:''}<small>${esc(r.note)}</small></td><td>${r.source.split(' ; ').map((s,i)=>`<a href="${esc(s)}" target="_blank" rel="noopener">Source ${i+1} ↗</a>`).join('<br>')}</td></tr>`).join('')}
 document.getElementById('atlas-search').oninput=render;render();
}catch(e){document.getElementById('atlas-count').textContent='Unable to load lookup table: '+e.message}
