export function cutoff(stats, amount, mode) {
  if (!stats.count) return Infinity;
  if (mode === 'absolute') return Math.max(0, Number(amount));
  const percent = Math.max(0, Math.min(100, Math.round(Number(amount))));
  return percent === 0 ? Infinity : stats.cutoffs[percent];
}
export function visible(value, positive, negative, sign = 'both') {
  return Number.isFinite(value) && ((value > 0 && sign !== 'negative' && value >= positive) || (value < 0 && sign !== 'positive' && -value >= negative));
}
export function maskValues(original, positive, negative, sign) {
  const masked = original.slice();
  for (let i = 0; i < original.length; i++) if (!visible(original[i], positive, negative, sign)) masked[i] = 0;
  return masked;
}
export function summarize(values) {
  const result = {};
  for (const [name, sign] of [['positive',1],['negative',-1]]) {
    const sorted = Array.from(values).filter(v=>Number.isFinite(v)&&v*sign>0).map(v=>v*sign).sort((a,b)=>a-b);
    const quantile = q=> {const n=(sorted.length-1)*q, i=Math.floor(n);return sorted[i]+(sorted[Math.min(i+1,sorted.length-1)]-sorted[i])*(n-i)};
    result[name]={count:sorted.length,max:sorted.at(-1)||0,cutoffs:Array.from({length:101},(_,p)=>sorted.length?quantile(1-p/100):0)};
  }
  return result;
}
export function matches(study, query, filters) {
  const text = JSON.stringify(study).normalize('NFD').replace(/[\u0300-\u036f]/g,'').toLowerCase();
  if (!query.normalize('NFD').replace(/[\u0300-\u036f]/g,'').toLowerCase().trim().split(/\s+/).every(w=>text.includes(w))) return false;
  return Object.entries(filters).every(([axis, selected])=>!selected.length||selected.some(tag=>(study[axis]||[]).includes(tag)|| (axis==='targets'&&study.maps.some(m=>m.targets.includes(tag)))));
}
export function parseDisplay(params) {
  const num=(key,fallback,max)=>{const v=Number(params.get(key));return params.has(key)&&Number.isFinite(v)?Math.max(0,Math.min(max,v)):fallback};
  const mode=params.get('mode')==='absolute'?'absolute':'percent';
  return {mode,positive:num('pos',35,mode==='percent'?100:1e12),negative:num('neg',35,mode==='percent'?100:1e12),sync:params.get('sync')!=='0',sign:['positive','negative'].includes(params.get('sign'))?params.get('sign'):'both',opacity:num('opacity',.8,1),anatomy:num('anatomy',1,1)};
}
