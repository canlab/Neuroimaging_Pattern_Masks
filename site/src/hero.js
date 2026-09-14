// A small 2D mask reveals the original image; the source PNG is never modified.
export function initHero(){
 const button=document.getElementById('hero-brain'),image=button.querySelector('img'),canvas=button.querySelector('canvas');
 const ctx=canvas.getContext('2d'),mask=document.createElement('canvas'),ink=mask.getContext('2d');
 const reduced=matchMedia('(prefers-reduced-motion: reduce)');let spots=[],frame=0,timer=0,visible=false,lastPoint=null;
 const lifetime=3600;
 function resize(){const box=button.getBoundingClientRect(),scale=Math.min(devicePixelRatio||1,2);canvas.width=mask.width=Math.max(1,Math.round(box.width*scale));canvas.height=mask.height=Math.max(1,Math.round(box.height*scale));if(spots.length)start()}
 function draw(now){
  frame=0;spots=spots.filter(s=>now-s.time<lifetime);ctx.clearRect(0,0,canvas.width,canvas.height);ink.clearRect(0,0,mask.width,mask.height);
  if(!spots.length)return;
  for(const spot of spots){
   const age=now-spot.time,fade=Math.min(1,(lifetime-age)/1800),appear=spot.ambient?Math.min(1,age/900):1;
   const radius=spot.radius*canvas.width,gradient=ink.createRadialGradient(spot.x*canvas.width,spot.y*canvas.height,0,spot.x*canvas.width,spot.y*canvas.height,radius);
   gradient.addColorStop(0,`rgba(255,255,255,${fade*appear})`);gradient.addColorStop(.5,`rgba(255,255,255,${fade*appear*.95})`);gradient.addColorStop(1,'rgba(255,255,255,0)');ink.fillStyle=gradient;ink.fillRect(spot.x*canvas.width-radius,spot.y*canvas.height-radius,radius*2,radius*2);
  }
  const scale=Math.min(canvas.width/image.naturalWidth,canvas.height/image.naturalHeight),w=image.naturalWidth*scale,h=image.naturalHeight*scale;
  ctx.globalCompositeOperation='source-over';ctx.drawImage(image,(canvas.width-w)/2,(canvas.height-h)/2,w,h);ctx.globalCompositeOperation='destination-in';ctx.drawImage(mask,0,0);ctx.globalCompositeOperation='source-over';frame=requestAnimationFrame(draw);
 }
 function start(){if(!frame&&image.complete&&image.naturalWidth)frame=requestAnimationFrame(draw)}
 function add(x,y,ambient=false){if(document.hidden||!visible)return;spots.push({x,y,radius:ambient?.15:.13,time:performance.now(),ambient});if(spots.length>48)spots.shift();start()}
 function pointer(e){const box=button.getBoundingClientRect(),x=(e.clientX-box.left)/box.width,y=(e.clientY-box.top)/box.height,now=performance.now();if(lastPoint&&now-lastPoint.time<45&&Math.hypot(x-lastPoint.x,y-lastPoint.y)<.035)return;lastPoint={x,y,time:now};add(x,y)}
 button.addEventListener('pointermove',pointer);button.addEventListener('pointerdown',pointer);button.addEventListener('pointerenter',pointer);button.addEventListener('pointerleave',()=>lastPoint=null);
 button.addEventListener('click',e=>{if(e.detail===0)add(.57,.55)});
 button.addEventListener('focus',()=>{if(button.matches(':focus-visible'))add(.57,.55)});
 function schedule(){clearTimeout(timer);if(!visible||document.hidden||reduced.matches)return;timer=setTimeout(()=>{if(visible&&!document.hidden&&!reduced.matches){const centers=[[.18,.32],[.57,.57],[.47,.68],[.3,.44]],point=centers[Math.floor(Math.random()*centers.length)];add(point[0],point[1],true)}schedule()},14000+Math.random()*8000)}
 function clear(){spots=[];cancelAnimationFrame(frame);frame=0;ctx.clearRect(0,0,canvas.width,canvas.height);lastPoint=null}
 new IntersectionObserver(entries=>{visible=entries[0].isIntersecting;if(!visible)clear();schedule()},{threshold:.25}).observe(button);
 new ResizeObserver(resize).observe(button);image.addEventListener('load',resize);resize();
 document.addEventListener('visibilitychange',()=>{if(document.hidden)clear();schedule()});reduced.addEventListener('change',()=>{spots=spots.filter(s=>!s.ambient);schedule()});
}
