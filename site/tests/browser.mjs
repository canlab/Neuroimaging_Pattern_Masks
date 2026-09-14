import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
const base=process.argv[2]||'http://127.0.0.1:8765/';
const browser=await chromium.launch({...(process.env.CI?{}:{channel:'chrome'}),headless:true,args:['--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
const page=await browser.newPage({viewport:{width:1440,height:1000}});const errors=[];
page.on('response',r=>{if(r.status()>=400)console.log('HTTP:',r.status(),r.url())});
page.on('pageerror',e=>errors.push(e.message));page.on('console',m=>{if(m.type()==='error')console.log('BROWSER:',m.text().slice(0,400))});
try{
 await page.goto(base);await page.locator('.card').first().waitFor();
 assert.equal(await page.locator('.card:visible').count(),25);assert.ok(await page.locator('#graph-panel').isHidden());
 assert.equal(await page.locator('#tile-view').getAttribute('aria-pressed'),'true');
 assert.equal(await page.locator('.intro h1').textContent(),'Neuromarkers: Population-level predictive brain models');
 const paragraphs=await page.locator('.intro-copy > p:not(.eyebrow):not(.inventory)').allTextContents();
 assert.deepEqual(paragraphs,[
  'Neuromarkers are population-level patterns of fMRI activity whose expression predicts a stimulus, behavior, mental state, or clinical outcome. In contrast to models individualized for a single person, the models shared here can be applied to data from new participants, allowing open-ended validation across populations, contexts, methodological variations, and outcomes.',
  'These models are in various stages of validation, but all have been successfully validated on at least one independent cohort after model training.',
  'This site focuses on linear models based on fMRI activity, which can be visualized and applied as pattern maps. It does not (yet) include connectivity-based models, or a broader set of neuromarkers based on pathology, electrophysiology, neurochemistry, or other measures.',
  'Browse, explore, and download these patterns for use in new studies.'
 ]);
 assert.equal(await page.locator('.review-link').count(),0);
 assert.deepEqual(await page.locator('header nav a').allTextContents(),['Search','Repository ↗','Reviews']);
 assert.equal(await page.locator('header nav a').first().getAttribute('href'),'./#browse');
 assert.equal(await page.locator('header nav a').last().getAttribute('href'),'reviews.html');
 const searchBox=await page.locator('#search').boundingBox();assert.ok(searchBox.y+searchBox.height<900,'Search is visible on a desktop landing viewport');
 const revealCount=()=>page.locator('.hero-reveal').evaluate(c=>{const d=c.getContext('2d').getImageData(0,0,c.width,c.height).data;let n=0;for(let i=3;i<d.length;i+=4)if(d[i])n++;return n});
 await page.emulateMedia({reducedMotion:'reduce'});
 assert.equal(await revealCount(),0);
 await page.locator('#hero-brain').scrollIntoViewIfNeeded();const hero=await page.locator('#hero-brain').boundingBox();await page.locator('#hero-brain').hover({position:{x:hero.width*.57,y:hero.height*.55}});
 await page.waitForFunction(()=>{const c=document.querySelector('.hero-reveal'),d=c.getContext('2d').getImageData(0,0,c.width,c.height).data;return d.some((v,i)=>i%4===3&&v>0)});
 const locality=await page.locator('.hero-reveal').evaluate(c=>{const d=c.getContext('2d').getImageData(0,0,c.width,c.height).data;let n=0,outside=0;for(let y=0;y<c.height;y++)for(let x=0;x<c.width;x++)if(d[(y*c.width+x)*4+3]){n++;if(Math.hypot(x-c.width*.57,y-c.height*.55)>c.width*.14)outside++}return {n,outside,total:c.width*c.height}});
 assert.ok(locality.n>100&&locality.n<locality.total*.15);assert.equal(locality.outside,0);
 await page.screenshot({path:'/tmp/neuromarker-hero-local.png'});await page.locator('header').hover();
 await page.waitForTimeout(700);assert.ok(await revealCount()>100);
 await page.waitForFunction(()=>{const c=document.querySelector('.hero-reveal');return !c.getContext('2d').getImageData(0,0,c.width,c.height).data.some((v,i)=>i%4===3&&v>0)},{},{timeout:6000});
 // An idle, visible hero occasionally blooms without pointer input.
 await page.emulateMedia({reducedMotion:'no-preference'});
 await page.waitForFunction(()=>{const c=document.querySelector('.hero-reveal');return c.getContext('2d').getImageData(0,0,c.width,c.height).data.some((v,i)=>i%4===3&&v>0)},{},{timeout:12000});
 await page.waitForTimeout(900);await page.screenshot({path:'/tmp/neuromarker-hero-ambient.png'});
 await page.emulateMedia({reducedMotion:'reduce'});
 await page.waitForFunction(()=>{const c=document.querySelector('.hero-reveal');return !c.getContext('2d').getImageData(0,0,c.width,c.height).data.some((v,i)=>i%4===3&&v>0)});
 // Keyboard activation reveals a local patch without a persistent toggle.
 await page.locator('#hero-brain').focus();await page.keyboard.press('Enter');await page.waitForFunction(()=>{const c=document.querySelector('.hero-reveal');return c.getContext('2d').getImageData(0,0,c.width,c.height).data.some((v,i)=>i%4===3&&v>0)});await page.locator('#hero-brain').blur();await page.locator('header').hover();
 assert.equal(await page.locator('link[rel="icon"]').getAttribute('href'),'brand/neuromarkers-ncs.png');
 const social=await page.locator('meta[property="og:image"]').getAttribute('content');assert.match(social,/brand\/neuromarkers-ncs.png$/);
 assert.match(await page.locator('footer').textContent(),/Built by Tor Wager with GPT Astra/);
 assert.equal(await page.locator('footer a').first().getAttribute('href'),'https://torwager.github.io/canlab/');
 await page.screenshot({path:'/tmp/neuromarker-home.png'});
 await page.locator('#graph-view').click();await page.locator('.graph-node').first().waitFor();
 await page.locator('.graph-node').first().scrollIntoViewIfNeeded();await page.locator('.graph-node').first().focus();await page.locator('#graph-popup').waitFor();await page.evaluate(()=>window.scrollBy(0,20));await page.locator('#graph-popup').waitFor();await page.keyboard.press('Escape');assert.ok(await page.locator('#graph-popup').count()===0);
 assert.equal(await page.locator('.graph-node').count(),25);assert.ok(await page.locator('#results').isHidden());
 await page.locator('.graph-node').first().hover();await page.locator('#graph-popup').waitFor();assert.ok((await page.locator('#graph-popup').textContent()).length>100);await page.keyboard.press('Escape');await page.locator('.graph-node').first().blur();await page.locator('.intro').hover();
 await page.screenshot({path:'/tmp/neuromarker-graph.png',fullPage:true});
 await page.locator('#tile-view').click();assert.equal(await page.locator('.card:visible').count(),25);await page.screenshot({path:'/tmp/neuromarker-tiles.png',fullPage:true});await page.locator('#graph-view').click();
 await page.locator('#search').fill('PINES');assert.equal(await page.locator('.card').count(),1);await page.locator('.graph-node').click();
 await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});
 console.log('Viewer status:',await page.locator('#viewer-status').textContent());console.log('Surface status:',await page.locator('#surface-status').textContent());
 assert.equal(await page.locator('#viewer-status').textContent(),'');assert.equal(await page.locator('#surface-status').textContent(),'');
 assert.ok((await page.locator('#canlab-niivue-controls select').count())>=2);
 assert.equal(await page.locator('#canlab-niivue-controls select').last().inputValue(),'off');assert.equal(await page.locator('#positive-number').inputValue(),'35');assert.equal(await page.locator('#negative-number').inputValue(),'35');
 const slices=await page.locator('#gl').boundingBox(),bar=await page.locator('#weight-colorbar').boundingBox();assert.ok(bar.y>=slices.y+slices.height+4);
 await page.screenshot({path:'/tmp/neuromarker-default-viewer.png',fullPage:true});
 await page.locator('#gl').click({position:{x:160,y:180}});
 await page.waitForFunction(()=>!document.querySelector('#canlab-niivue-readout').textContent.includes('value: ---'));
 console.log('Atlas readout:',await page.locator('#canlab-niivue-readout').textContent());
 assert.ok(await page.locator('#sync').isChecked());await page.locator('#canlab-niivue-controls select').last().selectOption('outline');await page.waitForFunction(()=>new URLSearchParams(location.search).get('atlas')==='outline');
 await page.locator('#sync').uncheck();await page.locator('#positive-number').fill('10');await page.locator('#negative-number').fill('35');
 await page.waitForFunction(()=>location.search.includes('pos=10')&&location.search.includes('neg=35'));
 assert.equal(await page.locator('#positive-number').inputValue(),'10');assert.equal(await page.locator('#negative-number').inputValue(),'35');
 assert.match(await page.locator('#canlab-niivue-readout').textContent(),/ — /);
 const share=page.url();await page.reload();await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});assert.equal(await page.locator('#negative-number').inputValue(),'35');assert.equal(await page.locator('#canlab-niivue-controls select').last().inputValue(),'outline');
 await page.locator('#sign').selectOption('positive');assert.match(await page.locator('#negative-cutoff').textContent(),/hidden/);
 await page.locator('#threshold-mode').selectOption('absolute');assert.equal(await page.locator('#threshold-mode').inputValue(),'absolute');
 await page.locator('#hemisphere').selectOption('L');
 const downloadPromise=page.waitForEvent('download');await page.locator('#download').click();const download=await downloadPromise;await download.saveAs('/tmp/neuromarker-download.nii');assert.ok((await fs.stat('/tmp/neuromarker-download.nii')).size>1000);
 const pngPromise=page.waitForEvent('download');await page.locator('#png').click();const png=await pngPromise;await png.saveAs('/tmp/neuromarker-view.png');assert.ok((await fs.stat('/tmp/neuromarker-view.png')).size>2000);
 await page.screenshot({path:'/tmp/neuromarker-desktop.png',fullPage:true});
 await page.setViewportSize({width:390,height:844});assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=window.innerWidth+1));await page.screenshot({path:'/tmp/neuromarker-mobile.png',fullPage:true});
 await page.goto(share+'&embed=1');await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});assert.ok(await page.locator('header').isHidden());
 // Exercise a different NIfTI container, 4D split, and local replacement.
 for(const query of ['NCS','Geuter','mentalizing']) {
  await page.goto(base);await page.locator('.card').first().waitFor();await page.locator('#search').fill(query);await page.locator('.card').first().click();
  await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});
  assert.equal(await page.locator('#viewer-status').textContent(),'');assert.equal(await page.locator('#surface-status').textContent(),'');
  const options=await page.locator('#map-select option').evaluateAll(os=>os.map(o=>o.value));
  if(options.length>1){await page.locator('#map-select').selectOption(options[1]);await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});assert.equal(await page.locator('#surface-status').textContent(),'');}
 }
 await page.locator('#local-file').setInputFiles('/tmp/neuromarker-download.nii');await page.waitForFunction(()=>document.querySelector('#map-role').textContent==='Local overlay');assert.ok(await page.locator('#share').isDisabled());assert.match(await page.locator('#surface-status').textContent(),/No precomputed/);
 await page.locator('#reset').click();await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});assert.ok(await page.locator('#share').isEnabled());
 await page.goto(base);await page.locator('.card').first().waitFor();await page.locator('#search').fill('nothing-matches-92837');assert.equal(await page.locator('.card').count(),0);
 await page.goto(base+'atlas.html');await page.locator('#atlas-rows tr').first().waitFor();assert.equal(await page.locator('#atlas-rows tr').count(),518);await page.locator('#atlas-search').fill('Ctx_p24_L');assert.equal(await page.locator('#atlas-rows tr').count(),1);assert.match(await page.locator('#atlas-rows').textContent(),/Cortex: Area posterior 24, Left/);
 await page.goto(base);await page.locator('.card').first().waitFor();assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));await page.screenshot({path:'/tmp/neuromarker-home-mobile.png'});await page.locator('.intro').screenshot({path:'/tmp/neuromarker-intro-mobile.png'});
 await page.locator('header nav a').last().click();await page.locator('.review-link').first().waitFor();
 assert.equal(await page.locator('.review-link').count(),5);
 assert.deepEqual(await page.locator('.review-link').evaluateAll(links=>links.map(a=>a.getAttribute('href'))),["https://www.nature.com/articles/nn.4478", "https://www.cell.com/neuron/fulltext/S0896-6273(18)30477-X", "https://www.sciencedirect.com/science/article/pii/S105381191600210X", "https://www.cell.com/neuron/fulltext/S0896-6273(14)00967-2", "https://www.sciencedirect.com/science/article/pii/S1053811908012263"]);
 assert.ok(await page.locator('.review-link img').evaluateAll(images=>images.every(img=>img.complete&&img.naturalWidth>0)));
 assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1));
 await page.screenshot({path:'/tmp/neuromarker-reviews-mobile.png',fullPage:true});
 await page.setViewportSize({width:1440,height:900});await page.screenshot({path:'/tmp/neuromarker-reviews.png',fullPage:true});
 await page.locator('header nav a').first().click();await page.locator('.card').first().waitFor();
 await page.locator('#search').waitFor();const jumpedSearch=await page.locator('#search').boundingBox();assert.ok(jumpedSearch.y>=0&&jumpedSearch.y<200);
 assert.deepEqual(errors,[]);console.log('PASS: catalog, search, volumes, cortical surfaces, independent/synced state, reload, sign controls, absolute mode, hemisphere, NIfTI, PNG, mobile, embed, empty state.');
}catch(error){await page.screenshot({path:'/tmp/neuromarker-test-failure.png',fullPage:true});throw error}finally{await browser.close()}
