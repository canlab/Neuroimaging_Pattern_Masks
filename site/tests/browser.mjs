import {chromium} from '@playwright/test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
const browser=await chromium.launch({...(process.env.CI?{}:{channel:'chrome'}),headless:true,args:['--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
const page=await browser.newPage({viewport:{width:1440,height:1000}});const errors=[];
page.on('pageerror',e=>errors.push(e.message));page.on('console',m=>{if(m.type()==='error')console.log('BROWSER:',m.text().slice(0,400))});
try{
 await page.goto('http://127.0.0.1:8765/');await page.locator('.card').first().waitFor();
 assert.equal(await page.locator('.card').count(),25);
 await page.locator('#search').fill('PINES');assert.equal(await page.locator('.card').count(),1);await page.locator('.card').click();
 await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});
 console.log('Viewer status:',await page.locator('#viewer-status').textContent());console.log('Surface status:',await page.locator('#surface-status').textContent());
 assert.equal(await page.locator('#viewer-status').textContent(),'');assert.equal(await page.locator('#surface-status').textContent(),'');
 assert.ok((await page.locator('#canlab-niivue-controls select').count())>=2);
 await page.locator('#gl').click({position:{x:160,y:180}});
 await page.waitForFunction(()=>!document.querySelector('#canlab-niivue-readout').textContent.includes('value: ---'));
 console.log('Atlas readout:',await page.locator('#canlab-niivue-readout').textContent());
 await page.locator('#sync').uncheck();await page.locator('#positive-number').fill('10');await page.locator('#negative-number').fill('35');
 await page.waitForFunction(()=>location.search.includes('pos=10')&&location.search.includes('neg=35'));
 assert.equal(await page.locator('#positive-number').inputValue(),'10');assert.equal(await page.locator('#negative-number').inputValue(),'35');
 const share=page.url();await page.reload();await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});assert.equal(await page.locator('#negative-number').inputValue(),'35');
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
  await page.goto('http://127.0.0.1:8765/');await page.locator('.card').first().waitFor();await page.locator('#search').fill(query);await page.locator('.card').first().click();
  await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});
  assert.equal(await page.locator('#viewer-status').textContent(),'');assert.equal(await page.locator('#surface-status').textContent(),'');
  const options=await page.locator('#map-select option').evaluateAll(os=>os.map(o=>o.value));
  if(options.length>1){await page.locator('#map-select').selectOption(options[1]);await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});assert.equal(await page.locator('#surface-status').textContent(),'');}
 }
 await page.locator('#local-file').setInputFiles('/tmp/neuromarker-download.nii');await page.waitForFunction(()=>document.querySelector('#map-role').textContent==='Local overlay');assert.ok(await page.locator('#share').isDisabled());assert.match(await page.locator('#surface-status').textContent(),/No precomputed/);
 await page.locator('#reset').click();await page.waitForFunction(()=>!document.querySelector('#map-select').disabled,{},{timeout:120000});assert.ok(await page.locator('#share').isEnabled());
 await page.goto('http://127.0.0.1:8765/');await page.locator('.card').first().waitFor();await page.locator('#search').fill('nothing-matches-92837');assert.equal(await page.locator('.card').count(),0);
 assert.deepEqual(errors,[]);console.log('PASS: catalog, search, volumes, cortical surfaces, independent/synced state, reload, sign controls, absolute mode, hemisphere, NIfTI, PNG, mobile, embed, empty state.');
}finally{await browser.close()}
