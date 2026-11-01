/* ScerTF 2012 database archive: self-contained browser search using local recommended.json.
 * Matrix alignment, log-odds scoring, and Bonferroni correction adapted from the original
 * Python 2 common.py, with defensive input validation and modern JavaScript.
 */
(() => {
  'use strict';
  const $ = (selector, scope=document) => scope.querySelector(selector);
  const root = document.body.dataset.root || './';
  const escapeHTML = str => String(str ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const bases=['A','C','G','T'];
  const bg=[.32,.18,.18,.32];
  let dataPromise;
  function getData(){ return dataPromise ||= fetch(root+'data/recommended.json').then(r => {if(!r.ok)throw Error('Could not load the archival matrix file.'); return r.json()}); }
  function panelError(message,info=false){const el=$('#search-message');if(el)el.innerHTML=`<div class="${info?'search-info':'search-error'}">${escapeHTML(message)}</div>`;}
  function clearResults(){const el=$('#search-results');if(el){el.hidden=true;el.innerHTML='';}const msg=$('#search-message');if(msg)msg.innerHTML='';}
  function fmtP(p){if(p===null || !Number.isFinite(p))return 'Not available';if(p===0)return '< 1e-16';return p<.0001?p.toExponential(2):p.toFixed(4);}
  const formulaMatrix=`A | 0.515 0.030 0.064 0.687 0.044 0.117 0.382\nC | 0.087 0.028 0.885 0.056 0.561 0.182 0.178\nG | 0.317 0.045 0.024 0.190 0.024 0.164 0.138\nT | 0.081 0.897 0.027 0.067 0.371 0.537 0.303`;
  function parseMatrix(text,allowCounts=true){
    const rows={};
    for(const line of String(text).trim().split(/\r?\n/)){
      if(!line.trim())continue;
      const m=line.match(/^\s*([ACGT])\s*\|?\s+(.+)$/i);
      if(!m)throw Error('Matrix rows must begin with A, C, G or T followed by values (optionally with |).');
      const letter=m[1].toUpperCase(); if(rows[letter])throw Error('Each nucleotide must have exactly one row.');
      const values=m[2].trim().split(/\s+/).map(Number);
      if(!values.length || values.some(v=>!Number.isFinite(v)))throw Error('Matrix contains invalid numeric values.');
      rows[letter]=values;
    }
    if(bases.some(b=>!rows[b]))throw Error('Matrix needs four nucleotide rows: A, C, G and T.');
    const len=rows.A.length;
    if(len<4||len>100||bases.some(b=>rows[b].length!==len))throw Error('Matrix must have 4–100 equally sized columns.');
    const matrix=bases.map(b=>rows[b].slice());
    if(!allowCounts) return matrix;
    const hasNegative=matrix.some(row=>row.some(v=>v<0));
    if(hasNegative){
      // Legacy PWM logs used distinct AT/CG background frequencies inferred per column.
      const inferred=[];
      for(let j=0;j<len;j++){
        const lat=Math.exp(matrix[0][j])+Math.exp(matrix[3][j]);
        const lcg=Math.exp(matrix[1][j])+Math.exp(matrix[2][j]);
        const den=lat-lcg;
        if(Math.abs(den)>1e-12){const p=(1-.5*lcg)/den;if(p>=0&&p<=.5)inferred.push(p);}
      }
      const at=inferred.length?inferred.reduce((s,n)=>s+n,0)/inferred.length:.32;
      for(let j=0;j<len;j++){
        let sum=0;
        for(let i=0;i<4;i++){matrix[i][j]=Math.exp(matrix[i][j])*(i===0||i===3?at:.5-at);sum+=matrix[i][j];}
        if(sum<=0)throw Error('Could not normalize the position weight matrix.');
        for(let i=0;i<4;i++)matrix[i][j]=Math.max(.001,matrix[i][j]/sum);
      }
    }else{
      for(let j=0;j<len;j++){
        const sum=matrix.reduce((s,row)=>s+row[j],0);
        if(sum<=0)throw Error('Every matrix column must have a positive sum.');
        for(let i=0;i<4;i++)matrix[i][j]=Math.max(.001,matrix[i][j]/sum);
      }
    }
    return matrix;
  }
  function matrixFromSequence(text){const seq=text.toUpperCase();return bases.map(base=>[...seq].map(char=>char===base?.85:.05));}
  function reverseComplementSeq(seq){return [...seq].reverse().map(n=>({A:'T',T:'A',C:'G',G:'C'})[n]).join('');}
  function reverseMatrix(matrix){return matrix.map((row,i)=>matrix[3-i].slice().reverse());}
  // Use ungapped local overlaps >=4 bp (same boundary-search objective as legacy align).
  function align(a,b){
    const m=a[0].length,n=b[0].length;
    let best=null;
    for(let offset=-(n-4);offset<=m-4;offset++){
      const iStart=Math.max(0,offset),jStart=Math.max(0,-offset);
      const length=Math.min(m-iStart,n-jStart);
      if(length<4)continue;
      let score=0;
      for(let k=0;k<length;k++){
        const i=iStart+k,j=jStart+k;
        for(let base=0;base<4;base++)score+=(a[base][i]*b[base][j]+b[base][j]*a[base][i])/(2*bg[base]);
      }
      if(!best || score>best.score)best={score,iStart,jStart,length};
    }
    if(!best)throw Error('The query matrix must share at least four positions with each reference.');
    return best;
  }
  // Standard-normal CDF formula used by the original Scertf common.py.
  function phi(z){const sign=z<0?-1:1;const x=Math.abs(z)/Math.sqrt(2),p=.3275911;
    const t=1/(1+p*x);
    const y=1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-.284496736)*t+.254829592)*t*Math.exp(-x*x);
    return .5*(1+sign*y);
  }
  function matrixPvalue(a,b,alignResult,sequence=null){
    const {iStart,jStart,length}=alignResult;
    let mean=0,variance=0,score=0;
    for(let j=0;j<length;j++){
      const logodds=bases.map((_,i)=>Math.log(Math.max(1e-9,b[i][jStart+j])/bg[i]));
      const ave=logodds.reduce((s,v)=>s+v,0)/4;
      mean+=ave;
      variance+=logodds.reduce((s,v)=>s+(v-ave)**2,0)/4;
      if(sequence)score+=logodds[bases.indexOf(sequence[iStart+j])];
      else for(let i=0;i<4;i++)score+=a[i][iStart+j]*logodds[i];
    }
    return variance>0?Math.max(0,Math.min(1,1-phi((score-mean)/Math.sqrt(variance)))):1;
  }
  function compareMotifs(query,seq,items){
    const results=[];
    for(const item of items){
      let ref;
      try{ref=parseMatrix(item.PFM,false);}catch(e){continue;}
      const fwd=align(query,ref),revQuery=reverseMatrix(query),rev=align(revQuery,ref);
      const reverse=rev.score>fwd.score;
      const chosen=reverse?rev:fwd;
      const testSeq=seq?(reverse?reverseComplementSeq(seq):seq):null;
      const p=matrixPvalue(reverse?revQuery:query,ref,chosen,testSeq);
      results.push({...item,probability:p*items.length,uncorrected:p,partial:chosen.length<ref[0].length});
    }
    results.sort((a,b)=>a.probability-b.probability);
    const significant=results.filter(x=>x.probability<=.01);
    return {results:significant.length?significant:results.slice(0,1).map(x=>({...x,probability:x.uncorrected,unadjusted:true})),found:significant.length>0};
  }
  function renderRows(results,title,annotation){
    const section=$('#search-results');
    if(!section)return;
    section.hidden=false;
    section.innerHTML=`<div class="results-top"><h3>${escapeHTML(title)}</h3><p class="helper muted">${escapeHTML(annotation||'Click a factor to inspect the original data and candidate matrices.')}</p></div>`+
       results.map(item=>`<div class="result-row"><div><h4><a href="${root}details/${encodeURIComponent(item.name)}/">${escapeHTML(item.name)} ↗</a></h4><div class="result-meta">${escapeHTML(item.source)}</div></div><div>${item.logofile?`<a href="${root}details/${encodeURIComponent(item.name)}/"><img class="motif" loading="lazy" src="${root}assets/${encodeURI(item.logofile)}" alt="Sequence logo for ${escapeHTML(item.name)}"></a>`:''}</div><div class="result-desc"><p>${escapeHTML((item.description||'Function unknown').slice(0,280))}</p>${item.probability!==undefined?`<div class="result-meta" style="margin-top:8px"><strong>${item.unadjusted?'Unadjusted P-value (no significant hits):':'Bonferroni-adjusted P-value:'}</strong> ${fmtP(item.probability)}${item.partial?' · Partial motif overlap':''}</div>`:`<div class="result-meta" style="margin-top:8px">Aliases: ${escapeHTML(item.alias||'None')}</div>`}</div></div>`).join('');
    section.scrollIntoView({behavior:'smooth',block:'start'});
  }
  async function submitName(query){
    clearResults();const normalized=query.replace(/\s+/g,'').toUpperCase();
    if(!normalized){panelError('Enter a TF name, systematic gene name or DNA sequence.');return;}
    if(normalized.length>20){panelError('Search input is limited to 20 characters.');return;}
    try{
      const items=await getData();
      let match=items.filter(x=>x.name.toUpperCase()===normalized);
      if(!match)match=[];
      if(!match.length)match=items.filter(x=>x.name.toUpperCase().includes(normalized));
      if(!match.length)match=items.filter(x=>x.alias.toUpperCase().includes(normalized));
      if(match.length){renderRows(match,`${match.length} TF ${match.length===1?'record':'records'} matching “${normalized}”`);return;}
      if(!/^[ACGT]+$/.test(normalized)){panelError('No transcription factor matches your query. For DNA search use only A, C, G and T.');return;}
      if(normalized.length<4){panelError('DNA sequence search requires at least four bases.');return;}
      const {results,found}=compareMotifs(matrixFromSequence(normalized),normalized,items);
      if(!results.length){panelError('No motif comparisons were available.');return;}
      if(!found)panelError('No matches passed the original 1% Bonferroni-adjusted significance threshold; displaying the highest-scoring candidate.',true);
      renderRows(results,`${results.length} ${found?'matching':'top-ranked'} DNA-binding ${results.length===1?'motif':'motifs'}`,`DNA sequence: ${normalized}. P-values reflect the restored 2012 normal-approximation test and 196-comparison Bonferroni correction; matches are not evidence of occupancy.`);
    }catch(error){panelError(error.message||'Unable to search archived data.');}
  }
  async function submitMatrix(source){
    clearResults();
    try{
      const matrix=parseMatrix(source);
      const consensus=Array.from({length:matrix[0].length},(_,j)=>bases[matrix.map(row=>row[j]).indexOf(Math.max(...matrix.map(row=>row[j])))]).join('');
      const items=await getData();
      const {results,found}=compareMotifs(matrix,null,items);
      if(!results.length){panelError('No motif comparisons were available.');return;}
      if(!found)panelError('No matches passed the original 1% Bonferroni-adjusted significance threshold; displaying the highest-scoring candidate.',true);
      renderRows(results,`${results.length} ${found?'matching':'top-ranked'} matrix ${results.length===1?'motif':'motifs'}`,`Query consensus: ${consensus}. Search uses historical aligned motif scores and normal-approximation P-values adjusted for multiple comparisons.`);
    }catch(error){panelError(error.message||'Invalid matrix.');}
  }
  function showTab(tab){const isMatrix=tab==='matrix';$('#tab-matrix').classList.toggle('active',isMatrix);$('#tab-sequence').classList.toggle('active',!isMatrix);$('#tab-matrix').setAttribute('aria-selected',String(isMatrix));$('#tab-sequence').setAttribute('aria-selected',String(!isMatrix));$('#panel-matrix').hidden=!isMatrix;$('#panel-sequence').hidden=isMatrix;clearResults();}
  function setupSearch(){
    if(!$('#name-search'))return;
    $('#tab-sequence').addEventListener('click',()=>showTab('sequence'));
    $('#tab-matrix').addEventListener('click',()=>showTab('matrix'));
    $('#example-seq').addEventListener('click',()=>{$('#q').value='TGACTC';submitName('TGACTC');});
    $('#example-matrix').addEventListener('click',()=>{$('#matrix-input').value=formulaMatrix;$('#matrix-input').focus();});
    $('#name-search').addEventListener('submit',ev=>{ev.preventDefault();const q=$('#q').value.trim();history.replaceState(null,'',location.pathname+(q?'?q='+encodeURIComponent(q):''));submitName(q);});
    $('#matrix-search').addEventListener('submit',ev=>{ev.preventDefault();submitMatrix($('#matrix-input').value);});
    const params=new URLSearchParams(location.search),q=params.get('q');
    if(q){$('#q').value=q;submitName(q);}
  }
  function setupBrowse(){const input=$('#browse-filter');if(!input)return;input.addEventListener('input',()=>{let shown=0;const term=input.value.toLowerCase().trim();document.querySelectorAll('.tf-card').forEach(card=>{const display=card.dataset.find.includes(term);card.hidden=!display;if(display)shown++;});$('#browse-count').textContent=`Showing ${shown} factors`;$('#browse-empty').hidden=shown!==0;});}
  function setupMatrices(){document.querySelectorAll('.copy-matrix').forEach(button=>button.addEventListener('click',async()=>{
    const pre=button.closest('.matrix-view').querySelector('pre');const content=pre.textContent;
    try{await navigator.clipboard.writeText(content);button.textContent='Copied';setTimeout(()=>{button.textContent='Copy';},1300);}catch(e){button.textContent='Select text to copy';}
  }));}
  // Expose a read-only helper only for local functional tests; never required by the UI.
  if(typeof window!=='undefined')window.ScerTFSearch={parseMatrix,matrixFromSequence,compareMotifs};
  setupSearch();setupBrowse();setupMatrices();
})();
