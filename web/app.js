const $ = s => document.querySelector(s);
let genes = [], selected = null, original = [], lastSvg = '', suggestTimer, appInfo = null, activeExperiment = 'default';
function decodeText(value){
  try { return decodeURIComponent(String(value || '')); }
  catch (_) { return String(value || ''); }
}

function displayName(g){
  const description = decodeText(g.description).trim();
  const identifiers = [g.name, g.alias, g.gene_id].map(decodeText).join('; ');
  const loc = identifiers.match(/\bLOC\d+\b/i)?.[0];
  if (description) return loc ? `${description}; ${loc}` : description;
  return loc || decodeText(g.name).trim() || String(g.gene_id || '');
}
const fmt = n => n == null ? '—' : Number(n).toLocaleString('es-CL', {maximumFractionDigits: 3});
const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
async function api(url){let r=await fetch(url);if(!r.ok)throw Error((await r.json()).error);return r.json()}
function sampleOrder(){return (appInfo?.samples || []).map(s => s.timepoint).filter((v,i,a)=>v && a.indexOf(v)===i)}
function experimentParam(){return '&experiment=' + encodeURIComponent(activeExperiment)}

async function loadInfo(){
  appInfo = await api('/api/info?experiment=' + encodeURIComponent(activeExperiment));
  activeExperiment = appInfo.active_experiment || activeExperiment;
  const selector = $('#experiment');
  selector.innerHTML = (appInfo.experiments || []).map(e => `<option value="${esc(e.experiment_id)}">${esc(e.name || e.experiment_id)}</option>`).join('');
  selector.value = activeExperiment;
  $('#status').textContent = `${fmt(appInfo.gene_count)} genes · ${appInfo.samples.length} muestras`;
}

loadInfo().catch(e => $('#status').textContent = e.message);
$('#experiment').onchange = async e => {activeExperiment = e.target.value; await loadInfo(); if(genes.length) await search()};
$('#search-button').onclick = search;
$('#query').onkeydown = e => {if(e.key==='Enter'){e.preventDefault();search()} if(e.key==='Escape') $('#results').classList.add('hidden')};
$('#query').oninput = () => {clearTimeout(suggestTimer); suggestTimer = setTimeout(suggest,180)};
document.addEventListener('click', e => {if(!e.target.closest('.hero')) $('#results').classList.add('hidden')});

async function suggest(){
  let input = $('#query').value, term = input.split(/[,;\n]/).pop().trim(), box = $('#results');
  if(term.length < 2){box.classList.add('hidden'); return}
  try{
    let rows = await api('/api/search?q=' + encodeURIComponent(term) + experimentParam());
    if(!rows.length){box.innerHTML = '<div class="result"><small>Sin coincidencias</small></div>'; box.classList.remove('hidden'); return}
    box.innerHTML = rows.map((r,i)=>`<div class="result" data-id="${esc(r.gene_id)}"><b>${i+1}. ${esc(displayName(r))}</b><small>${esc(r.gene_id)} · ${esc(displayName(r))}</small></div>`).join('');
    box.classList.remove('hidden');
    box.querySelectorAll('[data-id]').forEach(el => el.onclick = () => {
      let parts = input.split(/([,;\n])/); parts[parts.length-1] = ' ' + el.dataset.id;
      $('#query').value = parts.join('').trim(); box.classList.add('hidden'); if(!/[,;\n]/.test(input)) search();
    });
  }catch(e){box.classList.add('hidden')}
}

function stats(g){
  let by = {};
  sampleOrder().forEach(t => {
    let v = g.expression.filter(x => x.timepoint === t).map(x => +x.tpm);
    let mean = v.length ? v.reduce((a,b)=>a+b,0)/v.length : null;
    let sd = v.length > 1 ? Math.sqrt(v.reduce((a,b)=>a+(b-mean)**2,0)/(v.length-1)) : 0;
    by[t] = {values:v, mean, sd, n:v.length};
  });
  return by;
}

async function search(){
  let q = $('#query').value.trim(); if(!q) return;
  $('#search-button').disabled = true; $('#results').classList.add('hidden');
  try{
    genes = await api('/api/batch?q=' + encodeURIComponent(q) + experimentParam());
    original = [...genes];
    let requested = q.split(/[,;\n]+/).filter(Boolean).length;
    $('#not-found').classList.toggle('hidden', genes.length === requested);
    $('#not-found').textContent = genes.length ? `Se encontraron ${genes.length} de ${requested} términos. Revisa nombres, alias o descripciones no presentes.` : 'No se encontraron coincidencias.';
    if(genes.length){$('#empty').classList.add('hidden'); $('#detail').classList.remove('hidden'); render()}
  }catch(e){alert(e.message)} finally{$('#search-button').disabled = false}
}

function render(){
  renderSummary();
  $('#multi').classList.toggle('hidden', genes.length < 2);
  if(genes.length > 1){renderHeatmap(); renderMulti()}
  selectGene(selected && genes.find(g => g.gene_id === selected.gene_id) || genes[0]);
}

function renderSummary(){
  const times = sampleOrder();
  $('#summary-head').innerHTML = '<tr><th>Nombre</th><th>ID ITAG4</th><th>Descripción</th>' + times.map(t=>`<th>${esc(t)}</th>`).join('') + '</tr>';
  $('#summary-body').innerHTML = genes.map(g => {
    let s = stats(g);
    return `<tr data-id="${esc(g.gene_id)}"><td><b>${esc(displayName(g))}</b></td><td>${esc(g.gene_id)}</td><td class="desc">${esc(displayName(g))}</td>${times.map(t=>`<td class="num">${fmt(s[t]?.mean)}</td>`).join('')}</tr>`;
  }).join('');
  document.querySelectorAll('#summary-body tr').forEach(r => r.onclick = () => selectGene(genes.find(g => g.gene_id === r.dataset.id)));
}

function selectGene(g){
  selected = g; let s = stats(g), times = sampleOrder();
  $('#gene-id').textContent = g.gene_id;
  $('#gene-name').textContent = displayName(g);
  $('#description').textContent = g.description || 'Sin descripción funcional en el GFF';
  $('#location').textContent = `${g.seqid}:${fmt(g.start)}–${fmt(g.end)} (${g.strand})`;
  $('#max-tpm').textContent = fmt(Math.max(...times.map(t=>s[t]?.mean||0)));
  $('#tx-count').textContent = g.transcripts.length;
  $('#replicates').innerHTML = times.map(t => {
    let rows = g.expression.filter(x => x.timepoint === t);
    return rows.map((r,i)=>`<tr><td>${esc(t)}</td><td>${esc(r.label || r.sample)}</td><td>${esc(r.replicate || i+1)}</td><td class="num">${fmt(r.tpm)}</td><td class="num">${fmt(s[t].mean)}</td><td class="num">${fmt(s[t].sd)}</td><td class="srr">${esc(r.sample)}</td></tr>`).join('');
  }).join('');
  let official = [['Nombre',g.name],['Alias',g.alias],['ID ITAG4',g.gene_id],['GO',g.go],['InterPro',g.interpro],['Pfam',g.pfam],['Dbxref',g.dbxref],['Cromosoma',g.seqid],['Coordenadas',`${g.start}–${g.end}`],['Hebra',g.strand],['Exones',g.exon_count],['Longitud gen',g.gene_length],['Longitud mRNA',g.mrna_length],['Descripción',g.description],['Transcritos',g.transcripts.join(', ')||'No disponible']];
  let curated = (g.curated_annotations || []).map(a => [`${a.source} · ${a.field}`, a.value]);
  $('#annotation').innerHTML = official.concat(curated).filter(x=>x[1]).map(x=>`<div><b>${esc(x[0])}</b>${esc(x[1]||'No disponible')}</div>`).join('');
  lastSvg = lineSvg([{name:displayName(g), stats:s}], true);
  $('#chart').innerHTML = lastSvg;
}

function lineSvg(series, error=false){
  const times = sampleOrder(), W=850,H=350,p={l:70,r:35,t:36,b:62};
  let all = series.flatMap(x => times.map(t => (x.stats[t]?.mean||0)+(x.stats[t]?.sd||0)));
  let max = niceYAxisMax(Math.max(...all, 0)), mid = max/2;
  let x = i => times.length === 1 ? W/2 : p.l+i*(W-p.l-p.r)/(times.length-1);
  let y = v => H-p.b-v/max*(H-p.t-p.b);
  let colors = ['#0f7659','#c1522b','#4467a8','#8b5aa5','#aa8b19','#168a8a','#444'];
  let yGrid = [0, mid, max].map(v => `<line x1="${p.l}" y1="${y(v)}" x2="${W-p.r}" y2="${y(v)}" stroke="#e5e8e2"/><text x="${p.l-8}" y="${y(v)+4}" text-anchor="end">${fmt(v)}</text>`).join('');
  let xGrid = times.map((t,i)=>`<line x1="${x(i)}" y1="${p.t}" x2="${x(i)}" y2="${H-p.b}" stroke="#eef0ea"/><text x="${x(i)}" y="${H-25}" text-anchor="middle">${esc(t)}</text>`).join('');
  let paths = series.map((q,j)=>{
    let pts = times.map((t,i)=>`${x(i)},${y(q.stats[t]?.mean||0)}`).join(' '), c=colors[j%colors.length];
    let bars = error ? times.map((t,i)=>{let m=q.stats[t]?.mean||0, sd=q.stats[t]?.sd||0; return `<line x1="${x(i)}" y1="${y(Math.max(0,m-sd))}" x2="${x(i)}" y2="${y(m+sd)}" stroke="${c}"/><line x1="${x(i)-5}" y1="${y(Math.max(0,m-sd))}" x2="${x(i)+5}" y2="${y(Math.max(0,m-sd))}" stroke="${c}"/><line x1="${x(i)-5}" y1="${y(m+sd)}" x2="${x(i)+5}" y2="${y(m+sd)}" stroke="${c}"/>`}).join('') : '';
    let dots = times.map((t,i)=>(q.stats[t]?.values||[]).map((v,k)=>`<circle cx="${x(i)+(k-1)*5}" cy="${y(v)}" r="3" fill="${c}" opacity=".65"/>`).join('')).join('');
    return `<polyline points="${pts}" fill="none" stroke="${c}" stroke-width="3"/>${bars}${dots}<text x="${p.l+10}" y="${18+j*18}" fill="${c}">${esc(q.name)}</text>`;
  }).join('');
  return `<svg class="svg-chart" viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg"><style>text{font:12px Verdana;fill:#53605b}</style>${yGrid}<line x1="${p.l}" y1="${p.t}" x2="${p.l}" y2="${H-p.b}" stroke="#59645f"/>${xGrid}${paths}<text transform="translate(18 ${H/2}) rotate(-90)" text-anchor="middle">TPM</text></svg>`;
}

function niceYAxisMax(value){
  if (!(value > 0)) return 1;
  const target = value * 1.1;
  const exponent = Math.floor(Math.log10(target));
  const magnitude = 10 ** exponent;
  const fraction = target / magnitude;
  const niceFraction = fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10;
  return niceFraction * magnitude;
}

function renderMulti(){
  let s = genes.map(g => ({name:displayName(g), stats:stats(g)}));
  $('#multi-chart').innerHTML = lineSvg(s);
  $('#legend').innerHTML = s.map((g,i)=>`<span><i style="background:${['#0f7659','#c1522b','#4467a8','#8b5aa5','#aa8b19','#168a8a','#444'][i%7]}"></i>${esc(g.name)}</span>`).join('');
}

function renderHeatmap(){
  const times = sampleOrder();
  let maxima = genes.map(g => Math.max(...times.map(t=>stats(g)[t]?.mean||0), 1));
  $('#heatmap').innerHTML = `<div class="heat-cell"><span></span>${times.map(t=>`<span><b>${esc(t)}</b></span>`).join('')}</div>` + genes.map((g,j)=>{
    let s = stats(g);
    return `<div class="heat-cell"><span class="heat-gene">${esc(displayName(g))}</span>${times.map(t=>{let v=s[t]?.mean; let a=v==null?0:v/maxima[j]; return `<span title="${fmt(v)} TPM" style="background:${v==null?'#eee':`rgba(15,118,89,${.12+.88*a})`};color:${a>.55?'white':'#162621'}">${fmt(v)}</span>`}).join('')}</div>`;
  }).join('') + '<p class="heat-legend">Color más intenso = mayor TPM promedio dentro de ese gen. Columnas según metadatos del experimento activo.</p>';
}

$('#sort').onchange = e => {let mode=e.target.value, times=sampleOrder(); if(mode==='manual') genes=[...original]; if(mode==='max') genes.sort((a,b)=>Math.max(...times.map(t=>stats(b)[t]?.mean||0))-Math.max(...times.map(t=>stats(a)[t]?.mean||0))); if(mode==='fold') genes.sort((a,b)=>fold(b)-fold(a)); if(mode==='similarity') genes.sort((a,b)=>peak(a)-peak(b)); render()};
const fold = g => {let times=sampleOrder(), s=stats(g), a=s[times[0]]?.mean||.0001; return Math.max(...times.map(t=>s[t]?.mean||0))/a};
const peak = g => {let times=sampleOrder(), s=stats(g), values=times.map(t=>s[t]?.mean||0); return values.indexOf(Math.max(...values))};
function tableRows(){let times=sampleOrder(); return [['Nombre','ID ITAG4','Descripción',...times],...genes.map(g=>{let s=stats(g); return [displayName(g),g.gene_id,displayName(g),...times.map(t=>s[t]?.mean)]})]}
function blob(text,type,name){let a=document.createElement('a');a.href=URL.createObjectURL(new Blob([text],{type}));a.download=name;a.click();URL.revokeObjectURL(a.href)}
$('#csv').onclick=()=>blob('\ufeff'+tableRows().map(r=>r.map(v=>'"'+String(v??'').replaceAll('"','""')+'"').join(',')).join('\n'),'text/csv','expresion_transcriptomica.csv');
$('#excel').onclick=()=>blob(`<table>${tableRows().map(r=>'<tr>'+r.map(v=>`<td>${esc(v)}</td>`).join('')+'</tr>').join('')}</table>`,'application/vnd.ms-excel','expresion_transcriptomica.xls');
$('#svg').onclick=()=>blob(lastSvg,'image/svg+xml','perfil_expresion.svg');
$('#png').onclick=()=>{let img=new Image();img.onload=()=>{let c=document.createElement('canvas');c.width=1700;c.height=700;c.getContext('2d').drawImage(img,0,0,c.width,c.height);c.toBlob(b=>{let u=URL.createObjectURL(b),a=document.createElement('a');a.href=u;a.download='perfil_expresion.png';a.click();URL.revokeObjectURL(u)})};img.src='data:image/svg+xml;charset=utf-8,'+encodeURIComponent(lastSvg)};
$('#pdf').onclick=()=>window.print();
