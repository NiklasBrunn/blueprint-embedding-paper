/* Client-side recolor engine for the blueprint UMAP explorer.
   Ported from the standalone interactive-UMAP HTML. Recolors a single
   scattergl trace entirely in-browser (branch scrubber / single axis /
   metadata / gene expression). No server roundtrip except gene vectors,
   which are requested via a Dash clientside Store handshake. */

/* ---- module-scoped state (re-assigned inside initExplore) ---- */
let DATA = null;

let plot = null;
let modeSel = null;
let branchSel = null;
let paletteSel = null;
let timeSlider = null;
let stateLine = null;
let ticksDiv = null;
let sidebar = null;
let branchGroup = null;
let paletteGroup = null;
let timeGroup = null;
let catLegend = null;
let catTitle = null;
let footNote = null;

let currentMode = "branch";
let currentBranch = null;
let currentAxis = null;
let currentMetaCol = null;
let currentGene = null;
let currentPaletteKey = "viridis";
/* cache of the last fetched gene vector so palette changes / mode re-entry
   recolor locally without a redundant server round-trip. */
let lastGene = null;
let lastGeneValues = null;
/* the gene we're currently asking the plot to show (via the gene picker OR a
   clickable role chip in the axis/branch sidebar). A late server vector only
   applies if it still matches pendingGene; any axis/branch/metadata recolor
   clears it so a stale gene vector can't override the active coloring. */
let pendingGene = null;

/* Request a gene's expression vector from the server. We can't use
   window.dash_clientside.set_props here — it updates props WITHOUT triggering
   the callback graph, so the server _gene callback would never fire. Instead we
   drive the hidden dcc.Input through a real input event (native value setter so
   React/Dash registers the change), which fires _gene -> gene-vector store ->
   applyGene. */
function requestGene(name){
  const el=document.getElementById("exp-gene-select");
  if(!el)return;
  try{
    const d=Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype,"value");
    if(d&&d.set){d.set.call(el,name);}else{el.value=name;}
  }catch(e){el.value=name;}
  el.dispatchEvent(new Event("input",{bubbles:true}));
  el.dispatchEvent(new Event("change",{bubbles:true}));
}

/* Recolor by the current gene: reuse the cached vector if we already have it,
   else ask the server. */
function showGene(name){
  if(!name)return;
  pendingGene=name;
  if(lastGene===name&&lastGeneValues){
    applyGeneVector({gene:name,values:lastGeneValues});
  }else{
    requestGene(name);
  }
}

/* Click handler (delegated on the sidebar) for the role gene chips shown when
   an axis is selected: recolor the UMAP by that gene WITHOUT leaving the axis
   view, and highlight the clicked chip. */
function onSidebarClick(e){
  const chip=e.target&&e.target.closest?e.target.closest(".gene-chip[data-gene]"):null;
  if(!chip)return;
  const g=chip.getAttribute("data-gene");
  if(!g)return;
  sidebar.querySelectorAll(".gene-chip.gene-chip-on")
    .forEach(c=>c.classList.remove("gene-chip-on"));
  chip.classList.add("gene-chip-on");
  currentGene=g;
  showGene(g);
}

/* ---- small html helpers ---- */
function escapeHtml(s){return String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");}

function geneChipsHtml(moduleName){const genes=(DATA&&DATA.modules&&DATA.modules[moduleName])||[];if(genes.length===0)return '<div class="no-genes">(no gene list available)</div>';return '<div class="gene-chips">'+genes.map(g=>{const e=escapeHtml(g);const on=(g===pendingGene)?' gene-chip-on':'';return '<span class="gene-chip'+on+'" data-gene="'+e+'" title="color UMAP by log1p('+e+')">'+e+'</span>';}).join("")+'</div>';}

function roleBlockHtml(role,moduleName){if(!moduleName)return "";return '<div class="role '+role+'"><div class="role-heading">'+role+'</div><div class="role-mod">'+escapeHtml(moduleName)+'</div>'+geneChipsHtml(moduleName)+'</div>';}

function negativeBlockHtml(modules){if(!modules||modules.length===0)return "";const items=modules.map(nm=>'<div class="role-mod">'+escapeHtml(nm)+'</div>'+geneChipsHtml(nm)).join("");return '<div class="role negative"><div class="role-heading">negative</div>'+items+'</div>';}

function axisInfoHtml(axisName){const meta=(DATA&&DATA.axisMeta&&DATA.axisMeta[axisName])||{};const parts=[];parts.push('<h2>'+escapeHtml(axisName)+'</h2>');parts.push('<div class="kind-badge">'+escapeHtml(meta.kind||"axis")+'</div>');parts.push(roleBlockHtml("positive",meta.positive));if(meta.gate)parts.push(roleBlockHtml("gate",meta.gate));if(meta.negatives&&meta.negatives.length)parts.push(negativeBlockHtml(meta.negatives));return parts.join("");}

function renderBranchSidebar(axisName){sidebar.innerHTML=axisInfoHtml(axisName);}

/* generic searchable picker builder shared by axis/meta/gene sidebars */
function buildPicker(opts){
  /* opts: {title, items, selected, onPick, cap, infoHtml} */
  const wrap=document.createElement("div");
  const title=document.createElement("div");
  title.className="picker-title";
  title.textContent=opts.title;
  wrap.appendChild(title);

  const search=document.createElement("input");
  search.className="picker-search";
  search.type="text";
  search.placeholder="filter…";
  wrap.appendChild(search);

  const picker=document.createElement("div");
  picker.className="picker";
  wrap.appendChild(picker);

  const note=document.createElement("div");
  note.className="sidebar-note";
  note.style.display="none";
  wrap.appendChild(note);

  const cap=opts.cap||10000;
  function renderList(filter){
    picker.innerHTML="";
    const f=(filter||"").toLowerCase();
    const matched=opts.items.filter(name=>name.toLowerCase().indexOf(f)!==-1);
    const shown=matched.slice(0,cap);
    shown.forEach(name=>{
      const item=document.createElement("div");
      item.className="picker-item"+(name===opts.selected?" selected":"");
      item.textContent=name;
      item.title=name;
      item.addEventListener("click",()=>opts.onPick(name));
      picker.appendChild(item);
    });
    if(matched.length>shown.length){
      note.style.display="";
      note.textContent="Showing "+shown.length+" of "+matched.length+" — refine the filter to narrow.";
    } else {
      note.style.display="none";
      note.textContent="";
    }
  }
  search.addEventListener("input",()=>renderList(search.value));
  renderList("");

  sidebar.innerHTML="";
  sidebar.appendChild(wrap);

  if(opts.infoHtml){
    const divider=document.createElement("div");
    divider.className="picker-divider";
    sidebar.appendChild(divider);
    const info=document.createElement("div");
    info.innerHTML=opts.infoHtml;
    sidebar.appendChild(info);
  }
}

function renderAxisSidebar(selectedAxis){
  const axes=(DATA&&DATA.axisOrder)||[];
  buildPicker({
    title:"axes",
    items:axes,
    selected:selectedAxis,
    onPick:(name)=>{currentAxis=name;renderAxisSidebar(name);recolorAxis();},
    infoHtml:selectedAxis?axisInfoHtml(selectedAxis):""
  });
}

function renderMetaSidebar(selectedCol){
  const cols=(DATA&&DATA.metadata)?Object.keys(DATA.metadata):[];
  buildPicker({
    title:"metadata columns",
    items:cols,
    selected:selectedCol,
    onPick:(name)=>{currentMetaCol=name;renderMetaSidebar(name);recolorMeta(name);}
  });
}

function renderGeneSidebar(selected){
  const genes=(DATA&&DATA.geneList)||[];
  buildPicker({
    title:"genes",
    items:genes,
    selected:selected,
    cap:500,
    onPick:(name)=>{
      currentGene=name;
      renderGeneSidebar(name);
      showGene(name);
    }
  });
}

/* ---- palettes ---- */
const PALETTES={
  viridis:{colorscale:"Viridis",cmin:-1.0,cmax:2.5,barTitle:"activation (z)",tickAccent:"68, 1, 84"},
  divbo:{colorscale:[[0.0,"#2166ac"],[0.25,"#67a9cf"],[0.5,"#f7f7f7"],[0.75,"#fdb863"],[1.0,"#b35806"]],cmin:-2.5,cmax:2.5,barTitle:"activation (z, signed)",tickAccent:"179, 88, 6"}
};

function applyPalette(key){currentPaletteKey=key;const p=PALETTES[key];Plotly.restyle(plot,{"marker.colorscale":[p.colorscale],"marker.cmin":[p.cmin],"marker.cmax":[p.cmax],"marker.colorbar.title.text":[p.barTitle]},[0]);}

/* ---- modes ---- */
const MODES=[
  {key:"branch",label:"Branch (scrubber)"},
  {key:"axis",label:"Single axis activation"},
  {key:"meta",label:"Metadata"},
  {key:"gene",label:"Gene expression"}
];

const FOOTNOTES={
  branch:"Color = convex combination (1−α)·score(axis_i) + α·score(axis_{i+1}),  i = ⌊t⌋,  α = t−i.  Drag t (or click a timepoint) to scrub; the sidebar snaps to the dominant axis.",
  axis:"Single axis activation: search and click an axis on the left.",
  meta:"Metadata: categorical coloring of cells by the selected column.",
  gene:"Gene expression: log1p of the selected gene; search and click a gene on the left."
};

function shortName(axn){const m=axn.match(/^([A-Za-z][A-Za-z0-9]*)_(.*)$/);return m?m[2]:axn;}

function highlightTicks(i,alpha){const accent=PALETTES[currentPaletteKey].tickAccent;const chips=ticksDiv.querySelectorAll(".tick");chips.forEach((chip,k)=>{let w=0;if(k===i)w=1-alpha;else if(k===i+1)w=alpha;chip.style.backgroundColor=w===0?"#eeeeee":"rgba("+accent+", "+(0.20+w*0.80).toFixed(3)+")";chip.style.color=w>0.5?"white":"#333";chip.style.fontWeight=w>0.5?"600":"400";});}

function renderCategoricalLegend(col,meta){catTitle.textContent=col+"  ("+meta.kind+", "+meta.labels.length+" categories)";catTitle.classList.add("visible");catLegend.innerHTML="";meta.labels.forEach((lbl,i)=>{const item=document.createElement("span");item.className="cat-item";item.innerHTML='<span class="cat-swatch" style="background:'+meta.colors[i]+'"></span>'+escapeHtml(lbl);catLegend.appendChild(item);});catLegend.classList.add("visible");}

function hideCategoricalLegend(){catLegend.classList.remove("visible");catTitle.classList.remove("visible");}

/* ---- recolor ops ---- */
function recolorBranch(){if(!currentBranch)return;const t=parseFloat(timeSlider.value);const K=currentBranch.axes.length;const tClamp=Math.max(0,Math.min(K-1,t));const i=Math.min(K-2,Math.floor(tClamp));const alpha=tClamp-i;const lo=DATA.values[currentBranch.axes[i]];const hi=DATA.values[currentBranch.axes[i+1]];const n=lo.length;const interp=new Float32Array(n);for(let k=0;k<n;k++)interp[k]=(1-alpha)*lo[k]+alpha*hi[k];pendingGene=null;const pb=PALETTES[currentPaletteKey];Plotly.restyle(plot,{"marker.color":[Array.from(interp)],"marker.colorscale":[pb.colorscale],"marker.cmin":[pb.cmin],"marker.cmax":[pb.cmax],"marker.showscale":[true],"marker.colorbar.title.text":[pb.barTitle]},[0]);highlightTicks(i,alpha);const dominant=alpha<0.5?i:i+1;renderBranchSidebar(currentBranch.axes[dominant]);stateLine.textContent="t="+t.toFixed(2)+"   axis_i = "+currentBranch.axes[i]+"   axis_{i+1} = "+currentBranch.axes[i+1]+"   α = "+alpha.toFixed(2);}

function recolorAxis(){const arr=DATA.values[currentAxis];if(!arr)return;pendingGene=null;const p=PALETTES[currentPaletteKey];Plotly.restyle(plot,{"marker.color":[arr],"marker.colorscale":[p.colorscale],"marker.cmin":[p.cmin],"marker.cmax":[p.cmax],"marker.showscale":[true],"marker.colorbar.title.text":[p.barTitle]},[0]);stateLine.textContent="axis = "+currentAxis;}

function recolorMeta(col){const meta=DATA.metadata[col];if(!meta)return;pendingGene=null;const idx=meta.cellIdx;const colors=meta.colors;const n=idx.length;const out=new Array(n);for(let k=0;k<n;k++)out[k]=colors[idx[k]]||"#d0d0d0";Plotly.restyle(plot,{"marker.color":[out],"marker.showscale":[false]},[0]);renderCategoricalLegend(col,meta);stateLine.textContent="metadata = "+col;}

/* ---- ticks / branch selection ---- */
function renderTicks(){if(!currentBranch)return;const axes=currentBranch.axes;ticksDiv.innerHTML="";axes.forEach((axn,k)=>{if(k>0){const arrow=document.createElement("span");arrow.className="tick-arrow";arrow.textContent="→";ticksDiv.appendChild(arrow);}const chip=document.createElement("span");chip.className="tick";chip.textContent="t="+k+" · "+shortName(axn);chip.title=axn;chip.addEventListener("click",()=>{timeSlider.value=String(k);recolorBranch();});ticksDiv.appendChild(chip);});}

function selectBranch(idx){currentBranch=DATA.branches[idx];timeSlider.min="0";timeSlider.max=String(currentBranch.axes.length-1);timeSlider.step="0.01";timeSlider.value="0";renderTicks();recolorBranch();}

/* ---- mode switching ---- */
function applyMode(modeKey){
  currentMode=modeKey;
  const isBranch=modeKey==="branch";
  const isAxis=modeKey==="axis";
  const isMeta=modeKey==="meta";
  const isGene=modeKey==="gene";

  /* control groups */
  if(branchGroup)branchGroup.classList.toggle("hidden",!isBranch);
  if(timeGroup)timeGroup.classList.toggle("hidden",!isBranch);
  if(ticksDiv)ticksDiv.classList.toggle("hidden",!isBranch);
  /* palette applies to continuous colorings (branch / axis / gene) */
  if(paletteGroup)paletteGroup.classList.toggle("hidden",isMeta);

  if(isMeta){
    /* categorical legend shown only in meta mode */
  } else {
    hideCategoricalLegend();
  }

  if(footNote)footNote.textContent=FOOTNOTES[modeKey]||"";

  if(isBranch){
    if(currentBranch){renderTicks();recolorBranch();}
  } else if(isAxis){
    if(!currentAxis){const axes=(DATA&&DATA.axisOrder)||[];if(axes.length)currentAxis=axes[0];}
    renderAxisSidebar(currentAxis);
    if(currentAxis)recolorAxis();
  } else if(isMeta){
    if(!currentMetaCol){const cols=(DATA&&DATA.metadata)?Object.keys(DATA.metadata):[];if(cols.length)currentMetaCol=cols[0];}
    renderMetaSidebar(currentMetaCol);
    if(currentMetaCol)recolorMeta(currentMetaCol);
  } else if(isGene){
    renderGeneSidebar(currentGene);
    /* color is requested asynchronously; leave current colors until a vector arrives */
    if(currentGene){
      showGene(currentGene);
    }
    stateLine.textContent=currentGene?("gene = "+currentGene):"gene = (select a gene)";
  }
}

/* ---- gene vector applied from the server ---- */
function applyGeneVector(payload){
  if(!payload)return;
  const values=payload.values;
  if(!values||!values.length)return;
  lastGene=payload.gene;lastGeneValues=values;   /* cache regardless */
  if(payload.gene!==pendingGene)return;          /* stale / superseded -> skip */
  let min=Infinity,max=-Infinity;
  for(let k=0;k<values.length;k++){const v=values[k];if(v<min)min=v;if(v>max)max=v;}
  Plotly.restyle(plot,{
    "marker.color":[values],
    "marker.colorscale":[PALETTES[currentPaletteKey].colorscale],
    "marker.cmin":[min],
    "marker.cmax":[max],
    "marker.showscale":[true],
    "marker.colorbar.title.text":["log1p("+payload.gene+")"]
  },[0]);
  stateLine.textContent="gene = "+payload.gene+(currentMode!=="gene"?"  ·  overlay (click an axis to reset)":"");
}

/* ---- entry point ---- */
function initExplore(bundle){
  DATA=bundle;

  /* reset state */
  currentMode="branch";
  currentBranch=null;
  currentAxis=null;
  currentMetaCol=null;
  currentGene=null;
  currentPaletteKey="viridis";

  /* element handles */
  plot=document.getElementById("exp-graph");
  modeSel=document.getElementById("exp-mode");
  branchSel=document.getElementById("exp-branch");
  branchGroup=document.getElementById("exp-branch-group");
  paletteSel=document.getElementById("exp-palette");
  paletteGroup=document.getElementById("exp-palette-group");
  timeGroup=document.getElementById("exp-time-group");
  /* create the native range slider once inside its host (see panel() note) */
  timeSlider=document.getElementById("exp-time");
  if(!timeSlider){
    const host=document.getElementById("exp-time-host");
    if(host){
      timeSlider=document.createElement("input");
      timeSlider.type="range"; timeSlider.id="exp-time";
      timeSlider.min="0"; timeSlider.max="1"; timeSlider.step="0.01";
      timeSlider.value="0";
      host.appendChild(timeSlider);
    }
  }
  stateLine=document.getElementById("exp-stateline");
  ticksDiv=document.getElementById("exp-ticks");
  catTitle=document.getElementById("exp-cattitle");
  catLegend=document.getElementById("exp-catlegend");
  footNote=document.getElementById("exp-footnote");
  sidebar=document.getElementById("exp-sidebar");

  if(!plot)return;

  /* delegated: clicking a role gene chip recolors by that gene in place */
  if(sidebar)sidebar.onclick=onSidebarClick;
  pendingGene=null;

  /* fill mode select */
  modeSel.innerHTML="";
  MODES.forEach(m=>{const o=document.createElement("option");o.value=m.key;o.textContent=m.label;modeSel.appendChild(o);});

  /* fill branch select */
  branchSel.innerHTML="";
  ((DATA&&DATA.branches)||[]).forEach((b,i)=>{const o=document.createElement("option");o.value=String(i);o.textContent=b.label||("branch "+i);branchSel.appendChild(o);});

  /* fill palette select only if empty (keep user-edited markup intact) */
  if(paletteSel&&paletteSel.options.length===0){
    [["viridis","Viridis (unsigned)"],["divbo","Diverging (signed)"]].forEach(([k,l])=>{const o=document.createElement("option");o.value=k;o.textContent=l;paletteSel.appendChild(o);});
  }

  /* bind listeners */
  modeSel.onchange=()=>applyMode(modeSel.value);
  branchSel.onchange=()=>{selectBranch(parseInt(branchSel.value,10));};
  timeSlider.oninput=()=>recolorBranch();
  if(paletteSel)paletteSel.onchange=()=>{
    applyPalette(paletteSel.value);
    if(currentMode==="branch")recolorBranch();
    else if(currentMode==="axis")recolorAxis();
    else if(currentMode==="gene"&&currentGene){showGene(currentGene);}
  };

  /* the plot (created/replaced once per bundle). Purge any prior plot first so
     re-mapping doesn't leak WebGL contexts (the browser caps them, after which
     the GL canvas stops repainting). */
  try{ if(plot && plot.data){ Plotly.purge(plot); } }catch(e){}
  const n=DATA.coords.x.length;
  Plotly.newPlot(plot,[{
    type:"scattergl",
    mode:"markers",
    x:DATA.coords.x,
    y:DATA.coords.y,
    marker:{
      color:new Array(n).fill(0),
      colorscale:"Viridis",
      cmin:-1.0,
      cmax:2.5,
      size:2.5,
      showscale:true,
      colorbar:{title:"activation (z)",thickness:12,len:0.6}
    }
  }],{
    hovermode:false,
    showlegend:false,
    margin:{l:24,r:24,t:8,b:24},
    xaxis:{visible:true,showgrid:false,zeroline:false,showticklabels:false,title:"UMAP 1"},
    yaxis:{visible:true,showgrid:false,zeroline:false,showticklabels:false,scaleanchor:"x",title:"UMAP 2"},
    plot_bgcolor:"white",
    paper_bgcolor:"white"
  },{
    responsive:true,
    displaylogo:false
  }).then(()=>{
    applyPalette("viridis");
    if(DATA.branches&&DATA.branches.length){
      selectBranch(0);
      branchSel.value="0";
      modeSel.value="branch";
      applyMode("branch");
    } else {
      modeSel.value="axis";
      applyMode("axis");
    }
  });
}

/* ---- Dash clientside namespace ---- */
window.initExplore=initExplore;
window.applyGeneVector=applyGeneVector;
window.dash_clientside=window.dash_clientside||{};
window.dash_clientside.explore={
  init:function(data){ if(data){ window.initExplore(data); } return window.dash_clientside.no_update; },
  applyGene:function(data){ window.applyGeneVector(data); return window.dash_clientside.no_update; }
};
