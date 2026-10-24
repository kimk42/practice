import {parsePDB,contacts} from './model.mjs';
const $=id=>document.getElementById(id), canvas=$('scene'),ctx=canvas.getContext('2d');
let real=false,atoms=Array.from({length:150},(_,i)=>({x:(i-75)*2,y:Math.sin(i*.17)*10,z:Math.cos(i*.17)*10,chain:'demo',residue:'',id:String(i)}));
function draw(){
 const w=canvas.clientWidth,h=canvas.clientHeight,dpr=devicePixelRatio||1;canvas.width=w*dpr;canvas.height=h*dpr;ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,h);
 const angle=+$('angle').value*Math.PI/180,gap=+$('gap').value,cutoff=+$('cutoff').value;
 $('angleValue').textContent=$('angle').value+'°';$('gapValue').textContent=gap+' Å';$('cutoffValue').textContent=cutoff+' Å';
 const state=contacts(atoms,angle,gap,cutoff), span=Math.max(60,...state.atoms.map(a=>Math.abs(a.x)*2),...state.atoms.map(a=>Math.abs(a.z)*2),Math.abs(state.floor)*2+30),scale=Math.min(w*.85,h*.8)/span;
 const project=a=>[w/2+a.x*scale,h*.43-a.z*scale+a.y*scale*.25];
 const py=h*.43-state.floor*scale;const color={PTFE:'#64b9d6',PVDF:'#b8a0e6',PE:'#77c7a2'}[$('polymer').value];
 ctx.fillStyle=color+'25';ctx.fillRect(0,py,w,Math.max(0,h-py));ctx.strokeStyle=color;ctx.beginPath();ctx.moveTo(0,py);ctx.lineTo(w,py);ctx.stroke();
 ctx.fillStyle=color;ctx.font='12px system-ui';ctx.fillText($('polymer').value+' · 理想 평면'.replace('理想','이상적인'),14,Math.min(h-12,py+20));
 state.atoms.forEach((a,i)=>{const [x,y]=project(a),prev=state.atoms[i-1];if(prev&&prev.chain===a.chain&&Math.hypot(a.x-prev.x,a.y-prev.y,a.z-prev.z)<5){const p=project(prev);ctx.beginPath();ctx.strokeStyle='#7c9baa';ctx.moveTo(...p);ctx.lineTo(x,y);ctx.stroke();}ctx.beginPath();ctx.fillStyle=a.contact?'#ffd07c':'#84d6df';ctx.arc(x,y,a.contact?3.4:2.2,0,Math.PI*2);ctx.fill();});
 const near=state.atoms.filter(a=>a.contact);
 $('results').textContent=real?(near.length?near.map(a=>`${a.chain||'∅'}:${a.residue} ${a.id} — ${a.distance.toFixed(1)} Å`).join('\n'):'현재 거리 기준 안에 있는 Cα가 없습니다.'):'실제 PDB 파일을 열면 잔기별 근접 결과가 표시됩니다. 모식도에서는 연구용 수치를 제공하지 않습니다.';
}
$('file').addEventListener('change',async e=>{const f=e.target.files[0];if(!f)return;try{if(f.size>30*1024*1024)throw Error('30 MB 이하의 PDB 파일을 선택하세요.');const parsed=parsePDB(await f.text());atoms=parsed;real=true;$('badge').textContent='사용자 PDB · Cα 구조';$('status').textContent=`${f.name} · Cα ${atoms.length}개 · 단백질 종류와 구조 완전성은 원본 출처에서 확인하세요.`;draw();}catch(err){$('status').textContent='읽기 실패: '+err.message;}});
for(const id of ['angle','gap','cutoff','polymer'])$(id).addEventListener('input',draw);
$('reset').onclick=()=>{$('angle').value=0;$('gap').value=5;$('cutoff').value=8;draw();};
let last=null;canvas.onpointerdown=e=>{last=e.clientX;canvas.setPointerCapture(e.pointerId);};canvas.onpointermove=e=>{if(last===null)return;$('angle').value=Math.max(-180,Math.min(180,+$('angle').value+(e.clientX-last)));last=e.clientX;draw();};canvas.onpointerup=canvas.onpointercancel=()=>{last=null;};window.addEventListener('resize',draw);draw();
// Load the supplied, observed fibrinogen coordinates; no missing residues are invented.
fetch('data/3ghg-upload/assembly1-view.pdb').then(r=>{if(!r.ok)throw Error('structure unavailable');return r.text();}).then(text=>{
 if(real)return;
 atoms=parsePDB(text);real=true;
 $('badge').textContent='실제 3GHG · 조립체 1 · 관측 Cα 좌표';
 $('status').textContent=`사용자 제공 3GHG · 피브리노겐 6사슬 · Cα ${atoms.length}개. 누락 영역은 미복원이며, knob·당쇄·칼슘은 이 보기에서만 생략했습니다. 흡착 계산 결과가 아닙니다.`;
 draw();
}).catch(()=>{});
