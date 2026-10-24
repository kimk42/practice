const $=id=>document.getElementById(id);
function update(){
 const pvdf=$('polymer').value==='PVDF';
 $('phase').disabled=!pvdf;
 const alpha=pvdf&&$('phase').value==='alpha';
 $('polarization').disabled=alpha;
 if(alpha)$('polarization').value='nonpolar';
 $('phaseNote').textContent=pvdf?{alpha:'α: 대표적인 TGTG′ 사슬 / 비극성 결정 배열',beta:'β: all-trans 사슬 / 극성 결정. 도메인 배향은 별도 지정',gamma:'γ: 대표적인 T₃GT₃G′ 사슬 / 극성 결정. 도메인 배향은 별도 지정'}[$('phase').value]:'결정 구조와 노출면을 준비된 원자 모델에 명시해야 합니다.';
}
for(const id of ['polymer','phase'])$(id).addEventListener('change',update);
$('exportConfig').onclick=()=>{
 if(!$('mdForm').reportValidity())return;
 const config={schema_version:1,surface:{polymer:$('polymer').value,phase:$('polymer').value==='PVDF'?$('phase').value:'specified-in-model',polarization:$('polarization').value,face:$('face').value.trim()},temperature_K:Number($('temperature').value),timestep_fs:Number($('timestep').value),steps:Number($('steps').value),seed:Number($('seed').value)};
 const url=URL.createObjectURL(new Blob([JSON.stringify(config,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='simulation.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
 $('mdStatus').textContent='계산 조건을 저장했습니다. 실행 전 해당 결정상·표면·분극을 실제 좌표와 힘장에 반영한 모델 번들이 필요합니다.';
};
update();

fetch('artifacts/protein-pilot/result.json').then(r=>{if(!r.ok)throw Error('No completed run');return r.json();}).then(r=>{if(r.status!=='completed')return; $('pilotResults').hidden=false;$('pilotSummary').textContent=`업로드한 3GHG γ사슬 조각 · ${r.residues}잔기 · ${r.particles}원자 · ${r.steps*r.timestep_fs/1000} ps · 최소화 구조 대비 최종 Cα RMSD ${r.final_CA_RMSD_from_minimized_A.toFixed(3)} Å`;}).catch(()=>{});

fetch('artifacts/pe-pilot/result.json').then(r=>{if(!r.ok)throw Error('No completed mixed run');return r.json();}).then(r=>{if(r.status!=='completed')return;$('peResults').hidden=false;$('peSummary').textContent=`실제 3GHG γ사슬 조각 + PE ${r.model_manifest.polymer_chains}사슬 + 물·중화 이온 · ${r.atoms.toLocaleString()}원자 · ${r.steps*r.timestep_fs/1000} ps · 최종 최소 거리 ${r.final_min_distance_nm.toFixed(3)} nm · 근접 잔기 ${r.final_contact_residues.length}개`;}).catch(()=>{});

fetch('artifacts/core-model/summary.json').then(r=>{if(!r.ok)throw Error('No core model');return r.json();}).then(r=>{if(r.status!=='prepared_dry_input_not_adsorption_MD')return;$('coreModel').hidden=false;$('coreModelSummary').textContent=`관측 ${r.core.observed_residues}잔기 · 이황화 결합 ${r.core.source_disulfides}개 · 캡 12개 · 칼슘 ${r.core.calcium_ions}개 · PE ${r.mixed_model.polymer_chains}사슬 · 물 없는 입력 ${r.mixed_model.dry_atoms.toLocaleString()}원자`;}).catch(()=>{});

fetch('artifacts/staged-preparation/result.json').then(r=>{if(!r.ok)throw Error('No staged results');return r.json();}).then(r=>{if(r.status!=='staged_preparation_completed_equilibrium_not_established')return;$('stagedPreparation').hidden=false;$('stagedSummary').textContent=`${r.atoms.toLocaleString()}원자 · 새 준비 구간 ${r.time_ps} ps · 최종 온도 ${r.final_temperature_K.toFixed(1)} K · 전체 밀도 ${r.final_total_density_g_cm3.toFixed(3)} g/cm³ · 평형 확립 미완료`;}).catch(()=>{});

fetch('artifacts/large-hold/result.json').then(r=>{if(!r.ok)throw Error('No hold result');return r.json();}).then(r=>{if(r.status!=='100K_hold_MD_completed_numerical_checks_passed')return;const target=$('latestHoldSummary');if(target)target.textContent=`${r.atoms.toLocaleString()}원자 혼합계 · 추가 ${r.time_ps} ps 유지 · 큰 모델 MD 누적 ${r.final_context_time_ps.toFixed(1)} ps · 최종 실제 온도 ${r.final.temperature_K.toFixed(2)} K · 유지 시작 대비 Cα RMSD ${r.final.protein_CA_RMSD_from_hold_start_A.toFixed(3)} Å`;}).catch(()=>{});
