const state=document.getElementById('campaignState');
async function refresh(){try{const response=await fetch('artifacts/research-campaign/status.json',{cache:'no-store'});if(!response.ok)throw Error('No campaign record');const data=await response.json();const candidates=Object.entries(data.molecular_candidates||{}).map(([name,item])=>`${name}: ${item.status}`).join(' · ');state.textContent=`${data.status} · ${candidates||'분자 입력 준비 중'} · 후속 반복 MD와 PMF의 실행 여부는 보고서에서 확인하세요.`;}catch{state.textContent='아직 후속 실행 기록이 없거나 서버에 연결되지 않았습니다.';}}
await refresh();setInterval(refresh,10000);
