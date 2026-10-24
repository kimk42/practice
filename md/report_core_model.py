"""Publish preparation evidence for the six-chain observed-core/PE model; no MD plots invented."""
import argparse,base64,html,json,os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/fibrinogen-core-plot')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from openmm.app import AmberPrmtopFile,AmberInpcrdFile,element
from openmm import unit


def publish(core,bundle,output):
    core=Path(core);bundle=Path(bundle);out=Path(output);out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((core/'manifest.json').read_text());relax=json.loads((core/'local-relaxation-restrained/report.json').read_text());model=json.loads((bundle/'manifest.json').read_text())
    p=AmberPrmtopFile(str(bundle/'dry.prmtop'));c=AmberInpcrdFile(str(bundle/'dry.inpcrd'));xyz=np.array(c.positions.value_in_unit(unit.nanometer));mapping=manifest['residue_mapping']
    fig=plt.figure(figsize=(14,7));gs=fig.add_gridspec(2,1,height_ratios=[1.5,1])
    ax=fig.add_subplot(gs[0]);colors=['#cb7c25','#2897a8','#8756ad','#b8a52c','#389c70','#c76f94']
    for ch,color in zip('ABCDEF',colors):
        ca=[a.index for a in p.topology.atoms() if a.residue.index<len(mapping) and mapping[a.residue.index]['source_chain']==ch and mapping[a.residue.index]['source_label_seq'] is not None and a.name=='CA']
        ax.plot(xyz[ca,0],xyz[ca,2],lw=.8,color=color,label=f'{ch}: '+dict(A='alpha',B='beta',C='gamma',D='alpha',E='beta',F='gamma')[ch])
    pe=[a.index for a in p.topology.atoms() if a.residue.name=='PEO' and a.element==element.carbon]
    ax.scatter(xyz[pe,0],xyz[pe,2],s=.4,color='#76b88b',alpha=.18,label='PE carbons (prepared patch)')
    caps=[a.index for a in p.topology.atoms() if a.residue.index<len(mapping) and mapping[a.residue.index]['artificial_cap'] and a.element!=element.hydrogen]
    ax.scatter(xyz[caps,0],xyz[caps,2],s=8,color='#c33c32',marker='x',label='Artificial neutral caps')
    ax.set_xlabel('X (nm)');ax.set_ylabel('Z (nm)');ax.set_title('Side projection of actual prepared coordinates; no solvent or MD');ax.legend(loc='upper center',ncol=4,fontsize=8)
    ax.set_aspect('equal',adjustable='datalim');ax.grid(alpha=.15)
    bottom=fig.add_subplot(gs[1]);observed=[174,401,381,174,401,382];total=[562,461,411,562,461,411];missing=np.array(total)-observed
    bottom.bar(list('ABCDEF'),observed,color=colors,label='Observed residues retained')
    bottom.bar(list('ABCDEF'),missing,bottom=observed,color='#dbe0e4',label='Absent coordinates (not rebuilt)')
    bottom.set_ylabel('Residues in deposited sequence');bottom.legend(fontsize=8,ncol=2);bottom.set_ylim(0,650)
    fig.suptitle('Six-chain 3GHG observed core + size-matched PE: PREPARED INPUT MODEL',fontsize=15)
    fig.text(.5,.035,'1,913 observed residues | 12 neutral caps | 25 disulfides | 4 calcium ions | 690 PE oligomers | 81,275 dry atoms',ha='center',fontsize=10)
    fig.text(.5,.01,'Long alphaC/tails and glycans omitted. Not full-length fibrinogen; not equilibrated; no new adsorption MD or affinity result.',ha='center',fontsize=10,color='#a12d32')
    fig.tight_layout(rect=(0,.065,1,.93));fig.savefig(out/'model.png',dpi=160);fig.savefig(out/'model.svg')
    report={'status':'prepared_dry_input_not_adsorption_MD','core':{k:v for k,v in manifest.items() if k not in ('residue_mapping','sha256')},'local_preparation':relax,'mixed_model':{k:v for k,v in model.items() if k!='lower_layer_residue_ids_for_future_restraints'}}
    (out/'summary.json').write_text(json.dumps(report,indent=2))
    image=base64.b64encode((out/'model.png').read_bytes()).decode()
    text=f'''<!doctype html><html lang="ko"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>6사슬 피브리노겐 코어–PE 모델 준비 결과</title><style>body{{max-width:1100px;margin:35px auto;padding:0 20px;font:16px/1.7 system-ui;color:#172d39}}img{{width:100%}}table{{border-collapse:collapse;width:100%}}td,th{{padding:9px;border-bottom:1px solid #ddd;text-align:left}}.note{{background:#fff0dc;padding:18px;border-left:4px solid #b77729}}pre{{white-space:pre-wrap}}</style><h1>여섯 사슬 피브리노겐 관측 코어 + PE</h1><p><b>완료 상태: 힘장 파라미터를 적용한 건조 입력 모델. 새 모델의 흡착 MD는 미실행.</b></p><table><tr><td>관측 단백질</td><td>α·β·γ 각 2사슬, 1,913잔기</td></tr><tr><td>추가 처리</td><td>중성 ACE/NME 캡 각 6개 · 원본 칼슘 4개 보존</td></tr><tr><td>이황화 결합</td><td>원본 25개 보존, 그중 사슬 간 15개</td></tr><tr><td>관측 코어 국소 조정</td><td>중원자에 강한 위치 제약, 최대 이동 {relax['max_restrained_atom_displacement_nm']*10:.3f} Å</td></tr><tr><td>접촉 검사</td><td>완화 뒤 비결합 원자 쌍의 1.2 Å 미만 접촉 0개. 전체 물리적 검증을 뜻하지 않음.</td></tr><tr><td>PE와 전체 크기</td><td>PE 690사슬 · 물 없는 모델 81,275원자</td></tr><tr><td>용매화 규모 추정</td><td>상자 기준 약 69만 원자(용질 배제 전 추정). 실제 용매화·동역학은 미실행.</td></tr></table><p class="note">빠진 αC 영역과 말단을 임의로 복원하지 않았습니다. 중성 캡은 누락 영역을 대신하는 실제 구조가 아닙니다. 당쇄와 knob 펩타이드는 working model에서 생략했습니다. 전체 길이의 당화된 피브리노겐을 재현한 모델이 아니며, 앞선 작은 γ절편의 흡착 결과를 이 모델에 적용할 수 없습니다.</p><img src="data:image/png;base64,{image}" alt="여섯 사슬 코어와 PE 준비 좌표 및 사슬별 관측·미관측 잔기"><h2>다음 계산에 필요한 것</h2><p>물과 중화 이온 추가 → 용매와 계면의 제한적 최소화 → 압력·밀도·계면 평형화 → 여러 초기 배향과 시드로 반복 계산. 현재 물 없는 입력 모델에는 이러한 계산 결과가 없습니다.</p><details><summary>원본 검증 및 준비 메타데이터</summary><pre>{html.escape(json.dumps(report,ensure_ascii=False,indent=2))}</pre></details></html>'''
    (out/'report.html').write_text(text);print('Published preparation report:',out)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('core');p.add_argument('bundle');p.add_argument('output');a=p.parse_args();publish(a.core,a.bundle,a.output)
