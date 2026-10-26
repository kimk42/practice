# 사용자가 제공한 실제 실험 조건

| 항목 | 현재 확인 내용 |
|---|---|
| 합성 조성 | CBAA-1:BMA:BPAA=7:2:1 몰비. 실제 편입 조성 분석값은 미보고. |
| 코팅액 | 물에0.5wt% |
| 표면 접촉 |30분 |
| 후처리 | 드레인 후 건조 |
| UV | 회당9000mJ/cm²×3회,총27000mJ/cm²=27J/cm². 사용자 확인. |
| UV 파장 | 인용 논문을 따름. 논문의312nm를 참조하며 별도 장비 측정값은 미보고. |
| 단백질 온도·pH | 실온,pH7.4. 정확한°C는 미보고. SOP의37°C를 실제 조건으로 덮어쓰지 않음. |
| 버퍼 | cPBSzI1×. 업로드 SOP의 레시피 확인;config/experimental-buffer.json 참고. |
| 농도 기준 | 각 정제 단백질을 정상 혈장 농도의1%로 설정.1% 희석 혈장이 아님. 사용자 확인. |
| Fibrinogen 농도 | SOP 주 흐름의 최종1×0.03mg/mL=30µg/mL로 매핑. 뒤쪽 단독 Fbg 옵션의0.3mg/mL와 불일치는 별도 기록. |
| FXII 농도 | 이 SOP에는 수치가 없음. 문헌375nM 기준의 잠정1%값3.75nM을 별도 출처로 연결(하단 참고). |
| 단백질 접촉 |2시간=7200초 |

## cPBSzI1×,1L 기준

NaCl6.43g,NaI1.50g,citric acid monohydrate2.10g,sodium phosphate monobasic1.38g,sodium azide0.20g,pH7.4 조정. 대략 NaCl110mM,NaI10mM,total citrate10mM,azide3.08mM이다. Phosphate는 monohydrate라면10mM이며 문서에는 수화물 형태가 명시되지 않았다. 이온화 상태·Na+ 수지·검증된 MD 파라미터는 별도로 정한다.

## 모델과의 연결

기존300K·중화 이온 파일럿은 초기 수치 검사 기준으로 유지한다. 버퍼의 각 성분을 검증해 포함하기 전 이를 cPBSzI 조건의 계산으로 표시하지 않는다. 실온의 정확한 온도,막 두께·흡착량,분자량·분산도,기판 등급과 가교 분포는 계속 미정이다. 코팅액0.5wt%와 UV 조사량만으로 막 밀도나 접목률을 확정하지 않는다.

SOP는 albumin/fibrinogen 경쟁 흡착과 단독·순차 실험 옵션을 포함한다. 해당 혼합물을 현재 Fbg/FXII 개별 모델에 자동으로 추가하지 않는다. 방사성 표지·계수·SDS 잔류량 단계는 참고 자료이며 실행 지시나 현재 MD의 반응 경로가 아니다. 자세한 판독 근거는 research/protocols/RADIOLABEL_SOP_READOUT.md에 보존했다.

## FXII 농도의 문헌 기반 잠정 매핑

일차 연구 DOI10.3389/fcvm.2022.1008410은 FXII 결핍 마우스 전혈에 사람 유래 FXII를 보충하면서375nM을 생리적 기준으로 사용했다. 사용자1% 기준을 여기에 적용한 잠정값은3.75nM이다.80kDa 질량을 가정해 환산하면 약0.3µg/mL이지만,이는 실제 사용자 시료 농도 또는 현재 비당화 MD 단백질의 측정 질량이 아니다. 출처·적용 범위는 research/proteins/modeling-evidence/fxii-plasma-reference.json에 기록했다. 단일 단백질 주기 상자가 이 희석 농도를 직접 재현한다고 표시하지 않는다.
