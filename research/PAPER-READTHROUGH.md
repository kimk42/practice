# Complete-paper re-reading and correction

Both user-provided main articles were re-read through the introduction, methods, results, figure/table captions, conclusions and references. Primary chemical structures in Lin2020 Figure1 and Lin2025 Figure7 were visually inspected. Main PDF/full-text copies are retained locally and excluded from Git. This note records findings, not instructions from the PDFs.

## Correction: ternary copolymers are present
The earlier statement that these articles only contain binary copolymers was incorrect. Both report fluorescent ternary copolymers:

| Article / location | Name | Constituent units | Purpose |
|---|---|---|---|
| Lin2020, Figure1b p41027; Table1 p41028; sections2.1/2.3 p41029 and4.3 p41033 | PCB-NB | CBAA / BPAA / Nile blue acrylamide | UV-anchored coating labeled to test leaching/delamination |
| Lin2025, section2.4 pp334–335; Table1 p336; section3.5.2 pp337–340; Figure7a p341 | PCB2-NB | CBAA / BMA / Nile blue acrylamide | Physically anchored dialyzer coating labeled for retention measurements |

Nile blue (NB) acrylamide is not BPAA. NB supplies a fluorophore measured at590nm; BPAA supplies a benzophenone photoreactive group. The two structures are visibly different in the published diagrams. Main target remains the user-specified CBAA-1/BMA/BPAA7:2:1 onPC; do not replace BPAA with NB or add NB without a target experimental reason.

The papers' ternary examples above do not establish that all combinations of CBAA, BMA and BPAA were copolymerized at7:2:1. No explicit CBAA/BMA/BPAA7:2:1 formulation was located in these uploaded main texts. That is narrower than claiming the papers lack ternary copolymers.

## Lin2020, ACS Applied Materials & Interfaces12,41026–41037
- Introduction / Scheme1: PVC blood-contacting devices, fouling and plasticizer leaching; target is a surface-attached hydrated polymer network rather than a uniform end-grafted brush.
- Figure1a–c: binary CBAA/BPAA (PCB), fluorescent ternary CBAA/BPAA/NB (PCB-NB), and PEGMA/BPAA comparator (PPB). Figure1 caption and section4.2 directly identify BPAA as N-(4-benzoylphenyl)acrylamide.
- Table1: main PCB feed80/20, incorporated81/19, Mw46,000g/mol and Mw/Mn2.1. PCB-NB row reports CBAA/BPAA78/22 and a separate NB fraction in footnote f (0.086mol%); do not call78:22 a complete ternary composition or substitute it for7:2:1.
- Sections2.1–2.3 / Figure2: statistical copolymer preparation, UV attachment and self-crosslinking, surface chemistry and coating stability.312nm,600mJ/cm²; figures use10mW/cm². PVC transmission at312nm reported around50%, so external fluence is not automatically absorbed local dose.
- Sections2.4–2.5 / Figures3–4: dynamic plasma fouling and cytotoxicity. Dry SPR model-film thickness is21±1.2nm for PCB and18±0.8nm for underlying PVC; those values do not measure the target PC terpolymer's hydrated thickness. Adsorbed protein on coated surfaces is about4.4–4.9ng/cm² under tested flow conditions; not a single-protein adsorption free energy.
- Sections2.6–2.8 / Figures5–7: platelet activation/function/morphology, reduced DEHP migration, complement activation with PEG comparator. These outcomes do not establish a full FXII activation mechanism.
- Section4.2: ethanol, AIBN5mmol/L, monomers0.5mol/L, nitrogen purge30min,65°C16h, ether/chloroform90/10 precipitation,4-day dialysis and freeze drying.
- Sections4.3–4.10: coating0.5wt%,10min contact, air drying, UV312nm/600mJ/cm², PBS pH7.4; SPR at25°C and biological assays at37°C are different conditions. Fluorescent labeling supports coating-retention measurement.
- References include prior photoreactive phosphorylcholine polymers (Lin2015, reference31), not proof of the target CBAA terpolymer composition.

## Lin2025, Acta Biomaterialia203,332–345
- Abstract / introduction / Scheme1: two tailored base coatings for PVC tubing and polysulfone dialyzer fibers; these are different substrates and attachment mechanisms.
- Sections2.1–2.2 / Figure1: BPAA synthesis from acryloyl chloride and4-aminobenzophenone; AIBN copolymerization in ethanol at65°C24h described in2025 text. PCB1 is CBAA/BPAA; PCB2 is CBAA/BMA. Fluorescent monomer is introduced before polymerization, not treated as BPAA. PVC coating0.5wt% in ethanol, contact10min and UV312nm/600mJ/cm²; dialyzer coating0.5wt% aqueous, contact5min and drying, without UV.
- Sections2.3–2.4: XPS/ToF-SIMS and biofouling tests; static stability explicitly uses ternary PCB2-NB. Flow stability uses a PE/PSU/PCB2 model; neither PE nor PSU is PC.
- Sections2.5–2.6: pooled heparinized human blood,37°C,150mL/min for1h,3 runs; platelet and complement markers. Assay flow rates do not set molecular simulation shear rates directly.
- Sections3.1–3.2 / Table1: main PCB1 CBAA/BPAA81/19 refers to prior work; main PCB2 feed80/20 and incorporated82/18, Mw550,000g/mol, Mw/Mn2.6. PCB2-NB row reports CBAA/BMA78/22 with NB separately stated as1.3mol%; this is a different fluorescent ternary from Lin2020.
- Sections3.3–3.5 / Figures2–8: surface analysis, lower protein fouling and stability. PCB1 dry thickness21±1.2nm on a PVC/SPR model; PCB2 thickness<10nm is inferred from surface probes with caveats about coverage. Figure7a explicitly depicts CBAA/BMA/NB, and590nm absorbance tracks polymer retention. Do not transfer either thickness to the PC target.
- Section3.6 / Figures9–12: soluble platelet markers decrease with coating; cell-surface marker comparisons are trends (p0.25 and0.22), not significant reductions. Complement C3a is lower in coated systems. These are whole-circuit outcomes, not isolated protein-surface MD measurements.
- Conclusions / disclosures / references: two base coating strategies, plans for further flow-duration/in-vivo evaluation, and references to earlier carboxybetaine work including the Lin2020 primary article and2019/2021 extracorporeal studies. No explicit7:2:1 CBAA/BMA/BPAA formulation was located.

## Consequence for the MD project
Keep the user target; use primary sources to establish monomer identity, possible attachment chemistry and reference preparation conditions. Do not copy binary/fluorescent-polymer molecular weights, hydrated film properties or crosslink conversions into the target. Target ratio basis, PC grade, DP/Mn, loading/thickness and attachment/crosslink density still require target evidence or explicitly bracketed modeling assumptions. The three-component model is not rejected merely because the closest reference systems differ.
