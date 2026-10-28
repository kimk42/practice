# Current coating model: user target and verified literature

## User target
Polycarbonate substrate; UV-immobilized **CBAA-1-co-BMA-co-BPAA** coating; composition **7:2:1** in that order. Synthesis molar ratio7:2:1 is confirmed. PC grade,incorporated composition,molecular weight/DP,film loading/thickness and reacted fraction remain unreported. Actual coating:0.5wt% in water,30min contact,drain/dry,9000mJ/cm2 per exposure x3;312nm is the cited-paper reference. See USER_EXPERIMENT.md. Do not replace this target with a binary polymer or PVC substrate. Do not infer a periodic sequence, an end-grafted brush or crosslink conversion from composition.

## Primary papers now available
User-uploaded complete PDFs are preserved locally as lin2020-main.pdf and actbio2025-main.pdf, with extracted text (*-main.txt). PDFs/full texts are Git-ignored; uploaded-papers.json records source names, DOI and SHA-256. This resolves the earlier lack of main-text access, regardless of the ACS website's HTTP403.

1. Xiaojie Lin et al. **Photoreactive Carboxybetaine Copolymers Impart Biocompatibility and Inhibit Plasticizer Leaching on Polyvinyl Chloride**, ACS Applied Materials & Interfaces 12 (2020) 41026–41037. DOI https://doi.org/10.1021/acsami.0c09457. Xiaojie Lin first author, Shaoyi Jiang one of the corresponding authors.
2. Xiaojie Lin et al. **Carboxybetaine zwitterionic coatings on kidney dialyzer components: Surface characterization and blood interactions**, Acta Biomaterialia 203 (2025) 332–345. DOI https://doi.org/10.1016/j.actbio.2025.07.021. Buddy D. Ratner corresponding author; this is a distinct paper, not evidence of a Jiang-authored 7:2:1 terpolymer.

## Component structures
BPAA is **N-(4-benzoylphenyl)acrylamide**, directly confirmed by 2020 Figure1 (p41027) and section4.2 (p41033), and 2025 section2.1 (p333). It is made from acryloyl chloride and 4-aminobenzophenone. Monomer connectivity: `C=CC(=O)Nc1ccc(C(=O)c2ccccc2)cc1`. The linkage is an acrylamide, not a benzophenone acrylate ester or free acrylic acid. The existing BPAA candidate geometry matches this identity; its filename retains “candidate” because the actual target batch and production parameters are not established.

BMA is n-butyl methacrylate, `C=C(C)C(=O)OCCCC`, independently defined in Lin2018 SI (DOI10.1021/acs.langmuir.8b02540.s001). CBAA-1 connectivity remains as in chemistry/MODEL.md. All three polymerizable vinyl double bonds must be consumed when building saturated polymer backbones; BMA additionally has an alpha-methyl group. Post-UV benzophenone sites require different chemistry and parameters from unreacted BPAA.

## Verified synthesis and coating conditions — reference systems only

| Item | Lin2020 CBAA/BPAA (PCB) | Lin2025 CBAA/BMA (PCB2) |
|---|---|---|
| Feed / incorporated molar composition | 80:20 / 81:19 | 80:20 / 82:18 |
| Total monomer | 0.5 mol/L | 2.0 mol/L |
| AIBN | 5 mmol/L | 5 mmol/L |
| Polymerization | ethanol, N2 purge30min, 65°C16h | ethanol, N2 purge30min, 65°C24h |
| Mw / dispersity Mw/Mn | 46,000g/mol / 2.1 | 550,000g/mol / 2.6 |
| Substrate and attachment | PVC, UV covalent immobilization | polysulfone dialyzer, adsorptive attachment |

2020 source: Table1 p41028, sections4.2–4.3 p41033. 2025 source: sections2.1–2.2 pp333–334 and Table1 p336. The 2025 synthesis text describes65°C24h and refers back to the earlier PCB1 paper for characteristics; do not merge the different synthesis durations.

2020 coating procedure: polymer solution0.5wt%, contact on rotator10min, drain and dry with airflow at room temperature; **UV312nm and fluence600mJ/cm²**; rinse/prewet with sterile1×PBS,pH7.4. Section4.3 does not itself specify the coating solvent there. The 2025 PVC/PCB1 protocol explicitly uses0.5wt% in ethanol and the same312nm/600mJ/cm² exposure (section2.2.1). Figure2c in2020 reports312nm at10mW/cm² with spectra at0/10/30/60s. Numerically600mJ/cm² corresponds to60s at constant10mW/cm², but geometry, attenuation and actual target exposure cannot be inferred from this arithmetic.

These papers document binary base coatings **and fluorescent ternary copolymers**: CBAA/BPAA/Nile blue acrylamide (Lin2020 Figure1b,p41027) and CBAA/BMA/Nile blue acrylamide (Lin2025 Figure7a,p341). The earlier claim that they only contain binary copolymers was incorrect. No explicit CBAA-1/BMA/BPAA7:2:1 onPC formulation was located in the uploaded main texts. See PAPER-READTHROUGH.md for the full correction and page-specific reading notes. Their measured Mw, composition, coating thickness (if any) and adsorption mechanism must not be assigned to the target terpolymer. No measured target DP or crosslink fraction has been established.

## Model and validation gate
Fixed-topology GROMACS should model a post-UV network, not simulate photoexcitation. The SI mechanism depicts hydrogen abstraction, conversion of benzophenone carbonyl to an alcohol and formation of a new C–C attachment. Parameterize reacted/unreacted BPAA separately; constrain possible attachment to PC versus interchain links using target experimental evidence. UV fluence alone does not determine conversion.

Component geometry candidates pass connectivity/formal-charge/embedding/MMFF-geometry checks, but MMFF is not the production force field. DP100 shuffled 70/20/10 sequences are provisional composition examples, not measured chains or reaction kinetics. Next:evaluate finite-PC-chain,coating-length/coverage/attachment scenarios with explicit assumptions;validate physical properties of assigned GAFF2/AM1-BCC parameters and relaxation before comparative protein calculations. Unmeasured values can be sensitivity scenarios;they are not all mandatory new user inputs. Use separate fibrinogen and FXII simulations; 8OS5 is only FXII N-terminal domains.
