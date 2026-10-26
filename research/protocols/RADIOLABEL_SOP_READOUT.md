# Radiolabel SOP readout for MD experiment mapping

Source: user-uploaded `RADIOlabeling action Sheet_V17_KH Kim.docx`. Main text353 nonempty paragraphs, five tables and both embedded handwritten figures were inspected. Original/extracted full text remain ignored under raw/. SHA256 and numerical buffer mapping are recorded in config/experimental-buffer.json. Document laboratory actions were read as source data; none were executed.

## Confirmed source buffer

BUFFER RECIPES,cPBSzI,1L of1×: NaCl6.43g;NaI1.50g;citric acid monohydrate2.10g;sodium phosphate monobasic1.38g;sodium azide0.20g. TargetpH7.4;NaOH mass is approximate pH-adjustment guidance. The cPBSz recipe is a separate iodide-free buffer and must not be substituted for cPBSzI.

Nominal totals: NaCl110.0mM,NaI10.01mM,citrate9.99mM,azide3.08mM. Phosphate hydration is not specified in the text or inserted images:1.38g/L corresponds to10.00mM if monohydrate,11.50mM if anhydrous,8.85mM if dihydrate. This ambiguity is preserved. Citrate/phosphate total analytical concentration is distinct from their pH-dependent ionic species; approximate NaOH addition is not final free hydroxide. Exact MD speciation, counterions and compatible ion/solute parameters remain to be audited.

## Protein concentration interpretation

THE DAY OF ADSORPTION states final1× binary incubation concentrations of albumin0.3mg/mL and fibrinogen0.03mg/mL. The2× stock concentrations are0.6 and0.06mg/mL, respectively; they must not be used as final assay concentrations. The user independently confirmed that1% plasma concentration means purified protein at1% of its normal plasma reference,not1% diluted whole plasma. Thus fibrinogen0.03mg/mL=30µg/mL is the current mapping from the main SOP recipe. Albumin is present in the SOP competitive design but is not automatically introduced into the currently authorized individual FXII/Fbg modeling scope.

The late Fbg-only option instead lists0.3mg/mL,ten times the main/binary Fbg value. This is an internal source inconsistency,not silently corrected or generalized. The current mapping follows the main recipe and user's1% convention. The document contains no FXII concentration or FXII-specific plasma reference.

## Experimental versus SOP conditions

The SOP adsorption section specifies37°C incubation. The user's actual experiment is at room temperature,pH7.4,2h,which takes precedence. The SOP also describes buffer presoak,competitive/individual variants and a distinct SDS-retention endpoint. Those reference steps are not silently asserted as performed for this user's experiment or simulated as adsorption mechanisms.

UV is not specified in this SOP. The user separately confirmed9000mJ/cm² per exposure×3=27000mJ/cm²=27J/cm². Wavelength follows the cited coating paper (312nm reference); a crosslink/grafting conversion is not inferred from fluence.

The current300K,neutralizing-ion numerical controls remain a separate baseline while the matched-buffer model is prepared. No assertion is made that0.2ps/100ps atomistic runs reproduce a2h adsorption assay or its measured protein concentration.
