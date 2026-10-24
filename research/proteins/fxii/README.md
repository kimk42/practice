# Human Factor XII structure retrieval

Target reference: human F12, UniProtP00748. Raw coordinates/API records are retained under raw/ and excluded from Git; downloads.json records sources, sizes and SHA-256. These are downloaded references, not prepared MD inputs.

## RCSB downloads completed

| ID | Downloaded coordinate file | Scope / selection consequence |
|---|---|---|
| 8OS5 | raw/8OS5.cif | N-terminal domains1–5; construct295 residues, not full-length FXII. Useful experimental reference for the surface-binding region. Entity description says “Coagulation factor XII-Mie”; sequence comparison below flags a difference. |
| 7PRJ | raw/7PRJ.cif |71-residue FnII domain; domain-level reference, not full protein. |
| 6X0T | raw/6X0T.cif |259-residue protease construct; inhibitor-bound activated FXIIa. Do not substitute for intact zymogen surface adsorption. |

Organism and P00748 mapping were checked against RCSB polymer-entity metadata. Mapping/construct length is not equivalent to modeled-coordinate completeness. Biological assemblies, missing residues, mutations, glycosylation/disulfides and bound ligands must be checked before using any reference in MD. Do not blindly concatenate these fragments into a full protein.

## Full-length computed candidate located in RCSB

RCSB entry https://www.rcsb.org/structure/AF_AFP00748F1 is “Computed structure model of Coagulation factor XII”, source AlphaFoldDB (AF-P00748-F1). RCSB metadata reports615 residues, human origin and computational methodology; it is **not an experimental full-length structure**. The RCSB entity sequence is retained as P00748-reference.fasta. This is precursor sequence; mature-chain boundaries, confidence/domain orientation, glycans and disulfides still require preparation. A metadata disulfide count of0 is not evidence that physiological FXII has no disulfides.

The standard files.rcsb.org/download endpoint returned404 for the computed ID. Its actual model-server coordinate endpoint is https://models.rcsb.org/v1/AF_AFP00748F1/full?encoding=cif. This request and direct AlphaFold downloads were blocked by runtime network policy403. RCSB metadata was retrieved successfully, but the full computed coordinate file was not downloaded. Added search.rcsb.org,models.rcsb.org and alphafold.ebi.ac.uk to the environment network draft while preserving existing rules; saving a draft is not runtime connectivity proof.

Experimental discovery used the PDBe UniProtP00748 best-structure mapping API and verified selected entries against RCSB, because RCSB's separate search domain was blocked. The returned list did not identify a full-length experimental zymogen; this is not an exhaustive proof that no such structure exists.

## Sequence caution for8OS5

The295-residue8OS5 entity sequence differs from the RCSB P00748 computed-model reference over the mapped20–314 interval at one position: referenceA207 versus constructP207 (precursor numbering inferred from this mapping). RCSB's mutation-count metadata nevertheless reports0. Treat this as a sequence discrepancy requiring explicit verification, not a confirmed wild-type construct. Do not silently restore/introduce a residue based solely on this note; check UniProt and the construct's primary sequence documentation before selecting a wild-type model.

## Reproduce downloads

```bash
source scripts/activate.sh
python scripts/fetch_fxii.py --include-computed
```

The script reuses and validates existing coordinate files and records failures separately. Once model-server network settings apply, it can download the full computed candidate without asking the user to supply a protein file. No private credentials are required for these public structures.
