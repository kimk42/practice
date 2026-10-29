# Citrate numerical-control candidates

GAFF2/AM1-BCC candidates for isolated numerical controls only. CIT3 is citrate³⁻;
HCIT2 is one terminal-protonated hydrogen-citrate²⁻ microstate. No pH populations,
aqueous physical approval or full cPBSzI readiness are assigned.

Include `types.itp` once after compatible Amber defaults, then the required
molecule ITP. The `CBUF_` namespace prevents accidental host type redefinitions.
The atom charges match the preserved original candidates at exported precision;
no charge projection was applied. See `provenance.json` and the parent directory's
`refined-citrate-water-audit.json` for source hashes and completed2ps control scope.
