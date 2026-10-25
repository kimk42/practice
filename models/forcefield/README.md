# Exploratory force-field artifacts

Generated GAFF2/AM1-BCC parameters for PC4 and the single GRAFT molecule. provenance.json records original output sources, net charges and SHA256; files are preserved byte-for-byte including whitespace. Parameter gates and the combined surface EM/0.2ps MD passed, followed by a10ps fixed-volume numerical follow-up. These facts establish initial numerical usability, not physical validation of the coating, analog torsions or adsorption predictions.

The COP10 standalone parameter gate also passed; the actual surface uses the covalently connected GRAFT species rather than separately combining a COP10 and PC4 topology across a nonexistent bond. Molecule geometry, chemistry and assumptions are described elsewhere in models/. These force-field files do not include protein glycans, physiological salt or production sampling.
