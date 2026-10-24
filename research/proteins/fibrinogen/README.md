# Fibrinogen reference

`python scripts/fetch_fibrinogen.py` retrieves3GHG asymmetric-unit and biological-assembly1 mmCIF files plus entry metadata from RCSB. Raw files are ignored and hashes recorded in downloads.json. `python scripts/audit_protein_inputs.py` checks observed residues against each declared construct sequence without synthesizing missing regions.

Assembly1 has six protein chains; the observed1913 residues leave955 of the declared2868 construct residues unmodeled. This reference is not a ready full-native-fibrinogen MD input. See models/proteins/fibrinogen-audit.json and models/proteins/README.md for chain-specific gaps, geometry and modeling limitations.
