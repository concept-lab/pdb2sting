# pdb2sting

Part of the preprocessing pipeline used in the STING database (to be published)

1. (PypKa) Adds correct protonation state to .pdb files taken directly from the Protein Data Bank
2. (fixtermini) Fixes the termini nomenclature to be compatible with AMBER-derived

Tested on single and multi chain proteins, protein-DNA and protein-RNA complexes.

## Dependencies

It requires the instalation of PypKa. Check the [documentation](https://pypka.readthedocs.io/en/latest/installation.html) for details.

## Usage

Assign the variable input_file (on line 4) to the path of the desired input PDB .pdb file.

```
python3 pdb2sting.py
```