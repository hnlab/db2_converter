1. Install notes

You need to download, compile & install AMSOL7.1 here. The executable should be named
amsol7.1

2. Location

http://comp.chem.umn.edu/amsol/

***
- Here provides ***`amsol7.1_patched`***, which is actually in use. small patches on AMSOL7.1, see https://wiki.docking.org/index.php/AMSOL
- All scripts here are from DOCK src
- `calc_solvation.csh` is the workflow for mol2 -> amsol -> out
- scripts for amsol71 input:
  - `make_amsol71_input.pl` # perl version
  - `make_amsol71_input.2.py` # raw make_amsol71_input.py (needs OpenEye)
  - `make_amsol71_input.py` # modified by qyfu to substitute OpenEye by RDKit


- scripts for amsol71 output:
  - `process_amsol_mol2.py` # output mol2 and amsol results
  - `mol2amsol.py` # mol2amsol utils for processing amsol and mol2

3. LOG (2024-07-14)
Modify amsol7.1/include/SIZES.i and re-compile:
```bash
# from
PARAMETER (MAXHEV=60, MAXLIT=60, NCHAIN=60, NRELAX=25)
# to
PARAMETER (MAXHEV=80, MAXLIT=110, NCHAIN=60, NRELAX=25)
```