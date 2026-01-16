# MD simulations using AMOEBA Force Field for the study: "Using AMOEBA to Uncover Deformation and Hydration Anisotropy in Collagen Mimetic Peptides".

NOTE: Before the paper is published only Movie_S1.mp4 will be made available. The rest will be added afterwards.

### Structure files for the CHARMM MD simulations:

subdirectory `ppgn_charmm` with files as `ppg{n}_charmm.pdb` - in pdb file format, where `n` is the number of tripeptides (5, 12, 25). The same input files were used for CHARMM36+TIP3P, CHARMM36+TIP4P and CHARMM35mGP+TIP4P.

### Structure files for the AMBER MD simulations:

subdirectory `ppgn_amber` with files as `ppg{n}_amber.pdb` - in pdb file format, where `n` is the number of tripeptides (5, 12, 25). Used for AMBER(ff19sb)+TIP4P.

### Structure files for the AMOEBA MD simulations:

subdirectory `ppgn_amoeba` with files as `ppg{n}_{glc}mM_solvated_amoeba.pdb` - in Tinker XYZ file format, where `n` is the number of tripeptides (5, 12, 25) and `glc` is the concentration of glucose in mM (0, 500, 1000, 2000).

subdirectory `mutations` with files as `{mut}_{glc}.pdb` - in PDB file format, where `mut` is the name of mutation (ala_end, ala_mid, lys_end, lys_mid) and `glc` is the concentration of glucose in mM (0, 500, 1000, 2000).

### Parameter files for AMOEBA:

`parameters.prm` - common parameter file for all amoeba simulations.

`tinker.key` - sample `.key` file for amoeba, user needs to update box size.

### Python scripts for analysis:

subdirectory `python` with files:

`pyramidalization.py` - python functions to calculate condtions to satisfy $n \rightarrow \pi^*$ interactions.

`residence_time.py` - python functions to calculate residence time using thresholds defined in Impey's algorithm.

`waterdyn.py` - python code to use the quaternion method to quantify hydration dynamics.

### VMD scripts for analysis:

subdirectory `python` with files:

`vmd_lysglc.tcl` - glucose-lysine hbonds, water-lysine hbonds, glucose-lysine rdf, water-lysine rdf, glucoses around lysine, waters around lysine.

`vmd_pbwater_xxxsel_template` - Gly-Xaa(Pro) hbonds, CMP-water hbonds, waters around CMP with pbc, for a selection `xxxsel` (to be added by user).

`vmd_sharedwater.tcl` - get waters around protein, based on sections of CMP (`deformed` and `helical`), and check for waters shared between the section of collagen and glucose molecules.

`vmd_total` - various hbonds, backbone dihedrals, pyrroldine pucker angles, various rdfs, end-to-end distances, waters around protein (for residence time in `python/residence_time.py`), C $\alpha$ coordinates for deformations (in `python/HeliXplore.py`), coordinates of (N, C $\alpha$, C, O) atoms for $n \rightarrow \pi^*$ interactions (in `python/pyramidalization.py`), backbone RMSD, waters around glucose (for residence time in `python/residence_time.py`), glucoses around protein (for residence time in `python/residence_time.py`).

### SI Movie

Movie_S1.mp4 - Movie S1 of the SI.

### Citation

If you use this repository or scripts in your work, please cite the associated publications:

```Using AMOEBA to Uncover Deformation and Hydration Anisotropy in Collagen Mimetic Peptides, Ronnie Mondal and Valerie Vaissier Welborn, Journal Name, Year. [DOI/link to be added when available]```

See also: ```HeliXplore: A Python package for analyzing multi-strand helix deformations , Ronnie Mondal and Valerie Vaissier Welborn, Journal Name, Year. [DOI/link to be added when available]```
