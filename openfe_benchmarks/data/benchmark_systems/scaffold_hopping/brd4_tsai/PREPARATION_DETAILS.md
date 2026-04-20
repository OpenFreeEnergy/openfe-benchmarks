protein.pdb
- PDB ID 3mxf
- ligand and crystallographic solvent artifacts (EDO, DMS) were removed
- protein was prepared using pdbfixer (pdbfixer protein.pdb --add-atoms=all --add-residues --output=protein_fixed.pdb)

ligands.sdf
- PDB IDs 3MXF, 3SVG, 3U5J, 3U5L, 4E96, 4HBW, 4HXL, 4HXM, 4J0R, 4MR3, 4NUD, 4PCE, 4PCI, 4YH3, 5IGK, 5TI7, 5WUU were loaded into pymol, aligned onto ID 3MXF, and all ligands were extracted into a sdf.
- These PDBs are the ones that were used in Tsai et al. 10.1021/acs.jcim.5c02204
