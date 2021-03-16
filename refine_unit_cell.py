from ase.io import read, write
from ase import Atoms
import spglib as spg
from sys import argv

cif = read(argv[-1])
length_tolerance = input("What length tolerance would you like to use?? (in Angstroms): ")
angle_tolerance = input("What angle tolerance would you like to use?? (in degrees): ")
cif_cell = (cif.cell[:], cif.cell.scaled_positions(cif.positions), cif.get_atomic_numbers())
refined_cif = spg.refine_cell(cif_cell, symprec=float(length_tolerance), angle_tolerance=float(angle_tolerance))
print(spg.get_spacegroup(refined_cif))
proceed = input("Would you like to proceed? (Y/N) ")
if proceed == "Y":
	new_cif = Atoms(scaled_positions=refined_cif[1], numbers=refined_cif[2], cell=refined_cif[0])
	new_cif.write(str(argv[-1])[:-4] + "_refined.cif")
