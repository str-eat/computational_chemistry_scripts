import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from gpaw import GPAW
import sys

try:
	print(sys.argv[1])
	if sys.argv[1] is None:
		filename = input('Enter the wavefunction file (.gpw) here: ')
	else:
		filename = sys.argv[1]
except IndexError:
	filename = input('Enter the wavefunction file (.gpw) here: ')

calc = GPAW(filename, txt=None)
sz1, sz2 = calc.get_orbital_ldos(a=0, angular='d')
print(sz1)
print(len(sz1))
print(sz2)
print(len(sz2))

sz1 = len(sz1)
sz2 = len(sz2)


C = np.empty(shape=(sz1))
C_dos = np.empty(shape=(sz1))
O = np.empty(shape=(sz1))
O_dos = np.empty(shape=(sz1))
Zn = np.empty(shape=(sz1))
Zn_dos = np.empty(shape=(sz1))
Cu = np.empty(shape=(sz1))
Cu_dos = np.empty(shape=(sz1))

C_atoms = 0
O_atoms = 0
Zn_atoms = 0
Cu_atoms = 0

for an, atom in enumerate(calc.atoms):
	for orb in ['s', 'p', 'd']:
		energy, pdos = calc.get_orbital_ldos(a=atom.index, angular=orb)
		if atom.symbol == 'C':
			C_atoms = C_atoms + 1     
			C = C + energy
			C_dos = C_dos + pdos
		if atom.symbol == 'O':
			O_atoms = O_atoms + 1
			O = O + energy
			O_dos = O_dos + pdos
		if atom.symbol == 'Zn':
			Zn_atoms = Zn_atoms + 1
			Zn = Zn + energy
			Zn_dos = Zn_dos + pdos
		if atom.symbol == 'Cu':
			Cu_atoms = Cu_atoms + 1
			Cu = Cu + energy
			Cu_dos = Cu_dos + pdos

"""Cres = np.sum(C, 1)
Cdosres = np.sum(C_dos, 1)
Ores = np.sum(O, 1)
Odosres = np.sum(O_dos, 1)
Znres = np.sum(Zn, 1)
Zndosres = np.sum(Zn_dos, 1)
Cures = np.sum(Cu, 1)
Cudosres = np.sum(Cu_dos, 1)"""
C = C/C_atoms
O = O/O_atoms
Zn = Zn/Zn_atoms
Cu = Cu/Cu_atoms
C_dos = C_dos/C_atoms
O_dos = O_dos/O_atoms
Zn_dos = Zn_dos/Zn_atoms
Cu_dos = Cu_dos/Cu_atoms

plt.plot(C, C_dos, label='C orbitals')
plt.plot(O, O_dos, label='O orbitals')
plt.plot(Zn, Zn_dos, label='Zn orbitals')
plt.plot(Cu, Cu_dos, label='Cu orbitals')
plt.show()
plt.savefig('pdos.png')

np.savetxt("C_energies.csv", C, delimiter=",")
np.savetxt("O_energies.csv", O, delimiter=",")
np.savetxt("Zn_energies.csv", Zn, delimiter=",")
np.savetxt("Cu_energies.csv", Cu, delimiter=",")
np.savetxt("C_dos.csv", C_dos, delimiter=",")
np.savetxt("O_dos.csv", O_dos, delimiter=",")
np.savetxt("Zn_dos.csv", Zn_dos, delimiter=",")
np.savetxt("Cu_dos.csv", Cu_dos, delimiter=",")
