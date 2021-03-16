from ase.io import read

filename = input("Enter cif file: ")

cif = read(filename)

old_scaled_positions = cif.get_scaled_positions()

print("Old cell parameters are:\n\n"+str(cif.cell))

perc = float(input("How much would you like to scale your cell by (in %)? "))
print(perc/100)
new_cell = cif.cell * (perc/100)
print(str(new_cell))
cif.cell = new_cell
cif.set_scaled_positions(old_scaled_positions)
cif.cell[2][2]
cif.write('TAT_ethyl_EXPscaled.cif')

