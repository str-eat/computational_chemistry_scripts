from ase.io import read
from sys import argv

filename = input("Enter .cif file you would like to convert: ")
cif = read(filename)

cif.write(filename[:-4]+".gen")
