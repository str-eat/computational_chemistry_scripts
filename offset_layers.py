from ase.io import read

filename = input("Enter filename: ")
atoms = read(filename)
top = atoms[atoms.positions[:, 2] > sum(atoms.positions[:, 2])/len(atoms.positions[:, 2])]
bottom = atoms[atoms.positions[:, 2] < sum(atoms.positions[:, 2])/len(atoms.positions[:, 2])]

xlow = float(input("Enter the lower offset x value: "))
xhigh = float(input("Enter the upper offset x value: "))
xstep = float(input("Enter the offset range x step: "))
ylow = float(input("Enter the lower offset y value: "))
yhigh = float(input("Enter the upper offset y value: "))
ystep = float(input("Enter the offset range y step: "))

for xpoint in [x * 0.01 for x in range(int(xlow/0.01), int(xhigh/0.01)+1, int(xstep/0.01))]:
    for ypoint in [y * 0.01 for y in range(int(ylow/0.01), int(yhigh/0.01)+1, int(ystep/0.01))]: 
        top.positions[:, 0] = top.positions[:, 0] + xpoint
        top.positions[:, 1] = top.positions[:, 1] + ypoint
        new_atoms = top + bottom
        new_atoms.write(filename[:-4]+"_Xoffset_"+str(xpoint)+"_Yoffset_"+str(ypoint)+".xyz")
        top.positions[:, 0] = top.positions[:, 0] - xpoint
        top.positions[:, 1] = top.positions[:, 1] - ypoint


