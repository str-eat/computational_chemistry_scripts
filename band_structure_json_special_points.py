import pandas as pd
import numpy as np
filename = input("Enter filename: ")
d = pd.read_json(filename)
nrg = d['energies']['__ndarray__']
try:
 nrg2 = np.reshape(nrg[2], (nrg[0][1], nrg[0][2]))
except ValueError:
 nrg2 = np.reshape(nrg[2], (nrg[0][1], nrg[0][2]*nrg[0][0]))
 print("Spin Polarized")

sp = d['path']['special_points']
sp = [each['__ndarray__'][2] for each in list(sp.values())]
sp = np.reshape(sp, (len(sp), len(sp[0])))
kpts = np.reshape(d['path']['kpts']['__ndarray__'][2], (d['path']['kpts']['__ndarray__'][0][0], d['path']['kpts']['__ndarray__'][0][1]))
for each in sp:
 for key, each2 in enumerate(kpts):
  if (each == each2).all():
    print(key + 1)

print('Fermi Level:', d['reference'][0])
print('Path:', d['path']['labelseq'])
rng = np.arange(1, len(nrg2)+1)
try:
 nrg2 = np.column_stack((rng, nrg2))
except ValueError:
 rng = np.arange(1, len(nrg2[1])+1)
 nrg2 = np.column_stack((rng, nrg2[1]))

np.savetxt(filename[:-5] + filename[0:3] + '.csv', nrg2, delimiter=',')
