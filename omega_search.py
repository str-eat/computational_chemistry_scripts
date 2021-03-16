import glob
import numpy as np
files = glob.glob('*com')
omega = np.linspace(200, 300, 11, dtype=int)
print(omega)
con = input("If this is correct please press enter")
omega_lens = [len(str(o)) for o in omega]
omega_lefts = [10 - ln for ln in omega_lens]
for file in files:
 for k, o in enumerate(omega):
  with open(file, 'r') as f:
   with open(file[:-4] + "_omega0{}".format(o) + file[-4:], 'w') as w:
    for line in f:
     if "Toluene" in line:
      w.write(line + "iop(3/107=" + "0"*(5-omega_lens[k]) + "{}".format(o) + "0"*(10-omega_lens[k]-(5-omega_lens[k])) + ")"\
              + " iop(3/108=" + "0"*(5-omega_lens[k]) + "{}".format(o) + "0"*(10-omega_lens[k]-(5-omega_lens[k])) + ")\n")
     else:
      w.write(line)
