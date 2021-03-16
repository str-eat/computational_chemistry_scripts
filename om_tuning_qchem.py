import numpy as np
import os
import cclib
import math
import time
oms = np.linspace(100, 200, 9)

def parse_logfiles(n, c, a):
 nut = cclib.io.ccopen(n)
 cat = cclib.io.ccopen(c)
 ani = cclib.io.ccopen(a)
 nut = nut.parse()
 cat = cat.parse()
 ani = ani.parse()
 return nut, cat, ani

def get_homo_ens(nut, cat, ani):
 hn = nut.moenergies[0][nut.homos[0]]
 hc = cat.moenergies[0][cat.homos[0]]
 ha = ani.moenergies[0][ani.homos[0]]
 if ha > ani.moenergies[1][ani.homos[0]]:
  ha = ani.moenergies[1][ani.homos[0]]
 return hn, hc, ha
 
def get_scf_ens(nut, cat, ani):
 ne = nut.scfenergies[-1]
 ce = cat.scfenergies[-1]
 ae = ani.scfenergies[-1]
 return ne, ce, ae 

for n in ['N', 'P', 'M']:
 for o in oms:
  with open('{}.in'.format(n), 'r') as op:
   with open('{}_{}.in'.format(n, str(int(o))), 'w') as wr:
    for ln in op:
     if 'omega' in ln.lower() and 'exchange' not in ln.lower() and 'method' not in ln.lower():
      wr.write('omega     {}\n'.format(str(int(o))))
     else:
      wr.write(ln)
cwd = os.getcwd()
print(cwd)

ens = list()
for o in oms:
 o = int(o)
 n = 'N_{}.out'.format(o)
 p = 'P_{}.out'.format(o)
 m = 'M_{}.out'.format(o)
 om = ""
 with open(n[:-4]+'.in', 'r') as ff:
  for ln in ff:
   if 'omega' in ln:
    om = ln.split()[-1]
    break
 try:
  nl, pl, ml = parse_logfiles(n, p, m)
  homonut, homocat, homoani = get_homo_ens(nl, pl, ml)
  nut_energy, cat_energy, ani_energy = get_scf_ens(nl, pl, ml)
 except:
  rc = os.system("/mmfs1/home/2778streatd/bin/subqchem-optomega {}/N_{}.in & /mmfs1/home/2778streatd/bin/subqchem-optomega {}/P_{}.in & /mmfs1/home/2778streatd/bin/subqchem-optomega {}/M_{}.in".format(cwd, o, cwd, o, cwd, o))
  time.sleep(30)
  try:
   time.sleep(5)
   nl, pl, ml = parse_logfiles(n, p, m)
   homonut, homocat, homoani = get_homo_ens(nl, pl, ml)
  except AttributeError:
   time.sleep(15)
   nl, pl, ml = parse_logfiles(n, p, m)
   homonut, homocat, homoani = get_homo_ens(nl, pl, ml)
  except:
   time.sleep(30)
   homonut, homocat, homoani = get_homo_ens(nl, pl, ml)
  nut_energy, cat_energy, ani_energy = get_scf_ens(nl, pl, ml)
  
 ip_diff = homonut + (cat_energy - nut_energy)
 ea_diff = homoani + (nut_energy - ani_energy)
 j = math.sqrt((ip_diff**2)+(ea_diff**2))
 ens.append(j)

three_mins = sorted(set(ens))[0:3]
print(three_mins)
idx = [ens.index(t) for t in three_mins]
oms_idx = oms[idx]
min_om = min(oms_idx)
max_om = max(oms_idx)
new_oms = np.linspace(min_om, max_om, 11)
for n in ['N', 'P', 'M']:
 for o in new_oms:
  with open('{}.in'.format(n), 'r') as op:
   with open('{}_{}.in'.format(n, str(int(o))), 'w') as wr:
    for ln in op:
     if 'omega' in ln and 'exchange' not in ln:
      wr.write('omega     {}\n'.format(str(int(o))))
     else:
      wr.write(ln)
cwd = os.getcwd()
for o in new_oms:
 o = int(o)
 n = 'N_{}.out'.format(o)
 p = 'P_{}.out'.format(o)
 m = 'M_{}.out'.format(o)
 try:
  homonut, homocat, homoani = get_homo_ens(n, p, m)
 except:
  rc = os.system("/mmfs1/home/2778streatd/bin/subqchem-optomega {}/N_{}.in & /mmfs1/home/2778streatd/bin/subqchem-optomega {}/P_{}.in & /mmfs1/home/2778streatd/bin/subqchem-optomega {}/M_{}.in".format(cwd, o, cwd, o, cwd, o))

