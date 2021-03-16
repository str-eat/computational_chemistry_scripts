import cclib
import glob
import math
files = input("Enter qchem output glob here: ")
files = glob.glob(files)
files.sort(reverse=False)
d1 = int(len(files)/3)
d2 = 2*int(len(files)/3)
ani_files = files[0:d1]
nut_files = files[d1:d2]
cat_files = files[d2:]
homs = list()

def get_scf_energy(of):
    with open(of, 'r') as o:
        for line in o:
            if 'Solute Internal Energy' in line:
                oe = float(line.split()[-1])
                return oe
    return "Error"

if files:
    print("OMEGA HOMO_Anion HOMO_Neutral E_Anion E_Neutral E_Cation IP_DIFF EA_DIFF, J")
    for ag, nj, cf in zip(ani_files, nut_files, cat_files):
        om = ""
        with open(nj, 'r') as ff:
            for ln in ff:
                if 'omega' in ln:
                    om = ln.split()[-1]
        nut = cclib.io.ccopen(nj)
        cat = cclib.io.ccopen(cf)
        ani = cclib.io.ccopen(ag)
        nut = nut.parse()
        cat = cat.parse()
        ani = ani.parse()
        try:
            homonut = nut.moenergies[0][nut.homos[0]]
            homocat = cat.moenergies[0][cat.homos[0]]
            homoani = ani.moenergies[0][ani.homos[0]]
            if homoani > ani.moenergies[1][ani.homos[0]]:
                homoani = ani.moenergies[1][ani.homos[0]]
        except:
            pass
        nut_energy = get_scf_energy(nj)
        cat_energy = get_scf_energy(cf)
        ani_energy = get_scf_energy(ag)
        #nut_energy = nut.scfenergies[-1]
        #cat_energy  = cat.scfenergies[-1]
        #ani_energy  = ani.scfenergies[-1]
        ip_diff = homonut + (cat_energy - nut_energy)
        ea_diff = homoani + (nut_energy - ani_energy)
        j = math.sqrt((ip_diff**2)+(ea_diff**2))
        print(om, homoani, homonut, ani_energy, nut_energy, cat_energy, ip_diff, ea_diff, j)


