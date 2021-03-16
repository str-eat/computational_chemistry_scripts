import cclib
import glob
files = input("Enter qchem output glob here: ")
files = glob.glob(files)
files.sort(reverse=False)
nut_files = files[0:int(len(files)/2)]
cat_files = files[int(len(files)/2):]
homs = list()
if files:
    print("HOMO Neutral, E Neutral, E Cation")
    for nj, cf in zip(nut_files, cat_files):
        om = ""
        with open(nj, 'r') as ff:
            for ln in ff:
                if 'omega' in ln:
                    om = ln.split()[-1]
        nut = cclib.io.ccopen(nj)
        cat = cclib.io.ccopen(cf)
        nut = nut.parse()
        cat = cat.parse()
        try:
            homonut = nut.moenergies[0][nut.homos[0]]
            homocat = cat.moenergies[0][cat.homos[0]]
        except:
            continue
        homo_energy = nut.scfenergies[-1]
        cat_energy  = cat.scfenergies[-1]
        ip = (-1 * homonut) - (cat_energy - homo_energy)
        print(om, homonut, homo_energy, cat_energy, ip)


