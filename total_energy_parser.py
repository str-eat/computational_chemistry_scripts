def parse_z(single_file):
    zz = [0]
    z = 0
    for p in single_file.split('_'):
        try:
            try:
                if p[0] == 'z':
                    z = float(p[1:])
            except:
                if p[0] == 'z':
                    z = float(p[1:-4])
        except:
            pass
    zz[0] = z
    return zz

def parse_xy(single_file):
    xy = [0, 0]
    x = 0
    y = 0
    for p in single_file.split('_'):
        try:
            try:
                float(p[1:])
                if p[0] == 'x':
                    x = float(p[1:])
                elif p[0] == 'y':
                    y = float(p[1:])
            except:
                float(p[1:-4])
                if p[0] == 'x':
                    x = float(p[1:-4])
                elif p[0] == 'y':
                    y = float(p[1:-4])
        except:
            pass
    xy[0] = x
    xy[1] = y
    return xy

def find_tot_energy(dft_file, d3_file):
    with open(d3_file, 'r') as d:
        for line in reversed(list(d)):
            if "Edisp" in line:
                d3 = float(line.split()[-1])
                break
    with open(dft_file, 'r') as df:
        for line in reversed(list(df)):
            if "Extrapolated" in line:
                dft = float(line.split()[-1])
                break
    tot_energy = dft + d3
    #tot_energy = dft
    return float(tot_energy)

import glob

inclined = input("Enter regex for inclined dft files here: ")
incl_files = glob.glob(inclined)
if not incl_files:
    incl_files = glob.glob('*inclined*txt')
inc_output = list()
for inc in incl_files:
    try:
        xy = parse_xy(inc)
        x = xy[0]; y = xy[1]
        d3 = glob.glob('ase_dftd3*inclined*x{}*y{}*out'.format(x, y))[0]
        inc_output.append([x, y, find_tot_energy(inc, d3)])
    except:
        pass
inc_output.sort(key = lambda i: i[1])
inc_output.sort(key = lambda i: i[0])

staggered = input("Enter regex for staggered dft files here: ")
stag_files = glob.glob(staggered)
if not stag_files:
    stag_files = glob.glob('*staggered*txt')
stag_output = list()
for stag in stag_files:
    try:
        xy = parse_xy(stag)
        x = xy[0]; y = xy[1]
        d3 = glob.glob('ase_dftd3*stag*x{}*y{}*out'.format(x, y))[0]
        stag_output.append([x, y, find_tot_energy(stag, d3)])
    except:
        pass
stag_output.sort(key = lambda i: i[1])
stag_output.sort(key = lambda i: i[0])

interlayer = input("Enter regex for interlayer dft files here: ")
int_files = glob.glob(interlayer)
if not int_files:
    int_files = glob.glob('*_z[0-9]*txt')
int_output = list()
for intl in int_files:
    try:
        zz = parse_z(intl)
        z = zz[0]
        d3 = glob.glob('ase_dftd3*_z{}*out'.format(z))[0]
        int_output.append([z, find_tot_energy(intl, d3)])
    except:
        pass
int_output.sort(key = lambda i: i[0])

print("Inclined")
try:
    for inc in inc_output:
        print(inc[0], inc[1], inc[2])
except:
    pass
print("Staggered")
try:
    for stag in stag_output:
        print(stag[0], stag[1], stag[2])
except:
    pass
print("Interlayer")
try:
    for intl in int_output:
        print(intl[0], intl[1])
except:
    pass
