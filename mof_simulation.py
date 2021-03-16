from ase.io import read
from gpaw import GPAW, PW, Davidson
from ase.constraints import ExpCellFilter
from ase.optimize import BFGS
from ase.io import Trajectory
from ase.parallel import paropen, world, parprint
from ase.calculators.dftd3 import DFTD3
from ase.build import make_supercell
import numpy as np
import os.path
from os import path
import os
import itertools
import time

def setup_structure_calculation():
    """
    Edit these variables to change the conditions of the calculation

    Parameters: None

    Returns:
        starting_structure_file - the .cif file to read for starting the calculation
        pw_cutoff               - the planewave cutoff (in eV) to use for calculating the total energy of the systems
        xc_func                 - the exchange-correlation functional to be used in the structure optimization
        x_kpts, y_kpts, z_kpts  - the parameters to defined the k-points sampling of the Brillouin zone
                                   in the determination of the total energy of the systems
        txt_output              - the name of the file to which the output of the calculations is written
        opt_log                 - the name of the file to which the QuasiNewton optimization algorithm is written
        opt_traj                - the name of the trajectory file to which the QuasiNewton optimization algorithm
                                   is writted
        convergence             - the convergence criteria of the QuasiNewton optimization algorithm (in eV/A)
    """
    starting_structure_file = "Cu_S_Pyrene_tbutyl.cif" #index 0
    cif = read(starting_structure_file)
    pw_cutoff = 600 #index 1
    xc_func = "PBEsol" #index 2
    x_kpts = 1 #index 3
    y_kpts = 1 #index 4
    z_kpts = 3 #index 5
    txt_output = starting_structure_file[:-4] + ".txt" #index 6
    opt_log = starting_structure_file[:-4] + ".opt" #index 7
    opt_traj = starting_structure_file[:-4] + ".traj" #index 8
    convergence = 0.05
    min_x_coord = 0 #index 10
    max_x_coord = cif.cell[0][0]/2 #index 11
    num_x_coord = int(max_x_coord) #index 12
    min_y_coord = 0 #index 13
    max_y_coord = cif.cell[1][1]/2 #index 14
    num_y_coord = int(max_y_coord) #index 15
    min_z_coord = (2*cif.cell[2][2])/3 #index 16
    max_z_coord = cif.cell[2][2] #index 17
    num_z_coord = int((max_z_coord-min_z_coord)/0.5) #index 18
    return [starting_structure_file, pw_cutoff, xc_func, x_kpts, y_kpts, z_kpts, txt_output, opt_log, opt_traj, convergence,
            min_x_coord, max_x_coord, num_x_coord, min_y_coord, max_y_coord, num_y_coord, min_z_coord, max_z_coord, num_z_coord]

def write_disp_energy_file(ase_filename):
    with open('ase_dftd3.out', 'r') as disp:
        with open(ase_filename, 'w') as disp_out:
            for line in disp:
                disp_out.write(line)

def get_disp_energy_from_file(ase_filename, wait):
    if wait > 100:
        raise ValueError("Waited for 10 cycles")
    disp_energy = list()
    try:
        with paropen(ase_filename, 'r') as f:
            for line in reversed(list(f)):
                if 'Edisp' in line:
                    disp_energy.append(float(line.split()[-1]))
            if disp_energy[-1]: 
                return disp_energy[-1]
    except:
        return get_disp_energy_from_file(ase_filename, wait+1)

def get_disp_energy(structure, filename_base, ase_filename):
    cif = structure
    d3 = DFTD3(dft = calc, damping='bj',\
        command='/mmfs1/home/2778streatd/Software/dftd3/prg/dftd3',\
        abc=True)
    cif.set_calculator(d3)
    cif.get_potential_energy()
    write_disp_energy_file(ase_filename)
    disp_energy = get_disp_energy_from_file(ase_filename, 1)
    return disp_energy

def get_dft_energy_from_file(filename_base, wait):
    if wait > 100:
        raise ValueError("Waited for 100 cycles")
    dft_energy = list()
    try:
        with paropen(filename_base + '.txt', 'r') as f:
            for line in f:
                if 'Extrapolated' in line:
                    dft_energy.append(float(line.split()[-1]))
            if dft_energy[-1]:
                return dft_energy[-1]
    except:
        return get_dft_energy_from_file(filename_base, wait+1)

def get_dft_energy(structure, filename_base, ase_filename):
    cif = structure
    calc = GPAW(mode=PW(pw_cutoff, gammacentered=False, fftwflags=1, dedecut='estimate'),\
        basis='dzp',\
        eigensolver=Davidson(4),\
        kpts={'size': (x_kpts, y_kpts, z_kpts)},\
        xc=xc_func,\
        symmetry='off',\
        spinpol=True,\
        maxiter=500,\
        txt=filename_base + '.txt',
        parallel={'domain': world.size})
    d3 = DFTD3(dft = calc, damping='bj',\
        command='/mmfs1/home/2778streatd/Software/dftd3/prg/dftd3',\
        abc=True)
    cif.set_calculator(d3)
    dft = float(cif.get_potential_energy())
    write_disp_energy_file(ase_filename)
    return dft

def get_relaxed_dft_energy(structure, filename_base, ase_filename):
    from ase.constraints import FixCom
    cif = structure.copy()
    c = FixCom()
    cif.set_constraint(c)

    calc = GPAW(mode=PW(pw_cutoff, gammacentered=False, fftwflags=1, dedecut='estimate'),\
        basis='dzp',\
        eigensolver=Davidson(4),\
        kpts={'size': (x_kpts, y_kpts, z_kpts)},\
        xc=xc_func,\
        symmetry='off',\
        spinpol=True,\
        maxiter=500,\
        txt=filename_base + '.txt',
        parallel={'domain': world.size})
    cif.set_calculator(calc)
    d3 = DFTD3(dft=calc, damping='bj',\
        command='/mmfs1/home/2778streatd/Software/dftd3/prg/dftd3',\
        abc=True)
    cif.set_calculator(d3)
    qn = BFGS(cif, logfile=filename_base + '.log')
    traj = Trajectory(filename_base + ".traj", 'w', cif, master=True)
    qn.attach(traj)
    ol = paropen(filename_base + ".log", 'w')
    qn.run(convergence)
    cif.write(filename_base + '.cif')
    dft = float(cif.get_potential_energy())
    return dft

def get_nuclear_z_variance(structure):
    nuc_var = structure.positions[:,2].var()
    return nuc_var

def find_minimum_energy(xy_coords, tot_energies):
    min_energy_index = tot_energies.index(min(tot_energies))
    min_xy = xy_coords[min_energy_index]
    return min_xy

def run_layer_optimization(structure, layer_opts):
    cif = structure.copy()
 
    if layer_opts > 0:
        ecf = ExpCellFilter(cif, mask=[True, True, True, True, True, True])
        filename_base = txt_output[:-4]+'_final'
        ase_filename = 'ase_dftd3_final.out'
    else:
        if path.exists(txt_output[:-4]+"_layer_optimized.cif"):
            return read(txt_output[:-4]+"_layer_optimized.cif")
        ecf = ExpCellFilter(cif, mask=[True, True, False, False, False, False])
        filename_base = txt_output[:-4]+'_layer'
        ase_filename = 'ase_dftd3_layer.out'

    calc = GPAW(mode=PW(pw_cutoff, gammacentered=False, fftwflags=1, dedecut='estimate'),
        xc=xc_func, eigensolver=Davidson(4), kpts={'size': (x_kpts, y_kpts, z_kpts), 'gamma': False},
        symmetry='off', spinpol=True, maxiter=500, basis='dzp', txt=filename_base+'.txt')
    d3 = DFTD3(dft=calc, damping='bj',\
                        command='/mmfs1/home/2778streatd/Software/dftd3/prg/dftd3',\
                        abc=True)
    cif.set_calculator(d3)

    qn = BFGS(ecf, logfile=opt_log)
    traj = Trajectory(filename_base+'.traj', 'w', cif, master=True) 
    qn.attach(traj)
    ol = paropen(filename_base+'.log', 'w') 
    qn.run(convergence)
    write_disp_energy_file(ase_filename)

    try:
        if layer_opts > 0:
            dft = float(get_dft_energy_from_file(filename_base, 1))
            disp_energy = float(get_disp_energy_from_file('ase_dftd3_final.out', 1))
            opt_cif = Trajectory(opt_traj)[-1]
            opt_cif.write(txt_output[:-4]+"_layer_optimized_final.cif")
        else:
            dft = float(get_dft_energy_from_file(filename_base, 1))
            disp_energy = float(get_disp_energy_from_file('ase_dftd3_eclipsed.out', 1))
            opt_cif = Trajectory(opt_traj)[-1]
            opt_cif.write(txt_output[:-4]+"_layer_optimized.cif")
    except IndexError:
        opt_cif = cif
        opt_cif.write(txt_output[:-4]+"_layer_optimized.cif")
    layer_opts = layer_opts + 1
    return opt_cif

def run_il_optimization(structure, type='eclipsed'):
    scale = False

    if type != 'staggered':
        z_offsets = np.linspace(max_z_coord, min_z_coord, num_z_coord, endpoint=True)
        if type == 'inclined':
            tx = txt_output[:-4] + '_inclined'
            asf = 'ase_dftd3_inclined'
        else:
            tx = txt_output[:-4]
            asf = 'ase_dftd3'
    elif type == 'staggered':
        scale = True
        z_offsets = np.linspace(max_z_coord*2, min_z_coord*2, num_z_coord, endpoint=True)
        tx = txt_output[:-4] + '_staggered'
        asf = 'ase_dftd3_staggered'

    dft_energies = list()
    d3_energies = list()
    total_energies = list()

    for k, z in enumerate(z_offsets):
        filename_base = tx + '_z{}'.format(z)
        ase_filename = asf + '_z{}.out'.format(z)
        
        if k == 0:
            cif = structure.copy()
        else:
            cif = read(tx + '_z{}.cif'.format(z_offsets[k-1]))
        
        if type != 'staggered':
            nuc_var = get_nuclear_z_variance(cif)
        else:
            nuc_var = 0

        cif.set_cell([[cif.cell[0][0], cif.cell[0][1], cif.cell[0][2]],\
                [cif.cell[1][0], cif.cell[1][1], cif.cell[1][2]],\
                [cif.cell[2][0], cif.cell[2][1], z]], scale_atoms=scale)
        
        try:
            dft = float(get_dft_energy_from_file(filename_base, 99))
        except:            
            if nuc_var > 0.1:
                dft = get_relaxed_dft_energy(cif, filename_base, ase_filename)
                write_disp_energy_file(ase_filename)
            elif nuc_var <= 0.1:
                dft = float(get_dft_energy(cif, filename_base, ase_filename))

        dft_energies.append(dft)
        cif.write(filename_base + '.cif')

        try:
            disp_energy = float(get_disp_energy_from_file(ase_filename, 1))
        except:
            disp_energy = float(get_disp_energy(cif, filename_base, ase_filename))
        d3_energies.append(disp_energy)

    total_energies = [float(dft) + float(d3) for dft, d3 in zip(dft_energies, d3_energies)]
    min_idx = total_energies.index(min(total_energies))
    poly = np.polyfit(z_offsets[:min_idx-1-1:-1], total_energies[:min_idx-1-1:-1], 3)

    def fun(x):
        return poly[0]*x**3+poly[1]*x**2+poly[2]*x+poly[3]

    from scipy.optimize import fminbound
    x1 = z_offsets[-1]; x2 = z_offsets[min_idx-1-1]
    z_min = fminbound(fun, x1=x1, x2=x2)
    eclipsed_energy = fun(z_min)
    ts = z_offsets[min_idx]
    cif = read(tx + "_z{}.cif".format(ts))
    cif.set_cell([[cif.cell[0][0], cif.cell[0][1], cif.cell[0][2]], [cif.cell[1][0], cif.cell[1][1], cif.cell[1][2]], [0, 0, z_min]], scale_atoms=scale)
    cif.write(tx + "_z{}_interlayer_optimized.cif".format(z_min))
    return cif, eclipsed_energy

def run_inclined_scan(structure):

    x_y = list(); dft_energies = list(); d3_energies = list(); total_energies = list()
    
    alpha_offsets = np.linspace(min_x, max_x, num_x, endpoint=True)
    beta_offsets = np.linspace(min_y, max_y, num_y, endpoint=True)
    offsets = list(itertools.product(alpha_offsets, beta_offsets)) 
    for x, y in offsets:
        x_y.append([x, y])
        new_cif = structure.copy()
        new_cell = structure.cell
        new_cell[2][0] = x
        new_cell[2][1] = y
        new_cif.cell = new_cell
        filename_base = txt_output[:-4] + '_inclined_x{}_y{}'.format(x, y)
        ase_filename = 'ase_dftd3_inclined_x{}_y{}.out'.format(x, y)

        try:
            dft = float(get_dft_energy_from_file(filename_base, 99))
        except:
            nuc_var = get_nuclear_z_variance(new_cif)
            if nuc_var > 0.1:
                dft = get_relaxed_dft_energy(new_cif, filename_base, ase_filename)
                write_disp_energy_file(ase_filename)
            elif nuc_var <= 0.1:
                dft = float(get_dft_energy(new_cif, filename_base, ase_filename))
            new_cif.write(filename_base +'.cif')
        dft_energies.append(dft)
        try:
            disp_energy = float(get_disp_energy_from_file(ase_filename, 1))
            d3_energies.append(disp_energy)
        except:
            disp_energy = float(get_disp_energy(new_cif, filename_base, ase_filename))
            d3_energies.append(disp_energy)
    total_energies = [float(dft) + float(d3) for dft, d3 in zip(dft_energies, d3_energies)]
    
    min_xy = find_minimum_energy(x_y, total_energies)
    x = min_xy[0]
    y = min_xy[1]
    new_cell = structure.cell
    new_cell[2][0] = x
    new_cell[2][1] = y
    inclined_cif = structure.copy()
    inclined_cif.cell = new_cell
    inclined_cif.write(txt_output[:-4]+"_inclined_optimized.cif")
    return inclined_cif 

def run_staggered_scan(structure):
    x_y = list(); dft_energies = list(); d3_energies = list(); total_energies = list()

    x_offsets = np.linspace(min_x, max_x, int((num_x/2)+1), endpoint=True)
    y_offsets = np.linspace(min_y, max_y, int((num_y/2)+1), endpoint=True)
    offsets = list(itertools.product(x_offsets, y_offsets))

    for x, y in offsets:
        x_y.append([x, y])
        cif = structure.copy()
        cif = make_supercell(cif, [[1, 0, 0],[0, 1, 0],[0, 0, 2]])
        cif.center()
        layer1 = cif[cif.positions[:,2]>cif.cell[2][2]/2]
        nuc_var = get_nuclear_z_variance(layer1)
        layer2 = cif[cif.positions[:,2]<cif.cell[2][2]/2]
        layer1.positions[:,0] = layer1.positions[:,0] + x
        layer1.positions[:,1] = layer1.positions[:,1] + y
        new_cif = layer1 + layer2
        filename_base = txt_output[:-4] + '_staggered_x{}_y{}'.format(x, y)
        ase_filename = 'ase_dftd3_staggered_x{}_y{}.out'.format(x, y)
        try:
            dft = float(get_dft_energy_from_file(filename_base, 99))        
        except:
            if nuc_var > 0.1:
                dft = get_relaxed_dft_energy(new_cif, filename_base, ase_filename)
                write_disp_energy_file(ase_filename)
            elif nuc_var <= 0.1:
                dft = float(get_dft_energy(new_cif, filename_base, ase_filename))
            new_cif.write(filename_base + '.cif')
        dft_energies.append(dft)
        try:
            disp_energy = float(get_disp_energy_from_file(ase_filename, 1))
            d3_energies.append(disp_energy)
        except:
            disp_energy = float(get_disp_energy(new_cif, filename_base, ase_filename))
            d3_energies.append(disp_energy)
            
    total_energies = [float(dft) + float(d3) for dft, d3 in zip(dft_energies, d3_energies)]

    min_xy = find_minimum_energy(x_y, total_energies)
    x = min_xy[0]
    y = min_xy[1]
    cif = structure.copy()
    cif = make_supercell(cif, [[1, 0, 0],[0, 1, 0],[0, 0, 2]])
    cif.center()
    layer1 = cif[cif.positions[:,2]>cif.cell[2][2]/2]
    layer2 = cif[cif.positions[:,2]<cif.cell[2][2]/2]
    layer1.positions[:,0] = layer1.positions[:,0] + x
    layer1.positions[:,1] = layer1.positions[:,1] + y
    staggered_cif = layer1 + layer2
    min_energy_key = x_y.index([x, y])
    staggered_energy = total_energies[min_energy_key]
    staggered_cif.write(txt_output[:-4]+"_staggered_optimized.cif")
    return staggered_cif, staggered_energy

def run_inclined_2layer(structure):
    inclined_cif = structure.copy()
    inclined_2layer_cif = make_supercell(inclined_cif, [[1, 0, 0],[0, 1, 0],[0, 0, 2]])
    filename_base = txt_output[:-4] + '_inclined_2layer'
    ase_filename = 'ase_dftd3_inclined_2layer.out'
    try:
        dft = float(get_dft_energy_from_file(filename_base + '.txt', 99))
    except:
        dft = float(get_dft_energy(inclined_2layer_cif, filename_base, ase_filename))
    try:
        disp_energy = float(get_disp_energy_from_file(ase_filename, 1))
    except:
        disp_energy = float(get_disp_energy(inclined_2layer_cif, filename_base, ase_filename))

    return dft + disp_energy

def find_optimal_stacking(inclined_cif, staggered_cif, staggered_energy, eclipsed_energy):
    inclined_energy = run_inclined_2layer(inclined_cif)

    if staggered_energy < inclined_energy:
        parprint("Staggered more stable than inclined by {}eV".format(staggered_energy-inclined_energy))
        if staggered_energy < eclipsed_energy:
            parprint("Staggered more stable than eclipsed by {}eV".format(staggered_energy-eclipsed_energy))
            return staggered_cif
        else:
            parprint("Eclipsed more stable than staggered by {}eV".format(eclipsed_energy-staggered_energy))
            return eclipsed
    else:
        parprint("Inclined more stable than staggered by {}eV".format(inclined_energy-staggered_energy))
        if inclined_energy < eclipsed_energy:
            parprint("Inclined more stable than eclipsed by {}eV".format(inclined_energy-eclipsed_energy))
            return inclined_cif
        else:
            parprint("Eclipsed more stable than inclined by {}eV".format(eclipsed_energy-inclined_energy))

inputs = setup_structure_calculation()
print(inputs)
input_cif_filename = inputs[0]; pw_cutoff = inputs[1]; 
xc_func = inputs[2]; opt_log = inputs[7]; opt_traj = inputs[8]; txt_output = inputs[6]; convergence = inputs[9]
x_kpts = inputs[3]; y_kpts = inputs[4]; z_kpts = inputs[5]
min_x = inputs[10]; max_x = inputs[11]; num_x = inputs[12]; min_y = inputs[13]; max_y = inputs[14]; num_y = inputs[15]
min_z_coord = inputs[16]; max_z_coord = inputs[17]; num_z_coord = inputs[18]

input_cif = read(input_cif_filename)
layer_cif = run_layer_optimization(input_cif, 0)
cif, il_eclipsed_energy = run_il_optimization(layer_cif)
inclined_cif = run_inclined_scan(cif.copy())
staggered_cif, staggered_energy = run_staggered_scan(cif.copy())
il_staggered_cif, il_staggered_energy = run_il_optimization(staggered_cif, type='staggered')
il_inclined_cif, il_inclined_energy = run_il_optimization(inclined_cif, type='inclined')
cif = find_optimal_stacking(il_inclined_cif, il_staggered_cif, il_staggered_energy, il_eclipsed_energy)
run_layer_optimization(cif, 1)

