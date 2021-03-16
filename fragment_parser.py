try:
    from openbabel import openbabel, pybel
except ImportError:
    import openbabel, pybel
import matplotlib.pyplot as plt
from scipy.stats import circvar, circmean
import numpy as np
from math import sqrt


def open_file(filename):
    try:
        molecule = pybel.readfile(filename[-3:], filename).__next__().OBMol
    except AttributeError:
        obconversion = openbabel.OBConversion()
        obconversion.SetInAndOutFormats(filename[:-3], filename[:-3])
        molecule = openbabel.OBMol()
        obconversion.ReadFile(molecule, filename)
    return molecule

def array_of_separate_molecules(OBMolecule):
    molecule = OBMolecule
    separate_molecules = list(molecule.Separate())
    return separate_molecules

def array_of_imine_nitrogens(OBMolecule):
    molecule = OBMolecule
    imine_nitrogens = list()
    for atom in openbabel.OBMolAtomIter(molecule):
        atomic_num = atom.GetAtomicNum()
        atom_in_ring_size = atom.MemberOfRingSize()
        if atomic_num == 7:
            if atom_in_ring_size == 0:
                imine_nitrogens.append(atom)
    return imine_nitrogens

def array_of_imine_carbons(OBMolecule):
    molecule = OBMolecule
    imine_carbons = list()
    for atom in openbabel.OBMolAtomIter(molecule):
        atomic_num = atom.GetAtomicNum()
        atom_in_ring_size = atom.MemberOfRingSize()
        if atomic_num == 6:
            if atom_in_ring_size == 0:
                imine_carbons.append(atom)
    return imine_carbons

def array_of_imines(molecule):
    imines = list()
    imine_carbons = array_of_imine_carbons(molecule)
    for carbon in imine_carbons:
        imines.append(carbon)
        for atom in openbabel.OBAtomAtomIter(carbon):
            if atom.MemberOfRingSize() > 0:
                continue
            imines.append(atom)
    return imines

def array_of_imine_hydrogens(molecule):
    imine = array_of_imines(molecule)
    hydrogens = list()
    for atom in imine:
        if atom.GetAtomicNum() == 1:
            hydrogens.append(atom)
    return(hydrogens)

def array_of_separate_layer_imines(molecule):
    separate_molecules = array_of_separate_molecules(molecule)
    separate_layer_imines = list()
    for molecule in separate_molecules:
        separate_layer_imines.append(array_of_imines(molecule))
    return separate_layer_imines

def array_of_others(molecule):
    
    def get_valence(valence_atom):
        count = 0
        for each in openbabel.OBAtomAtomIter(valence_atom):
            count = count + 1
        return count

    def return_fragment(atom, fragment_array):
        """
            start by adding the current atom to the fragment array
            then step through all the bonds of the atom adding terminal atoms to fragment array
            only step in the direction of an atom which has not yet been added to fragment array
            when all options are exhausted return fragment array
        """

        safe = False
        non_terminal = list()
        
        for valence_atom in openbabel.OBAtomAtomIter(atom):
            if get_valence(valence_atom) > 1:
                if valence_atom not in fragment_array and\
                    valence_atom not in non_terminal and valence_atom not in imine_atoms:
                    fragment_array.append(valence_atom)
                    non_terminal.append(valence_atom)
                    safe = True
            else:
                if valence_atom not in fragment_array and\
                    valence_atom not in imine_hydrogens:
                    fragment_array.append(valence_atom)

        if safe:
            for each in non_terminal:
                return_fragment(each, fragment_array)
        else:
            return True       
    
    bridge_array = list()
    node_array = list()

    for layer in array_of_separate_molecules(molecule):
        imine_carbons = array_of_imine_carbons(layer)
        imine_nitrogens = array_of_imine_nitrogens(layer)
        imine_hydrogens = array_of_imine_hydrogens(layer)
        imine_atoms = array_of_imines(layer)
        run1 = list()

        for carbon in imine_carbons:
            fragment_array = list()
            return_fragment(carbon, fragment_array)
            if run1:
                run2 = set(fragment_array)
            else:
                run1 = set(fragment_array)
        
        bridge = ([each.GetId() for each in run1].sort() == [each2.GetId() for each2 in run2].sort())
        
        run3 = list()

        for nitrogen in imine_nitrogens:
            fragment_array = list()
            return_fragment(nitrogen, fragment_array)
            if run3:
                run4 = set(fragment_array)
            else:
                run3 = set(fragment_array)

        if bridge:
            for each in run2:
                bridge_array.append(each)
            for each1, each2 in zip(run3, run4):
                node_array.append(each1)
                node_array.append(each2)
        else:
            for each in run2:
                node_array.append(each)
            for each1, each2 in zip(run3, run4):
                bridge_array.append(each1)
                bridge_array.append(each2)

    return bridge_array, node_array
        
def array_of_nodes(molecule):
     
    not_node = list()
    node = list()
        
    def get_valence(valence_atom):
        count = 0        
        for bonded_atom in openbabel.OBAtomAtomIter(valence_atom):                         
            if bonded_atom.GetAtomicNum() != 1:    
                count = count + 1
        return count
        
    def get_hydrogens(in_atom):
        for valence in openbabel.OBAtomAtomIter(in_atom):
            if valence.GetAtomicNum() == 1:
                hydrogen = valence
                return hydrogen

    for layer in array_of_separate_molecules(molecule):
        for ring in openbabel.OBMolRingIter(molecule):            
            bridge_array, node_array = array_of_others(layer)    
            atom_array = []
            first = False
            second = False 
            for atomi in node_array:
                atom = molecule.GetAtomById(atomi.GetId())
                if ring.IsMember(atom):
                    if get_valence(atom) == 3 and first:
                        second = True                    
                        atom_array.append(atom)
                    elif get_valence(atom) == 3:
                        first = True
                        atom_array.append(atom)
                    else:
                        atom_array.append(atom)

            if second and atom_array:                
                for each in atom_array:
                    hydrogen = get_hydrogens(each)
                    if hydrogen:
                        atom_array.append(hydrogen)
                for a in atom_array:
                    not_node.append(a)

            elif first and atom_array:
                for each in atom_array:
                    hydrogen = get_hydrogens(each)
                    if hydrogen:
                        atom_array.append(hydrogen)
                for a in atom_array:
                    node.append(a)
    return node, not_node

def separate_fragment_by_layer(array, molecule):
    layer_fragment = list()
    layers = array_of_separate_molecules(molecule)
    for layer in layers:
        fragment = list()
        for atom in array:
            if layer.GetAtomById(atom.GetId()):                
                fragment.append(molecule.GetAtomById(atom.GetId()))                            
    return layer_fragment
