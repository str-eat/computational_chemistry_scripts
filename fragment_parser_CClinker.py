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
        obconversion.SetInAndOutFormats(filename[-3:], filename[-3:])
        molecule = openbabel.OBMol()
        obconversion.ReadFile(molecule, filename)
    return molecule

def array_of_linker_carbons(OBMolecule):
    molecule = OBMolecule
    linker_carbons = list()
    for atom in openbabel.OBMolAtomIter(molecule):
        atomic_num = atom.GetAtomicNum()
        atom_in_ring_size = atom.MemberOfRingSize()
        if atomic_num == 6:
            if atom_in_ring_size == 0:
                linker_carbons.append(atom)
    return linker_carbons

def array_of_linkers(molecule):
    linkers = list()
    linker_carbons = array_of_linker_carbons(molecule)
    for carbon in linker_carbons:
        if carbon.GetId() not in [each.GetId() for each in linkers]:
            linkers.append(carbon)
        for atom in openbabel.OBAtomAtomIter(carbon):
            if atom.MemberOfRingSize() > 0:
                continue
            if atom.GetId() not in [each.GetId() for each in linkers]:
                linkers.append(atom)
    return linkers

def array_of_linker_hydrogens(molecule):
    linker = array_of_linkers(molecule)
    hydrogens = list()
    for atom in linker:
        if atom.GetAtomicNum() == 1:
            hydrogens.append(atom)
    return(hydrogens)

def array_of_others(molecule):

    array1 = list()
    array2 = list()
    array_catch_all = list()

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
                    valence_atom not in non_terminal and valence_atom not in linker_atoms:
                    fragment_array.append(valence_atom)
                    non_terminal.append(valence_atom)
                    safe = True
            else:
                if valence_atom not in fragment_array and\
                    valence_atom not in linker_atoms:
                    fragment_array.append(valence_atom)

        if safe:
            for each in non_terminal:
                return_fragment(each, fragment_array)
        else:
            return True       
    

    linker_atoms = array_of_linkers(molecule)
    fragment1_size = 0
    fragment2_size = 0
    for carbon in linker_atoms:
        fragment_array = list()
        return_fragment(carbon, fragment_array)
        if fragment1_size == 0:
            fragment1_size = len(set(fragment_array))
        elif fragment2_size == 0:
            fragment2_size = len(set(fragment_array))
        if fragment1_size == len(set(fragment_array)):
            for each in set(fragment_array):
                if each not in array1:
                    array1.append(each)
        elif fragment2_size == len(set(fragment_array)):
            for each in set(fragment_array):
                if each not in array2:
                    array2.append(each)
        else:
            for each in set(fragment_array):
                array_catch_all.append(each)
            
    return array1, array2, array_catch_all

def separate_fragment_by_layer(array, molecule):
    layer_fragment = list()
    layers = array_of_separate_molecules(molecule)
    for layer in layers:
        fragment = list()
        for atom in array:
            if layer.GetAtomById(atom.GetId()):                
                fragment.append(molecule.GetAtomById(atom.GetId()))                            
    return layer_fragment
