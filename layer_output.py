try:
    import fragment_parser as utils
except ImportError:
    from . import fragment_parser as utils

import math
import glob

filename = input("Enter filename: ")
output_filename = input("Enter output filename: ")
filenames = glob.glob(filename)

with open(output_filename, 'w+') as f:
    for filename in filenames:
        molecule = utils.open_file(filename)
        print('\n', filename, '\n')
        f.write('layer_' + filename + '\n')
        all_info = list()
        for layer in utils.array_of_separate_molecules(molecule):
            all_info.append([layer.GetAtom(each).GetId()+1 for each in range(1, layer.NumAtoms()+1)])
           
        f.write(str(all_info) + '\n') 
    

"""
    separate_layer_imine_array = utils.separate_fragment_by_layer(imine_array, molecule)
    for each in separate_layer_imine_array:
        print("Separated Imine: ", [each2.GetId() for each2 in each])
    separate_layer_bridge_array = utils.separate_fragment_by_layer(bridge_array, molecule)
    for each in separate_layer_bridge_array:
        for each2 in each:
            print("Separated Bridge: ", each2.GetId())
    separate_layer_node_array = utils.separate_fragment_by_layer(node_array, molecule)
    for each in separate_layer_node_array:
        print("Separated Node: ", [each2.GetId() for each2 in each])
    separate_layer_other_array = utils.separate_fragment_by_layer(other_array, molecule)
    for each in separate_layer_other_array:
        print("Separated Phenyl: ", [each2.GetId() for each2 in each])
        """

        
