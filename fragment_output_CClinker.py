try:
    import fragment_parser_CClinker as utils
except ImportError:
    from . import fragment_parser_CClinker as utils

import math
import glob

filename = input("Enter filename: ")
output_filename = input("Enter output filename: ")
filenames = glob.glob(filename)

with open(output_filename, 'w+') as f:
    for filename in filenames:
        molecule = utils.open_file(filename)
        print('\n', filename, '\n')
        f.write(filename + '\n')

        linker_array = utils.array_of_linkers(molecule)
        print("Linker: \n" + " ".join([str(each.GetId() + 1) for each in linker_array]) + "\n")
        f.write("Linker: \n" + " ".join([str(each.GetId() + 1) for each in linker_array]) + "\n")
        
        bridge_array, node_array, rest_array = utils.array_of_others(molecule)
        print("Bridge: \n" + " ".join([str(each.GetId() + 1) for each in bridge_array]) + "\n")
        f.write("Bridge: \n" + " ".join([str(each.GetId() + 1) for each in bridge_array]) + "\n")
        print("Node: \n" + " ".join([str(each.GetId() + 1) for each in node_array]) + "\n")
        f.write("Node: \n" + " ".join([str(each.GetId() + 1) for each in node_array]) + "\n")
        print("Leftovers: \n" + " ".join([str(each.GetId() + 1) for each in rest_array]) + "\n")

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

        
