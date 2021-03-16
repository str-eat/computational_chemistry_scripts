import sys
data = ''
try:
    if sys.argv[1]:
        arg1 = sys.argv[1]
    if sys.argv[2]:
        arg2 = sys.argv[2]
except IndexError:
     print("Note! TheoDORE input file needs to have same filename as the xyz coordinate files used to generate the fragments file!")
     arg1 = input("Please enter a TheoDORE input file: ")
     arg2 = input("Please enter a fragments file: ")

multiple = input("Does the fragments file have multiple xyz files? (Y/N)")
with open(arg1, 'r+') as f:
    right_file = False
    linker = False
    bridge = False
    node = False
    phenyl = False
    at_lists = list()
    for line in f:
        if line.split('=')[0] == 'at_lists':
            with open(arg2, 'r') as f2:
                for line2 in f2:
                    lists = list()
                    
                    if line2.split('.')[0] == arg1.split('.')[0]:
                        right_file = True
                    elif line2.split('.')[-1] == 'xyz\n':
                        right_file = False

                    if multiple=='N' or right_file:
                        if linker:
                            lists = [int(each) for each in line2.split()]
                        if bridge:
                            lists = [int(each) for each in line2.split()]
                        if node:
                            lists = [int(each) for each in line2.split()]
                        if phenyl:
                            lists = [int(each) for each in line2.split()]
 
                        if "Imine" in line2:
                            linker = True
                            bridge = False
                            node = False
                            phenyl = False
                        elif "Bridge" in line2:
                            bridge = True
                            linker = False
                            node = False
                            phenyl = False
                        elif "Node" in line2:
                            linker = False
                            bridge = False
                            node = True
                            phenyl = False
                        elif "Phenyl" in line2:
                            linker = False
                            bridge = False
                            node = False
                            phenyl = True
                        else:
                            linker = False
                            bridge = False
                            node = False
                            phenyl = False
    
                        if lists:
                            at_lists.append(lists) 
                newline = "at_lists="+str(at_lists)
                data = data + newline + '\n'
                continue
        data = data + line
with open(arg1, 'w+') as f:
    f.write(data)
                
