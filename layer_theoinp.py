import glob
# blank input files go here
filename = input("Enter input filename: ")
fragments = input("Enter fragment filename: ")
filenames = glob.glob(filename)
for sdfs in filenames:
    data = ''
    with open(sdfs, 'r+') as f:
        right_file = False
        for line in f:
            if line.split('=')[0] == 'at_lists':
                with open(fragments, 'r') as f2:
                    lists = ''
                    for line2 in f2:
                        # xyz filenames printed in fragments file must equal input filenames
                        if line2.split('.')[0] == sdfs.split('.')[0]:
                            right_file = True
                        elif line2.split('.')[-1] == 'xyz\n':
                            right_file = False
                        elif right_file == True:
                            lists = lists + line2
                    if lists:                        
                        newline = "at_lists="+lists
                        data = data + newline 
                    continue
            data = data + line
    with open(sdfs, 'w+') as f:
        f.write(data)
                
