import os
import sys
from glob import glob
convert = input("Do you want to convert the xyz files to input files? (Y/N) ")

user_input = sys.argv[1:]
#user_input = [glob(each) for each in user_input]
#user_input = user_input[0]
print(user_input)
if convert == 'Y':
    for input_file in user_input:
        #with open('./' + input_file[0:-4] + '_' + 'SP_GD3BJ_TD.com', 'w+') as calc_file:
        with open('./' + input_file[0:-4] + '_TD.com', 'w+') as calc_file:
            """
            calc_file.write('--Link1--\n\n'
                            'NProc=3\n' + 
                            'Mem=3GB\n' +
                            '%chk=' + input_file[0:-4] + '_' + 'Opt_GD3BJ.chk\n' + 
                            '#P B3LYP/6-311G** Opt EmpiricalDispersion=GD3BJ\n' +
                            '\n' + 
                            'Title Required\n' + 
                            '\n' + 
                            '0 1\n')
            with open(input_file, 'r+') as geom_file:
                for line in geom_file:
                    try:
                        int(line)
                    except ValueError:
                        #print(bool(line))
                        if line != '\n':
                            calc_file.write(line)
                calc_file.write('\n')
            """
            calc_file.write('--Link1--\n\n'
                            '%NProc=24\n' + 
                            '%Mem=96GB\n' +
                            #'%oldchk=' + input_file[0:-4] + '_' + 'Opt_GD3BJ.chk\n' +
                            '%chk=' + input_file[0:-4] + '_' + 'TD.chk\n' + 
                            '#P CAM-B3LYP/6-311G** TD(NStates=10) EmpiricalDispersion=GD3BJ Pop=Full GFINPUT iop(9/40=3)\n'                            + '\n' + 
                            'Title Required\n' + 
                            '\n' + 
                            '0 1\n')
            #calc_file.write('\n')

            with open(input_file, 'r+') as geom_file:
                for line in geom_file:
                    try:
                        int(line)
                        continue
                    except ValueError:
                        if line != '\n':
                            calc_file.write(line)
                calc_file.write('\n')



