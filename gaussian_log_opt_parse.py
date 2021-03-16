import glob
import cclib
filenames = input("Enter regex for gaussian optimization log files here: ")
if not filenames:
    filenames = "*log"
files = glob.glob(filenames)
for file in files:
    parser = cclib.io.ccopen(file)
    data=parser.parse()
    data.metadata['comments'] = data.metadata['methods']
    cclib.io.ccwrite(data, outputtype='xyz', outputdest=file[:-4]+'_optimized.xyz', indices=[-1])
