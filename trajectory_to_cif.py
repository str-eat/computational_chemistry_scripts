from ase.io import Trajectory
from sys import argv
infile = input("Enter trajectory file: ")
step = input("Enter the trajectory step you would like to output: ")
traj = Trajectory(infile)
print("There are {} steps in the trajectory file.".format(len(traj)))
traj[int(step)].write('traj_output.cif')
