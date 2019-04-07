import os



# the input of parameters starts here

# this is the simulation parameters


Lx = 16 					# system size in x direction
Ly = 16	 				# system size in y direction
Lz = 1 					# system size in z direction
BCx = 1 					# 1: PBC, 0: PBC in x direction
BCy = 1 					# 1: PBC, 0: PBC in y direction
BCz = 0 					# 1: PBC, 0: PBC in z direction
T = 0.1              # Temperature of the system
Ntherm = 10000       # number of sweeps for thermalization
Nbins = 100          # number of bins to collect
Nsweeps = 1000 	   # number of sweeps per bin
Jx = 1.              # hopping in x-direction
Jy = 1.              # hopping in y-direction
Delta = 0. 				# anisotropy
V = 1.						# off_diagonal disorder
epsilon = 0.5			# to ensure positive weights
h = 1. 					# diagonal disorder

run_time = 96*3600

number_of_realizations = 1  # number of disorder realizations for the disorder average

# make a executable out of the src files

os.system('g++ -std=c++0x -O3 -o SSE main.cpp')

if os.path.exists("qc_obs") == False:
  os.mkdir('qc_obs')

# write a new parameter file
ofFile = open('params', 'w')
ofFile.write(str(Lx)+'\n')
ofFile.write(str(Ly)+'\n')
ofFile.write(str(Lz)+'\n')
ofFile.write(str(BCx)+'\n')
ofFile.write(str(BCy)+'\n')
ofFile.write(str(BCz)+'\n')
ofFile.write(str(T)+'\n')
ofFile.write(str(Ntherm)+'\n')
ofFile.write(str(Nbins)+'\n')
ofFile.write(str(Nsweeps)+'\n')
ofFile.write(str(Jx)+'\n')
ofFile.write(str(Jy)+'\n')
ofFile.write(str(Delta)+'\n')
ofFile.write(str(V)+'\n')
ofFile.write(str(epsilon)+'\n')
ofFile.write(str(h)+'\n')
ofFile.close()

cmd_line_args = str(Lx) + " " + str(Jx) + " " + str(V) + " " + str(V) + " " + str(number_of_realizations)
os.system("python disorder_average.py " + cmd_line_args )
