import os



# the input of parameters starts here

# this is the simulation parameters

params = {}
params['Lx'] = 16 				# system size in x direction
params['Ly'] = 16	 				# system size in y direction
params['Lz'] = 1 					# system size in z direction
params['BCx'] = 1 			   # 1: PBC, 0: PBC in x direction
params['BCy'] = 1 				# 1: PBC, 0: PBC in y direction
params['BCz'] = 0 				# 1: PBC, 0: PBC in z direction
params['T'] = 0.1             # Temperature of the system
params['Ntherm'] = 10000      # number of sweeps for thermalization
params['Nbins'] = 100         # number of bins to collect
params['Nsweeps'] = 1000 	   # number of sweeps per bin
params['Jx'] = 1.             # hopping in x-direction
params['Jy'] = 1.             # hopping in y-direction
params['Delta'] = 0. 			# anisotropy
params['V'] = 1.					# off_diagonal disorder
params['epsilon'] = 0.5			# to ensure positive weights
params['h'] = 1. 					# diagonal disorder

run_time = 96*3600
number_of_realizations = 1  # number of disorder realizations for the disorder average

# make a executable out of the src files

os.system('g++ -std=c++11 -O3 -o sse_XXZ.o -c sse_XXZ.cpp')
os.system('g++ -std=c++11 -O3 -o main.o -c main.cpp')
os.system('g++ -std=c++11 -o SSE main.o sse_XXZ.o')

if os.path.exists("qc_obs") == False:
  os.mkdir('qc_obs')

# write a new parameter file
ofFile = open('params', 'w')
for key in params:
  ofFile.write(str(params[key])+'\n')
ofFile.close()

cmd_line_args = str(params['Lx']) + " " + str(params['Jx']) + " " + str(params['V']) + " " + str(params['h']) + " " + str(number_of_realizations)
os.system("python disorder_average.py " + cmd_line_args )
