import subprocess, sys, random
import numpy as np

counter = 0
for arg in sys.argv:
  if counter == 1:
    L = int(arg)
  elif counter == 2:
    J = float(arg)
  elif counter == 3:
    V = float(arg)
  elif counter == 4:
    h = float(arg)
  elif counter == 5:
    number_of_realizations = int(arg)
  counter += 1

for i in range(number_of_realizations):
  random_hopping = np.zeros((L*L,2))
  random_fields = np.zeros((L*L))

  for x in range(L):
    for y in range(L):
      if h != 0:
        random_hopping[x + L*y][0] = J * np.random.power(1./h)
        random_hopping[x + L*y][1] = J * np.random.power(1./h)
      else:
        random_hopping[x + L*y][0] = J
        random_hopping[x + L*y][1] = J
      random_fields[x + L*y] = np.random.uniform(-h,h)

  np.savetxt('./qc_obs/t_hop_b', random_hopping)
  np.savetxt('./qc_obs/h_sites', random_fields)

  p = subprocess.call('./SSE', shell=True)
