The repository provides a minimal implementation for the stochastic series
expansion algorithm -- a highly efficient Quantum Monte Carlo algorithm for
the simulation of the Spin-1/2 XXZ model or hard-core boson model in
thermal equilibrium. The detailed description of the algorithm can be found at
https://arxiv.org/abs/cond-mat/0202316. The focus of the implementation is
two-fold: first, I provide an efficient implementation without additional
coding overhead in order to allow beginners to easily modify the Quantum Monte
Carlo simulation for their purposes as fast as possible. Second, the
implementation adds support for diagonal and off-diagonal disorder, i.e., it
extents the algorithm explicitly discussed in the above publication to the case
of non-translational invariant Hamiltonians. Therefore, no symmetries between
possible vertex configurations exists. For all possible vertex configurations
consult Fig. 8 in https://arxiv.org/abs/cond-mat/0202316. To my knowledge this
is the only open source implementation that provides support for the general,
non-translational invariant case. In this general case the directed loop
equations have to be solved for each site in real space separately. Depending on
the random field and random hopping at the site, bounces in the directed loop
construction can be either avoided or the probability can be minimized by
bouncing on the vertex configuration with the highest weight.


General Usage
-------------

To run the Quantum Monte Carlo simulation specify the system parameters in the
python script run_SSE.py and execute the script. The specific disorder
distribution for the disorder average can be modified in the python script
disorder_average.py.

License
-------

Copyright © 2019  Tobias Pfeffer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available in the
file [LICENSE.txt](LICENSE.txt).