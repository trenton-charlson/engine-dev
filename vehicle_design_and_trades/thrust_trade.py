# Python file to loosely trade total thrust vs vehicle sizing

import numpy as np
from _1D_rocket_traj import _1D_rocket_traj

mwet = 81.5#kg
mdry = 69.5 #kg
mdot = 1.22 #kg/s
thrust = 3000.0 #Newton
D = 6.17 # inches - 6" tube
A_cs_i = np.pi*(D/2)**2
A_cs = A_cs_i * 0.00064516 # in**2 -> m**2

print(A_cs)

traj = _1D_rocket_traj(mwet, mdry, mdot, thrust, A_cs)
print(f'Max Altitude: {np.max(traj["x"])}')
