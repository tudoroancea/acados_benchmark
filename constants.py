import numpy as np
from casadi import MX

sym_t = MX
constant_controls = False
nx = 6
nu = 2
Nf = 50  # number of discretization steps
dt = 0.02  # sampling time
Tsim = 10.0  # maximum simulation time[s]
Nsim = int(Tsim / dt) + 1
sref_N = 3.0  # reference for final reference progress
ny = nx + nu
ny_e = nx

nh = 4
nsbx = 2
nsh = nh
ns = nsh + nsbx

## Race car parameters
m = 0.043
l_R = 1 / 31
l_F = 1 / 31
L = 0.1
W = 0.05
C1 = l_R / (l_R + l_F)
C2 = 1 / (l_R + l_F)
Cm1 = 0.28
Cm2 = 0.05
Cr0 = 0.011
Cr2 = 0.006
Cr3 = 5.0
tau_T = 1e-3
tau_delta = 0.02

# model bounds
v_min = 0.0
v_max = 3.0
alpha_min = -np.pi / 2
alpha_max = np.pi / 2
n_min = -0.12
n_max = 0.12
T_min = -1.0
T_max = 1.0
delta_min = -0.4
delta_max = 0.4
Tdot_min = -10.0
Tdot_max = 10.0
deltadot_min = -2.0
deltadot_max = 2.0
alat_min = -4.0
alat_max = 4.0
along_min = -4.0
along_max = 4.0


# initial state
x0 = np.array([-1.7, 0, 0, 0, 0, 0])

# costs
# Q = np.diag([1e-1, 1e-8, 1e-8, 1e-8, 1e-3, 5e-3])
# R = np.diag([1e-3, 5e-3])
# Qe = np.diag([5e0, 1e-8, 1e-8, 1e-8, 5e-3, 2e-3])
# x = (s, n, alpha, v, T, delta), u = (Tdot, deltadot)
Q = np.diag([1e-1, 5e-7, 5e-7, 5e-7, 5e-2, 2.5e-2])
R = np.diag([5e-1, 1e-1])
Qe = np.diag([1e-1, 2e-10, 2e-10, 2e-10, 1e-4, 4e-5])
zl = 100 * np.ones((ns,))
zu = 100 * np.ones((ns,))
Zl = 100 * np.ones((ns,))
Zu = 100 * np.ones((ns,))
