from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
from icecream import ic
from tracks.readDataFcn import load_track
from utils import *

from sing_free_mpc import gen_solver

track_file = "LMS_Track.txt"
[s_ref, X_ref, Y_ref, psi_ref, kappa_ref] = load_track(track_file)

# load model
t = perf_counter()
solver, alat_fun, along_fun = gen_solver("sing_free_mpc", s_ref, kappa_ref)
elapsed = perf_counter() - t
print(f"Generating solver took {elapsed:.3f} s")

# initialize data structs
simX = np.ndarray((Nsim, nx))
simU = np.ndarray((Nsim, nu))
runtimes = []

# simulate
start_sim = perf_counter()
for i in range(Nsim):
    # update reference
    solver.set(0, "lbx", x0)
    solver.set(0, "ubx", x0)
    s0 = x0[0]
    for j in range(Nf):
        # yref = np.array([0, 0, 0, 0, 0, 0, 0, 0])
        yref = np.array([s0 + sref_N * j / Nf, 0, 0, 0, 0, 0, 0, 0])
        solver.set(j, "yref", yref)
    yref_N = np.array([s0 + sref_N, 0, 0, 0, 0, 0])
    solver.set(Nf, "yref", yref_N)

    # solve ocp
    t = perf_counter()
    status = solver.solve()
    if status != 0:
        print(
            "acados returned status {} in closed loop iteration {}.".format(status, i)
        )
    elapsed = perf_counter() - t
    runtimes.append(1000*elapsed)

    # get solution
    x0 = solver.get(0, "x")
    u0 = solver.get(0, "u")
    simX[i] = x0
    simU[i] = u0

    # update initial condition
    x0 = solver.get(1, "x")

    # check if one lap is done and break and remove entries beyond
    if x0[0] > s_ref[-1] + 1.0:
        # find where vehicle first crosses start line
        simX = simX[:i, :]
        simU = simU[:i, :]
        Nsim = i
        break
stop_sim = perf_counter()
print(f"Simulation took {stop_sim - start_sim:.3f} s\n")

# find where there are big jumps in s
idx = np.argwhere(np.abs(np.diff(np.mod(simX[:, 0], s_ref[-1]))) > 0.1).ravel()
start_lap = idx[0]
end_lap = idx[1]
# Nsim = idx[1] - idx[0]
# simX = simX[idx[0] : idx[1], :]
# simU = simU[idx[0] : idx[1], :]

# Print some stats
runtimes = runtimes[10:]
print(f"Lap time: {dt*(end_lap - start_lap):.3f} s")
print(f"Average computation time: {np.mean(runtimes):.3f} ms")
print(f"Maximum computation time: {np.max(runtimes):.3f} ms")
print(f"Average speed: {np.mean(simX[start_lap-1:end_lap, 3]):.3f} m/s")

# Plot Results
t = np.linspace(0.0, Nsim * dt, Nsim)
plot_cl(simX, simU, t)
plot_traj(s_ref, X_ref, Y_ref, psi_ref, simX)
plot_accels(t, simX, simU, alat_fun, along_fun)
plt.figure()
plt.hist(runtimes, bins=50)
plt.title("Computation time histogram")
plt.xlabel("Computation time [ms]")
plt.ylabel("Number of occurences")
plt.tight_layout()
plt.savefig('hist.png', dpi=300)
plt.show()
