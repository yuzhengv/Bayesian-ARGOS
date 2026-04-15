# %%
import os
import sys
import numpy as np
import pandas as pd
import math
from scipy import signal
from scipy.integrate import odeint

import matplotlib.pyplot as plt

# workspace_path = os.path.dirname(os.getcwd())
# local_package_path = os.path.join(workspace_path, "DataGeneration")
# sys.path.append(local_package_path)

import ode_auto

# %%
# ! -------------------------- Compare different solver ------------------------
# ~ -------------------------- solve odes with solve_ivp -----------------------
# ^ RK45
# %%
num_init = 100
num_state_var = 2

lotka_volterra_coeff = [[1, -1], [-1, 1]]
lotka_volterra_names = [["x1", "x1x2"], ["x2", "x1x2"]]


n = 5000
dt = 0.01
snr = 49
t = np.arange(0, float(n) * dt, dt)
out = ode_auto.solve_ode_ivp(
    lotka_volterra_coeff, lotka_volterra_names, [1, 2], t, method="RK45"
)

plt.plot(out[:, 0], out[:, 1])

# ^ LSODA
# %%
out_LSODA = ode_auto.solve_ode_ivp(
    lotka_volterra_coeff, lotka_volterra_names, [1, 2], t, method="LSODA"
)
plt.plot(out_LSODA[:, 0], out_LSODA[:, 1])

# ~ -------------------------- solve odes with ode_int __-----------------------
# %%
out_odeint = ode_auto.solve_ode_odeint(
    lotka_volterra_coeff, lotka_volterra_names, [1, 2], t
)
plt.plot(out_odeint[:, 0], out_odeint[:, 1])

# ^ ------------------------ Generate the Rossler system -----------------------
# %%
rossler_coeff = [[-1, -1], [1, 0.2], [0.2, 1, -5.7]]
rossler_names = [["x2", "x3"], ["x1", "x2"], ["", "x1x3", "x3"]]

n = 5000
dt = 0.01
snr = 49
t = np.arange(0, float(n) * dt, dt)

out_odeint = ode_auto.solve_ode_odeint(rossler_coeff, rossler_names, [1, 2, 3], t)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot(out_odeint[:, 0], out_odeint[:, 1], out_odeint[:, 2], lw=0.5)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("rossler")

# %%
# ^ ------------------------ Generate the vdp system -----------------------
vdp_coeff = [[1], [1.2, -1.2, -1]]
vdp_names = [["x2"], ["x2", "x2x1^2", "x1"]]

n = 5000
dt = 0.01
snr = 49
t = np.arange(0, float(n) * dt, dt)
out_odeint = ode_auto.solve_ode_odeint(vdp_coeff, vdp_names, [1, 2], t)

plt.plot(out[:, 0], out[:, 1])

# ^ ------------------------ Generate the Rossler system -----------------------
# %%
halvorsen_coeff = [[-1.89, -4, -4, -1], [-1.89, -4, -4, -1], [-1.89, -4, -4, -1]]
halvorsen_names = [
    ["x1", "x2", "x3", "x2^2"],
    ["x2", "x3", "x1", "x3^2"],
    ["x3", "x1", "x2", "x1^2"],
]

n = 5000
dt = 0.01
snr = 49
t = np.arange(0, float(n) * dt, dt)

out_odeint = ode_auto.solve_ode_odeint(halvorsen_coeff, halvorsen_names, [1, 2, 3], t)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot(out_odeint[:, 0], out_odeint[:, 1], out_odeint[:, 2], lw=0.5)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("halvorsen")


# %%
# ^ ------------------------ Generate the dadras system -----------------------
# %%
dadras_coeff = [[1, -3, 2.7], [1.7, -1, 1], [2, -9]]
dadras_names = [["x2", "x1", "x2x3"], ["x2", "x1x3", "x3"], ["x1x2", "x3"]]


n = 5000
dt = 0.01
snr = 49
t = np.arange(0, float(n) * dt, dt)
x0 = [1, 2, 3]
x0 = [2.675467, -0.09765433, -1.56746]
out_odeint = ode_auto.solve_ode_odeint(dadras_coeff, dadras_names, x0, t)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot(out_odeint[:, 0], out_odeint[:, 1], out_odeint[:, 2], lw=0.5)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("dadras " + str(x0))


# %%
