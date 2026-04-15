import os
import sys

util_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(util_dir)

import ode_auto

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# from mpl_toolkits.mplot3d import Axes3D


def generate_initial_value_df(
    seed,
    num_init_samples,
    x1_range,
    x2_range,
    x3_range=None,
    x4_range=None,
    num_columns=3,
):
    if num_columns not in [2, 3, 4]:
        raise ValueError("num_columns must be 2, 3, or 4")
    if num_columns == 3 and x3_range is None:
        raise ValueError("x3_range must be provided when num_columns is 3")
    if num_columns == 4 and (x3_range is None or x4_range is None):
        raise ValueError("x3_range and x4_range must be provided when num_columns is 4")

    np.random.seed(seed)

    ranges = [x1_range, x2_range]
    if num_columns == 3:
        ranges.append(x3_range)
    if num_columns == 4:
        ranges.extend([x3_range, x4_range])
    initial_values = [np.random.uniform(r[0], r[1], num_init_samples) for r in ranges]

    initial_value_df = pd.DataFrame(
        np.vstack(initial_values).T, columns=[f"x{i+1}" for i in range(num_columns)]
    )

    return initial_value_df


def plot_systems_trajectories(
    seed,
    num_init_samples,
    x1_range,
    x2_range,
    x3_range=None,
    x4_range=None,
    num_columns=3,
    plot_2d=None,
    system_names=None,
    system_coeff=None,
    n_obs=5000,
    dt=0.01,
    #   snr = 49,
    num_plot=100,
):
    initial_value_df = generate_initial_value_df(
        seed, num_init_samples, x1_range, x2_range, x3_range, x4_range, num_columns
    )
    t = np.arange(0, float(n_obs) * dt, dt)

    if num_columns == 3 and plot_2d == None:
        fig = plt.figure(figsize=(40, 40))
        for i in range(num_plot):
            init_condition = initial_value_df.iloc[i, :]
            traj_data = ode_auto.solve_ode_odeint(
                system_coeff, system_names, init_condition, t
            )
            ax = fig.add_subplot(10, 10, i + 1, projection="3d")
            ax.plot3D(
                traj_data[:, 0], traj_data[:, 1], traj_data[:, 2], color="#481b51"
            )
            ax.scatter(
                traj_data[n_obs - 1, 0],
                traj_data[n_obs - 1, 1],
                traj_data[n_obs - 1, 2],
                color="red",
            )
            ax.text2D(
                0.05,
                0.95,
                np.round(initial_value_df.iloc[i, :].values.tolist(), 4),
                transform=ax.transAxes,
            )

    elif num_columns == 3 and plot_2d is True:
        fig = plt.figure(figsize=(40, 40))
        for i in range(num_plot):
            init_condition = initial_value_df.iloc[i, :]
            traj_data = ode_auto.solve_ode_odeint(
                system_coeff, system_names, init_condition, t
            )
            traj_data_with_t = np.column_stack((traj_data, t))
            ax = fig.add_subplot(10, 10, i + 1)
            ax.plot(t, traj_data_with_t[:, 0], label="x1", color="#ca3e3d")
            ax.plot(t, traj_data_with_t[:, 1], label="x2", color="#0a73b0")
            ax.plot(t, traj_data_with_t[:, 2], label="x3", color="#f0c05a")
            ax.scatter(t[-1], traj_data_with_t[-1, 0], color="#ca3e3d", label="x1 end")
            ax.scatter(t[-1], traj_data_with_t[-1, 1], color="#0a73b0", label="x2 end")
            ax.scatter(t[-1], traj_data_with_t[-1, 2], color="#f0c05a", label="x3 end")
            ax.text(
                0.05,
                0.95,
                np.round(initial_value_df.iloc[i, :].values.tolist(), 4),
                transform=ax.transAxes,
            )

    elif num_columns == 4:
        fig = plt.figure(figsize=(40, 40))
        for i in range(num_plot):
            init_condition = initial_value_df.iloc[i, :]
            traj_data = ode_auto.solve_ode_odeint(
                system_coeff, system_names, init_condition, t
            )
            traj_data_with_t = np.column_stack((traj_data, t))
            ax = fig.add_subplot(10, 10, i + 1)
            ax.plot(t, traj_data_with_t[:, 0], label="x1", color="#ca3e3d")
            ax.plot(t, traj_data_with_t[:, 1], label="x2", color="#0a73b0")
            ax.plot(t, traj_data_with_t[:, 2], label="x3", color="#f0c05a")
            ax.plot(t, traj_data_with_t[:, 3], label="x4", color="#1d7139")
            ax.scatter(t[-1], traj_data_with_t[-1, 0], color="#ca3e3d", label="x1 end")
            ax.scatter(t[-1], traj_data_with_t[-1, 1], color="#0a73b0", label="x2 end")
            ax.scatter(t[-1], traj_data_with_t[-1, 2], color="#f0c05a", label="x3 end")
            ax.scatter(t[-1], traj_data_with_t[-1, 3], color="#1d7139", label="x4 end")
            ax.text(
                0.05,
                0.95,
                np.round(initial_value_df.iloc[i, :].values.tolist(), 4),
                transform=ax.transAxes,
            )
    else:
        raise ValueError("Currently only support 3D and 4D systems")
    plt.show()
