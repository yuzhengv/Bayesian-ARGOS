# Resources related to the pysindy package
# https://github.com/dynamicslab/pysindy/issues/551
# %%
import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))

src_path = os.path.abspath(os.path.join(current_path, "..", "src"))
sys.path.insert(0, src_path)


# %%
import ode_auto as ode_auto
import warnings
from contextlib import contextmanager
from copy import copy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import LinAlgWarning
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import Lasso

import pysindy as ps
from pysindy.utils import enzyme
from pysindy.utils import lorenz
from pysindy.utils import lorenz_control


# ! -------------------- Test the basic pysindy package --------------------
# %%

integrator_keywords = {}
integrator_keywords["rtol"] = 1e-12
integrator_keywords["method"] = "LSODA"
integrator_keywords["atol"] = 1e-12
# %%
# basic usage
t_end_train = 0.04
t_end_test = 0.04

# %%
# Train the model
# Generate measurement data
dt = 0.002

t_train = np.arange(0, t_end_train, dt)
x0_train = [-8, 8, 27]
t_train_span = (t_train[0], t_train[-1])
x_train = solve_ivp(
    lorenz, t_train_span, x0_train, t_eval=t_train, **integrator_keywords
).y.T
# %%
# Instantiate and fit the SINDy model
model = ps.SINDy()
model.fit(x_train, t=dt)
model.print()
model.get_feature_names()

# ! -------- Use pysindy package and ode_auto package together ----------------
# - Test the lorenz system
# %%
lorenz_coeff = [[-10, 10], [28, -1, -1], [1, -8 / 3]]
lorenz_names = [["x1", "x2"], ["x1", "x2", "x1x3"], ["x1x2", "x3"]]

n = 5000
dt = 0.001
# snr = 49
t = np.arange(0, float(n) * dt, dt)
init_condition = [2, 2, 20]

out_odeint = ode_auto.solve_ode_odeint(lorenz_coeff, lorenz_names, init_condition, t)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot(out_odeint[:, 0], out_odeint[:, 1], out_odeint[:, 2], lw=0.5)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("lorenz")
# %%
model = ps.SINDy()
model.fit(out_odeint, t=dt)
model.print()
# %%
model.coefficients()[0]
# %%
model.get_feature_names()
# %%
model.coefficients()[1]

# %%
# Define the Thomas system with Fourier terms
def thomas(t, x):
    return [
        np.sin(x[1]) - 0.208186 * x[0],
        np.sin(x[2]) - 0.208186 * x[1],
        np.sin(x[0]) - 0.208186 * x[2],
    ]


# Generate measurement data for the Thomas system
t_end_train = 100
dt = 0.01
t_train = np.arange(0, t_end_train, dt)
x0_train = [1, 1, 1]
t_train_span = (t_train[0], t_train[-1])
x_train = solve_ivp(
    thomas, t_train_span, x0_train, t_eval=t_train, **integrator_keywords
).y.T

# Instantiate and fit the SINDy model with Fourier library
fourier_library = ps.FourierLibrary(n_frequencies=6)
model = ps.SINDy(feature_library=fourier_library)
model.fit(x_train, t=dt)
model.print()
# %%
type(model.print())


# Test the function to obtain the results produced by the pysindy algorithm
# %%
def perform_STLS_3D(
    threshold_sequence,
    poly_order,
    t_span,
    max_iter_num,
    training_data,
    validation_init_condition_matrix_new,
    validation_data,
):
    """Perform SINDy with AIC for 3D systems"""
    # generate data
    training_data = training_data
    validation_data = validation_data
    # poly_order
    polyorder = poly_order
    # time span
    t_total = t_span
    # load training data
    x_dot_smoothed, x_smoothed, dt = training_data
    # load validation data
    (
        nterms_sindy_AIC,
        validation_theta_list,
        x_validation_list,
        x_dot_validation_smoothed_list,
    ) = validation_data
    # perform STLS on original data
    models = []
    models_coeff = []
    for j in np.arange(0, len(threshold_sequence)):
        model = ps.SINDy(
            optimizer=ps.STLSQ(
                threshold=threshold_sequence[j], alpha=0, max_iter=int(max_iter_num)
            ),
            feature_library=ps.PolynomialLibrary(degree=poly_order),
        )
        model.fit(
            x=x_smoothed,
            t=dt,
            x_dot=x_dot_smoothed,
            quiet=True,
        )
        models.append(model)
        models_coeff.append(model.coefficients())
    # put models in dataframe
    xdot_eta0_100_pe_df = pd.DataFrame()
    ydot_eta0_100_pe_df = pd.DataFrame()
    zdot_eta0_100_pe_df = pd.DataFrame()
    for i in range(0, len(models)):
        xdot_eta0_100_pe_df[i,] = models[i].coefficients()[0]
        ydot_eta0_100_pe_df[i,] = models[i].coefficients()[1]
        zdot_eta0_100_pe_df[i,] = models[i].coefficients()[2]
    # Create combinations of STLS
    sindy_combinations_list = []
    for i in np.arange(0, xdot_eta0_100_pe_df.shape[1]):
        x_pe_df = xdot_eta0_100_pe_df.iloc[:, i]
        for j in np.arange(0, ydot_eta0_100_pe_df.shape[1]):
            y_pe_df = ydot_eta0_100_pe_df.iloc[:, j]
            for k in np.arange(0, zdot_eta0_100_pe_df.shape[1]):
                z_pe_df = zdot_eta0_100_pe_df.iloc[:, k]
                sindy_combinations = np.vstack([x_pe_df, y_pe_df, z_pe_df]).T
                sindy_combinations_list.append(sindy_combinations)
    sindy_combinations_list = np.unique(sindy_combinations_list, axis=0)
    sindy_combinations_list_new = [
        element for element in sindy_combinations_list if np.count_nonzero(element) != 0
    ]
    aicc_min_list = []
    for i in np.arange(0, len(sindy_combinations_list_new)):
        current_sindy_model = sindy_combinations_list_new[i]
        x_y_sindy_abs_error = np.zeros(len(validation_theta_list))
        for j in np.arange(0, len(validation_theta_list)):
            x_hat_sindy = np.matmul(validation_theta_list[j], current_sindy_model[:, 0])
            y_hat_sindy = np.matmul(validation_theta_list[j], current_sindy_model[:, 1])
            z_hat_sindy = np.matmul(validation_theta_list[j], current_sindy_model[:, 2])
            x_int_sindy = (
                integrate.cumtrapz(y=x_hat_sindy, dx=dt, initial=0)
                + validation_init_condition_matrix_new[j][0]
            )
            y_int_sindy = (
                integrate.cumtrapz(y=y_hat_sindy, dx=dt, initial=0)
                + validation_init_condition_matrix_new[j][1]
            )
            z_int_sindy = (
                integrate.cumtrapz(y=z_hat_sindy, dx=dt, initial=0)
                + validation_init_condition_matrix_new[j][2]
            )
            length_data = len(validation_theta_list[j])
            nterms_sindy_abserror = current_sindy_model.shape[1]
            x_y_sindy_abs_error[j] = np.sum(
                np.sum(
                    (
                        np.absolute(x_validation_list[j][:, 0] - x_int_sindy)
                        + np.absolute(x_validation_list[j][:, 1] - y_int_sindy)
                        + np.absolute(x_validation_list[j][:, 2] - z_int_sindy)
                    )
                )
                / length_data
                / nterms_sindy_abserror
            )
        log_L = (
            (-1 * nterms_sindy_AIC)
            * np.log(
                np.sum(x_y_sindy_abs_error * x_y_sindy_abs_error) / nterms_sindy_AIC
            )
            / 2
        )
        num_nonzero_coeff = np.count_nonzero(current_sindy_model)
        if (num_nonzero_coeff > 0) & ((length_data - num_nonzero_coeff - 1) != 0):
            aic = -2 * log_L + 2 * num_nonzero_coeff
            aicc_min = aic + 2 * num_nonzero_coeff * (num_nonzero_coeff + 1) / (
                length_data - num_nonzero_coeff - 1
            )
        else:
            aicc_min = -2 * log_L
        aicc_min_list.append(aicc_min)
    try:
        sindy_final_model = sindy_combinations_list_new[np.argmin(aicc_min_list)]
    except:
        sindy_final_model = np.zeros(shape=(len(x_pe_df), 3))
    return sindy_final_model
