# %%
import sys
from pathlib import Path

src_path = Path.cwd().parent / "src"
sys.path.append(str(src_path))

import numpy as np
import ode_auto as ode_auto
import pysindy as ps

# %%
# : -------------------- Use lorenz system as a demo --------------------
# Lorenz system using pysindy class and ode_auto function
# Paramters settings for lorenz system
lorenz_coeff = [[-10, 10], [28, -1, -1], [1, -8 / 3]]
lorenz_names = [["x1", "x2"], ["x1", "x2", "x1x3"], ["x1x2", "x3"]]

n = 5000
dt = 0.001
# snr = 49
t = np.arange(0, float(n) * dt, dt)
init_condition = [2, 2, 20]

# %%
# Generate the training dataset
out_odeint = ode_auto.solve_ode_odeint(
    lorenz_coeff, lorenz_names, init_condition, t
)

# %%
model = ps.SINDy()
model.fit(out_odeint, t=dt)
model.print()


# %%
# $ -------------------- Function dedvelopment --------------------
# ! -------------------- Extract experiment output --------------------
# - Dealing with output of the model to get the coefficients and results
feature_names = model.get_feature_names()
coefficients = model.coefficients()

# Pair feature names with each row of coefficients
paired_features = np.array(
    [
        (feature_names, coefficients[i, :].tolist())
        for i in range(coefficients.shape[0])
    ]
)
paired_features.shape

# %%
# ! Pair feature names with each coefficient for each list
feature_names = model.get_feature_names()
coefficients = model.coefficients()

paired_features = [
    [(feature_names[j], coefficients[i, j]) for j in range(len(feature_names))]
    for i in range(coefficients.shape[0])
]

paired_features_cleaned = [
    [
        (feature_names[j], coefficients[i, j])
        for j in range(len(feature_names))
        if coefficients[i, j] != 0
    ]
    for i in range(coefficients.shape[0])
]


# %%
# ! -------------------- Post processing of the extracted results --------------------
# Extract the feature names of the three functions as tuples
feature_names_list = []
for i, feature_list in enumerate(paired_features_cleaned):
    feature_names_only = tuple(feature for feature, coeff in feature_list)
    feature_names_list.append(feature_names_only)

# %%
# ! -------------------- Function Integration --------------------


def post_processing_pysindy(model: ps.SINDy()):
    feature_names = model.get_feature_names()
    coefficients = model.coefficients()

    paired_features_cleaned = [
        [
            (feature_names[j], coefficients[i, j])
            for j in range(len(feature_names))
            if coefficients[i, j] != 0
        ]
        for i in range(coefficients.shape[0])
    ]

    return paired_features_cleaned


# %%
# : -------------------- Check the function including fourier terms --------------------
thomas_coeff = [[1, -0.208186], [1, -0.208186], [1, -0.208186]]
thomas_names = [["sin(x2)", "x1"], ["sin(x3)", "x2"], ["sin(x1)", "x3"]]

n = 5000
dt = 0.01
# snr = 49
t = np.arange(0, float(n) * dt, dt)
init_condition = [-0.5, 0.5, 0.5]

thomas_data = ode_auto.solve_ode_odeint(
    thomas_coeff, thomas_names, init_condition, t
)

# %%
fourier_library = ps.FourierLibrary(n_frequencies=6)
model = ps.SINDy(feature_library=fourier_library)
model.fit(thomas_data, t=dt)
model.print()
# %%
post_processing_pysindy(model)

# %%
# : -------------------- Check 2D function --------------------
lv_coeff = [[1, -1], [-1, 1]]
lv_names = [["x1", "x1x2"], ["x2", "x1x2"]]

n = 5000
dt = 0.01
# snr = 49
t = np.arange(0, float(n) * dt, dt)
init_condition = [3, 5]

lv_data = ode_auto.solve_ode_odeint(lv_coeff, lv_names, init_condition, t)

# %%
model = ps.SINDy()
model.fit(lv_data, t=dt)
model.print()
# %%
post_processing_pysindy(model)
# %%
np <- import("numpy")
smoothed_data_np <- np$array(as.matrix(smoothed_data))
smoothed_data_dot_np <- np$array(as.matrix(smoothed_data_dot))

# Fit the libraries first
polynomial_library$fit(smoothed_data_np)
fourier_library$fit(smoothed_data_np)

# Create combined library and fit it
combined_library <- ps$feature_library$GeneralizedLibrary(
    libraries = list(polynomial_library, fourier_library)
)
combined_library$fit(smoothed_data_np)