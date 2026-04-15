# %%
import sys
from pathlib import Path

src_path = Path.cwd().parent / "src"
sys.path.append(str(src_path))

import numpy as np

# import ode_auto as ode_auto
import pysindy as ps


# %%
def post_processing_pysindy(model : ps.SINDy()):
    feature_names = model.get_feature_names()
    coefficients = model.coefficients()
    
    paired_features_cleaned = [
        [(feature_names[j], coefficients[i, j]) for j in range(len(feature_names)) if coefficients[i, j] != 0]
        for i in range(coefficients.shape[0])
        ]

    return paired_features_cleaned