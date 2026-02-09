import os
import sys
import pathlib
import pyerrors as pe
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import toml
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from utils.h5_utils import(
    import_observable,
    import_physical_parameter,
)

# <m>
def compute_average_magnetisation(directory:pathlib.Path) -> pe.Obs:
    magnetisation_array = import_observable(directory, "magnetisation")
    obs = pe.Obs([magnetisation_array], ["magnetisation"])
    obs.gamma_method()
    return obs

# <E> or <U> or <H>
def compute_average_energy(directory:pathlib.Path) -> pe.Obs:
    energy_array = import_observable(directory, "energy")
    obs = pe.Obs([energy_array], ["energy"])
    obs.gamma_method()
    return obs

def compute_normalised_energy(directory:pathlib.Path) -> pe.Obs:
    energy_array = import_observable(directory, "energy")
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    J = import_physical_parameter(directory, "J")
    obs = pe.Obs([energy_array], ["energy"])
    norm_obs = obs / (J * N)

    norm_obs.gamma_method()
    return norm_obs


# kT^2C_v = <H^2> - <H>^2 
def compute_specific_heat(directory:pathlib.Path) -> pe.Obs:
    energy_array = import_observable(directory, "energy") 
    temp = import_physical_parameter(directory, "temperature")
    obs = pe.Obs([energy_array], ["energy"])
    obs2 = pe.Obs([energy_array**2], ["energy2"])
    cvk_obs = (obs2 - obs ** 2) / (temp**2)
    cvk_obs.gamma_method()

    return cvk_obs

# kTX = <M^2> - <M>^2 
def compute_susceptibility(directory:pathlib.Path) -> pe.Obs:
    magnetisation_array = import_observable(directory, "magnetisation") 
    temp = import_physical_parameter(directory, "temperature")
    obs = pe.Obs([magnetisation_array], ["magnetisation"])
    obs2 = pe.Obs([magnetisation_array**2], ["magnetisation2"])
    X_obs = (obs2 - obs ** 2) / temp
    X_obs.gamma_method()
    return X_obs

def compute_binder_cumulant(directory:pathlib.Path) -> pe.Obs:
    order_parameter_array = import_observable(directory, "magnetisation")
    obs2 = pe.Obs([order_parameter_array**2], ["magnetisation2"])
    obs4 = pe.Obs([order_parameter_array**4], ["magnetisation4"])
    binder_obs = 1 - (obs4 / (2 * obs2**2))
    binder_obs.gamma_method()
    return binder_obs

# C_H \sim |\epsilon|^{-\alpha}
def compute_alpha_critical_exponent(directory:pathlib.Path, critical_temp:float) -> float:
    params = [x for x in directory.iterdir() if x.is_dir()]
    epsilon_array = np.array([])
    specific_heat_array = np.array([])

    for direc in params:
        temperature = import_physical_parameter(direc, "temperature")
        epsilon = np.abs(temperature - critical_temp)/critical_temp
        if temperature == critical_temp:
            continue
        epsilon_array = np.append(epsilon_array, epsilon)
        specific_heat_array = np.append(specific_heat_array, compute_specific_heat(direc).value)

    p_fitted = np.polynomial.Polynomial.fit(np.log(epsilon_array), np.log(specific_heat_array), deg=1)
    log_fit = p_fitted.convert().coef
    return -log_fit[1]

# M \sim |\epsilon|^{\beta}, only for T < T_c
def compute_beta_critical_exponent(directory:pathlib.Path, critical_temp:float) ->  float:
    params = [x for x in directory.iterdir() if x.is_dir()]
    epsilon_array = np.array([])
    magnetisation_array = np.array([])

    for direc in params:
        temperature = import_physical_parameter(direc, "temperature")
        if temperature < critical_temp:
            epsilon = np.abs(temperature - critical_temp)/critical_temp
        elif temperature == critical_temp:
            continue
        else:
            epsilon = 1e12
        epsilon_array = np.append(epsilon_array, epsilon)
        magnetisation_array = np.append(magnetisation_array, compute_average_magnetisation(direc).value)

    p_fitted = np.polynomial.Polynomial.fit(np.log(epsilon_array), np.log(magnetisation_array), deg=1)
    log_fit = p_fitted.convert().coef
    return log_fit[1]

# X_T \sim |\epsilon| ^{-\gamma}
def compute_gamma_critical_exponent(directory:pathlib.Path, critical_temp:float) ->  float:
    params = [x for x in directory.iterdir() if x.is_dir()]
    epsilon_array = np.array([])
    susceptibility_array = np.array([])

    for direc in params:
        temperature = import_physical_parameter(direc, "temperature")
        if temperature == critical_temp:
            continue
        epsilon = np.abs(temperature - critical_temp)/critical_temp
        epsilon_array = np.append(epsilon_array, epsilon)
        susceptibility_array = np.append(susceptibility_array, compute_susceptibility(direc).value)

    p_fitted = np.polynomial.Polynomial.fit(np.log(epsilon_array), np.log(susceptibility_array), deg=1)
    log_fit = p_fitted.convert().coef
    return -log_fit[1]

