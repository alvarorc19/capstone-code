import os
import sys
import pathlib
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
def compute_average_magnetisation(directory:pathlib.Path) -> float:
    magnetisation_array = import_observable(directory, "magnetisation")
    return np.average(magnetisation_array)

# <E> or <U> or <H>
def compute_average_energy(directory:pathlib.Path) -> float:
    magnetisation_array = import_observable(directory, "energy")
    return np.average(magnetisation_array)

# kT^2C_v = <H^2> - <H>^2 
def compute_specific_heat(directory:pathlib.Path) -> float:
    energy_array = import_observable(directory, "energy") 
    kT2c_v = np.average(energy_array**2) - (np.average(energy_array))**2
    temp = import_physical_parameter(directory, "temperature")
    return kT2c_v / temp**2

# kTX = <M^2> - <M>^2 
def compute_susceptibility(directory:pathlib.Path) -> float:
    magnetisation_array = import_observable(directory, "magnetisation") 
    kTX = np.average(magnetisation_array**2) - (np.average(magnetisation_array))**2
    temp = import_physical_parameter(directory, "temperature")
    return kTX / temp

def compute_binder_cumulant(directory:pathlib.Path) -> float:
    order_parameter_array = import_observable(directory, "magnetisation")
    forth_avg = np.average(order_parameter_array**4)
    square_avg = np.average(order_parameter_array**2)
    binder_cumulant = 1 - (forth_avg / (2 * square_avg**2))
    return binder_cumulant

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
        specific_heat_array = np.append(specific_heat_array, compute_specific_heat(direc))

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
        magnetisation_array = np.append(magnetisation_array, compute_average_magnetisation(direc))

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
        susceptibility_array = np.append(susceptibility_array, compute_susceptibility(direc))

    p_fitted = np.polynomial.Polynomial.fit(np.log(epsilon_array), np.log(susceptibility_array), deg=1)
    log_fit = p_fitted.convert().coef
    return -log_fit[1]
