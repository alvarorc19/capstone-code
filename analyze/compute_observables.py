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
from utils.h5_utils import(
    import_observable,
    import_physical_parameter,
)

k = 1.380649 * 10**(-23)

# <m>
def compute_average_magnetisation(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    magnetisation = import_observable(directory ,"magnetisation")
    x_obs = pe.Obs([[i[0] for i in magnetisation[start:]]], ["x_magnetisation"])
    y_obs = pe.Obs([[i[1] for i in magnetisation[start:]]], ["y_magnetisation"])
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    mag_obs = np.sqrt(x_obs**2 + y_obs**2) / N
    mag_obs.gamma_method()
    return mag_obs

# <E> or <U> or <H>
def compute_average_energy(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    energy_array = import_observable(directory, "energy")
    obs = pe.Obs([energy_array[start:]], ["energy"])
    obs.gamma_method()
    return obs

def compute_normalised_energy(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    energy_array = import_observable(directory, "energy")
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    J = import_physical_parameter(directory, "J")
    obs = pe.Obs([energy_array[start:]], ["energy"])
    norm_obs = obs / (J * N)

    norm_obs.gamma_method()
    return norm_obs


def compute_reduced_temperature(temp_array:np.ndarray, t_c:float) -> np.ndarray:
    t = ( temp_array - t_c) / t_c
    return t

def compute_susceptibility_scaling_function(susc_obs:pe.Obs, length:int, gamma:float, nu:float) -> pe.Obs:
    reduced_susc_obs = length ** (- gamma / nu) * susc_obs
    reduced_susc_obs.gamma_method()
    return reduced_susc_obs


# kT^2C_v = <H^2> - <H>^2 from Landau and Binder, is total one, not per spin 
# so its different. According to Newman & Barkema, this can be
# c = (1/kT^2 N) (<E^2> - <E>^2), susceptibility per spin
# We shall calculate kc
def compute_specific_heat_per_spin(directory:pathlib.Path, start:int = 0) -> pe.Obs:

    energy_array = import_observable(directory, "energy") 
    temp = import_physical_parameter(directory, "temperature")
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim

    obs = pe.Obs([energy_array[start:]], ["energy"])
    obs2 = pe.Obs([energy_array[start:]**2], ["energy2"])
    cvk_obs = (obs2 - obs ** 2) / (temp**2 * N)
    cvk_obs.gamma_method()

    return cvk_obs


def compute_specific_heat(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    energy_array = import_observable(directory, "energy") 
    temp = import_physical_parameter(directory, "temperature")
    obs = pe.Obs([energy_array[start:]], ["energy"])
    obs2 = pe.Obs([energy_array[start:]**2], ["energy2"])
    cvk_obs = (obs2 - obs ** 2) / (temp**2)
    cvk_obs.gamma_method()

    return cvk_obs



def compute_susceptibility(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    magnetisation = import_observable(directory ,"magnetisation")
    x_magnetisation = np.array([i[0] for i in magnetisation[start:]])
    y_magnetisation = np.array([i[1] for i in magnetisation[start:]])
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    temp = import_physical_parameter(directory, "temperature")
    magnetisation_array = np.sqrt(x_magnetisation**2 + y_magnetisation**2) / N

    obs = pe.Obs([magnetisation_array], ["magnetisation"])
    obs2 = pe.Obs([magnetisation_array**2], ["magnetisation2"])
    X_obs = ((obs2 - obs ** 2) * N**2) / temp
    X_obs.gamma_method()
    return X_obs

# kTX = <M^2> - <M>^2 This is definition with big M, I calculate little m
# so its different. According to Newman & Barkema, this can be
# X = (N/kT) (<m^2> - <m>^2), susceptibility per spin
# We shall calculate kX

def compute_susceptibility_per_spin(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    magnetisation = import_observable(directory ,"magnetisation")
    x_magnetisation = np.array([i[0] for i in magnetisation[start:]])
    y_magnetisation = np.array([i[1] for i in magnetisation[start:]])
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    temp = import_physical_parameter(directory, "temperature")
    magnetisation_array = np.sqrt(x_magnetisation**2 + y_magnetisation**2) / N

    obs = pe.Obs([magnetisation_array], ["magnetisation"])
    obs2 = pe.Obs([magnetisation_array**2], ["magnetisation2"])
    X_obs = ((obs2 - obs ** 2) * N) / temp
    X_obs.gamma_method()
    return X_obs

def compute_binder_cumulant(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    order_parameter_array = import_observable(directory, "magnetisation")
    obs2 = pe.Obs([order_parameter_array[start:]**2], ["magnetisation2"])
    obs4 = pe.Obs([order_parameter_array[start:]**4], ["magnetisation4"])
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


