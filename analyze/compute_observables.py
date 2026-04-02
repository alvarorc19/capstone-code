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
    import_renormalised_magnetisation,
    import_renormalised_energy,
)

k = 1.380649 * 10**(-23)

# <m>
def compute_average_magnetisation(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    magnetisation = import_observable(directory ,"magnetisation")
    x_mag = magnetisation[:,0]
    y_mag = magnetisation[:,1]
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    magnetisation = np.sqrt((x_mag)**2 + (y_mag)**2) / N
    mag_obs = pe.Obs([magnetisation[start:]], ["ens"])
    mag_obs.gamma_method()
    return mag_obs

def compute_renormalised_magnetisation(directory:pathlib.Path, start:int = 0, b:int = 1) -> pe.Obs:
    if b == 1:
        magnetisation = import_observable(directory ,"magnetisation")
    else:
        magnetisation = import_renormalised_magnetisation(directory , b)

    x_mag = magnetisation[:,0]
    y_mag = magnetisation[:,1]
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = (length / b) ** dim
    magnetisation = np.sqrt(x_mag**2 + y_mag**2) / N
    mag_obs = pe.Obs([magnetisation[start:]], ["ens"])
    mag_obs.gamma_method()
    return mag_obs

# <E> or <U> or <H>
def compute_average_energy(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    energy_array = import_observable(directory, "energy")
    obs = pe.Obs([energy_array[start:]], ["ens"])
    obs.gamma_method()
    return obs

def compute_renormalised_energy(directory:pathlib.Path, start:int = 0, b:int = 1) -> pe.Obs:
    if b == 1:
        energy_array = import_observable(directory, "energy")
    else:
        energy_array = import_renormalised_energy(directory, b)

    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = (length/b) ** dim
    J = import_physical_parameter(directory, "J")
    obs = pe.Obs([energy_array[start:]], ["ens"])
    norm_obs = obs / (J * N)
    norm_obs.gamma_method()
    return norm_obs

def compute_normalised_energy(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    energy_array = import_observable(directory, "energy")
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    J = import_physical_parameter(directory, "J")
    obs = pe.Obs([energy_array[start:]], ["ens"])
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

    obs = pe.Obs([energy_array[start:]], ["ens"])
    obs2 = pe.Obs([energy_array[start:]**2], ["ens"])
    cvk_obs = (obs2 - obs ** 2) / (temp**2 * N)
    cvk_obs.gamma_method()

    return cvk_obs


def compute_specific_heat(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    energy_array = import_observable(directory, "energy") 
    temp = import_physical_parameter(directory, "temperature")
    obs = pe.Obs([energy_array[start:]], ["ens"])
    obs2 = pe.Obs([energy_array[start:]**2], ["ens"])
    cvk_obs = (obs2 - obs ** 2) / (temp**2)
    cvk_obs.gamma_method()

    return cvk_obs

def compute_cluster_size(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    n_array = import_observable(directory, "average_cluster_size")
    obs = pe.Obs([n_array[start:]], ["ens"])
    obs.gamma_method()
    return obs

def compute_cluster_size_per_spin(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    n_array = import_observable(directory, "average_cluster_size")
    obs = pe.Obs([n_array[start:]], ["ens"])
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    obs = obs / N
    obs.gamma_method()
    return obs

def compute_cluster_susceptibility(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    n_array = import_observable(directory, "average_cluster_size")
    obs = pe.Obs([n_array[start:]], ["ens"])
    temp = import_physical_parameter(directory, "temperature")
    x_obs = 1 / temp * obs
    x_obs.gamma_method()
    return x_obs

def compute_cluster_susceptibility_per_spin(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    n_array = import_observable(directory, "average_cluster_size")
    obs = pe.Obs([n_array[start:]], ["ens"])
    temp = import_physical_parameter(directory, "temperature")
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    x_obs = ((1 / temp) * obs) / N
    x_obs.gamma_method()
    return x_obs


def compute_susceptibility(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    magnetisation = import_observable(directory ,"magnetisation")
    x_magnetisation = np.array(magnetisation[:,0])
    y_magnetisation = np.array(magnetisation[:,1])
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    temp = import_physical_parameter(directory, "temperature")
    magnetisation_array = np.sqrt(x_magnetisation**2 + y_magnetisation**2) / N

    obs = pe.Obs([magnetisation_array[start:]], ["ens"])
    obs2 = pe.Obs([magnetisation_array[start:]**2], ["ens"])
    X_obs = ((obs2 - obs ** 2) * N**2) / temp
    X_obs.gamma_method()
    return X_obs

# kTX = <M^2> - <M>^2 This is definition with big M, I calculate little m
# so its different. According to Newman & Barkema, this can be
# X = (N/kT) (<m^2> - <m>^2), susceptibility per spin
# We shall calculate kX

def compute_susceptibility_per_spin(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    magnetisation = import_observable(directory ,"magnetisation")
    x_magnetisation = np.array(magnetisation[:,0])
    y_magnetisation = np.array(magnetisation[:,1])
    length = import_physical_parameter(directory, "L")
    dim = import_physical_parameter(directory, "dimension")
    N = length ** dim
    temp = import_physical_parameter(directory, "temperature")
    magnetisation_array = np.sqrt(x_magnetisation**2 + y_magnetisation**2) / N

    obs = pe.Obs([magnetisation_array[start:]], ["ens"])
    obs2 = pe.Obs([magnetisation_array[start:]**2], ["ens"])
    X_obs = ((obs2 - obs**2) * N) / temp
    X_obs.gamma_method()
    return X_obs

def compute_binder_cumulant(directory:pathlib.Path, start:int = 0) -> pe.Obs:
    order_parameter_array = import_observable(directory, "magnetisation")
    obs2 = pe.Obs([order_parameter_array[start:]**2], ["ens"])
    obs4 = pe.Obs([order_parameter_array[start:]**4], ["ens"])
    binder_obs = 1 - (obs4 / (3 * obs2**2))
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


