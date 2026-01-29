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
from utils.h5_utils import import_observable

# <m>
def compute_average_magnetisation(directory:pathlib.Path) -> float:
    magnetisation_array = import_observable(directory, "magnetisation")
    return np.average(magnetisation_array)

# <E> or <U> or <H>
def compute_average_energy(directory:pathlib.Path) -> float:
    magnetisation_array = import_observable(directory, "energy")
    return np.average(magnetisation_array)

# kT^2C_v = <H^2> - <H>^2 
def compute_specific_heat(directory:pathlib.Path, config:dict) -> float:
    energy_array = import_observable(directory, "energy") 
    kT2c_v = np.average(energy_array**2) - (np.average(energy_array))**2
    temp = config["physical_settings"]["temperature"]
    return kT2c_v / temp**2

# kTX = <M^2> - <M>^2 
def compute_susceptibility(directory:pathlib.Path, config:dict) -> float:
    magnetisation_array = import_observable(directory, "magnetisation") 
    kTX = np.average(magnetisation_array**2) - (np.average(magnetisation_array))**2
    temp = config["physical_settings"]["temperature"]
    return kTX / temp

def compute_binder_cumulant(directory:pathlib.Path) -> float:
    order_parameter_array = import_observable(directory, "magnetisation")
    forth_avg = np.average(order_parameter_array**4)
    square_avg = np.average(order_parameter_array**2)
    binder_cumulant = 1 - (forth_avg / (2 * square_avg**2))
    return binder_cumulant


