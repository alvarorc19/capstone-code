import os
import sys
import pathlib
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import toml
import numpy as np
import matplotlib.pyplot as plt
from compute_observables import (
        compute_average_magnetisation,
        compute_average_energy,
        compute_susceptibility,
        compute_specific_heat,
        compute_binder_cumulant,
    )

def magnetisation_plot(directory:pathlib.Path, axs:plt.axes, xaxis: str = "temperature") -> plt.axes:
    params = [x for x in directory.iterdir() if x.is_dir()]
    x_array = np.array([])
    magnetisation_array = np.array([])

    for direc in params:
        config = toml.load(direc / "config.toml")
        print("direc", direc)
        x_array = np.append(x_array, config["physical_settings"][xaxis])
        magnetisation_array = np.append(magnetisation_array, compute_average_magnetisation(direc))

    axs.scatter(x_array, magnetisation_array)

    return axs

def energy_plot(directory:pathlib.Path, axs:plt.axes, xaxis:str = "temperature") -> plt.axes:
    params = [x for x in directory.iterdir() if x.is_dir()]
    x_array = np.array([])
    energy_array = np.array([])

    i = 0
    for direc in params:
        config = toml.load(direc / "config.toml")
        print("direc", direc)
        x_array = np.append(x_array, config["physical_settings"][xaxis])
        energy_array = np.append(energy_array, compute_average_energy(direc))
        i+=1

    axs.scatter(x_array, energy_array)

    return axs

def susceptibility_plot(directory:pathlib.Path, axs:plt.axes, xaxis:str = "temperature") -> plt.axes:
    params = [x for x in directory.iterdir() if x.is_dir()]
    x_array = np.array([])
    susceptibility_array = np.array([])

    i = 0
    for direc in params:
        config = toml.load(direc / "config.toml")
        print("direc", direc)
        x_array = np.append(x_array, config["physical_settings"][xaxis])
        susceptibility_array = np.append(susceptibility_array, compute_susceptibility(direc, config))
        i+=1

    axs.scatter(x_array, susceptibility_array)

    return axs


def specific_heat_plot(directory:pathlib.Path, axs:plt.axes, xaxis:str = "temperature") -> plt.axes:
    params = [x for x in directory.iterdir() if x.is_dir()]
    x_array = np.array([])
    specific_heat_array = np.array([])

    i = 0
    for direc in params:
        config = toml.load(direc / "config.toml")
        print("direc", direc)
        x_array = np.append(x_array, config["physical_settings"][xaxis])
        specific_heat_array = np.append(specific_heat_array, compute_specific_heat(direc, config))
        i+=1

    axs.scatter(x_array, specific_heat_array)

    return axs

def binder_cumulant_plot(directory:pathlib.Path, ax:plt.axes, xaxis: str = "temperature") -> plt.axes:
    params = [x for x in directory.iterdir() if x.is_dir()]
    binder_cumulant_array = np.array([])
    x_array = np.array([])

    for direc in params:
        config = toml.load(direc / "config.toml")
        print("direc", direc)
        x_array = np.append(x_array, config["physical_settings"][xaxis])
        binder_cumulant_array = np.append(binder_cumulant_array, compute_binder_cumulant(direc))

    ax.scatter(x_array, binder_cumulant_array)

    return ax
