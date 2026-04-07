import os
import sys
import pathlib

from utils.h5_utils import import_observable
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import toml
import json
import pandas as pd
import pyerrors as pe
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import itertools
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
from matplotlib.widgets import Slider

from compute_observables import (
    compute_renormalised_energy,
    compute_renormalised_magnetisation,
)
from utils.h5_utils import (
    import_observable,
    import_physical_parameter,
)

from observables_plots import (
    get_observables_csv,
    _add_scatter_data,
    _add_format_plot,
)

def do_renormalisation_plot(directory:pathlib.Path, is_deep:bool = False, start:int = 0):
    saving_path = directory.parent.parent / "analyze" / "output"/"img_dump"
    # plt.tight_layout()
    csv_file = directory / "ensemble_observables.csv"
    if csv_file.exists():
        df = pd.read_csv(csv_file)
    else:
        get_observables_csv(directory, is_deep, start, False)
        df = pd.read_csv(csv_file)
    df = df[df["energy_value"] != 0.0].reset_index(drop=True)

    fig_mag, ax1 = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(11,6)
    )

    fig_energy, ax2 = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(11,6)
    )

    b_list = [2,4]
    cmap = plt.cm.tab20
    colors = cmap(np.arange(10))
    # markers = ['.-',':.','-.','.-',':.','-.','.-',':.','-.']
    for l, group in df.groupby("L"):
        group = group.sort_values("temperature")
        xaxis = group["temperature"].to_numpy()
        yaxis = group[f"energy_per_spin_value"].to_numpy()
        yerr = group[f"energy_per_spin_error"].to_numpy()

        ax1 = _add_scatter_data(
            axs = ax1,
            xaxis = xaxis,
            yaxis = yaxis,
            yerr = yerr,
            data_label = f"$b = 1, L = {l}$",
            main_color = colors[0],
            secondary_color = colors[1],
            marker = ".-",
        )

        xaxis = group["temperature"].to_numpy()
        yaxis = group[f"magnetisation_value"].to_numpy()
        yerr = group[f"magnetisation_error"].to_numpy()
        ax2 = _add_scatter_data(
            axs = ax2,
            xaxis = xaxis,
            yaxis = yaxis,
            yerr = yerr,
            data_label = f"$b = 1, L = {l}$",
            main_color = colors[0],
            secondary_color = colors[1],
            marker = ".-",
        )
    for i,b in enumerate(b_list):
        for l, group in df.groupby("L"):
            group = group.sort_values("temperature")
            xaxis = group["temperature"].to_numpy()
            yaxis = group[f"energy_b_{b}_value"].to_numpy()
            yerr = group[f"energy_b_{b}_error"].to_numpy()

            ax1 = _add_scatter_data(
                axs = ax1,
                xaxis = xaxis,
                yaxis = yaxis,
                yerr = yerr,
                data_label = f"$b = {b}, L = {l}$",
                main_color = colors[i+1],
                secondary_color = colors[i+1+1],
                marker = ".-",
            )

            xaxis = group["temperature"].to_numpy()
            yaxis = group[f"magnetisation_b_{b}_value"].to_numpy()
            yerr = group[f"magnetisation_b_{b}_error"].to_numpy()
            ax2 = _add_scatter_data(
                axs = ax2,
                xaxis = xaxis,
                yaxis = yaxis,
                yerr = yerr,
                data_label = f"$b = {b}, L = {l}$",
                main_color = colors[i+1],
                secondary_color = colors[i+2],
                marker = ".-",
            )

    
    ax1 = _add_format_plot(
        axs = ax1,
        xlabel="Temperature, $k_BT$",
        ylabel="Energy E",
        title = f"Energy vs Temperature renorm",
    )
    ax2 = _add_format_plot(
        axs = ax2,
        xlabel="Temperature, $k_BT$",
        ylabel="Magnetisation $m$",
        title = f"Magnetisation vs Temperature renorm",
    )

    print("rg save path ", saving_path)
    fig_mag.savefig(saving_path / f"{directory.name}_renormalised_magnetisation.pdf", bbox_inches = "tight")
    fig_energy.savefig(saving_path / f"{directory.name}_renormalised_energy.pdf", bbox_inches = "tight")

def do_biggest_L_renormalisation_plot(directory:pathlib.Path, is_deep:bool = False, start:int = 0):
    plt.tight_layout()

    saving_path = directory.parent.parent / "analyze" / "output"/f"renormalisation_{directory.name}"
    saving_path.mkdir(parents = True, exist_ok = True)
    print("save path : ", saving_path)
    if is_deep:
        sub_dir = [x for x in directory.iterdir() if x.is_dir()]
        params = []
        for dir in sub_dir:
            if dir.is_dir():
                for direc in dir.iterdir():
                    if direc.is_dir():
                        params.append(direc)
    else:
        params = [x for x in directory.iterdir() if x.is_dir()]

    fig_mag, ax1 = plt.subplots(
        ncols=1,
        nrows=1
    )

    fig_energy, ax2 = plt.subplots(
        ncols=1,
        nrows=1
    )

    b_list = [1,2,4,16]
    cmap = plt.cm.tab20
    colors = cmap(np.arange(10))
    markers = ['.-',':.','-.','.-',':.','-.','.-',':.','-.']
    for i,b in enumerate(b_list):
        temp_array = np.array([])
        length_array = np.array([])
        magnetisation_array = np.array([])
        magnetisation_error = np.array([])
        magnetisation_obs_array = np.array([])
        energy_array = np.array([])
        energy_error = np.array([])
        energy_obs_array = np.array([])

        with open(directory / "global_parameters.json", "r") as f:
            global_config = json.load(f)
        unique_lengths = global_config["physical_settings"]["L"]
        unique_temp = global_config["physical_settings"]["temperature"]

        start = 1000
        for direc in params:
            temp_array = np.append(temp_array, import_physical_parameter(direc, "temperature"))
            length_array = np.append(length_array, import_physical_parameter(direc, "L")) 
            magnetisation_obs = compute_renormalised_magnetisation(direc, start, b)
            magnetisation_obs_array = np.append(magnetisation_obs_array, magnetisation_obs)
            magnetisation_array = np.append(magnetisation_array, magnetisation_obs.value)
            magnetisation_error = np.append(magnetisation_error, magnetisation_obs.dvalue)

            energy_obs = compute_renormalised_energy(direc, start, b)
            energy_obs_array = np.append(energy_obs_array, energy_obs)
            energy_array = np.append(energy_array, energy_obs.value)
            energy_error = np.append(energy_error, energy_obs.dvalue)


        # temp_array, magnetisation_array, magnetisation_error, energy_array, energy_error = zip(*sorted(zip(temp_array, magnetisation_array, magnetisation_error, energy_array, energy_error)))
        length_array, temp_array, magnetisation_array, magnetisation_error, energy_array, energy_error = zip(*sorted(zip(length_array,temp_array, magnetisation_array, magnetisation_error, energy_array, energy_error)))

        if type(unique_lengths) == int:
            unique_lengths = [unique_lengths]
        for j, l in enumerate(unique_lengths):
            if ((b == 1) or (l==256)):
                xaxis = np.array(temp_array[j*len(unique_temp): (j+1)*len(unique_temp)])
                yaxis = np.array(magnetisation_array[j * len(unique_temp):(j+1)*len(unique_temp)])
                yerr = np.array(magnetisation_error[j * len(unique_temp):(j+1)*len(unique_temp)])
                ax1.errorbar(xaxis, yaxis, yerr, label = f"b = {b}, L = {l}", color = colors[i], fmt = markers[j])

                xaxis = np.array(temp_array[j*len(unique_temp): (j+1)*len(unique_temp)])
                yaxis = np.array(energy_array[j * len(unique_temp):(j+1)*len(unique_temp)])
                yerr = np.array(energy_error[j * len(unique_temp):(j+1)*len(unique_temp)])
                ax2.errorbar(xaxis, yaxis, yerr, label = f"b = {b}, L = {l}", color = colors[i], fmt =markers[j])
        
    ax1.set_xlabel(r"Temperature $T$")
    ax1.set_ylabel(r"Magnetisation $m$")
    ax1.set_title(f"RG for $L = ${max(unique_lengths)} and its fractions $b = ${b_list}")
    ax1.legend()

    ax2.set_xlabel(r"Temperature $T$")
    ax2.set_ylabel(r"Energy per spin $e$")
    ax2.set_title(f"RG for $L = ${max(unique_lengths)} and its fractions $b = ${b_list}")
    ax2.legend()

    fig_mag.savefig(saving_path / f"{direc.name}_biggest_rg_magnetisation.pdf", bbox_inches = "tight")
    fig_energy.savefig(saving_path / f"{direc.name}_biggest_rg_energy.pdf", bbox_inches = "tight")
