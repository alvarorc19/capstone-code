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
    compute_average_magnetisation,
    compute_average_energy,
    compute_susceptibility,
    compute_susceptibility_per_spin,
    compute_specific_heat,
    compute_binder_cumulant,
    compute_normalised_energy,
    compute_specific_heat_per_spin,
    compute_renormalised_energy,
    compute_renormalised_magnetisation,
    compute_cluster_susceptibility,
    compute_cluster_susceptibility_per_spin,
    compute_cluster_size,
    compute_cluster_size_per_spin,
)
from utils.h5_utils import (
    import_observable,
    import_physical_parameter,
)


def _find_observable_function(observable:str):
    observable_functions = {
        "magnetisation":compute_average_magnetisation,
        "energy":compute_average_energy,
        "susceptibility":compute_susceptibility,
        "specific_heat":compute_specific_heat,
        "binder_cumulant":compute_binder_cumulant,
        "energy_per_spin":compute_normalised_energy,
        "susceptibility_per_spin":compute_susceptibility_per_spin,
        "specific_heat_per_spin":compute_specific_heat_per_spin,
        "cluster_susceptibility":compute_cluster_susceptibility,
        "cluster_susceptibility_per_spin":compute_cluster_susceptibility_per_spin,
        "correlation_length":compute_cluster_size,
        "correlation_length_per_spin":compute_cluster_size_per_spin,
    }
    return observable_functions[observable]

def get_observables_csv(
        directory:pathlib.Path, 
        is_deep:bool = False,
        start:int = 0,
        rg:bool = False
    ):

    saving_path = directory
    df = pd.DataFrame()
    observables = ["magnetisation", "energy", "susceptibility", "specific_heat", "energy_per_spin", "susceptibility_per_spin", "specific_heat_per_spin", "cluster_susceptibility"]

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

    temp_array = np.array([])
    length_array = np.array([])

    for direc in params:
        temp_array = np.append(temp_array, import_physical_parameter(direc, "temperature"))
        length_array = np.append(length_array, import_physical_parameter(direc, "L")) 
    df["temperature"] = temp_array
    df["L"] = length_array


    for observable in observables:
        observable_array = np.array([])
        observable_error = np.array([])
        observable_tauint = np.array([])
        compute_observable = _find_observable_function(observable)
        for direc in params:
            try:
                observable_obs = compute_observable(direc, start)
                observable_array = np.append(observable_array, observable_obs.value)
                observable_error = np.append(observable_error, observable_obs.dvalue)
                observable_tauint = np.append(observable_tauint, observable_obs.e_tauint["ens"])

            except Exception:
                observable_array = np.append(observable_array, 0)
                observable_error = np.append(observable_error, 0)
                observable_tauint = np.append(observable_tauint, 0)

        df[observable + "_value"] = observable_array
        df[observable + "_error"] = observable_error
        df[observable + "_tau_int"] = observable_tauint

    if rg:
        dim = import_physical_parameter(direc, "dimension")
        if dim == 2:
            blist = [2,4,16]

        elif dim ==3:
            blist = [2,4]


        for b in blist:
            magnetisation_array = np.array([])
            magnetisation_error = np.array([])
            magnetisation_tauint = np.array([])
            energy_array = np.array([])
            energy_error = np.array([])
            energy_tauint = np.array([])
            for direc in params:
                magnetisation_obs = compute_renormalised_magnetisation(direc, start, b)
                magnetisation_array = np.append(magnetisation_array, magnetisation_obs.value)
                magnetisation_error = np.append(magnetisation_error, magnetisation_obs.dvalue)
                magnetisation_tauint = np.append(magnetisation_tauint, magnetisation_obs.e_tauint["ens"])

                energy_obs = compute_renormalised_energy(direc, start, b)
                energy_array = np.append(energy_array, energy_obs.value)
                energy_error = np.append(energy_error, energy_obs.dvalue)
                energy_tauint = np.append(energy_tauint, energy_obs.e_tauint["ens"])
            df[f"magnetisation_b_{b}"] = magnetisation_array
            df[f"magnetisation_b_{b}_error"] = magnetisation_error
            df[f"magnetisation_b_{b}_tauint"] = magnetisation_tauint
            df[f"energy_b_{b}"] = energy_array
            df[f"energy_b_{b}_error"] = energy_error
            df[f"energy_b_{b}_tauint"] = energy_tauint
        

    df = df.sort_values(["L", "temperature"]).reset_index(drop = True)
    df.to_csv(saving_path / "ensemble_observables.csv", index = False)
    print("Finished writing to ", saving_path / "ensemble_observables.csv")



def do_observable_plot(
        observable:str,
        observable_title:str,
        directory:pathlib.Path, 
        is_deep:bool = False,
        x_data:str = "temperature",
        log_plot:bool = False,
        log_fit:bool = False,
        linear_fit:bool = False,
        start:int = 0
    ):

    saving_path = directory.parent.parent / "analyze" / "output"/"img_dump"
    plt.tight_layout()

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

    temp_array = np.array([])
    length_array = np.array([])
    observable_array = np.array([])
    observable_error = np.array([])
    observable_obs_array = np.array([])

    compute_observable = _find_observable_function(observable)

    with open(directory / "global_parameters.json", "r") as f:
        global_config = json.load(f)
    unique_lengths = global_config["physical_settings"]["L"]
    unique_temp = global_config["physical_settings"]["temperature"]


    for direc in params:
        temp_array = np.append(temp_array, import_physical_parameter(direc, "temperature"))
        length_array = np.append(length_array, import_physical_parameter(direc, "L")) 
        observable_obs = compute_observable(direc, start)
        observable_obs_array = np.append(observable_obs_array, observable_obs)
        observable_array = np.append(observable_array, observable_obs.value)
        observable_error = np.append(observable_error, observable_obs.dvalue)


    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize = (9,6)
    )

    cmap = plt.cm.tab20
    if x_data == "temperature":

        if isinstance(unique_lengths, list):
            length_array,temp_array, observable_array, observable_error = zip(*sorted(zip(length_array, temp_array, observable_array, observable_error)))
            colors = cmap(np.arange(len(unique_lengths)*2))
            for i,l in enumerate(unique_lengths):
                ax = _add_scatter_data(
                    axs = ax,
                    xaxis = temp_array[i*len(unique_temp): (i+1)*len(unique_temp)],
                    yaxis = observable_array[i * len(unique_temp):(i+1)*len(unique_temp)],
                    yerr = observable_error[i * len(unique_temp):(i+1)*len(unique_temp)],
                    data_label = f"$L = {l}$",
                    linear_fit = linear_fit,
                    log_fit = log_fit,
                    main_color = colors[i],
                    secondary_color = colors[i+1],
                )

            ax = _add_format_plot(
                axs = ax,
                xlabel="Temperature, $T$",
                ylabel=observable_title,
                title = f"{observable.capitalize()} vs Temperature",
                logscale = log_plot,
                linear_fit = linear_fit,
            )
        elif isinstance(unique_lengths, int):
            l = unique_lengths
            i = 0
            colors = cmap([0,1])
            ax = _add_scatter_data(
                axs = ax,
                xaxis = temp_array,
                yaxis = observable_array,
                yerr = observable_error,
                data_label = f"$L = {l}$",
                linear_fit = linear_fit,
                log_fit = log_fit,
                main_color = colors[i],
                secondary_color = colors[i+1],
            )

            ax = _add_format_plot(
                axs = ax,
                xlabel="Temperature, $T$",
                ylabel=observable_title,
                title = f"{observable.capitalize()} vs Temperature",
                logscale = log_plot,
                linear_fit = linear_fit,
            )

    elif x_data == "length":
        length_array,temp_array, observable_array, observable_error = zip(*sorted(zip(length_array, temp_array, observable_array, observable_error)))
        colors = cmap(np.arange(len(unique_temp)*2))

        for i,t in enumerate(unique_temp):
            ax = _add_scatter_data(
                axs = ax,
                xaxis = length_array[len(unique_lengths) * i : (i+1) * len(unique_lengths)],
                yaxis = observable_array[i * len(unique_lengths):(i+1)*len(unique_lengths)],
                yerr = observable_error[i * len(unique_lengths):(i+1)*len(unique_lengths)],
                data_label = f"$T = {t}$",
                linear_fit = linear_fit,
                log_fit = log_fit,
                main_color = colors[i],
                secondary_color = colors[i+1],
            )

        ax = _add_format_plot(
            axs = ax,
            xlabel="Length, $L$",
            ylabel=observable_title,
            title = f"{observable.capitalize()} vs Length",
            logscale = log_plot,
            linear_fit = linear_fit,
        )
    else:
        print("This is not incorporated, say temperature or length")


    saving_path.mkdir(parents = True, exist_ok = True)
    print("saving path ", saving_path)
    fig.savefig(saving_path / f"{directory.name}_{observable}.pdf",bbox_inches = "tight")
    plt.close(fig)
    print(f"finished {observable} plots")


def _add_scatter_data(
        axs:plt.axes,
        xaxis:np.ndarray,
        yaxis:np.ndarray,
        yerr:np.ndarray | float,
        data_label:str = "",
        linear_fit:bool = False,
        log_fit:bool = False,
        main_color: str = 'k',
        secondary_color: str = 'coral',
        marker: str = ".",
    ) -> plt.axes:

    assert len(xaxis) == len(yaxis), "X and Y data need to be the same length"
    xaxis = np.array(xaxis)
    yaxis = np.array(yaxis)
    yerr = np.array(yerr)

    # Sorting
    idx = np.argsort(xaxis)
    xaxis = np.array(xaxis)
    yaxis = np.array(yaxis)
    yerr = np.array(yerr)
    xaxis = xaxis[idx]
    yaxis = yaxis[idx]
    yerr = yerr[idx]
    
    if data_label != "":
        axs.errorbar(xaxis, yaxis, yerr = yerr, label = data_label, color = main_color, fmt = marker)
    else:
        axs.errorbar(xaxis, yaxis, yerr = yerr,color = main_color, fmt = marker)


    if linear_fit:
        p_fitted = np.polynomial.Polynomial.fit(xaxis, yaxis, deg=1)
        lin_fit = p_fitted.convert().coef
        axs.plot(xaxis, xaxis * lin_fit[1] + lin_fit[0], label="Linear fit", alpha = 0.7, ls="--", color =secondary_color)
    
    if log_fit:
        p_fitted = np.polynomial.Polynomial.fit(np.log(xaxis), np.log(yaxis), deg=1)
        log_fit = p_fitted.convert().coef
        axs.plot(xaxis, np.exp(log_fit[0]) * (xaxis ** log_fit[1]), label="Log fit", alpha = 0.5, ls="--", color = secondary_color)

    return axs

def _add_format_plot(
        axs:plt.axes,
        xlabel:str,
        ylabel:str,
        title:str,
        logscale:bool = False,
        linear_fit:bool = False,
    ) -> plt.axes:

    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_title(title)

    if (logscale and not linear_fit):
        axs.set_yscale('log')
        axs.set_xscale('log')

    # Make latex formatter
    latex_formatter = ticker.FuncFormatter(lambda x, pos: f'${x:g}$')

    # Apply formatting
    axs.xaxis.set_major_formatter(latex_formatter)
    axs.yaxis.set_major_formatter(latex_formatter)
    # axs.yaxis.set_minor_formatter(latex_formatter)
    # axs.xaxis.set_minor_formatter(latex_formatter)
    handles, labels = axs.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    # axs.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.1,1.05), fancybox = True, shadow = True)
    box = axs.get_position()
    axs.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])

# Put a legend below current axis
    axs.legend(by_label.values(), by_label.keys(),loc='upper center', bbox_to_anchor=(0.5, -0.1),
              fancybox=True, shadow=True, ncol=5)

    axs.grid()
    axs.minorticks_on()
    axs.tick_params(axis='both', which='major', direction='in', length=7, top=True, right = True)
    axs.tick_params(axis='both', which='minor', direction='in', length=4, top=True, right = True)

    return axs

def do_order_parameter_plot(directory:pathlib.Path, is_deep:bool = False, start:int = 0):
    # Add Obs crap

    plt.tight_layout()
    saving_path = directory.parent.parent / "analyze" / "output"/f"thermalisation_{directory.name}"

    saving_path.mkdir(parents = True, exist_ok = True)
    if is_deep:
        directory = directory / f"{directory.name}_0"
        params = [x for x in directory.iterdir() if x.is_dir()]
        # sub_dir = [x for x in directory.iterdir() if x.is_dir()]
        # params = []
        # for dir in sub_dir:
        #     if dir.is_dir():
        #         for direc in dir.iterdir():
        #             if direc.is_dir():
        #                 params.append(direc)
    else:
        params = [x for x in directory.iterdir() if x.is_dir()]


    for direc in params:
        magnetisation = import_observable(direc, "magnetisation")
        x_magnetisation = np.array([i[0] for i in magnetisation[start:]])
        y_magnetisation = np.array([i[1] for i in magnetisation[start:]])
        length = import_physical_parameter(direc, "L")
        dim = import_physical_parameter(direc, "dimension")
        N = length ** dim
        magnetisation_array = np.sqrt(x_magnetisation**2 + y_magnetisation**2) / N
        energy = import_observable(direc, "energy")
        cluster_size = import_observable(direc, "average_cluster_size")
        temperature = import_physical_parameter(direc, "temperature")
        susceptibility = ( 1 / temperature) * cluster_size

        fig_mag, ax1 = plt.subplots(
            ncols=1,
            nrows=1,
            figsize = (9,6)
        )

        ax1.plot(magnetisation_array, color = "orangered", alpha = 0.5)
        ax1.set_ylabel(r"Magnetisation $m$")
        ax1.set_xlabel(r"Time $t$ in MCS")
        fig_mag.savefig(saving_path / f"{direc.name}_magnetisation_time.pdf", bbox_inches = "tight")
        plt.close(fig_mag)


        fig_energy, ax2 = plt.subplots(
            ncols=1,
            nrows=1,
            figsize = (9,6)
        )
        fig_energy.set_facecolor("#ECEFF4")
        ax2.set_facecolor("#ECEFF4")

        ax2.plot(energy[start:], color = "chartreuse", alpha = 0.5)
        ax2.set_ylabel(r"Energy $E$")
        ax2.set_xlabel(r"Time $t$ in MCS")

        fig_energy.savefig(saving_path / f"{direc.name}_energy_time.pdf", bbox_inches = "tight")
        plt.close(fig_energy)

        fig_cluster_size, ax3 = plt.subplots(
            ncols=1,
            nrows=1,
            figsize = (9,6)
        )

        ax3.plot(cluster_size[start:], color = "navy", alpha = 0.5)
        ax3.set_ylabel(r"Average cluster size $\langle n \rangle$")
        ax3.set_xlabel(r"Time $t$ in MCS")

        fig_cluster_size.savefig(saving_path / f"{direc.name}_cluster_size_time.pdf", bbox_inches = "tight")
        plt.close(fig_cluster_size)

        fig_susc, ax3 = plt.subplots(
            ncols=1,
            nrows=1,
            figsize = (9,6)
        )

        ax3.plot(susceptibility[start:], color = "red", alpha = 0.5)
        ax3.set_ylabel(r"Susceptibility $\chi$")
        ax3.set_xlabel(r"Time $t$ in MCS")
        fig_susc.savefig(saving_path / f"{direc.name}_susceptibility_time.pdf", bbox_inches = "tight")
        plt.close(fig_susc)
        print("finished order parameter plots")


def do_magnetisation_inflection_plot(
        directory:pathlib.Path, 
        is_deep:bool = False,
        start:int = 0
    ):

    saving_path = directory.parent.parent / "analyze" / "output"/"img_dump"
    plt.tight_layout()

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

    temp_array = np.array([])
    length_array = np.array([])
    magnetisation_array = np.array([])
    magnetisation_error = np.array([])

    with open(directory / "global_parameters.json", "r") as f:
        global_config = json.load(f)
    unique_lengths = global_config["physical_settings"]["L"]
    unique_temp = global_config["physical_settings"]["temperature"]


    for direc in params:
        # print("direc", direc)
        temp_array = np.append(temp_array, import_physical_parameter(direc, "temperature"))
        length_array = np.append(length_array, import_physical_parameter(direc, "L")) 
        magnetisation_obs = compute_average_magnetisation(direc, start)
        magnetisation_array = np.append(magnetisation_array, magnetisation_obs.value)
        magnetisation_error = np.append(magnetisation_error, magnetisation_obs.dvalue)


    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize = (9,6)
    )

    cmap = plt.cm.tab20
    if isinstance(unique_lengths, list):
        length_array,temp_array, magnetisation_array, magnetisation_error = zip(*sorted(zip(length_array, temp_array, magnetisation_array, magnetisation_error)))
        colors = cmap(np.arange(len(unique_lengths)*2))
        for i,l in enumerate(unique_lengths):
            xaxis = np.array(temp_array[i*len(unique_temp): (i+1)*len(unique_temp)])
            yaxis = np.array(magnetisation_array[i * len(unique_temp):(i+1)*len(unique_temp)])
            yerr = np.array(magnetisation_error[i * len(unique_temp):(i+1)*len(unique_temp)])
            # Sorting
            idx = np.argsort(xaxis)
            xaxis = xaxis[idx]
            yaxis = yaxis[idx]
            yerr = yerr[idx]
            x_infl, y_infl, infl_err = _find_inflection_points(xaxis, yaxis)
            ax.errorbar(x_infl, y_infl, infl_err, label="Inflection points", color = "r", fmt = ".")
            ax = _add_scatter_data(
                axs = ax,
                xaxis = xaxis,
                yaxis = yaxis,
                yerr = yerr,
                data_label = f"$L = {l}$",
                main_color = colors[i],
                secondary_color = colors[i+1],
                marker = ".-",
            )

        ax = _add_format_plot(
            axs = ax,
            xlabel="Temperature, $T$",
            ylabel=r"Magnetisation $\langle |\mathbf{m}| \rangle$",
            title = f"Magnetisation vs Temperature, inflection points",
        )
    elif isinstance(unique_lengths, int):
        l = unique_lengths
        i = 0
        colors = cmap([0,1])
        xaxis = np.array(temp_array)
        yaxis = np.array(magnetisation_array)
        yerr = np.array(magnetisation_error)
        # Sorting
        idx = np.argsort(xaxis)
        xaxis = xaxis[idx]
        yaxis = yaxis[idx]
        yerr = yerr[idx]
        x_infl, y_infl, x_infl_err = _find_inflection_points(xaxis, yaxis)
        print("x_infl ",x_infl) 
        print("y_infl ",y_infl) 
        print("ingl_Er ",x_infl_err) 
        ax = _add_scatter_data(
            axs = ax,
            xaxis = temp_array,
            yaxis = magnetisation_array,
            yerr = magnetisation_error,
            data_label = f"$L = {l}$",
            main_color = colors[i],
            secondary_color = colors[i+1],
            marker = ".-",
        )
        ax.errorbar(x_infl, y_infl, xerr = x_infl_err, label=f"Inflection $T_c = {x_infl[5]:.4f}\pm{x_infl_err[5]:.4f}$", color = "r", fmt=".")

        ax = _add_format_plot(
            axs = ax,
            xlabel="Temperature, $T$",
            ylabel=r"Magnetisation $\langle \mathbf{m} \rangle$",
            title = f"Magnetisation vs Temperature, inflection points",
        )

    saving_path.mkdir(parents = True, exist_ok = True)
    fig.savefig(saving_path / f"{directory.name}_magnetisation_inflection_point.pdf")
    print("finished magnetisation inflection plots")

def do_inflection_vs_length_plot(
        directory:pathlib.Path, 
        is_deep:bool = False,
        start:int = 0
    ):

    if is_deep:
        saving_path = directory.parent.parent.parent / "analyze" / "output"/"img_dump"
    else:
        saving_path = directory.parent.parent / "analyze" / "output"/"img_dump"
    plt.tight_layout()

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

    temp_array = np.array([])
    length_array = np.array([])
    magnetisation_array = np.array([])
    magnetisation_error = np.array([])

    with open(directory / "global_parameters.json", "r") as f:
        global_config = json.load(f)
    unique_lengths = global_config["physical_settings"]["L"]
    unique_temp = global_config["physical_settings"]["temperature"]


    for direc in params:
        # print("direc", direc)
        temp_array = np.append(temp_array, import_physical_parameter(direc, "temperature"))
        length_array = np.append(length_array, import_physical_parameter(direc, "L")) 
        magnetisation_obs = compute_average_magnetisation(direc, start)
        magnetisation_array = np.append(magnetisation_array, magnetisation_obs.value)
        magnetisation_error = np.append(magnetisation_error, magnetisation_obs.dvalue)


    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )

    cmap = plt.cm.tab20
    length_array,temp_array, magnetisation_array, magnetisation_error = zip(*sorted(zip(length_array, temp_array, magnetisation_array, magnetisation_error)))
    colors = cmap(np.arange(len(unique_lengths)*2))
    inflection_array = []
    inflection_err_array = []
    for i,l in enumerate(unique_lengths):
        xaxis = np.array(temp_array[i*len(unique_temp): (i+1)*len(unique_temp)])
        yaxis = np.array(magnetisation_array[i * len(unique_temp):(i+1)*len(unique_temp)])
        yerr = np.array(magnetisation_error[i * len(unique_temp):(i+1)*len(unique_temp)])
        # Sorting
        idx = np.argsort(xaxis)
        xaxis = xaxis[idx]
        yaxis = yaxis[idx]
        yerr = yerr[idx]
        x_infl, y_infl, infl_err = _find_inflection_points(xaxis, yaxis)
        for i,x in enumerate(x_infl):
            if ( x < 1.3 and x > 0.7):
                ax.errorbar(1/l, x, yerr = infl_err[i], color = colors[0], fmt = ".")
                inflection_array.append(x)
                inflection_err_array.append(infl_err[i])

    ax = _add_format_plot(
        axs = ax,
        xlabel=r"Critical temperature $T_c$",
        ylabel="Length $L$",
        title = f"Inflection point vs length",
        logscale=False,
    )
    # ax.set_xlim(left=0.8)
    ax.grid(which="both")
    ax.set_xscale('log')
    print("infl_arr = ", inflection_array)

    saving_path.mkdir(parents = True, exist_ok = True)
    fig.savefig(saving_path / f"{directory.name}_t_vs_l.pdf")
    print("finished t vs l inflection plots")

def _find_inflection_points(x:np.ndarray, y:np.ndarray):
    assert len(y) == len(x), "magnetisation and temp should be same length" 

    dy = np.gradient(y, x)
    d2y = np.gradient(dy, x)
    d3y = np.gradient(d2y, x)

    refined_x = []
    refined_y = []
    x_errors = []

    for i in range(len(d2y) - 1):
        if d2y[i] == 0:
            refined_x.append(x[i])
            errors.append(0)
        elif d2y[i] * d2y[i+1] < 0:
            x0, x1 = x[i], x[i+1]
            y0, y1 = d2y[i], d2y[i+1]

            x_zero = x0 - y0 * (x1 - x0) / (y1 - y0)

            y0, y1 = y[i], y[i + 1]
            y_zero = y0 + (y1 - y0) * (x_zero - x0) / (x1 - x0)


            refined_x.append(x_zero)
            refined_y.append(y_zero)

            hd = x[i+1] - x[i]
            hs = x[i] - x[i-1] if i > 0 else hd

            d3 = abs(d3y[i])

            error = (hd * hs) / (max(d3, 1e-12))
            x_errors.append(error)
    return np.array(refined_x), np.array(refined_y), np.array(x_errors)
