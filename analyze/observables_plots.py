import os
import sys
import pathlib

from utils.h5_utils import import_observable
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import toml
import json
import pyerrors as pe
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from compute_observables import (
    compute_average_magnetisation,
    compute_average_energy,
    compute_susceptibility,
    compute_specific_heat,
    compute_binder_cumulant,
    compute_normalised_energy,
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
        "normalised_energy":compute_normalised_energy,
        }
    return observable_functions[observable]

def do_observable_plot(
        observable:str,
        observable_title:str,
        directory:pathlib.Path, 
        x_data:str = "temperature",
        x2_data:str | None = None,
        log_plot:bool = False,
        log_fit:bool = False,
        linear_fit:bool = False,
    ):

    saving_path = directory.parent.parent / "analyze" / "output"/"img_dump"
    params = [x for x in directory.iterdir() if x.is_dir()]
    temp_array = np.array([])
    length_array = np.array([])
    tauint_array = np.array([])
    tauint_error = np.array([])
    observable_array = np.array([])
    observable_error = np.array([])

    compute_observable = _find_observable_function(observable)

    with open(directory / "global_parameters.json", "r") as f:
        global_config = json.load(f)
    unique_lengths = global_config["physical_settings"]["L"]
    unique_temp = global_config["physical_settings"]["temperature"]


    for direc in params:
        config = toml.load(direc / "config.toml")
        print("direc", direc)
        temp_array = np.append(temp_array, import_physical_parameter(direc, "temperature"))
        length_array = np.append(length_array, import_physical_parameter(direc, "L")) 
        observable_obs = compute_observable(direc)
        keys = list(observable_obs.e_tauint.keys())
        print("tauint keys", keys)
        tauint_array = np.append(tauint_array, observable_obs.e_tauint[keys[0]])
        tauint_error = np.append(tauint_error, observable_obs.e_dtauint[keys[0]])
        observable_array = np.append(observable_array, observable_obs.value)
        observable_error = np.append(observable_error, observable_obs.dvalue)

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )

    figtau, axtau = plt.subplots(
        ncols=1,
        nrows=1
    )

    if x_data == "temperature":
        _,sorted_temp, observable_array, observable_error = zip(*sorted(zip(length_array, temp_array, observable_array, observable_error)))
        _,sorted_temp, tauint_array, tauint_error = zip(*sorted(zip(length_array, temp_array, tauint_array, tauint_error)))
        cmap = plt.cm.tab20
        colors = cmap(np.arange(len(unique_lengths)*2))

        for i,l in enumerate(unique_lengths):
            ax = _add_scatter_data(
                axs = ax,
                xaxis = sorted_temp[i*len(unique_temp): (i+1)*len(unique_temp)],
                yaxis = observable_array[i * len(unique_temp):(i+1)*len(unique_temp)],
                yerr = observable_error[i * len(unique_temp):(i+1)*len(unique_temp)],
                data_label = f"$L = {l}$",
                linear_fit = linear_fit,
                log_fit = log_fit,
                main_color = colors[i],
                secondary_color = colors[i+1],
            )

            axtau = _add_scatter_data(
                axs = axtau,
                xaxis = sorted_temp[i*len(unique_temp): (i+1)*len(unique_temp)],
                yaxis = tauint_array[i * len(unique_temp):(i+1)*len(unique_temp)],
                yerr = tauint_error[i * len(unique_temp):(i+1)*len(unique_temp)],
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

        axtau = _add_format_plot(
            axs = axtau,
            xlabel="Temperature, $T$",
            ylabel=r"$\tau_{\textrm{int}}$",
            title = f" Autocorrelated time for {observable.capitalize()}",
            logscale = log_plot,
            linear_fit = linear_fit,
        )

    elif x_data == "length":
        _,sorted_length = zip(*sorted(zip(temp_array, length_array)))
        cmap = plt.cm.tab20
        colors = cmap(np.arange(len(unique_temp)*2))

        for i,t in enumerate(unique_temp):
            ax = _add_scatter_data(
                axs = ax,
                xaxis = length_array[len(unique_lengths) * i : (i+1) * len(unique_lengths)],
                yaxis = observable_array,
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
    fig.savefig(saving_path / f"{directory.name}_{observable}.pdf")
    figtau.savefig(saving_path / f"{directory.name}_tauint_{observable}.pdf")
    plt.close(fig)


def _add_scatter_data(
        axs:plt.axes,
        xaxis:np.ndarray,
        yaxis:np.ndarray,
        yerr:np.ndarray | float,
        data_label:str = "",
        linear_fit:bool = False,
        log_fit:bool = False,
        main_color: str = 'k',
        secondary_color: str = 'coral'
    ) -> plt.axes:

    assert len(xaxis) == len(yaxis), "X and Y data need to be the same length"
    
    if data_label != "":
        axs.errorbar(xaxis, yaxis, yerr = yerr, label = data_label, color = main_color, fmt = ".")
    else:
        axs.errorbar(xaxis, yaxis, yerr = yerr,color = main_color, fmt = ".")


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
    axs.yaxis.set_minor_formatter(latex_formatter)
    axs.xaxis.set_minor_formatter(latex_formatter)

    axs.legend()

    return axs


def _generate_scatter_plot(
        axs:plt.axes,
        xaxis:np.ndarray,
        xlabel:str,
        yaxis:np.ndarray,
        ylabel:str,
        title:str,
        data_label:str = "",
        logscale:bool = False,
        linear_fit:bool = False,
        log_fit:bool = False,
    ) -> plt.axes:

    assert len(xaxis) == len(yaxis), "X and Y data need to be the same length"
    main_color = 'k'
    secondary_color = 'coral'
    
    if data_label != "":
        axs.scatter(xaxis, yaxis, label = data_label, color = main_color)
    else:
        axs.scatter(xaxis, yaxis, color = main_color)

    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_title(title)


    if linear_fit:
        p_fitted = np.polynomial.Polynomial.fit(xaxis, yaxis, deg=1)
        lin_fit = p_fitted.convert().coef
        axs.plot(xaxis, xaxis * lin_fit[1] + lin_fit[0], label="Linear fit", alpha = 0.7, ls="--", color =secondary_color)
    
    if log_fit:
        p_fitted = np.polynomial.Polynomial.fit(np.log(xaxis), np.log(yaxis), deg=1)
        log_fit = p_fitted.convert().coef
        axs.plot(xaxis, np.exp(log_fit[0]) * (xaxis ** log_fit[1]), label="Log fit", alpha = 0.5, ls="--", color = secondary_color)


    if (logscale and not linear_fit):
        axs.set_yscale('log')
        axs.set_xscale('log')

    # Make latex formatter
    latex_formatter = ticker.FuncFormatter(lambda x, pos: f'${x:g}$')

    # Apply formatting
    axs.xaxis.set_major_formatter(latex_formatter)
    axs.yaxis.set_major_formatter(latex_formatter)
    axs.yaxis.set_minor_formatter(latex_formatter)
    axs.xaxis.set_minor_formatter(latex_formatter)

    axs.legend()

    return axs
