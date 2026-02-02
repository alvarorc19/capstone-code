import os
import sys
import pathlib

from utils.h5_utils import import_observable
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import toml
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from compute_observables import (
        compute_average_magnetisation,
        compute_average_energy,
        compute_susceptibility,
        compute_specific_heat,
        compute_binder_cumulant,
    )
from utils.h5_utils import import_observable


def _find_observable_function(observable:str):
    observable_functions = {
        "magnetisation":compute_average_magnetisation,
        "energy":compute_average_energy,
        "susceptibility":compute_susceptibility,
        "specific_heat":compute_specific_heat,
        "binder_cumulant":compute_binder_cumulant,
        }
    return observable_functions[observable]

def do_observable_plot(
        observable:str,
        observable_title:str,
        directory:pathlib.Path, 
        x_data:str = "temperature",
        log_plot:bool = False,
        log_fit:bool = False,
        linear_fit:bool = False,
    ):

    saving_path = directory.parent.parent / "analyze" / "img_dump"
    params = [x for x in directory.iterdir() if x.is_dir()]
    x_array = np.array([])
    observable_array = np.array([])

    compute_observable = _find_observable_function(observable)

    for direc in params:
        config = toml.load(direc / "config.toml")
        print("direc", direc)
        x_array = np.append(x_array, import_physical_parameter(direc, x_data))
        observable_array = np.append(observable_array, compute_observable(direc, config))

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )
    if x_data == "temperature":
        ax = _generate_scatter_plot(
            axs = ax,
            xaxis = x_array,
            xlabel="Temperature, $T$",
            yaxis = observable_array,
            ylabel=observable_title,
            title = f"{observable} vs temperature",
            data_label = "Data",
            logscale = log_plot,
            linear_fit = linear_fit,
            log_fit = log_fit,
        )

    elif x_data == "length":
        ax = _generate_scatter_plot(
            axs = ax,
            xaxis = x_array,
            xlabel="Length, $L$",
            yaxis = observable_array,
            ylabel= observable_title,
            title = f"{observable} vs temperature",
            data_label = "Data",
            logscale = log_plot,
            linear_fit = linear_fit,
            log_fit = log_fit,
        )
    else:
        print("This is not incorporated, say temperature or length")


    saving_path.mkdir(parents = True, exist_ok = True)
    fig.savefig(saving_path / f"{directory.name}_magnetisation.pdf")


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
