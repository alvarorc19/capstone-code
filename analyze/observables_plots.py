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
from scipy.optimize import curve_fit
from tqdm import tqdm
from matplotlib.widgets import Slider

from compute_observables import (
    compute_average_magnetisation,
    compute_average_energy,
    compute_susceptibility,
    compute_susceptibility_per_spin,
    compute_specific_heat,
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
    observables = ["magnetisation", "energy", "susceptibility", "specific_heat", "energy_per_spin", "susceptibility_per_spin", "specific_heat_per_spin", "cluster_susceptibility","cluster_susceptibility_per_spin", "correlation_length", "correlation_length_per_spin"]

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
    length_array = np.array([], dtype = int)

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
            df[f"magnetisation_b_{b}_value"] = magnetisation_array
            df[f"magnetisation_b_{b}_error"] = magnetisation_error
            df[f"magnetisation_b_{b}_tauint"] = magnetisation_tauint
            df[f"energy_b_{b}_value"] = energy_array
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
    # plt.tight_layout()
    csv_file = directory / "ensemble_observables.csv"
    if csv_file.exists():
        df = pd.read_csv(csv_file)
    else:
        get_observables_csv(directory, is_deep, start, False)
        df = pd.read_csv(csv_file)
    df = df[df["energy_value"] != 0.0].reset_index(drop=True)

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize = (9,6)
    )

    cmap = plt.cm.tab20
    colors = cmap(np.arange(30))
    if x_data == "temperature":
        i = 0
        for l, group in df.groupby("L"):
            group = group.sort_values("temperature")
            xaxis = group["temperature"].to_numpy()
            yaxis = group[f"{observable}_value"].to_numpy()
            yerr = group[f"{observable}_error"].to_numpy()
            ax = _add_scatter_data(
                axs = ax,
                xaxis = xaxis,
                yaxis = yaxis,
                yerr = yerr,
                data_label = f"$L/a = {l}$",
                main_color = colors[i],
                secondary_color = colors[i+1],
                marker = ".-",
            )

            title = observable.replace("_", " ")
            ax = _add_format_plot(
                axs = ax,
                xlabel="Temperature, $k_BT$",
                ylabel=observable_title,
                # title = f"{title.capitalize()} vs Temperature",
                logscale = log_plot,
                linear_fit = linear_fit,
            )
            i+=1

    elif x_data == "length":
        j = 0
        for t, group in df.groupby("temperature"):
            group = group.sort_values("L")
            xaxis = group["L"].to_numpy()
            yaxis = group[f"{observable}_value"].to_numpy()
            yerr = group[f"{observable}_error"].to_numpy()
            ax = _add_scatter_data(
                axs = ax,
                xaxis = xaxis,
                yaxis = yaxis,
                yerr = yerr,
                data_label = f"$k_BT = {t}$",
                linear_fit = linear_fit,
                log_fit = log_fit,
                main_color = colors[i],
                secondary_color = colors[i+1],
            )

            title = observable.replace("_", " ")
            ax = _add_format_plot(
                axs = ax,
                xlabel="Length, $L/a$",
                ylabel=observable_title,
                # title = f"{title.capitalize()} vs Length",
                logscale = log_plot,
                linear_fit = linear_fit,
            )
            j+=1
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
        title:str | None = None,
        logscale:bool = False,
        linear_fit:bool = False,
    ) -> plt.axes:

    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    if title is not None:
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

    # plt.tight_layout()
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
    csv_file = directory / "ensemble_observables.csv"
    if csv_file.exists():
        df = pd.read_csv(csv_file)
    else:
        get_observables_csv(directory, is_deep, start, False)
        df = pd.read_csv(csv_file)
    # plt.tight_layout()
    df = df[df["energy_value"] != 0.0].reset_index(drop=True)


    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize = (9,6)
    )

    cmap = plt.cm.tab20
    colors = cmap(np.arange(20))
    i = 0
    for l, group in df.groupby("L"):
        group = group.sort_values("temperature")
        xaxis = group["temperature"].to_numpy()
        yaxis = group[f"magnetisation_value"].to_numpy()
        yerr = group[f"magnetisation_error"].to_numpy()

        idx = np.argsort(xaxis)
        xaxis = xaxis[idx]
        yaxis = yaxis[idx]
        yerr = yerr[idx]

        params, pcov = _fit_mag_to_tanh(xaxis, yaxis, yerr)
        temperature_fit = np.linspace(np.min(xaxis), np.max(xaxis), 300)
        ax.plot(temperature_fit, _tanh(temperature_fit, *params),"--", color = colors[i])

        ax = _add_scatter_data(
            axs = ax,
            xaxis = xaxis,
            yaxis = yaxis,
            yerr = yerr,
            data_label = f"$L/a = {int(l)}$",
            main_color = colors[i],
            secondary_color = colors[i+1],
            marker = "."
        )

        ax = _add_format_plot(
            axs = ax,
            xlabel="Temperature, $k_BT$",
            ylabel=r"Magnetisation $\langle |\mathbf{m}| \rangle$",
            # title = f"Magnetisation vs Temperature, inflection points",
        )
        i+=1

    saving_path.mkdir(parents = True, exist_ok = True)
    fig.savefig(saving_path / f"{directory.name}_magnetisation_inflection_point.pdf", bbox_inches="tight")
    print("finished magnetisation inflection plots")

def do_inflection_vs_length_plot(
        directory:pathlib.Path, 
        is_deep:bool = False,
        start:int = 0,
        omit_last:int = 0,
    ):

    saving_path = directory.parent.parent / "analyze" / "output"/"img_dump"
    csv_file = directory / "ensemble_observables.csv"
    if csv_file.exists():
        df = pd.read_csv(csv_file)
    else:
        get_observables_csv(directory, is_deep, start, False)
        df = pd.read_csv(csv_file)

    df = df[df["energy_value"] != 0.0].reset_index(drop=True)
    if is_deep:
        dim = import_physical_parameter(directory / f"{directory.name}_0" / "parameter-config-0", "dimension")
    else:
        dim = import_physical_parameter(directory / "parameter-config-0", "dimension")

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )

    i = 0
    critical_temp = []
    critical_temp_err = []
    l_values = []
    for l, group in df.groupby("L"):
        group = group.sort_values("temperature")
        xaxis = group["temperature"].to_numpy()
        yaxis = group[f"magnetisation_value"].to_numpy()
        yerr = group[f"magnetisation_error"].to_numpy()

        idx = np.argsort(xaxis)
        xaxis = xaxis[idx]
        yaxis = yaxis[idx]
        yerr = yerr[idx]

        params, pcov = _fit_mag_to_tanh(xaxis, yaxis, yerr)
        perr = np.sqrt(np.diag(pcov))

        l_values.append(l)
        critical_temp.append(params[1])
        critical_temp_err.append(perr[1])



    p0 = [1.0, 0.9]  # Initial guess for slope and intercept
    t_c_array = critical_temp[omit_last:]
    print("t_c values for fit:", t_c_array)
    values_to_fit = 1 / np.log(np.array(l_values))**dim
    popt, pcov = curve_fit(_linear_model, values_to_fit[omit_last:], t_c_array, p0=p0)
    x_fit = np.linspace(-1, max(values_to_fit), 100)
    y_fit = _linear_model(x_fit, *popt)
    ax.errorbar(values_to_fit, critical_temp, yerr=critical_temp_err, fmt=".-", capsize=5)
    ax.plot(x_fit, y_fit, "--", label=f"Linear fit, $k_BT_{{BKT}} = {popt[1]:.3f}({pcov[1,1]:.1g})$")
    ax = _add_format_plot(
        axs = ax,
        xlabel=f"$ 1 / \ln(a / L) ^{dim}$",
        ylabel=r"$k_BT_c(L)$",
        # title = f"Critial temperature",
    )
    ax.set_xlim(-0.01,max(values_to_fit)*1.1)

    with open(saving_path / f"{directory.name}_critical_temps.csv", "w") as f:
        f.write("l_values,critical_temp,critical_temp_err")
        for i, t, terr, l in enumerate(zip(critical_temp,critical_temp_err,l_values):
            f.write(f"{l},{t},{terr}")

    with open(saving_path / f"{directory.name}_t_bkt.txt", "w") as f:
        f.write(f"T_BKT = {popt[1]} with error {pcov[1,1]}")

    saving_path.mkdir(parents = True, exist_ok = True)
    fig.savefig(saving_path / f"{directory.name}_t_vs_l.pdf", bbox_inches="tight")
    print("finished t vs l inflection plots")


def _tanh(t, a, x0, b, c):
    t = np.asarray(t, dtype=float)
    return a * np.tanh((t - x0) / b) + c

def _fit_mag_to_tanh(x, y, yerr = None):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.shape != y.shape:
        raise ValueError("L and y must have the same shape")
    if x.size < 4:
        raise ValueError("Need at least 4 points for a 4-parameter tanh fit")

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    amplitude_guess = 0.5 * (np.max(y) - np.min(y))
    if np.isclose(amplitude_guess, 0.0):
        amplitude_guess = 1.0

    centre_guess = np.median(x)
    scale_guess = max((np.max(x) - np.min(x)) / 4.0, 1e-6)
    offset_guess = np.mean(y)
    initial_guess = [amplitude_guess, centre_guess, scale_guess, offset_guess]

    fit_kwargs = {"p0": initial_guess, "maxfev": 10000}
    if yerr is not None:
        yerr = np.asarray(yerr, dtype=float)
        if yerr.shape != y.shape:
            raise ValueError("yerr must have the same shape as y")
        yerr = yerr[order]
        fit_kwargs["sigma"] = yerr
        fit_kwargs["absolute_sigma"] = True

    params, covariance = curve_fit(_tanh, x, y, **fit_kwargs)
    return params, covariance


def _linear_model(x, m,n):
    return m * x + n
