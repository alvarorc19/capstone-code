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
from scipy.optimize import curve_fit

from compute_observables import (
    compute_reduced_temperature,
    compute_susceptibility_scaling_function,
    compute_binder_cumulant,
)
from utils.h5_utils import (
    import_physical_parameter,
)
from observables_plots import (
    get_observables_csv,
    _add_format_plot,
    _add_scatter_data,
    _linear_model,
    _obtain_numbers_format,
)

def do_finite_size_analysis_susceptibility(directory:pathlib.Path, is_deep:bool = False, start: int = 0):
    observable = "susceptibility"
    saving_path = directory.parent.parent / "analyze" / "output"/"vid_dump"
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
    observable_obs_array = np.array([])
    observable_array = np.array([])
    observable_error = np.array([])

    plt.tight_layout()

    compute_observable = _find_observable_function(observable)

    with open(directory / "global_parameters.json", "r") as f:
        global_config = json.load(f)
    unique_lengths = global_config["physical_settings"]["L"]
    unique_temp = global_config["physical_settings"]["temperature"]


    for direc in params:
        print("direc", direc)
        temp_array = np.append(temp_array, import_physical_parameter(direc, "temperature"))
        length_array = np.append(length_array, import_physical_parameter(direc, "L")) 
        observable_obs = compute_observable(direc, start)
        observable_obs_array = np.append(observable_obs_array, observable_obs)
        observable_array = np.append(observable_array, observable_obs.value)
        observable_error = np.append(observable_error, observable_obs.dvalue)


    length_array,temp_array, observable_obs_array = map(np.array,zip(*sorted(zip(length_array, temp_array, observable_obs_array))))

    cmap = plt.cm.tab20
    colors = cmap(np.arange(len(unique_lengths)*2))

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )
    fig.subplots_adjust(bottom=0.25)
    
    # plot it first
    initial_tc = 0.8
    initial_gamma = 1.
    initial_nu = 1.

    for i,(obs,L) in enumerate(zip(observable_obs_array, length_array)):
        reduced_susc_obs = compute_susceptibility_scaling_function(obs, L, initial_gamma, initial_nu)
        observable_array[i] = reduced_susc_obs.value
        observable_error[i] = reduced_susc_obs.dvalue

    reduced_temp = compute_reduced_temperature(temp_array, initial_tc)
    
    scatter_array, ax = _make_finite_size_plot(ax, reduced_temp, unique_lengths, unique_temp, observable_array, observable_error, colors)

    ax = _add_format_plot(
        axs = ax,
        xlabel="reduced temperature $t$",
        ylabel=r"Susceptibility scaling function $\tilde{\chi}$",
        title = r"Scaled susceptibilty $\tilde{\chi}$ vs reduced temperature $t$",
    )

    # create sliders
    ax_slider1 = plt.axes([0.2, 0.12, 0.6, 0.03])
    gamma_slider = Slider(ax_slider1, r'$\gamma$ exponent', 0.1, 5.0, valinit=0.1, valstep=0.1)
    ax_slider2 = plt.axes([0.2, 0.08, 0.6, 0.03])
    nu_slider = Slider(ax_slider2, r'$\nu$ exponent', 0.1, 5.0, valinit=0.1, valstep=0.1)
    ax_slider3 = plt.axes([0.2, 0.04, 0.6, 0.03])
    critical_temp_slider = Slider(ax_slider3, 'Critical temperature', 0.1, 1.0, valinit=0.1, valstep=0.1)

    # create update function
    def update(val):
        for i,(obs,L) in enumerate(zip(observable_obs_array, length_array)):
            reduced_susc_obs = compute_susceptibility_scaling_function(obs, L, gamma_slider.val, nu_slider.val)
            observable_array[i] = reduced_susc_obs.value
            observable_error[i] = reduced_susc_obs.dvalue

        reduced_temp = compute_reduced_temperature(temp_array, critical_temp_slider.val)

        for i,plot in enumerate(scatter_array):
            plotline, caplines, barlinecols = plot

            x_data = reduced_temp[i * len(unique_temp): (i+1) * len(unique_temp)]
            y_data =observable_array[i * len(unique_temp): (i+1)*len(unique_temp)] 
            y_err = observable_error[i * len(unique_temp): (i+1)*len(unique_temp)]

            # value to update
            plotline.set_ydata(y_data)
            plotline.set_xdata(x_data)

            # Update error bar segments
            segments = [[[xi, yi - ei], [xi, yi + ei]] for xi, yi, ei in zip(x_data, y_data, y_err)]
            barlinecols[0].set_segments(segments)

            if caplines:
                caplines[0].set_xdata(x_data)
                caplines[0].set_ydata(y_data - y_err)
                caplines[1].set_xdata(x_data)
                caplines[1].set_ydata(y_data + y_err)

        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw_idle()

    gamma_slider.on_changed(update)
    nu_slider.on_changed(update)
    critical_temp_slider.on_changed(update)

    saving_path.mkdir(parents = True, exist_ok = True)

    # Iterate over everything
    gamma_values = np.arange(0.7, 3.0, 0.1)
    nu_values = np.arange(0.7, 2.4, 0.1)
    tc_values = np.arange(0.4, 1.0, 0.1)

    # If you want all combinations, use itertools
    frames = list(itertools.product(gamma_values, nu_values, tc_values))

    def animate(frame):
        gamma, nu, tc = frame
        gamma_slider.set_val(gamma)
        nu_slider.set_val(nu)
        critical_temp_slider.set_val(tc)
        update(None)

    ani = FuncAnimation(fig, animate, frames=frames, interval=100)
    # ani.save(saving_path / f"{directory.name}_finite_size_susceptibility.mp4", writer='ffmpeg', fps=10)

    with tqdm(total=len(frames), desc="Saving animation") as pbar:
        ani.save(
            saving_path / f"{directory.name}_finite_size_susceptibility.mp4",
            writer='ffmpeg',
            fps=10,
            extra_args=['-vcodec', 'libx264', '-preset', 'fast'],
            progress_callback=lambda i, n: pbar.update(1)
        )

    # plt.show()

def _make_finite_size_plot(
        ax,
        reduced_temp,
        unique_lengths,
        unique_temp,
        scaled_susceptibility,
        scaled_susceptibility_error,
        colors,
    ):
    scatter_array=[]

    for i,l in enumerate(unique_lengths):
        xaxis = reduced_temp[i*len(unique_temp): (i+1)*len(unique_temp)]
        yaxis = scaled_susceptibility[i * len(unique_temp):(i+1)*len(unique_temp)]

        assert len(xaxis) == len(yaxis), "X and Y data need to be the same length"

        yerr = scaled_susceptibility_error[i * len(unique_temp):(i+1)*len(unique_temp)]
        error_plot = ax.errorbar(xaxis, yaxis, yerr = yerr, label = f"$L = {l}$", color = colors[2*i], fmt = ".")
        scatter_array.append(error_plot)

    return scatter_array, ax

def compute_critical_temp_binder(
        directory:pathlib.Path, 
        is_deep:bool = False,
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
    sigmoids_params = []
    sigmoids_err = []

    i = 0
    for l, group in df.groupby("L"):
        group = group.sort_values("temperature")
        xaxis = group["temperature"].to_numpy()
        yaxis = group[f"binder_cumulant_value"].to_numpy()
        yerr = group[f"binder_cumulant_error"].to_numpy()

        idx = np.argsort(xaxis)
        xaxis = xaxis[idx]
        yaxis = yaxis[idx]
        yerr = yerr[idx]

        params, pcov = _fit_to_sigmoid(xaxis, yaxis, yerr)
        temperature_fit = np.linspace(np.min(xaxis), np.max(xaxis), 300)
        ax.plot(temperature_fit, _sigmoid(temperature_fit, *params),"--", color = colors[i])
        sigmoids_params.append(params)
        sigmoids_err.append(pcov)

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
        i+=1

    ax = _add_format_plot(
        axs = ax,
        xlabel="Temperature, $k_BT$",
        ylabel=r"Binder cumulant $U_L$",
        # title = f"Magnetisation vs Temperature, inflection points",
    )

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True,ncol=5)

    crossing_points = _extract_crossing_points(xaxis, sigmoids_params)
    print("crossing_points = ", crossing_points)
    print("t_bkt = ", np.average(crossing_points))
    print("err t_bkt = ", np.std(crossing_points))

    with open(saving_path / f"{directory.name}_cumulant.txt", "w") as f:
        f.write("Crossing points\n[")
        for point in crossing_points:
            f.write(f"{point}")
        f.write("]\n\n\n")
        f.write(f"T_BKT = {np.average(crossing_points)} pm {np.std(crossing_points)}")

    saving_path.mkdir(parents = True, exist_ok = True)
    fig.savefig(saving_path / f"{directory.name}_cumulant_crossing.pdf", bbox_inches="tight")
    print("finished magnetisation inflection plots")

def _fit_to_sigmoid(x, y, yerr = None):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.shape != y.shape:
        raise ValueError("L and y must have the same shape")
    if x.size < 4:
        raise ValueError("Need at least 4 points for a 4-parameter tanh fit")

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    amplitude_guess = np.max(y) - np.min(y)
    centre_guess = np.median(x)
    slope_guess = 1.0
    offset_guess = np.min(y)
    initial_guess = [amplitude_guess, centre_guess, slope_guess, offset_guess]

    fit_kwargs = {"p0": initial_guess, "maxfev": 100000}
    if yerr is not None:
        yerr = np.asarray(yerr, dtype=float)
        if yerr.shape != y.shape:
            raise ValueError("yerr must have the same shape as y")
        yerr = yerr[order]
        fit_kwargs["sigma"] = yerr
        fit_kwargs["absolute_sigma"] = True

    params, covariance = curve_fit(_sigmoid, x, y, **fit_kwargs)
    return params, covariance

def _sigmoid(x, a, b, c, d):
    denom = 1 + np.exp(-c*(x-b))
    return a / denom + d

def _extract_crossing_points(x_values,y_params):
    y_params = np.asarray(y_params)
    x_values = np.linspace(np.min(x_values), np.max(x_values), 500)
    crossings = []

    for i in range(len(y_params) -1):
        params1 = y_params[i]
        params2 = y_params[i+1]

        curve1 = _sigmoid(x_values, *params1)
        curve2 = _sigmoid(x_values, *params2)

        diff = curve2 - curve1
        crossing_idx = np.where(np.diff(np.sign(diff)) != 0)[0]
        for idx in crossing_idx:
            crossings.append(x_values[idx])

    return crossings

def do_susceptibility_vs_length_plot(
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
        dim = import_physical_parameter(directory / f"{directory.name}_last" / "parameter-config-0", "dimension")
    else:
        dim = import_physical_parameter(directory / "parameter-config-0", "dimension")

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )
    cmap = plt.cm.tab20
    colors = cmap(np.arange(20))

    xaxis = df["L"].to_numpy()
    yaxis = df[f"cluster_susceptibility_value"].to_numpy()
    yerr = df[f"cluster_susceptibility_error"].to_numpy()

    idx = np.argsort(xaxis)
    xaxis = xaxis[idx]
    yaxis = yaxis[idx]
    yerr = yerr[idx]

    p0 = [1.0, 0.9]  # Initial guess for slope and intercept
    log_x = np.log(xaxis)
    log_y = np.log(yaxis)
    log_y_err = yerr / yaxis
    popt, pcov = curve_fit(_linear_model, log_x, log_y,sigma = log_y_err, p0=p0)
    perr = np.sqrt(np.diag(pcov))
    error = perr[0]
    f, err = _obtain_numbers_format(error)
    x_fit = np.linspace(min(xaxis), max(xaxis), 300)
    y_fit = np.exp(popt[1]) * x_fit**popt[0]
    ax.plot(x_fit, y_fit, "--", label=f"Linear fit, $\gamma / \\nu = {popt[0]:.{int(f)}f}({int(err)})$", color = colors[1])

    ax = _add_scatter_data(
        axs = ax,
        xaxis = xaxis,
        yaxis = yaxis,
        yerr = yerr,
        data_label = f"$k_B T_c = {df.at[0,'temperature']}$",
        main_color = colors[0],
        secondary_color = colors[1],
        marker = "."
    )
    ax.loglog()

    ax = _add_format_plot(
        axs = ax,
        xlabel = "Length $L/a$",
        ylabel = "Susceptibility $\\chi$",
    )
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True)
    saving_path.mkdir(parents = True, exist_ok = True)
    fig.savefig(saving_path / f"{directory.name}_gamma_nu_critical.pdf", bbox_inches="tight")
    print("finished susceptibility for critical exponents plot")

def do_correlation_length_vs_length_plot(
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
        dim = import_physical_parameter(directory / f"{directory.name}_last" / "parameter-config-0", "dimension")
    else:
        dim = import_physical_parameter(directory / "parameter-config-0", "dimension")

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )
    cmap = plt.cm.tab20
    colors = cmap(np.arange(20))

    xaxis = df["L"].to_numpy()
    yaxis = df[f"correlation_length_per_spin_value"].to_numpy()
    yerr = df[f"correlation_length_per_spin_error"].to_numpy()

    idx = np.argsort(xaxis)
    xaxis = xaxis[idx]
    yaxis = yaxis[idx]
    yerr = yerr[idx]

    p0 = [1.0, 0.9]  # Initial guess for slope and intercept
    log_x = np.log(xaxis)
    log_y = np.log(yaxis)
    log_y_err = np.log(yerr)
    popt, pcov = curve_fit(_linear_model, log_x, log_y, p0=p0)
    perr = np.sqrt(np.diag(pcov))
    error = perr[0]
    f, err = _obtain_numbers_format(error)
    x_fit = np.linspace(min(xaxis), max(xaxis), 300)
    y_fit = np.exp(popt[1]) * x_fit**popt[0]
    ax.plot(x_fit, y_fit, "--", label=f"Linear fit, $\\nu = {-popt[0]:.{int(f)}f}({int(err)})$", color = colors[1])

    ax = _add_scatter_data(
        axs = ax,
        xaxis = xaxis,
        yaxis = yaxis,
        yerr = yerr,
        data_label = f"$k_B T_c = {df.at[0,'temperature']}$",
        main_color = colors[0],
        secondary_color = colors[1],
        marker = "."
    )
    ax.loglog()

    ax = _add_format_plot(
        axs = ax,
        xlabel = "Length $L/a$",
        ylabel = "Correlation length $\\xi$"
    )
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True)
    saving_path.mkdir(parents = True, exist_ok = True)
    fig.savefig(saving_path / f"{directory.name}_nu_critical.pdf", bbox_inches="tight")
    print("finished correlation length for critical exponents plot")

def do_specific_heat_vs_length_plot(
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
        dim = import_physical_parameter(directory / f"{directory.name}_last" / "parameter-config-0", "dimension")
    else:
        dim = import_physical_parameter(directory / "parameter-config-0", "dimension")

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )
    cmap = plt.cm.tab20
    colors = cmap(np.arange(20))

    xaxis = df["L"].to_numpy()
    yaxis = df[f"specific_heat_per_spin_value"].to_numpy()
    yerr = df[f"specific_heat_per_spin_error"].to_numpy()

    idx = np.argsort(xaxis)
    xaxis = xaxis[idx]
    yaxis = yaxis[idx]
    yerr = yerr[idx]

    p0 = [1.0, 0.9]  # Initial guess for slope and intercept
    log_x = np.log(xaxis)
    log_y = np.log(yaxis)
    log_y_err = np.log(yerr)
    popt, pcov = curve_fit(_linear_model, log_x, log_y, p0=p0)
    perr = np.sqrt(np.diag(pcov))
    error = perr[0]
    f, err = _obtain_numbers_format(error)
    x_fit = np.linspace(min(xaxis), max(xaxis), 300)
    y_fit = np.exp(popt[1]) * x_fit**popt[0]
    ax.plot(x_fit, y_fit, "--", label=f"Linear fit, $\\alpha / \\nu = {popt[0]:.{int(f)}f}({int(err)})$", color = colors[1])

    ax = _add_scatter_data(
        axs = ax,
        xaxis = xaxis,
        yaxis = yaxis,
        yerr = yerr,
        data_label = f"$k_B T_c = {df.at[0,'temperature']}$",
        main_color = colors[0],
        secondary_color = colors[1],
        marker = "."
    )
    ax.loglog()

    ax = _add_format_plot(
        axs = ax,
        xlabel = "Length $L/a$",
        ylabel = "Specific heat $C$",
    )
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True)
    saving_path.mkdir(parents = True, exist_ok = True)
    fig.savefig(saving_path / f"{directory.name}_alpha_nu_critical.pdf", bbox_inches="tight")
    print("finished specific heat for critical exponents plot")

def do_finite_size_analysis_observable(
        directory:pathlib.Path, 
        is_deep:bool = False,
        start:int = 0,
        observable: str = "susceptibility",
        nu:float = 0.5,
        t_c:float = 0.895,
        exponent:float = 1,
    ):
    saving_path = directory.parent.parent / "analyze" / "output"/f"{directory.name}_finite_size"
    csv_file = directory / "ensemble_observables.csv"
    if csv_file.exists():
        df = pd.read_csv(csv_file)
    else:
        get_observables_csv(directory, is_deep, start, False)
        df = pd.read_csv(csv_file)

    df = df[df["energy_value"] != 0.0].reset_index(drop=True)
    if is_deep:
        dim = import_physical_parameter(directory / f"{directory.name}_last" / "parameter-config-0", "dimension")
    else:
        dim = import_physical_parameter(directory / "parameter-config-0", "dimension")

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )
    cmap = plt.cm.tab20
    colors = cmap(np.arange(20))


    i = 0
    for l, group in df.groupby("L"):
        group = group.sort_values("temperature")
        xaxis = group["temperature"].to_numpy()
        obs_array = group[f"{observable}_value"].to_numpy()
        obs_err = group[f"{observable}_error"].to_numpy()

        xaxis = ((xaxis-t_c) / t_c) * l**(1 / nu)
        if observable == "magnetisation":
            yaxis = obs_array * l ** (exponent / nu)
            yerr = obs_err* l ** (exponent / nu)
        else:
            yaxis = obs_array * l ** (-exponent / nu)
            yerr = obs_err * l ** (-exponent / nu)

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

        i+=1


    if observable == "magnetisation":
        exponentstr = r"\beta"
        exponentletter = "beta"
        exponenttitle = r"$\beta$"
        obs_letter = "m"
    elif observable == "specific_heat":
        exponentstr = r"-\alpha"
        exponentletter = "alpha"
        exponenttitle = r"$\alpha$"
        obs_letter = "C"
    elif observable == "specific_heat_per_spin":
        exponentstr = r"-\alpha"
        exponentletter = "alpha"
        exponenttitle = r"$\alpha$"
        obs_letter = "c"
    elif observable == "cluster_susceptibility":
        exponentstr = r"-\gamma"
        exponentletter = "gamma"
        exponenttitle = r"$-\gamma$"
        obs_letter = "\\chi"

    title = f"$\\nu = {nu:.5f}$, {exponenttitle} $= {exponent:.5f}$"
    ax = _add_format_plot(
        axs = ax,
        xlabel=r"$tL^{1 / \nu}$",
        ylabel=f"${obs_letter} L^{{{exponentstr} / \\nu}}$",
        title = title,
    )

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol = 4)
    nu = f"{nu:.4f}".replace(".", "_")
    saving_path = saving_path / observable / f"nu_{nu}"
    saving_path.mkdir(parents = True, exist_ok = True)
    expontent = f"{exponent:.4f}".replace(".", "_")
    fig.savefig(saving_path / f"{exponentletter}_{exponent}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"finished finite size for {observable} for critical exponents plot")

def do_finite_size_analysis_nu(
        directory:pathlib.Path, 
        is_deep:bool = False,
        start:int = 0,
        observable: str = "cluster_susceptibility",
        nu:float = 0.5,
        t_c:float = 0.895,
        gamma_nu:float = 2.0,
    ):
    saving_path = directory.parent.parent / "analyze" / "output"/f"{directory.name}_finite_size"
    csv_file = directory / "ensemble_observables.csv"
    if csv_file.exists():
        df = pd.read_csv(csv_file)
    else:
        get_observables_csv(directory, is_deep, start, False)
        df = pd.read_csv(csv_file)

    df = df[df["energy_value"] != 0.0].reset_index(drop=True)
    if is_deep:
        dim = import_physical_parameter(directory / f"{directory.name}_last" / "parameter-config-0", "dimension")
    else:
        dim = import_physical_parameter(directory / "parameter-config-0", "dimension")
    exponent = gamma_nu * nu

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )
    cmap = plt.cm.tab20
    colors = cmap(np.arange(20))


    i = 0
    for l, group in df.groupby("L"):
        group = group.sort_values("temperature")
        xaxis = group["temperature"].to_numpy()
        obs_array = group[f"{observable}_value"].to_numpy()
        obs_err = group[f"{observable}_error"].to_numpy()

        xaxis = ((xaxis-t_c) / t_c) * l**(1 / nu)
        if observable == "magnetisation":
            yaxis = obs_array * l ** (exponent / nu)
            yerr = obs_err* l ** (exponent / nu)
        else:
            yaxis = obs_array * l ** (-gamma_nu)
            yerr = obs_err * l ** (-gamma_nu)

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

        i+=1


    if observable == "cluster_susceptibility":
        exponentstr = r"-\gamma"
        exponentletter = "gamma"
        exponenttitle = r"$\gamma$"
        obs_letter = "\\chi"

    title = f"$\\nu = {nu:.5f}$, {exponenttitle} $= {exponent:.5f}$"
    ax = _add_format_plot(
        axs = ax,
        xlabel=r"$tL^{1 / \nu}$",
        ylabel=f"${obs_letter} L^{{{exponentstr} / \\nu}}$",
        title = title,
    )

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol = 4)
    nu = f"{nu:.4f}".replace(".", "_")
    saving_path = saving_path / observable 
    saving_path.mkdir(parents = True, exist_ok = True)
    expontent = f"{exponent:.4f}".replace(".", "_")
    fig.savefig(saving_path / f"{exponentletter}_{exponent}_nu_{nu}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"finished finite size for {observable} for critical exponents plot")
