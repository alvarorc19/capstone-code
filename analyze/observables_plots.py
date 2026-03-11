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
    compute_reduced_temperature,
    compute_susceptibility_scaling_function,
    compute_specific_heat_per_spin,
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
    }
    return observable_functions[observable]

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
        print("direc", direc)
        temp_array = np.append(temp_array, import_physical_parameter(direc, "temperature"))
        length_array = np.append(length_array, import_physical_parameter(direc, "L")) 
        observable_obs = compute_observable(direc, start)
        observable_obs_array = np.append(observable_obs_array, observable_obs)
        observable_array = np.append(observable_array, observable_obs.value)
        observable_error = np.append(observable_error, observable_obs.dvalue)

    fig, ax = plt.subplots(
        ncols=1,
        nrows=1
    )

    if x_data == "temperature":
        cmap = plt.cm.tab20

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
                log_fit = log_fit,
                main_color = colors[i],
                secondary_color = colors[i+1],
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

def do_order_parameter_plot(directory:pathlib.Path, is_deep:bool = False, start:int = 0):
    # Add Obs crap

    plt.tight_layout()

    saving_path = directory.parent.parent / "analyze" / "output"/f"thermalisation_{directory.name}"
    saving_path.mkdir(parents = True, exist_ok = True)
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


    for direc in params:
        magnetisation = import_observable(direc, "magnetisation")
        x_magnetisation = np.array([i[0] for i in magnetisation[start:]])
        y_magnetisation = np.array([i[1] for i in magnetisation[start:]])
        length = import_physical_parameter(direc, "L")
        dim = import_physical_parameter(direc, "dimension")
        N = length ** dim
        magnetisation_array = np.sqrt(x_magnetisation**2 + y_magnetisation**2) / N
        energy = import_observable(direc, "energy")

        fig_mag, ax1 = plt.subplots(
            ncols=1,
            nrows=1
        )

        ax1.plot(magnetisation_array, color = "orangered", alpha = 0.5)
        ax1.set_ylabel(r"Magnetisation $m$")
        ax1.set_xlabel(r"Time $t$ in MCS")


        fig_energy, ax2 = plt.subplots(
            ncols=1,
            nrows=1
        )

        ax2.plot(energy[start:], color = "chartreuse", alpha = 0.5)
        ax2.set_ylabel(r"Energy $E$")
        ax2.set_xlabel(r"Time $t$ in MCS")
        plt.show()

        fig_mag.savefig(saving_path / f"{direc.name}_magnetisation_time.pdf", bbox_inches = "tight")
        plt.close()
        fig_energy.savefig(saving_path / f"{direc.name}_energy_time.pdf", bbox_inches = "tight")
        plt.close()
