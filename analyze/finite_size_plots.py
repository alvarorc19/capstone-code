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
    compute_reduced_temperature,
    compute_susceptibility_scaling_function,
)
from utils.h5_utils import (
    import_physical_parameter,
)
from observables_plots import _add_format_plot

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
