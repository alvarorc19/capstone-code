import os
import sys
import pathlib
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import toml
import numpy as np
import matplotlib.pyplot as plt
from observables_plots import (
        magnetisation_plot,
        energy_plot,
        susceptibility_plot,
        specific_heat_plot,
        binder_cumulant_plot,
    )

def main():
    # project_name = "temperature50_0-3_1-5_l128_dim2_10-3sweeps"
    project_name = "temperature30_0-3_1-5_l100_dim2_long"
    # parameter_combination = 2
    project_path = pathlib.Path("/home/alvaro/Documents/trinity/year4/capstone/capstone-code/projects") / project_name / "parameter-config-0"
    config = toml.load(project_path / "config.toml")

    do_all_plots(project_path.parent, config, project_name)

def do_all_plots(directory:pathlib.Path, config:dict, project_name:str):
    saving_path = directory.parent.parent / "analyze" / "img_dump"
    print("saving_path", saving_path)
    sim_settings = config["simulation_settings"]
    
    x_axis = "temperature"

    # Magnetisation plot
    fig, ax = plt.subplots(
        nrows = 1,
        ncols = 1,
        sharey = all,
        # figsize=(15,6),
    )
    ax = magnetisation_plot(directory, ax, x_axis)

    ax.set_title("<m> vs T")

    ax.set_ylabel("<m>")
    ax.set_xlabel("T")

    fig.savefig(saving_path / f"{project_name}_magnetisation.pdf")

    # Energy plot
    fig, ax = plt.subplots(
        nrows = 1,
        ncols = 1,
        sharey = all,
        # figsize=(15,6),
    )
    ax = energy_plot(directory, ax, x_axis)

    ax.set_title("<E> vs T")

    ax.set_ylabel("<E>")
    ax.set_xlabel("T")

    fig.savefig(saving_path / f"{project_name}_energy.pdf")
    # Susceptibility plot
    fig, ax = plt.subplots(
        nrows = 1,
        ncols = 1,
        sharey = all,
        # figsize=(15,6),
    )
    ax = susceptibility_plot(directory, ax,x_axis)

    ax.set_title("X vs T")

    ax.set_ylabel("X")
    ax.set_xlabel("T")

    fig.savefig(saving_path / f"{project_name}_susceptibility.pdf")
    # specific heat plot
    fig, ax = plt.subplots(
        nrows = 1,
        ncols = 1,
        sharey = all,
        # figsize=(15,6),
    )
    ax = specific_heat_plot(directory, ax, x_axis)

    ax.set_title("c_v vs t")

    ax.set_ylabel("c_v")
    ax.set_xlabel("t")

    fig.savefig(saving_path / f"{project_name}_specific_heat.pdf")

    # Binder cumulant plot
    fig, ax = plt.subplots(
        nrows = 1,
        ncols = 1,
        sharey = all,
        # figsize=(15,6),
    )
    ax = binder_cumulant_plot(directory, ax, x_axis)

    ax.set_title("Binder cumulant vs T")

    ax.set_ylabel("U_infty")
    ax.set_xlabel("T")

    fig.savefig(saving_path / f"{project_name}_binder_cumulant.pdf")




if __name__=="__main__":
    main()

