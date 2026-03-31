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
    do_observable_plot,
    do_order_parameter_plot,
    get_observables_csv,
    do_magnetisation_inflection_plot,
)
from finite_size_plots import do_finite_size_analysis_susceptibility
from rg_plots import(
    do_renormalisation_plot,
    do_biggest_L_renormalisation_plot,
)

from generate_critical_exponents import (
    generate_critical_exponents_file,
)

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size':14, 'figure.autolayout':True})
# plt.rcParams.update({'font.size':14})

def main():
    # project_name = "temperature50_0-3_1-5_l128_dim2_10-3sweeps"
    project_name = "20260330_30t0_8-1_2_4l16-40_dim2_10-3sweeps"
    is_deep = True
    rg = False
    start_step = 200
    # parameter_combination = 2
    # project_root = pathlib.Path("/home/alvaro/Documents/trinity/year4/capstone/capstone-code/projects")
    project_root = pathlib.Path("/home/users/romeroca/capstone-code/projects")
    project_path = project_root / project_name
    # config = toml.load(project_path / "config.toml")
    observables = [
        "magnetisation",
        "energy",
        "specific_heat",
        "susceptibility",
        # "binder_cumulant",
        "energy_per_spin",
        "susceptibility_per_spin",
        "specific_heat_per_spin",
        "cluster_susceptibility",
        "correlation_length",
    ]

    observables_titles = [
        r"Magnetisation $\langle |m| \rangle$",
        r"Energy $\langle E \rangle$",
        "Specific Heat $C_H$",
        r"Magnetic Susceptibility $\chi$",
        # r"Binder cumulant $U_\infty$",
        r"Energy per spin $e$",
        r"Magnetic Susceptibility per spin $\chi / N$",
        r"Specific Heat per spin $c_H$",
        r"Susceptibility $\chi$",
        r"Correlation length $\xi = \langle n \rangle$"
    ]

    do_order_parameter_plot(project_path, is_deep,0)
    # do_magnetisation_inflection_plot(project_path, is_deep, start_step)
    # do_renormalisation_plot(project_path, is_deep, start_step)
    # get_observables_csv(project_path, is_deep, start_step, rg)
    # do_biggest_L_renormalisation_plot(project_path, is_deep, start_step)

    # Create plots
    # for observable, observables_title in zip(observables, observables_titles):
    #     do_observable_plot(
    #         observable = observable,
    #         observable_title = observables_title,
    #         directory = project_path, 
    #         is_deep = is_deep,
    #         x_data = "temperature",
    #         log_plot = False,
    #         log_fit = False,
    #         linear_fit = False,
    #         start = start_step,
    #     )


    # do_finite_size_analysis_susceptibility(project_path, is_deep, start_step)

    # # Critical exponents
    # saving_path = project_path.parent.parent / "analyze" /"output"/ "critical_exponents"
    #
    # saving_path.mkdir(parents = True, exist_ok = True)
    # critical_temp = 0.89
    # generate_critical_exponents_file(project_path, saving_path, critical_temp)

if __name__=="__main__":
    main()

