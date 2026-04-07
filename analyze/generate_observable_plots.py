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
    do_inflection_vs_length_plot,
)
from finite_size_plots import (
    do_finite_size_analysis_susceptibility,
    compute_critical_temp_binder,
    do_susceptibility_vs_length_plot,
    do_specific_heat_vs_length_plot,
    do_correlation_length_vs_length_plot,
    do_finite_size_analysis_observable,
    do_finite_size_analysis_nu,
)
from rg_plots import(
    do_renormalisation_plot,
    do_biggest_L_renormalisation_plot,
)

from generate_critical_exponents import (
    generate_critical_exponents_file,
)
plt.rcParams.update({
    "font.size": 16,        # default text size
    "axes.titlesize": 19,   # title size
    "axes.labelsize": 18,   # axis labels
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 14,
    "figure.autolayout":True,
    "text.usetex":True
})
# plt.rcParams.update({'font.size':14})

def main():
    is_deep = False
    rg = False
    start_step = 200
    # parameter_combination = 2
    project_root = pathlib.Path("/home/alvaro/Documents/trinity/year4/capstone/capstone-code/projects")
    # project_root = pathlib.Path("/home/users/romeroca/capstone-code/projects") # hpc
    project_paths = [
        # HPC
        # project_root / "20260326_15t0_8-1_2_13l32-1024_dim3_10-4sweeps", # it's 2d_run1
        # project_root / "20260330_30t0_8-1_2_4l16-40_dim2_10-3sweeps", # 2d_run2
        # project_root / "20260403_t2-201_10l16-100_dim3_10-4sweeps", # 3d_critical
        # project_root / "20260403_t0-895_10l16-1024_dim2_10-4sweeps" # 3d_critical
        # project_root / "20260401_30t1_9-2_5_5l16-64_dim3_10-3sweeps"
        # # Thinkpad
        # project_root / "20260401_30t2-0_2-5_4l8_40_dim3_rg_10-3sweeps",
        project_root / "20260402_30t2-0_2-2_3l8_40_dim3_20-3sweeps", # 3d_run1
        # project_root / "rg_test2"
    ]

    # config = toml.load(project_path / "config.toml")
    observables = [
        "magnetisation",
        "energy",
        "specific_heat",
        "susceptibility",
        "energy_per_spin",
        "susceptibility_per_spin",
        "specific_heat_per_spin",
        "cluster_susceptibility",
        "correlation_length",
        "correlation_length_per_spin",
        "cluster_susceptibility_per_spin",
        "binder_cumulant",
        # "binder_cumulant_slope",
    ]

    observables_titles = [
        r"Magnetisation $\langle |\textbf{m}| \rangle$",
        r"Energy $\langle E \rangle$",
        "Specific Heat $C/k_B$",
        r"Magnetic Susceptibility $\chi$",
        r"Energy per spin $e$",
        r"Magnetic Susceptibility per spin $\chi / N$",
        r"Specific Heat per spin $c/k_B$",
        r"Susceptibility $\chi$",
        r"Correlation length $\xi$",
        r"Correlation length per spin $\xi /N$",
        r"Susceptibility per spin $\chi / N$",
        r"Binder cumulant $U_L$",
        # r"Binder cumulant slope $\partial U_L / \partial \beta|_{\beta_c}$",
    ]

    omit_values = [0,1,0,0,0,0,0,0,0,0]
    i = 0
    nu = 0.70204
    nus = np.linspace(0.6,0.8,50)
    alphas = np.linspace(0.16,0.23,40)
    alphas2 = np.linspace(-10,10,40)
    betas = np.linspace(0.3,0.4,40)
    gamma_nu = 2.0
    gamma_nu = 1.976
    t_c = 2.201
    for project_path in project_paths:

        # get_observables_csv(project_path, is_deep, start_step, rg)
        # do_susceptibility_vs_length_plot(project_path, is_deep, start_step)
        # do_specific_heat_vs_length_plot(project_path, is_deep, start_step)
        # do_correlation_length_vs_length_plot(project_path, is_deep, start_step)
        for nu in nus:
            do_finite_size_analysis_nu(project_path, is_deep, start_step, "cluster_susceptibility", nu,t_c,gamma_nu)
        # nu = 0.5
        nu = 0.67
        for i in range(len(alphas)):
            do_finite_size_analysis_observable(project_path, is_deep, start_step, "magnetisation", nu, t_c,betas[i])
            do_finite_size_analysis_observable(project_path, is_deep, start_step, "specific_heat_per_spin", nu, t_c,alphas[i])
            do_finite_size_analysis_observable(project_path, is_deep, start_step, "specific_heat", nu, t_c,alphas2[i])
        # do_order_parameter_plot(project_path, is_deep,0)
        # do_magnetisation_inflection_plot(project_path, is_deep, start_step)
        # do_inflection_vs_length_plot(project_path, is_deep, start_step, omit_values[i])
        if rg:
            do_renormalisation_plot(project_path, is_deep, start_step)
        if rg:
            # do_biggest_L_renormalisation_plot(project_path, is_deep, start_step)
            do_renormalisation_plot(project_path, is_deep, start_step)

        # compute_critical_temp_binder(project_path, is_deep, start_step)


        # Create plots
        for observable, observables_title in zip(observables, observables_titles):
            do_observable_plot(
                observable = observable,
                observable_title = observables_title,
                directory = project_path, 
                is_deep = is_deep,
                x_data = "temperature",
                log_plot = False,
                log_fit = False,
                linear_fit = False,
                start = start_step,
            )
        i+=1
    # project_path = project_root / "20260326_15t0_8-1_2_13l32-1024_dim3_10-4sweeps"
    # do_order_parameter_plot(project_path, True,0)


    # do_finite_size_analysis_susceptibility(project_path, is_deep, start_step)

    # # Critical exponents
    # saving_path = project_path.parent.parent / "analyze" /"output"/ "critical_exponents"
    #
    # saving_path.mkdir(parents = True, exist_ok = True)
    # critical_temp = 0.89
    # generate_critical_exponents_file(project_path, saving_path, critical_temp)

if __name__=="__main__":
    main()

