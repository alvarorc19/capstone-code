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
)

from generate_critical_exponents import (
    generate_critical_exponents_file,
)

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size':14, 'figure.autolayout':True})

def main():
    # project_name = "temperature50_0-3_1-5_l128_dim2_10-3sweeps"
    project_name = "temp20_l4_dim2_20-3sweeps"
    # parameter_combination = 2
    # project_root = pathlib.Path("/home/alvaro/Documents/trinity/year4/capstone/capstone-code/projects")
    project_root = pathlib.Path("/home/users/romeroca/capstone-code/projects")
    project_path = project_root / project_name
    # config = toml.load(project_path / "config.toml")
    observables = ["magnetisation", "energy", "specific_heat", "susceptibility", "binder_cumulant","normalised_energy"]
    observables_titles = [r"Magnetisation $\langle |m| \rangle$", r"Energy $\langle E \rangle$", "Specific Heat $C_H$", r"Magnetic Susceptibility $\chi$", r"Binder cumulant $U_\infty$", r"Energy $E/J L^2$"]

    # Create plots
    for observable, observables_title in zip(observables, observables_titles):
        do_observable_plot(
            observable = observable,
            observable_title = observables_title,
            directory = project_path, 
            x_data = "temperature",
            log_plot = False,
            log_fit = False,
            linear_fit = False,
        )

    # Critical exponents
    saving_path = project_path.parent.parent / "analyze" /"output"/ "critical_exponents"

    saving_path.mkdir(parents = True, exist_ok = True)
    critical_temp = 0.89
    generate_critical_exponents_file(project_path, saving_path, critical_temp)

if __name__=="__main__":
    main()

