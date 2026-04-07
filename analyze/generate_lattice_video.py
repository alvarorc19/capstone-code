import matplotlib
matplotlib.use('Agg')

import os
import sys
import pathlib
import matplotlib.pyplot as plt
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

print("path: ", os.getcwd())

import matplotlib as mpl
from analyze.lattice_plots import (
    do_lattice_smooth_plot,
    do_lattice_arrow_plot,
    do_lattice_temp_plot_smooth,
    do_lattice_temp_plot_arrows,
)

mpl.rcParams['animation.embed_limit'] = 200
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size':14, 'figure.autolayout':True})


def main():
    project_name = "highish_temp"
    project_root = pathlib.Path("/home/alvaro/Documents/trinity/year4/capstone/capstone-code/projects") 
    # project_root = pathlib.Path("/home/users/romeroca/capstone-code/projects") 
    parameter_combination = 0
    project_path = project_root/ project_name / f"parameter-config-{parameter_combination}"
    fps = 5
    frames_per_iter = 5
    # do_lattice_smooth_plot(project_path,project_name, parameter_combination)
    do_lattice_arrow_plot(project_path,project_name, parameter_combination,5)
    # do_lattice_temp_plot_smooth(project_root, project_name, fps, frames_per_iter)
    # do_lattice_temp_plot_arrows(project_root, project_name, fps, frames_per_iter)


if __name__ == "__main__":
    main()
