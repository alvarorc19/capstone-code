import matplotlib
matplotlib.use('Agg')

import os
import sys
import pathlib
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

print("path: ", os.getcwd())

import matplotlib as mpl
from analyze.lattice_plots import (
    do_lattice_smooth_plot,
    do_lattice_arrow_plot,
    do_lattice_temp_plot,
)

mpl.rcParams['animation.embed_limit'] = 200


def main():
    project_name = "mid_temp_test"
    project_root = pathlib.Path("/home/alvaro/Documents/trinity/year4/capstone/capstone-code/projects") 
    # project_root = pathlib.Path("/home/users/romeroca/capstone-code/projects") 
    parameter_combination = 0
    project_path = project_root/ project_name / f"parameter-config-{parameter_combination}"
    # do_lattice_smooth_plot(project_path,project_name, parameter_combination)
    do_lattice_arrow_plot(project_path,project_name, parameter_combination,5)
    do_lattice_temp_plot(project_path, project_name,5)


if __name__ == "__main__":
    main()
