import os
import sys
import pathlib
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import toml
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json

from analyze.compute_observables import (
    compute_alpha_critical_exponent,
    compute_beta_critical_exponent,
    compute_gamma_critical_exponent,
)

def generate_critical_exponents_file(directory:pathlib.Path, saving_path:pathlib.Path, critical_temp:float):
    alpha = compute_alpha_critical_exponent(directory, critical_temp)
    beta = compute_beta_critical_exponent(directory, critical_temp)
    gamma = compute_gamma_critical_exponent(directory, critical_temp)
    with open(directory / "global_parameters.json", "r") as f:
        global_config = json.load(f)


    with open(saving_path / f"{directory.name}_critical_exponents.txt", "w") as f:
        f.write(f"Critical exponents for the project {directory.name} ")
        f.write("\n" * 3)
        f.write("#" + "="*30 + "\n\n")
        f.write("Characteristics of system:\n")

        temps_import = global_config["physical_settings"]["temperature"]
        temps = temps_import if isinstance(temps_import, list) else [temps_import]
        formatted_temps = ", ".join([f"{t:.2f}" for t in temps])
        f.write(f"Temperature(s) = {formatted_temps}\n")
        
        lengths_import = global_config["physical_settings"]["L"]
        lengths = lengths_import if isinstance(lengths_import, list) else [lengths_import]
        formatted_lengths = ", ".join([f"{L:.2f}" for L in lengths])
        f.write(f"Length(s) = {formatted_lengths}\n")
        f.write(f"Dimensions = {global_config['physical_settings']['dimension']}\n")
        f.write(f"J = {global_config['physical_settings']['J']}\n\n")


        f.write(f"Total sweeps = {global_config["simulation_settings"]["total_sweeps"]}")
        f.write(f"Recording sweeps = {global_config["simulation_settings"]["recording_sweeps"]}")
        f.write(f"Git hash = {global_config['git_hash']}\n")
        f.write(f"Seed = {global_config['seed']}\n")

        f.write("\n" * 2)
        f.write("#" + "="*30 + "\n\n")
        f.write("Critical exponents :\n")
        f.write(f"Alpha = {alpha}\n")
        f.write(f"Beta = {beta}\n")
        f.write(f"Gamma = {gamma}\n")


