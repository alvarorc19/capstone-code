from pathlib import Path
import copy
import numpy as np
import toml
import subprocess
import json

project_folder = "./projects"
project_name = "test"

total_sweeps = 1000
L = 100
dim = 2
N = L**dim
temp = 0.2

def get_git_hash():
    git_hash = subprocess.check_output(
        ["git", "describe", "--dirty", "--always", "--tags"],
    ).decode('utf-8').strip()
    
    return git_hash

parameters = {
    "physical_settings" : {
        "temperature" : temp,
        "L":L,
        "dimension" : dim,
        "J" : 1,
        # "potts_q" : 3,
        # "H" : 1.0,
        # "vec_H" : [0.0,0.0,2.0],
    },
    "simulation_settings" : {
        "total_sweeps" : total_sweeps,
        # "record_correlation_length" : False,
        # "record_correlation_function" : False,
        "save_last_state": False,
    },
    "git_hash":get_git_hash(),
    "seed":42
}
base = Path(project_folder) / project_name
base.mkdir(parents=True, exist_ok=False)
print("Base directory: ", base)

seed_number = 0

# Various temperatures
# temp_array = np.linspace(0.3,1.5, 20).tolist()
temp_array = []

# Varius sizes
# length_array = np.linspace(10, 1000, 5, dtype=int).tolist()
length_array = []

global_parameters = copy.deepcopy(parameters)

if (temp_array and not length_array):
    global_parameters["physical_settings"]["temperature"] = temp_array
    for i,t in enumerate(temp_array):
        parameters["physical_settings"]["temperature"] = t

        config_path = base / f"parameter-config-{i}"
        config_path.mkdir(exist_ok=True)
        with open(config_path/ "config.toml","w") as f:
            toml.dump(parameters, f)

elif (not temp_array and length_array):
    global_parameters["physical_settings"]["L"] =length_array
    for i,l in enumerate(length_array):
        parameters["physical_settings"]["L"] = l

        config_path = base / f"parameter-config-{i}"
        config_path.mkdir(exist_ok=True)
        with open(config_path/ "config.toml","w") as f:
            toml.dump(parameters, f)

elif (temp_array and length_array):
    global_parameters["physical_settings"]["L"] =length_array
    global_parameters["physical_settings"]["temperature"] = temp_array
    for i,t in enumerate(temp_array):
        for j, l in enumerate(length_array):
            parameters["physical_settings"]["temperature"] = t
            parameters["physical_settings"]["L"] = l

            config_path = base / f"parameter-config-{i*len(temp_array) + j}"
            config_path.mkdir(exist_ok=True)
            with open(config_path/ "config.toml","w") as f:
                toml.dump(parameters, f)

else:
    config_path = base / f"parameter-config-0"
    config_path.mkdir(exist_ok=True)
    with open(config_path / "config.toml", "w") as f:
        toml.dump(parameters, f)

print("parameters: ",global_parameters)
with open(base / "global_parameters.json", "w") as f:
    json.dump(global_parameters, f, indent = 4)
