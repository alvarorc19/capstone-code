from pathlib import Path
import copy
import numpy as np
import toml
import subprocess
import json
from itertools import product

project_folder = "./projects"
project_name = "direc_test"

total_sweeps = 1000
L = 100
dim = 2
N = L**dim
temp = 0.2
seed_number = 0

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
        "rg_method": False,
    },
    "git_hash":get_git_hash(),
    "seed":seed_number
}
base = Path(project_folder) / project_name
base.mkdir(parents=True, exist_ok=False)
print("Base directory: ", base)


# Various temperatures
temp_array = np.linspace(0.3,1.5, 20).tolist()
# temp_array = []

# Varius sizes
length_array = np.linspace(10, 1000, 5, dtype=int).tolist()
# length_array = []

global_parameters = copy.deepcopy(parameters)


if (temp_array and not length_array):
    global_parameters["physical_settings"]["temperature"] = temp_array
    total_parameters = len(temp_array)
    num_sub_direc = int(total_parameters / 16)
    remaining = total_parameters % 16

    for i in range(num_sub_direc):
        direc = base / f"{project_name}_{i}"
        direc.mkdir(parents = True, exist_ok = False)
        for j in range(16):
            t = temp_array[16* i + j].item()
            parameters["physical_settings"]["temperature"] = t

            config_path = direc / f"parameter-config-{j}"
            config_path.mkdir(exist_ok=True)
            with open(config_path/ "config.toml","w") as f:
                toml.dump(parameters, f)

    direc = base / f"{project_name}_last"
    direc.mkdir(parents = True, exist_ok = False)
    for k in range(remaining):
        t = temp_array[num_sub_direc + k].item()
        parameters["physical_settings"]["temperature"] = t

        config_path = direc / f"parameter-config-{k}"
        config_path.mkdir(exist_ok=True)
        with open(config_path/ "config.toml","w") as f:
            toml.dump(parameters, f)

elif (not temp_array and length_array):
    global_parameters["physical_settings"]["L"] =length_array
    total_parameters = len(length_array)
    num_sub_direc = int(total_parameters / 16)
    remaining = total_parameters % 16


    for i in range(num_sub_direc):
        direc = base / f"{project_name}_{i}"
        direc.mkdir(parents = True, exist_ok = False)
        for j in range(16):
            l = length_array[16* i + j].item()
            parameters["physical_settings"]["L"] = l

            config_path = direc / f"parameter-config-{j}"
            config_path.mkdir(exist_ok=True)
            with open(config_path/ "config.toml","w") as f:
                toml.dump(parameters, f)

    direc = base / f"{project_name}_last"
    direc.mkdir(parents = True, exist_ok = False)
    for k in range(remaining):
        l = length_array[num_sub_direc + k].item()
        parameters["physical_settings"]["L"] = l

        config_path = direc / f"parameter-config-{k}"
        config_path.mkdir(exist_ok=True)
        with open(config_path/ "config.toml","w") as f:
            toml.dump(parameters, f)


elif (temp_array and length_array):
    total_parameters = len(length_array) * len(temp_array)
    num_sub_direc = int(total_parameters / 16)
    remaining = total_parameters % 16
    combinations = np.array(list(product(temp_array, length_array)))

    global_parameters["physical_settings"]["L"] =length_array
    global_parameters["physical_settings"]["temperature"] = temp_array

    for i in range(num_sub_direc):
        direc = base / f"{project_name}_{i}"
        direc.mkdir(parents = True, exist_ok = False)
        for j in range(16):
            t, l = combinations[16 * i + j]
            parameters["physical_settings"]["L"] = l.item()
            parameters["physical_settings"]["temperature"] = t.item()

            config_path = direc / f"parameter-config-{j}"
            config_path.mkdir(exist_ok=True)
            with open(config_path/ "config.toml","w") as f:
                toml.dump(parameters, f)

    direc = base / f"{project_name}_last"
    direc.mkdir(parents = True, exist_ok = False)
    for k in range(remaining):
        t, l = combinations[num_sub_direc + k]
        parameters["physical_settings"]["L"] = l.item()
        parameters["physical_settings"]["temperature"] = t.item()

        config_path = direc / f"parameter-config-{k}"
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
