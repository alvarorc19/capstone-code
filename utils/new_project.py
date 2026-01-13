from pathlib import Path
import toml
import subprocess

project_folder = "./projects"
project_name = "demo2"

def get_git_hash():
    git_hash = subprocess.check_output(
        ["git", "describe", "--dirty", "--always", "--tags"],
    ).decode('utf-8').strip()
    
    return git_hash

parameters = {
    "physical_settings" : {
        "temperature" : 1,
        "L":10,
        "dimension" : 2,
        "J" : 1,
        "potts_q" : 3,
        "H" : 1.0,
        "vec_H" : [0.0,1.0,2.0],
    },
    "simulation_settings" : {
        "time_steps" : 10000,
        "recording_steps" : 1000,
        "record_magnetisation" : True,
        "record_energy" : True,
        "record_susceptibility" : False,
        "record_specific_heat" : False,
        "record_correlation_length" : False,
        "record_correlation_function" : False
    },
    "git_hash":get_git_hash(),
    "seed":42
}
base = Path(project_folder) / project_name
base.mkdir(parents=True, exist_ok=True)
print("Base directory: ", base)

seed_number = 0

# Various temperatures
temp_array = [2,3,4]

# Varius sizes
length_array = [10,100,1000]

if temp_array and not length_array:
    for i,t in enumerate(temp_array):
        parameters["physical_settings"]["temperature"] = t

        config_path = base / f"parameter-config-{i}"
        config_path.mkdir(exist_ok=True)
        with open(config_path/ "config.toml","w") as f:
            toml.dump(parameters, f)

elif not temp_array and length_array:
    for i,t in enumerate(temp_array):
        parameters["physical_settings"]["temperature"] = t

        config_path = base / f"parameter-config-{i}"
        config_path.mkdir(exist_ok=True)
        with open(config_path/ "config.toml","w") as f:
            toml.dump(parameters, f)

elif temp_array and length_array:
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

