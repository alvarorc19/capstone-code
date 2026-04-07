import h5py
import numpy as np
import pandas as pd
import toml
import pathlib

# Import observable arrays
def import_observable(directory: pathlib.Path, observable_name: str)->np.ndarray:
    f = h5py.File(directory / "results.h5", "r")    
    observable = np.array(f["observables"][observable_name])
    f.close()
    return observable

def import_physical_parameter(directory: pathlib.Path, physical_parameter:str) -> int | float:
    config = toml.load(directory / "config.toml")
    parameter = config["physical_settings"][physical_parameter]
    return parameter

# Import lattice matrix
def import_lattice(directory: pathlib.Path, time_step: int) -> pd.DataFrame:
    f = h5py.File(directory / "results.h5", "r")
    lattice = pd.DataFrame(f["lattice"][time_step])
    f.close()
    return lattice

def import_lattice_size(directory: pathlib.Path) -> int:
    f = h5py.File(directory / "results.h5", "r")
    size = f["lattice"].shape[1]
    f.close()
    return size

def import_renormalised_magnetisation(directory: pathlib.Path, b: int = 4) -> np.ndarray:
    f = h5py.File(directory / "results.h5", "r")
    rg_keys = list(f["renormalisation"].keys())
    array = np.array([])
    for key in rg_keys:
        observable, scale, b_file = key.split("_")
        if scale != "b":
            exit()
        if (int(b_file) == b and observable == "magnetisation"):
            array = np.array(list(f["renormalisation"][key]))
    f.close()

    return array

def import_renormalised_energy(directory: pathlib.Path, b: int = 4) -> np.ndarray:
    f = h5py.File(directory / "results.h5", "r")
    rg_keys = list(f["renormalisation"].keys())
    for key in rg_keys:
        observable, scale, b_file = key.split("_")
        if scale != "b":
            exit()
        if (int(b_file) == b and observable == "energy"):
            array = np.array(list(f["renormalisation"][key]))
    f.close()

    return array
