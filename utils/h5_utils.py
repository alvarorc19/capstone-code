import h5py
import numpy as np
import pandas as pd
import toml
import pathlib

# Import observable arrays
def import_observable(directory: pathlib.Path, observable_name: str)->int | float:
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
    size = f["lattice"].shape[0]
    f.close()
    return size
