# Characterisation of the XY model

## How to use

### Installation NixOS
If you are using a NixOS systems you just need to execute
```nix
nix develop .#exec -c zsh
```
This will drop you in a nix shell with all the dependencies necessary to develop and run the code!
You can skip


### Installation non-NixOS

#### Requirements

To run the C++ simulation you will need the following requirements:
- gcc or clang with 
- cmake
- hdf5

To generate the configuration you will need:
- doxygen
- Graphviz

To generate mp4 videos you will need:
- ffpmeg

In order to generate results Python is used and the whole system is built with the **pixi** system in mind, however if not in availability, the python packages to be installed using pip are:
- Jupyter
- Numpy
- Matplotlib
- Pandas
- toml
- pyerrors
- opencv

If you are using a conda environment the hdf5 library can be installed from there similarly to the case done with pixi.

### Running the programme
Once you have all the dependencies installed you shall create a system configuration. To do this you copy and rename `new_project_sample.py`to `new_project.py` and fill all of your desired settings. Running the programme,
```bash
$ python3 utils/new_project.py
```
creates a folder under `projects` with the name given and the parameters specified.

#### Build and run
create and go to the build directory, or run this for short:
```bash
$ mkdir very_long_folder_name && cd $_
```
Then you shall follow the usual cmake commands
```bash
$ cmake .. && make
```

Once you have created the executable, you run it with the following command:
```bash
$ ./run_simulation -p PROJECT_FOLDER -m MODEL -j JOBS
```
The parameters that one can specify are:
- `-p`: This is the project where the programme is going to get the information. It gets as path the root directory so you just need to specify the project folder and the name of project. For example: `projects/demo`.
- `-m`: This gets the model that the system is going to use to perform the Monte Carlo simulations. The options implemented are:
    - `xy`: The XY model
    - `potts`: The Potts model
    - `ising`: The Ising model
- `-j`: Indicates the number of cores to be used in the simulation, managed by OpenMP.

## Development guide(TODO)
