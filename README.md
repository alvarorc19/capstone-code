# Characterisation of the XY model

## How to use

### Installation NixOS
If you are using a NixOS systems you just need to execute
```zsh
$ nix develop .#exec -c zsh
```
This will drop you in a nix shell with all the dependencies necessary to develop and run the code!
You can skip to [running the simulation](#running-the-programme).


### Installation non-NixOS

#### Requirements

To run the C++ simulation you will need the following requirements:
- gcc or clang with 
- cmake
- hdf5

To generate the documentation you will need:
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
Once you have all the dependencies installed you shall create a system configuration. To do this you copy and rename `new_project_sample.py`to `new_project.py`
```zsh
$ cp utils/new_project_sample.py utils/new_project.py && mkdir projects
```
 and fill all of your desired settings.
```zsh
$ python3 utils/new_project.py
```
creates a folder under `projects` with the name given and the parameters specified in the file.

#### Build and run
create and go to the build directory, or run this for short:
```zsh
$ mkdir build && cd $_
```
Then you shall follow the usual cmake commands
```zsh
$ cmake .. && make
```

Once you have created the executable, you run it with the following command:
```zsh
$ ./run_simulation -p PROJECT_FOLDER -m MODEL -j JOBS
```
The parameters that one can specify are:
- `-p`: This is the project where the programme is going to get the information. It gets as path the root directory so you just need to specify the project folder and the name of project. For example: `projects/demo`.
- `-m`: This gets the model that the system is going to use to perform the Monte Carlo simulations. The options implemented are:
    - `xy`: The XY model
    - `potts`: The Potts model (Might not work)
    - `ising`: The Ising model (Might not work)
- `-j`: Indicates the number of cores to be used in the simulation, managed by OpenMP.

## Analysing the data
Once you have the run done, you can analyse the data using the Python files in the `analyze/` folder. You just need to fill in the data for your project name and other parameters such as your preferred start for the observables.