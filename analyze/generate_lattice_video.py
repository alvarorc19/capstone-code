import matplotlib
matplotlib.use('Agg')

import os
import sys
import pathlib
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

print("path: ", os.getcwd())

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
import toml
from utils.h5_utils import (
    import_lattice,
    import_lattice_size,
)


def generate_ising_grid(time_step:int, df:pd.DataFrame)->np.array:
    L = int(np.sqrt(df.shape[1]-1))
    ground_1d_grid = df.loc[time_step].to_numpy()[1:]
    ground_grid = np.array([ground_1d_grid[:L-1]])
    for i in range(1,L):
        a = ground_1d_grid[i*L: i*L + (L-1)]
        ground_grid = np.vstack((ground_grid,a))
    return ground_grid

def generate_lattice(time_step:int, df:pd.DataFrame, config:dict) -> np.array:
    L = config["physical_settings"]["L"]
    dim = config["physical_settings"]["dimension"]
    flat_grid = df.to_numpy().flatten()
    shape = tuple([L] * dim)
    full_dim_grid = flat_grid.reshape(shape)
    return full_dim_grid

def main(project_name: str, parameter_combination: int):
    project_root = pathlib.Path("/home/alvaro/Documents/trinity/year4/capstone/capstone-code/projects") 
    # project_root = pathlib.Path("/home/users/romeroca/capstone-code/projects") 
    project_path = project_root/ project_name / f"parameter-config-{parameter_combination}"
    config = toml.load(project_path / "config.toml")

    fig, ax = plt.subplots()
    ax.set_title(f"XY Model with J = 1 and T = {config['physical_settings']['temperature']:.2f}")
    # # Good cmap for the wrapping of angles
    # im = ax.imshow(generate_lattice(0, import_lattice(project_path,-1), config), cmap="twilight_shifted", vmin = 0, vmax = 2*np.pi)

    lattice_size = import_lattice_size(project_path)
    # Good cmap for visualising vortices
    im = ax.imshow(generate_lattice(0, import_lattice(project_path,0), config), cmap="hsv", vmin = 0, vmax = 2*np.pi, interpolation='lanczos')

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, label='Spin Angle (radians)') 
    cbar.set_ticks([0, np.pi/2, np.pi, (3*np.pi) / 2, 2*np.pi])
    cbar.set_ticklabels(['0', 'π/2', 'π','3π/2',  '2π'])

    # # get indices to not show every change
    # stride = 100
    # frame_indices = range(0, 10000, stride)

    frames_ps = 10
    print("lattice size", lattice_size)

    def update(frame):
        # actual_frame = frame_indices[frame]
        print(f"Time = {frame}")
        im.set_array(generate_lattice(frame, import_lattice(project_path,frame), config))
        return [im]
    ani = FuncAnimation(fig, update, frames = lattice_size, interval = frames_ps, blit = True)
    # html = ani.to_jshtml()
    # with open("analyze/vid_dump/temperature50.html", "w") as f:
    #     f.write(html)
    save_path = project_path.parent.parent.parent / "analyze" / "vid_dump" / f"{project_name}_par_{parameter_combination}_lattice.mp4"
    save_path.parent.mkdir(parents = True, exist_ok = True)
    print("save path", save_path)
    ani.save(save_path, writer = "ffmpeg", fps =frames_ps)


if __name__ == "__main__":
    # project_name = "temperature50_0-3_1-5_l128_dim2_10-3sweeps"
    project_name = "test1000"
    parameter_combination = 0
    main(project_name, parameter_combination)
    # for i in range(30):
    #     main(project_name, i)
