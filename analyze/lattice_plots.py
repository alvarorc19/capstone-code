import pathlib
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from matplotlib import colors
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
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

def _generate_lattice(time_step:int, df:pd.DataFrame, config:dict) -> np.array:
    L = config["physical_settings"]["L"]
    dim = config["physical_settings"]["dimension"]
    flat_grid = df.to_numpy().flatten()
    shape = tuple([int(L)] * dim)
    full_dim_grid = flat_grid.reshape(shape)
    return full_dim_grid

def _get_arrow_data(frame, lattice, config):
    angles = _generate_lattice(frame,lattice, config)
    y,x = np.indices(angles.shape)

    u = np.cos(angles)
    v = np.sin(angles)

    return x,y,u,v, angles


def do_lattice_smooth_plot(project_path: pathlib.Path, project_name:str, parameter_combination: int, fps:int = 5):
    config = toml.load(project_path / "config.toml")

    fig, ax = plt.subplots()
    ax.set_title(f"XY Model with $J = 1$ and $T = {config['physical_settings']['temperature']:.2f}/k_B$")
    # # Good cmap for the wrapping of angles
    # im = ax.imshow(generate_lattice(0, import_lattice(project_path,-1), config), cmap="twilight_shifted", vmin = 0, vmax = 2*np.pi)

    lattice_size = import_lattice_size(project_path)
    # Good cmap for visualising vortices
    im = ax.imshow(_generate_lattice(0, import_lattice(project_path,0), config), cmap="hsv", vmin = 0, vmax = 2*np.pi, interpolation='bilinear', origin = "lower")

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, label='Spin Angle (radians)') 
    cbar.set_ticks([0, np.pi/2, np.pi, (3*np.pi) / 2, 2*np.pi])
    cbar.set_ticklabels(['0', r'$\pi$/2', r'$\pi$',r'3$\pi$/2',  r'2$\pi$'])

    fig.set_facecolor("#ECEFF4")  
    ax.set_facecolor("#ECEFF4")    

    # # get indices to not show every change
    # stride = 100
    # frame_indices = range(0, 10000, stride)

    nframes = min(100, lattice_size)
    print("lattice size", lattice_size)

    # save_path = project_path.parent.parent.parent / "analyze" /"output"/ "vid_dump" / f"{project_name}_par_{parameter_combination}_lattice.mp4"
    save_path = project_root.parent / "analyze" /"output"/ "vid_dump" 
    save_path.parent.mkdir(parents = True, exist_ok = True)
    save_path_frame = save_path / "frames"
    save_path_frame.mkdir(parents = True, exist_ok = True)
    print("save path", save_path)

    def update(frame):
        # actual_frame = frame_indices[frame]
        # print(f"Time = {frame}")
        im.set_array(_generate_lattice(frame, import_lattice(project_path,frame), config))
        # fig.savefig(save_path_frame / f"frame_{frame:03d}.png", dpi=300)
        return [im]
    ani = FuncAnimation(fig, update, frames = tqdm(range(nframes), desc="Rendering"), interval = fps, blit = True, cache_frame_data = False)
    # html = ani.to_jshtml()
    # with open(save_path / f"{project_name}_par{parameter_combination}_lattice.html", "w") as f:
    #     f.write(html)
    ani.save(save_path / f"{project_name}_par{parameter_combination}_lattice.mp4", writer = "ffmpeg", fps =fps, dpi = 300)


def do_lattice_arrow_plot(project_path: pathlib.Path, project_name:str, parameter_combination: int, fps:int):
    config = toml.load(project_path / "config.toml")

    fig, ax = plt.subplots()
    ax.set_title(f"XY Model with $J = 1$ and $T = {config['physical_settings']['temperature']:.2f}/k_B$")
    # # Good cmap for the wrapping of angles
    # im = ax.imshow(generate_lattice(0, import_lattice(project_path,-1), config), cmap="twilight_shifted", vmin = 0, vmax = 2*np.pi)

    x,y,u,v, angles = _get_arrow_data(0, import_lattice(project_path,0), config)
    step = 1
    norm = colors.Normalize(vmin=0, vmax=2*np.pi)
    q= ax.quiver(
        x[::step, ::step],
        y[::step, ::step],
        u[::step, ::step],
        v[::step, ::step],
        angles[::step, ::step],
        cmap="hsv",
        norm=norm,
        pivot="mid",
        scale=25,
        width=0.003,
        headwidth=3,
        headlength=4,
        headaxislength=3.5,
    )
    # Add colorbar
    cbar = fig.colorbar(q, ax=ax, label='Spin Angle (radians)') 
    cbar.set_ticks([0, np.pi/2, np.pi, (3*np.pi) / 2, 2*np.pi])
    cbar.set_ticklabels(['0', r'$\pi$/2', r'$\pi$',r'3$\pi$/2',  r'2$\pi$'])
    ax.set_aspect('equal')

    fig.set_facecolor("#ECEFF4")  
    ax.set_facecolor("#ECEFF4")    

    # # get indices to not show every change
    # stride = 100
    # frame_indices = range(0, 10000, stride)

    nframes = 100

    # save_path = project_path.parent.parent.parent / "analyze" /"output"/ "vid_dump" / f"{project_name}_par_{parameter_combination}_lattice.mp4"
    save_path = project_path.parent.parent.parent / "analyze" /"output"/ "vid_dump" 
    save_path.parent.mkdir(parents = True, exist_ok = True)
    save_path_frame = save_path / "frames"
    save_path_frame.mkdir(parents = True, exist_ok = True)
    print("save path", save_path)
    plt.tight_layout()

    def update(frame):
        # actual_frame = frame_indices[frame]
        # print(f"Time = {frame}")
        x,y,u,v, angles = _get_arrow_data(frame, import_lattice(project_path,frame), config) 
        q.set_UVC(u[::step, ::step],v[::step, ::step],angles[::step, ::step])
        fig.savefig(save_path_frame / f"frame_{frame:03d}.png", dpi=400, bbox_inches='tight')
        return [q]
    ani = FuncAnimation(fig, update, frames = tqdm(range(nframes), desc="Rendering"), interval = fps, blit = True, cache_frame_data = False)
    # with open(save_path / f"{project_name}_par{parameter_combination}_lattice.html", "w") as f:
    #     f.write(html)
    ani.save(save_path / f"{project_name}_arrow_video.mp4", writer = "ffmpeg", fps =fps, dpi = 300)


def do_lattice_temp_plot_smooth(project_root:pathlib.Path, project_name:str, fps:int = 5, frames_per_iter:int = 10):
    project_path = project_root/ project_name
    lattice = import_lattice(project_path / "parameter-config-0", 0)
    config = toml.load(project_path / "parameter-config-0" / "config.toml")
    frames_per_iter = 4
    num_dir = 0
    for dir in project_path.iterdir():
        if dir.is_dir():
            num_dir += 1

    fig, ax = plt.subplots()

    im = ax.imshow(
        _generate_lattice(0, lattice, config),
          cmap="hsv",
          vmin = 0,
          vmax = 2*np.pi,
          interpolation='bilinear',
          origin = "lower")

    cbar = fig.colorbar(im, ax=ax, label="Spin angle (radians)")
    cbar.set_ticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    cbar.set_ticklabels(["0", r"$\pi$/2", r"$\pi$", r"3$\pi$/2", r"2$\pi$"])

    ax.set_aspect('equal')
    fig.set_facecolor("#ECEFF4")  
    ax.set_facecolor("#ECEFF4")    

    nframes = frames_per_iter * num_dir

    save_path = project_root.parent / "analyze" /"output"/ "vid_dump" 
    save_path.parent.mkdir(parents = True, exist_ok = True)
    save_path_frame = save_path / "frames"
    save_path_frame.mkdir(parents = True, exist_ok = True)
    print("save path", save_path)

    def update(frame):
        current_project = project_path / f"parameter-config-{frame // frames_per_iter}"
        config = toml.load(current_project / "config.toml")
        ax.set_title(f"XY model for $J = 1$, $T = {config['physical_settings']['temperature']} / k_B$")
        im.set_array(_generate_lattice(frame % frames_per_iter, import_lattice(current_project,frame % frames_per_iter), config))
        return [im]
    ani = FuncAnimation(fig, update, frames = tqdm(range(nframes), desc="Rendering"), interval = fps, blit = True, cache_frame_data = False)
    ani.save(save_path / f"{project_name}_temp_video.mp4", writer = "ffmpeg", fps =fps, dpi = 300)


def do_lattice_temp_plot_arrows(project_root:pathlib.Path, project_name:str, fps:int = 5, frames_per_iter:int = 10):
    project_path = project_root/ project_name
    lattice = import_lattice(project_path / "parameter-config-0", 0)
    config = toml.load(project_path / "parameter-config-0" / "config.toml")
    frames_per_iter = 4
    num_dir = 0
    for dir in project_path.iterdir():
        if dir.is_dir():
            num_dir += 1

    fig, ax = plt.subplots()

    x,y,u,v, angles = _get_arrow_data(0, lattice, config)
    step = 6
    norm = colors.Normalize(vmin=0, vmax=2*np.pi)
    q= ax.quiver(
        x[::step, ::step],
        y[::step, ::step],
        u[::step, ::step],
        v[::step, ::step],
        angles[::step, ::step],
        cmap="hsv",
        norm=norm,
        pivot="mid",
        scale=15,
        width=0.003,
        headwidth=5,
        headlength=10,
        headaxislength=3.5,
    )

    cbar = fig.colorbar(q, ax=ax, label="Spin angle (radians)")
    cbar.set_ticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    cbar.set_ticklabels(["0", r"$\pi$/2", r"$\pi$", r"3$\pi$/2", r"2$\pi$"])

    ax.set_aspect('equal')
    # fig.set_facecolor("#ECEFF4")  
    # ax.set_facecolor("#ECEFF4")    

    nframes = frames_per_iter * num_dir

    save_path = project_root.parent / "analyze" /"output"/ "vid_dump" 
    save_path.parent.mkdir(parents = True, exist_ok = True)
    save_path_frame = save_path / "frames"
    save_path_frame.mkdir(parents = True, exist_ok = True)
    print("save path", save_path)

    def update(frame):
        current_project = project_path / f"parameter-config-{frame // frames_per_iter}"
        config = toml.load(current_project / "config.toml")
        ax.set_title(f"XY model for $J = 1$, $T = {config['physical_settings']['temperature']} / k_B$")
        x,y,u,v, angles = _get_arrow_data(frame, import_lattice(current_project,frame % frames_per_iter), config) 
        q.set_UVC(u[::step, ::step],v[::step, ::step],angles[::step, ::step])
        # fig.savefig(save_path_frame / f"frame_{frame:03d}.png", dpi=150)
        return [q]

    ani = FuncAnimation(fig, update, frames = tqdm(range(nframes), desc="Rendering"), interval = fps, blit = True, cache_frame_data = False)
    ani.save(save_path / f"{project_name}_arrow_temp_video.mp4", writer = "ffmpeg", fps =fps, dpi = 300)

