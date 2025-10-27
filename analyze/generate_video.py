import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation
from IPython.display import HTML


def generate_ising_grid(time_step:int, df:pd.DataFrame)->np.array:
    L = int(np.sqrt(df.shape[1]-1))
    ground_1d_grid = df.loc[time_step].to_numpy()[1:]
    ground_grid = np.array([ground_1d_grid[:L-1]])
    for i in range(1,L):
        a = ground_1d_grid[i*L: i*L + (L-1)]
        ground_grid = np.vstack((ground_grid,a))
    return ground_grid

def main():
    df = pd.read_csv("../build/results.csv", skiprows=1)

    fig, ax = plt.subplots()
    im = ax.imshow(generate_ising_grid(0, df))

    def update(frame):
        im.set_array(generate_ising_grid(frame, df))
        return [im]
    ani = FuncAnimation(fig, update, frames = 300, interval = 50, blit = True)
    html = ani.to_jshtml()
    with open("ising_animation.html", "w") as f:
        f.write(html)
    ani.save("ising.mp4", writer = "ffmpeg", fps = 10)


if __name__ == "__main__":
    main()
