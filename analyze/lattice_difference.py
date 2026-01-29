import os
import sys
import pathlib
# Set path
project_root = pathlib.Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

import cv2
import numpy as np
import matplotlib.pyplot as plt

def main():
    project_path = pathlib.Path("/home/alvaro/Documents/trinity/year4/capstone/capstone-code/analyze/vid_dump") 
    parameter_combination_1 = 24
    parameter_combination_2 = 25
    # Load videos
    v1 = cv2.VideoCapture(project_path / f"temperature50_0-3_1-5_l128_dim2_10-3sweeps_par_{parameter_combination_1}_lattice.mp4")
    v2 = cv2.VideoCapture(project_path / f"temperature50_0-3_1-5_l128_dim2_10-3sweeps_par_{parameter_combination_2}_lattice.mp4")

    # Read first frame (or loop over frames)
    ret1, f1 = v1.read()
    ret2, f2 = v2.read()

    # Convert to grayscale (optional but cleaner)
    g1 = cv2.cvtColor(f1, cv2.COLOR_BGR2GRAY)
    g2 = cv2.cvtColor(f2, cv2.COLOR_BGR2GRAY)

    # Absolute difference
    diff = cv2.absdiff(g1, g2)

    # Display
    plt.imshow(diff, cmap="hot")
    plt.colorbar(label="Pixel difference")
    plt.title(f"Difference Frame between parameter combination {parameter_combination_1} and {parameter_combination_2}")
    plt.show()

if __name__ == "__main__":
    main()
