import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import brentq

plt.rcParams.update({
    "figure.facecolor" : "white",
    "axes.facecolor"   : "#FEFFFF",
    "axes.edgecolor"   : "gray",
    "axes.spines.top"  : False,
    "axes.spines.right": False,
    "axes.grid"        : True,
    "grid.color"       : "white",
    "grid.linewidth"   : 0.3,
    "font.family"      : "sans-serif",
    "axes.titlesize"   : 14,
    "axes.labelsize"   : 12,
})

def find_front(x, u, thres=0.5):
    """
    Find x where u(x) = level using Brentq (bracketing root finder).
    Assumes u is monotone (e.g. a travelling wave front).
    """
    # Find index where u crosses the level
    for i in range(len(u) - 1):
        if (u[i] - thres) * (u[i+1] - thres) < 0:  # sign change
            # Interpolate within this interval
            u_interp = lambda x_: np.interp(x_, x, u) - thres
            x_front = brentq(u_interp, x[i], x[i+1])
            return x_front

    return None  # no crossing found

def front_position(fileX, times, thres=0.5, save_as="fkpp_front_pos.png"):
    fig, ax = plt.subplots(figsize=(8, 5))
    x = pd.read_csv(fileX, header=None).values.flatten()
    front_positions = []
    for i, t in enumerate(times):  
        fileU = f"../bin/solution_{t:.3f}.txt"
        u = pd.read_csv(fileU, header=None).values.flatten()
        front_pos = find_front(x, u, thres)
        front_positions.append(front_pos)
    
    ax.plot(times, front_positions, '-o', color='blue', linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel("Front Position")
    ax.set_title("Front Position over Time")
    plt.tight_layout()
    if save_as:
        plt.savefig(f"figures/{save_as}", dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    times_to_plot = np.arange(0.1, 15, 0.01)  # Plot at t = 0, 20, 40, 60, 80, 100
    front_position("../bin/grid_x.txt", times_to_plot)