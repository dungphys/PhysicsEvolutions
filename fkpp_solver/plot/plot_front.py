import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

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

def plot_solution(fileX, times, cmap_name="viridis", save_as="fkpp_solution.png"):
    n = len(times)
    cmap = mpl.colormaps[cmap_name].resampled(n)
    fig, ax = plt.subplots(figsize=(8, 5))
    x = pd.read_csv(fileX, header=None).values.flatten()
    u = [None for i in range(len(times))]
    for i, t in enumerate(times):  
        fileU = f"../bin/solution_{t:.3f}.txt"
        u[i] = pd.read_csv(fileU, header=None).values.flatten()
        ax.plot(x, u[i], color=cmap(i / (n - 1)), label=f"t = {t:.2f}", linewidth=2)
    
    # Colorbar to show time progression
    sm = mpl.cm.ScalarMappable(cmap=cmap,
                           norm=plt.Normalize(vmin=times[0], vmax=times[-1]))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Time")

    ax.set_xlabel("x")
    ax.set_xlim(-40, 80)
    ax.set_xscale("linear")
    ax.set_ylabel("u")
    ax.set_yscale("linear")
    ax.set_title("FKPP Solution over Time")
    plt.tight_layout()
    if save_as:
        plt.savefig(f"figures/{save_as}", dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    times_to_plot = np.arange(0., 50., 0.5)  #
    plot_solution("../bin/grid_x.txt", times_to_plot)
 