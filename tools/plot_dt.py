import os
import math
import numpy as np
import matplotlib.pyplot as plt


def logDynamicTimestep(K, t_N, Nsteps, n):
    # Handle K ~ 0 (uniform steps)
    if abs(K) < 1e-12:
        return t_N / Nsteps
    num = Nsteps + (n + 1) * (math.pow(10.0, K) - 1.0)
    den = Nsteps + n * (math.pow(10.0, K) - 1.0)
    return (t_N / K) * math.log10(num / den)


def make_plots(outpath):
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    t_N = 1.0
    Ks = [-1.0, -0.5, -0.1, 0.0, 0.1, 0.5, 1.0]
    Nsteps_list = [10, 50]

    fig, axes = plt.subplots(1, len(Nsteps_list), figsize=(10, 4), sharey=True)

    for ax, Nsteps in zip(axes, Nsteps_list):
        n = np.arange(0, Nsteps)
        for K in Ks:
            dt_vals = [logDynamicTimestep(K, t_N, Nsteps, int(i)) for i in n]
            ax.plot(n, dt_vals, marker='o', label=f'K={K}')
        ax.set_title(f'Nsteps={Nsteps}')
        ax.set_xlabel('Step index n')
        ax.grid(True)
        ax.legend(fontsize='small')

    axes[0].set_ylabel('dt (same units as t_N)')
    fig.suptitle('Dynamic timestep dt vs step index for various K and Nsteps')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(outpath, dpi=150)
    print('Saved plot to', outpath)


if __name__ == '__main__':
    make_plots('plots/dt_examples.png')
