import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches
import glob
import os


L = 26.0
a = L / 10.0
barrier_center = L / 2.0
barrier_start = barrier_center - a / 2.0
barrier_end = barrier_center + a / 2.0

file_paths = glob.glob('./outdir/*.dat')
files = sorted(file_paths, key=lambda x: int(os.path.basename(x).split('.')[0]))


fig, ax = plt.subplots(figsize=(10, 5))
ax.set_xlim(0, L)
ax.set_ylim(0, 0.5) 
ax.set_xlabel('Położenie x')
ax.set_ylabel(r'Gęstość prawdopodobieństwa $|\Psi(x, t)|^2$')

barrier = patches.Rectangle(
    (barrier_start, 0), a, 0.4, 
    color='gray', alpha=0.4, label='Bariera'
)
ax.add_patch(barrier)



line, = ax.plot([], [], lw=2, color='blue', label=r'$|\Psi|^2$')
ax.legend(loc='upper right')


def init():
    line.set_data([], [])
    return line,

def update(frame):
    data = np.loadtxt(files[frame])
    x = data[:, 0]
    abs_psi = data[:, 3]
    prob_density = abs_psi**2
    
    line.set_data(x, prob_density)
    

    step = os.path.basename(files[frame]).split(".")[0]
    ax.set_title(f'Ewolucja czasowa paczki falowej | Krok: {step}')
    return line,


ani = FuncAnimation(fig, update, frames=len(files), init_func=init, blit=True)
ani.save('wave_function.gif', writer='pillow', fps=30)
