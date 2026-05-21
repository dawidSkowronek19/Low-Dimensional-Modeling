import numpy as np
import matplotlib.pyplot as plt

Nx = 500
Ny = 200
Lx = 5.0
Ly = 2.0


center_x = Lx / 2.0
mid_x_idx = Nx // 2

for i in range(1, 10):
    filename = f"outdir_wire/SOR_qW_iterr{i}"

    Vd = -1.0 - i 
    data = np.loadtxt(filename)
    
    X = data[:, 0].reshape((Ny, Nx))
    Y = data[:, 1].reshape((Ny, Nx))
    V = data[:, 2].reshape((Ny, Nx))

    y_profile = Y[:, mid_x_idx]
    v_profile = V[:, mid_x_idx]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    #mapa
    heatmap = ax1.pcolormesh(X, Y, V, shading='auto', cmap='viridis')
    ax1.set_title(f"Mapa potencjału V(x,y) dla drutu, V_d = {Vd}")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_aspect('equal', 'box')
    
    ax1.axvline(x=center_x, color='blue', linestyle='--', linewidth=1.5, label=f'Przekrój x={center_x}')
    ax1.legend(loc='upper right')
    fig.colorbar(heatmap, ax=ax1, label="Potencjał V")

    # profil
    ax2.plot(y_profile, v_profile, color='blue', linewidth=2)
    ax2.set_title(f"Profil potencjału w poprzek drutu (x = {center_x})")
    ax2.set_xlabel("y")
    ax2.set_ylabel("V")
    ax2.grid(True, linestyle=':', alpha=0.7)

    plt.tight_layout()
    plot_filename = f"outdir/plot_qW_iterr{i}.png"
    plt.savefig(plot_filename, dpi=300)
 
    plt.close()