import numpy as np
import matplotlib.pyplot as plt



Nx = 500
Ny = 500
Lx = 5.0
Ly = 5.0


center_y = Ly / 2.0
mid_y_idx = Ny // 2  


V_MIN = -4.0
V_MAX = 1.5
# --------------------------------------------

for i in range(1, 8):
    filename = f"outdir/SOR_quantumDot_iterr{i}"

    Vg = -i * 0.5

    data = np.loadtxt(filename)
    
    X = data[:, 0].reshape((Ny, Nx))
    Y = data[:, 1].reshape((Ny, Nx))
    V = data[:, 2].reshape((Ny, Nx))
    
    x_profile = X[mid_y_idx, :]
    v_profile = V[mid_y_idx, :]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    #mapa
    heatmap = ax1.pcolormesh(X, Y, V, shading='auto', cmap='magma', vmin=V_MIN, vmax=V_MAX)
    ax1.set_title(f"Mapa potencjału V(x,y) dla kropki, Vg = {Vg}")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_aspect('equal', 'box') 
    
    ax1.axhline(y=center_y, color='cyan', linestyle='--', linewidth=1.5, label=f'Przekrój y={center_y}')
    ax1.legend(loc='upper right')
    fig.colorbar(heatmap, ax=ax1, label="Potencjał V")

    #profil
    ax2.plot(x_profile, v_profile, color='purple', linewidth=2)
    ax2.set_title(f"Profil potencjału V(x) dla y = {center_y}")
    ax2.set_xlabel("x")
    ax2.set_ylabel("V")
    

    ax2.set_ylim([V_MIN, V_MAX])
    ax2.grid(True, linestyle=':', alpha=0.7)

    plt.tight_layout()
    plot_filename = f"outdir/plot_qD_Vg_iterr{i}.png"
    plt.savefig(plot_filename, dpi=300)
    
    plt.close()