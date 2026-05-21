import numpy as np
import matplotlib.pyplot as plt
import os

# Parametry siatki dla kropki kwantowej zgodne z najnowszym kodem w C++
Nx = 500
Ny = 500
Lx = 5.0
Ly = 5.0

# Położenie przekroju y = Ly/2
center_y = Ly / 2.0
mid_y_idx = Ny // 2  # Indeks środkowego wiersza w osi Y

# --- USTAWIENIE STAŁYCH GRANIC POTENCJAŁU ---
# Zmień te wartości w zależności od tego, jak wysoko rośnie potencjał w Twojej symulacji
V_MIN = -1.5
V_MAX = 2.0
# --------------------------------------------

# Pętla po iteracjach (od 1 do 7)
for i in range(1, 8):
    filename = f"outdir/SOR_quantumDot_iterr{i}"
    
    if not os.path.exists(filename):
        print(f"Ostrzeżenie: Nie znaleziono pliku {filename}. Pomijam.")
        continue

    # Obliczenie wartości parametru rho_cst zdefiniowanego w pętli C++
    rho_cst = i * 0.5

    # Wczytanie danych z pliku (x, y, V)
    data = np.loadtxt(filename)
    
    # Przekształcenie płaskich kolumn w macierze 2D (wiersze: Y, kolumny: X)
    X = data[:, 0].reshape((Ny, Nx))
    Y = data[:, 1].reshape((Ny, Nx))
    V = data[:, 2].reshape((Ny, Nx))
    
    # Wyciągnięcie danych do wykresu 1D (wzdłuż osi X dla stałego Y)
    x_profile = X[mid_y_idx, :]
    v_profile = V[mid_y_idx, :]

    # Tworzenie panelu z dwoma wykresami
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # --- WYKRES 1: Mapa cieplna 2D potencjału ---
    # Dodano parametry vmin i vmax, aby zamrozić skalę kolorów dla wszystkich klatek
    heatmap = ax1.pcolormesh(X, Y, V, shading='auto', cmap='magma', vmin=V_MIN, vmax=V_MAX)
    ax1.set_title(f"Mapa potencjału V(x,y) dla kropki, rho_cst = {rho_cst}")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_aspect('equal', 'box') 
    
    # Zaznaczenie poziomą, przerywaną linią miejsca cięcia (y = Ly/2)
    ax1.axhline(y=center_y, color='cyan', linestyle='--', linewidth=1.5, label=f'Przekrój y={center_y}')
    ax1.legend(loc='upper right')
    fig.colorbar(heatmap, ax=ax1, label="Potencjał V")

    # --- WYKRES 2: Profil 1D przez środek (V vs x dla y = Ly/2) ---
    ax2.plot(x_profile, v_profile, color='purple', linewidth=2)
    ax2.set_title(f"Profil potencjału V(x) dla y = {center_y}")
    ax2.set_xlabel("x")
    ax2.set_ylabel("V")
    
    # ZAMROŻENIE OSI Y NA WYKRESIE 1D
    ax2.set_ylim([V_MIN, V_MAX])
    
    ax2.grid(True, linestyle=':', alpha=0.7)

    # Optymalizacja odstępów i zapis do pliku
    plt.tight_layout()
    plot_filename = f"outdir/plot_qD_V_iterr{i}.png"
    plt.savefig(plot_filename, dpi=300)
    print(f"Zapisano wykres: {plot_filename}")
    
    plt.close()