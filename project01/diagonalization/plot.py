import numpy as np
import matplotlib.pyplot as plt
import os
import re


input_dir = "./L_test"
output_dir = "./Drabinki_L_Wartosci_Pelne"


Y_LIMIT = 600 

print(f"Tworzenie wykresów w '{output_dir}'...")
os.makedirs(output_dir, exist_ok=True)

if not os.path.exists(input_dir):
    print(f"Błąd: Folder '{input_dir}' nie istnieje.")
else:
    count = 0
    for filename in os.listdir(input_dir):
        match = re.search(r"L=([0-9.]+)\.txt", filename)
        if match:
            L = float(match.group(1))
            filepath = os.path.join(input_dir, filename)
            
            try:
                energies = np.loadtxt(filepath)
                if energies.size == 0: continue
            except:
                continue

            fig, ax = plt.subplots(figsize=(11, 7), layout="constrained")
            
            max_states = min(len(energies) if isinstance(energies, np.ndarray) else 1, 5)
            if max_states == 1 and not isinstance(energies, np.ndarray):
                energies = [energies]

            colors = ['blue', 'orange', 'green', 'red', 'purple']
            
            for i in range(max_states):
                E = energies[i]
                ax.hlines(y=E, xmin=0, xmax=L, color=colors[i], linewidth=2.5, label=f"Stan {i}")

                y_offset = Y_LIMIT * 0.015
                ax.text(L * 0.02, E + y_offset, fr"$E_{{{i}}} = {E}$", 
                        color=colors[i], fontsize=10, fontweight='bold',
                        bbox=dict(facecolor='white', edgecolor='none', alpha=0.6, pad=0.5))
            
            ax.set_title(f"Poziomy energetyczne, $L = {L:.3f}$", fontsize=15)
            ax.set_xlabel("Położenie $x$", fontsize=13)
            ax.set_ylabel("Energia $E$", fontsize=13)
            
            ax.set_ylim(0, Y_LIMIT)
            ax.set_xlim(0, L)
            
            ax.axvline(x=0, color='black', linewidth=1.5)
            ax.axvline(x=L, color='black', linewidth=1.5)
            ax.grid(True, linestyle=':', alpha=0.5)

            ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0)
            
            save_path = os.path.join(output_dir, f"drabinka_L_{L:.3f}.png")
            plt.savefig(save_path, dpi=200)
            plt.close(fig)
            
            count += 1

    print(f"Gotowe! Wygenerowano {count} wykresów.")