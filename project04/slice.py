import matplotlib.pyplot as plt
import numpy as np
import os

def plot_by_frame():
    outdir = './outdir'
    
    frame_input = input("Podaj numer klatki ")
    filename = os.path.join(outdir, f"{frame_input}.dat")

    try:
        data = np.loadtxt(filename)
        x = data[:, 0]
        abs_psi = data[:, 3]

        prob_density = abs_psi**2

        plt.figure(figsize=(10, 6))
        plt.plot(x, prob_density, color='crimson', linewidth=2, label=f'Klatka: {frame_input}')
        
        plt.title(f"Gęstość prawdopodobieństwa - Plik: {frame_input}.dat", fontsize=14)
        plt.xlabel("Pozycja (x)", fontsize=12)
        plt.ylabel(r"$|\Psi(x)|^2$", fontsize=12)
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.legend()

        plt.ylim(0, max(prob_density) * 1.1)
        plt.show()

    except Exception as e:
        print(f"Wystąpił błąd podczas wczytywania danych: {e}")

if __name__ == "__main__":
    plot_by_frame()