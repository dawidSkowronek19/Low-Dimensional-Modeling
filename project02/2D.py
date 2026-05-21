import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.sparse import diags, eye
from scipy.sparse.linalg import eigsh

L=1
N=100
delta=1e-6

V0=1
m=1
hbar=1


dx = L/(N-1)
dt = 0.5 * m*dx*dx/(2*hbar)
sigma=L/8

Psi = np.zeros(N*N)



# =========================== Rho, SAVE & PLOT =================

def rho(const_val, Psi):
    cst_idx = int(const_val/dx)
    rho_x=np.zeros(N)
    rho_y=np.zeros(N)
    for idx in range(N):
        l_x = idx + cst_idx*N
        rho_x[idx] = Psi[l_x]**2

        l_y=cst_idx+idx*N
        rho_y[idx]=Psi[l_y]**2
    
    return rho_x, rho_y



def save_and_plot(eigenvector, eigenvalue, potential_type):

    output_dir = f"2D_results_{potential_type}"
    os.makedirs(output_dir, exist_ok=True)
    
    cut_val = L/2
    print(f"--- Plotting rho(x,{cut_val}): {potential_type}")
    rho_x, rho_y = rho(cut_val, eigenvector)
    
    x_grid = np.linspace(0, L, N)

    plt.figure(figsize=(8,6))
    plt.title(f"rho(x, y={cut_val}) - {potential_type}")
    plt.plot(x_grid, rho_x)
    plt.xlabel("x")
    plt.ylabel("rho")
    plt.grid(True)
    plt.savefig(f"{output_dir}/rho_x.png", dpi=300)
    plt.close()

    plt.figure(figsize=(8,6))
    plt.title(f"rho(x={cut_val}, y) - {potential_type}")
    plt.plot(x_grid, rho_y)
    plt.xlabel("y")
    plt.ylabel("rho")
    plt.grid(True)
    plt.savefig(f"{output_dir}/rho_y.png", dpi=300)
    plt.close()
    print(f"--- Plotting RGB Heatmap: {potential_type} ---")
    psi_matrix = eigenvector.reshape((N, N))
    prob_density = psi_matrix**2

    plt.figure(figsize=(8, 6))
    img = plt.imshow(prob_density, extent=[0, L, 0, L], origin='lower', cmap='jet')
    plt.colorbar(img, label=r'$|\Psi(x,y)|^2$')
    plt.title(f"rho(x,y) ({potential_type})\nE = {eigenvalue:.5f}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig(f"{output_dir}/heatmap_rgb.png", dpi=300)
    plt.close()

    
# ======================== Potential ======================
def V_init(potential_type, N, L):
    V = np.zeros(N*N)
    if (potential_type == "std_well"):
        pass
    elif (potential_type == "exp"):
        for idx_y in range(N):
            for idx_x in range(N):
                l = idx_x + idx_y * N
                x_coord = idx_x * dx
                y_coord = idx_y * dx
                V[l] = -V0 * np.exp(-((x_coord - L/2)**2 + (y_coord - L/2)**2) / (sigma**2))
    return V

   

# ============================== W.P. ====================
def Psi_init(Psi, N, L):
    for idx_x in range(N):
        for idx_y in range(N):
            l = idx_x+idx_y*N
            x = idx_x*dx
            y = idx_y*dx
            Psi[l]=x*(x-L)*y*(y-L)
# =========================================================

# ======================== Engine =========================
def Energy(Psi, V, L, N):
    dx = L/(N-1)
    E=0
    a=hbar*hbar/(2*m*dx*dx)
    for idx_y in range(1,N-1):
        for idx_x in range(1,N-1):
            l=idx_x+idx_y*N
            E+=Psi[l]*(-a*(Psi[l-1]-4.0*Psi[l]+Psi[l+1]+Psi[l-N]+Psi[l+N])+V[l]*Psi[l])*dx*dx

    return E

def normalization(Psi, L, N):
    dx=L/(N-1)
    sum=0
    for idx in range(N*N):
        sum+=Psi[idx]*Psi[idx]*dx*dx
    Psi*=1.0/np.sqrt(sum)


def solver_tau(Psi_old, V, L, N):
    #print(V)
    normalization(Psi_old, L, N)
    epsilon=1
    dx = L/(N-1)
    a=hbar*dt/(2*m*dx*dx)
    E_old=Energy(Psi_old, V, L, N)
    it=0
    Psi_new=np.zeros(N*N)
    while(epsilon>delta):
        it+=1
        #print(it)
        for idx_y in range(1,N-1):
            for idx_x in range(1,N-1):
                l = idx_x+idx_y*N
                Psi_new[l]=a*(Psi_old[l-1]-4.0*Psi_old[l]+Psi_old[l+1] +Psi_old[l-N]+Psi_old[l+N]) - (V[l]*dt/hbar -1.0)*Psi_old[l]
        normalization(Psi_new, L, N)
        E_new=Energy(Psi_new, V,L, N)
        epsilon=np.abs(E_new-E_old)

        E_old=E_new
        Psi_old=Psi_new

    return E_new, Psi_new, it

def solver_diag(V):    
    main_diag = np.full(N*N, -4.0)
    side_diag = np.full(N*N - 1, 1.0)
    
    side_diag[np.arange(1, N*N) % N == 0] = 0.0
    
    far_diag = np.full(N*N - N, 1.0)

    laplacian = diags(
        [main_diag, side_diag, side_diag, far_diag, far_diag],
        [0, -1, 1, -N, N],
        format='csc'
    )
    coeff = - (hbar**2) / (2 * m * dx**2)
    T = coeff * laplacian
    V_sparse = diags([V], [0], format='csc')

    H = T + V_sparse

    energies, wavefunctions = eigsh(H, 1, which='SA')
    normalization(wavefunctions[:, 0], L, N)


    return energies[0], wavefunctions[:,0]
    


print("Solving standard well (Imaginary Time)...")
Psi_init(Psi,N, L)
V = V_init("std_well", N, L)
E, Psi, it_nb = solver_tau(Psi, V, L, N)
print(f"STD well solved in {it_nb}, Energy = {E}")
save_and_plot(Psi, E, "std_well")


print("Solving exp_pot well (Imaginary Time)...")
Psi_init(Psi,N, L)
V = V_init("exp", N, L)
E, Psi, it_nb = solver_tau(Psi, V, L, N)
print(f"STD well solved in {it_nb}, Energy = {E}")
save_and_plot(Psi, E, "exp_well")

""""
print("Solving exp_pot well (Diagonalization)...")
Psi_init(Psi,N, L)
V = V_init("exp", N, L)
E, Psi = solver_diag(V)
print(f"STD well solved, Energy = {E}")
save_and_plot(Psi, E, "exp_well")
"""