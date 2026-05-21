import numpy as np
import os
import matplotlib.pyplot as plt

N=500
delta=1e-8

m=1
hbar=1
V0=1
a=2

#dx = L/(N-1)



Psi = np.zeros(N)


# =========================== SAVE & PLOT =================
def save_and_plot(eigenvector, eigenvalue, potential_type, acc_L):
    print("="*30)
    print("Saving & Plotting")
    dx=acc_L/(N-1)
    
    x_grid=np.zeros(N)
    psi_analitical=np.zeros(N)
    for i in range(N):
        x_grid[i]=dx*i
        psi_analitical[i]=np.sqrt(2.0/acc_L)*np.sin(np.pi * x_grid[i]/acc_L)
    E_analitical = ((np.pi*hbar/acc_L)**2)/(2*m)
    

    if (potential_type=="std_well"):
        eigenVec_dir="Psi_std_well"
    elif (potential_type=="fin_well"):
        eigenVal_filename="E_fin_well.txt"
        eigenVec_dir="Psi_fin_well"
        V_plot = np.zeros(N)
        mask= (x_grid > (acc_L/2)-a) & (x_grid < (acc_L/2)+a)
        V_plot[mask] = -1 #better plot
    
    os.makedirs(eigenVec_dir, exist_ok=True)
    #np.savetxt(eigenVal_filename, eigenvalue, fmt='%10.5f')
    file_name = f"./{eigenVec_dir}/eigen_vec_{acc_L}.txt"
    plot_file =f"./{eigenVec_dir}/eigen_vec_{acc_L}.png"


    plt.figure(figsize=(8, 5))
    data = np.column_stack((x_grid, eigenvector))
    np.savetxt(file_name, data, fmt='%10.5f', header="x_grid \t psi")
    plt.plot(x_grid, eigenvector**2, label=fr"Numeryczne $|\Psi(x)|^2$", color='blue', linewidth=2)
    if (potential_type=="std_well"):
        plt.plot(x_grid, psi_analitical**2, label=fr"Analityczne $|\Psi_1(x)|^2$", color='red', linewidth=2, linestyle='--')
        plt.title(fr"$E^{{\mathrm{{num}}}}_1 = {eigenvalue:.6f}$" + "\n" + fr"$E^{{\mathrm{{analit}}}} = {E_analitical:.6f}$", fontsize=14)
    if potential_type == "fin_well":
        plt.title(fr"$E^{{\mathrm{{num}}}} = {eigenvalue:.6f}$", fontsize=14)
        plt.plot(x_grid, V_plot, label=r"$V(x)$", color='black', linewidth=1.5)
    plt.xlabel(r"$x$", fontsize=13)
    plt.ylabel(r"$|\Psi(x)|^2$", fontsize=13)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left', borderaxespad=0., frameon=True, shadow=True)
    plt.tight_layout()

    plt.savefig(plot_file, dpi=300)
    plt.close()
        
    
    
        
    print("="*30)


# ============================== W.P. ====================
def Psi_init(Psi, N, acc_L):
    dx=acc_L/(N-1)
    for idx in range(N):
        x = idx*dx
        Psi[idx]=x*(x-acc_L)
# =========================================================


# ======================== Potential ======================

def V_init(potential_type, N, acc_L):
    V=np.zeros(N)
    if (potential_type=="std_well"):
        pass
    elif (potential_type=="fin_well"):
        dx=acc_L/(N-1)
        for idx in range(N):
            #print(idx*dx)
            x = idx*dx
            if (x>acc_L/2-a and x<acc_L/2+a):
                V[idx]=-V0
    return V


# ======================== Engine =========================
def Energy(Psi, V, acc_L, N):
    dx = acc_L/(N-1)
    E=0
    a=hbar*hbar/(2*m*dx*dx)
    for idx in range(1,N-1):
        E+=Psi[idx]*(-a*(Psi[idx-1]-2.0*Psi[idx]+Psi[idx+1])+V[idx]*Psi[idx])*dx

    return E

def normalization(Psi, acc_L, N):
    dx=acc_L/(N-1)
    sum=0
    for idx in range(N):
        sum+=Psi[idx]*Psi[idx]*dx
    Psi*=1.0/np.sqrt(sum)


def solver(Psi_old, V, acc_L, N):
    #print(V)
    normalization(Psi_old, acc_L, N)
    epsilon=1
    dx = acc_L/(N-1)
    dt = 0.5 * m*dx*dx/hbar
    c=hbar*dt/(2*m*dx*dx)
    E_old=Energy(Psi_old, V, acc_L, N)
    it=0
    Psi_new=np.zeros(N)
    while(epsilon>delta):
        it+=1
        for idx in range(1,N-1):
            Psi_new[idx]=c*(Psi_old[idx-1]-2.0*Psi_old[idx]+Psi_old[idx+1]) - (V[idx]*dt/hbar -1.0)*Psi_old[idx]
        normalization(Psi_new, acc_L, N)
        E_new=Energy(Psi_new, V,acc_L, N)
        epsilon=np.abs(E_new-E_old)

        E_old=E_new
        Psi_old=Psi_new

    return E_new, Psi_new, it



L=1
print("Solving standard well...")
Psi_init(Psi,N, L)
V = V_init("std_well", N, L)
E, Psi, it_nb = solver(Psi, V, L, N)
print(f"STD well solved in {it_nb}, Energy = {E}")
save_and_plot(Psi, E, "std_well", L)

print(f"Solving finite well, a = {a}...")
L_min=2*a+1
L_step=1
L_max=8*a
E_tab=[]

for acc_L in range(L_min, L_max, L_step):
    V = V_init("fin_well", N, acc_L)
    #print(V)
    Psi_init(Psi, N, acc_L)
    E, Psi, it_nb = solver(Psi, V, acc_L, N)
    print(f"FIN well solved in {it_nb}, L = {acc_L}, Energy = {E}")
    save_and_plot(Psi, E, "fin_well", acc_L)
    E_tab.append(E)

plt.figure(figsize=(8, 5))

plt.plot(range(L_min, L_max, L_step), E_tab, marker='o', linestyle='none', label=fr"$E(L)$", color='blue')
plt.xlabel("L", fontsize=13)
plt.ylabel("E", fontsize=13)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left', borderaxespad=0., frameon=True, shadow=True)
plt.tight_layout()
plt.savefig("E(L).png", dpi=300)
plt.close()