import numpy as np
import matplotlib.pyplot as plt
import os


import time

start_time=time.time()

potential_type="std_well" #finite_well, parabolic_well

hbar=1
m=1
V0=1000
find_optimal_N=True
L=1
N=5000
dx=L/(N-1)
dx_analitical=0.001
N_analitical=int(1+L/dx_analitical)
k=1 #state nb
state_list = [1, 2, 3, 4, 5, 6]
epsilon=1e-8

# ============================ ERROR SENSITIVITY ===================================

def ErrorAccepted(threshold, E_num, E_a):
    if (np.abs(E_num-E_a)<=threshold): return True
    else : return False


# ============================ Analitical solutions =========================
def Analitical_Psi(state_nb):
    
    x_grid=np.empty(N_analitical)
    psi=np.empty(N_analitical)

    for i in range(0,N_analitical):
        x=i*dx_analitical
        x_grid[i]=x
        psi[i]=np.sqrt(2.0/L)*np.sin((state_nb)*np.pi*x/L)

    return x_grid, psi
def Analitical_E(state_nb):
    E=(hbar*hbar)/(2*m)*((state_nb)*np.pi/L)**2
    return E



# ========================== POTENTIAL + NORMALIZATION =====================
def normalization(psi, dx):
    #print("="*30)
    #print("Normalization process")
    A=1.0/np.sqrt(np.inner(psi,psi)*dx)
    psi*=A
    #print("="*30)

def V_init(V,potential_type, L):
    N=len(V)
    dx=L/(N-1)
    if (potential_type=="std_well"):
        pass
    elif (potential_type=="parabolic_well"):
        for i in range(0, N):
            x=i*dx
            V[i]=V0*(x-L/2)**2
    elif (potential_type=="finite_well"):
        for i in range(0, N):
            x=i*dx
            if(x>L/3 and x<2*L/3):
                V[i]=-V0
# ========================== SAVE + PLOT ==========================

def save_and_plot(eigenvectors, eigenvalues, potential_type):
    print("="*30)
    print("Saving & Plotting")
    if (potential_type=="std_well"):
        eigenVal_filename="E_std_well.txt"
        eigenVec_dir="Psi_std_well"
    elif (potential_type=="parabolic_well"):
        eigenVal_filename="E_par_well.txt"
        eigenVec_dir="Psi_par_well"
        V_plot = np.zeros(N)
        V_plot = (4/L)*(x_grid - L/2)**2 #better plot
    elif (potential_type=="finite_well"):
        eigenVal_filename="E_fin_well.txt"
        eigenVec_dir="Psi_fin_well"
        V_plot = np.zeros(N)
        mask= (x_grid > L/3) & (x_grid < 2*L/3)
        V_plot[mask] = -1 #better plot
    
        
    np.savetxt(eigenVal_filename, eigenvalues, fmt='%10.5f')
    

    for idx, it in enumerate(state_list):
        os.makedirs(eigenVec_dir, exist_ok=True)
        file_name = f"./{eigenVec_dir}/eigen_vec{it}.txt"
        plot_file =f"./{eigenVec_dir}/eigen_vec{it}.png"
        psi = eigenvectors[:, idx]

        data = np.column_stack((x_grid, psi))
        np.savetxt(file_name, data, fmt='%10.5f', header="x_grid \t psi")
        x_grid_analitcal, psi_analitical = Analitical_Psi(it)
        E_analitical=Analitical_E(it)

        plt.figure(figsize=(8, 5))
        plt.plot(x_grid, psi**2, label=fr"Numeryczne $|\Psi_{{{it}}}(x)|^2$", color='blue', linewidth=2)
        if (potential_type=="std_well"):
            plt.plot(x_grid_analitcal, psi_analitical**2, label=fr"Analityczne $|\Psi_{{{it}}}(x)|^2$", color='red', linewidth=2, linestyle='--')
            plt.title(fr"$E^{{\mathrm{{num}}}}_{{{it}}} = {eigenvalues[idx]:.6f}$" + "\n" + fr"$E^{{\mathrm{{analit}}}} = {E_analitical:.6f}$", fontsize=14)
        elif potential_type == "parabolic_well":
            plt.title(fr"$E^{{\mathrm{{num}}}}_{{{it}}} = {eigenvalues[idx]:.6f}$", fontsize=14)
            plt.plot(x_grid, V_plot, label=r"$V(x)$", color='black', linewidth=1.5)
        elif potential_type == "finite_well":
            plt.title(fr"$E^{{\mathrm{{num}}}}_{{{it}}} = {eigenvalues[idx]:.6f}$", fontsize=14)
            plt.plot(x_grid, V_plot, label=r"$V(x)$", color='black', linewidth=1.5)
    
        plt.xlabel(r"$x$", fontsize=13)
        plt.ylabel(r"$|\Psi(x)|^2$", fontsize=13)
        plt.grid(True, linestyle=':', alpha=0.7)
        plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left', borderaxespad=0., frameon=True, shadow=True)
        plt.tight_layout()

        plt.savefig(plot_file, dpi=300)
        plt.close()
    print("="*30)

    
# ================================= SOLVER ==================================

def shooting(N, L, potential_type, E_L, E_P):
    dx = L/(N-1)
    Psi=np.zeros(N)
    Psi_deriv=np.zeros(N)
    V=np.zeros(N)
    V_init(V, potential_type, L)
    E_mid=E_L
    eigenvalues=[]
    eigenvectors=np.zeros((N, len(state_list)))
    Psi[N-1]=1000
    iter_tmp=0
    for state in state_list:
        E_L=E_mid
        E_P=1000
        Psi[N-1]=1000
        while(np.abs(Psi[N-1])>epsilon and E_P-E_L>epsilon):
            nodes=0
            Psi[0]=0
            Psi_deriv[0]=1
            E_mid=(E_P+E_L)/2.0
            k2 = (2.0 * m / hbar**2) * (E_mid - V)
            dx2_12 = (dx**2) / 12.0

            Psi[0] = 0.0
            Psi[1] = 1e-5

            for i in range(1, N - 1):
                term_i_minus = (1.0 + dx2_12 * k2[i-1]) * Psi[i-1]
                term_i = 2.0 * (1.0 - 5.0 * dx2_12 * k2[i]) * Psi[i]
                denom = 1.0 + dx2_12 * k2[i+1]
    
                Psi[i+1] = (term_i - term_i_minus) / denom
    
                if (Psi[i] * Psi[i+1]) < 0:
                    nodes += 1
            if (nodes<=state - 1): E_L=E_mid
            else: E_P=E_mid
    
        normalization(Psi, dx)
        eigenvalues.append(E_mid)
        eigenvectors[:,iter_tmp]=Psi
        iter_tmp+=1

    return eigenvalues, eigenvectors

x_grid=np.zeros(N)
for it in range(len(x_grid)):
    x_grid[it]=it*dx


# ============================= ENGINE INFITNITE WELL/ N_OPTIM ==================
solv_start_time=time.time()
if (find_optimal_N==False):
    print(f"Computing solutions for N = {N}...")
    eigenvalues, eigenvectors = shooting(N, L, potential_type, 0, 5000)
    print(f"DONE\n")
    print("*"*40)  
else:
    N_min = N
    N_max = 20*N

    print("Checking max N...\n")
    while(True):
        eigenvalues, eigenvectors = shooting(N_max, L, potential_type, 0, 5000)
        N_accept=True
        for idx, state in enumerate(state_list):
            N_accept=ErrorAccepted(0.0001, eigenvalues[idx], Analitical_E(state))
            if(not N_accept): break
        
        if (not N_accept):
            N_max=int(N_max*1.5)
            print(f"N_max not accepted\nNew N_max is {N_max}")
        else:
            print(f"N_max accepted\nN_max is {N_max}")
            print("="*40)
            print("\n")
            break
    print(f"Checking minimal effective N...")
    iterration=0
    while(N_max-N_min>10):
        iterration+=1
        N_mid=int((N_max+N_min)/2)
        #print(f"acc N = {N_mid}\n N_min = {N_min}\n N_max={N_max}\n")
        
        if (N_mid - 2 < k):
            N_min = N_mid
            continue

        eigenvalues, eigenvectors = shooting(N_mid,L, potential_type, 0, 5000)
        N_accept=True
        for idx, state in enumerate(state_list):
            N_accept=ErrorAccepted(0.0001, eigenvalues[idx], Analitical_E(state))
            if(not N_accept): break
        
        if(not N_accept):
            N_min=N_mid
            #print("Error not accepted\n")
        else:
            N_max=N_mid
            #print("Error accepted\n")
        
    print("\nLinear optimalization...")
    N= N_max
    
    for test_N in range(N_min, N_max + 1):
        iterration+=1
        eigenvalues, eigenvectors = shooting(test_N, L, potential_type, 0, 5000)
        
        N_accept = True
        for idx, state in enumerate(state_list):
            if (not ErrorAccepted(0.0001, eigenvalues[idx], Analitical_E(state))):
                N_accept = False
                break
                
        if N_accept:
            N = test_N
            break

    print("*"*40)
    print(f"Optimal N found, N = {N}\n\n Error for each state:\n")
    for idx, state in enumerate(state_list):
        print(f"{state}\t{np.abs(eigenvalues[idx] - Analitical_E(state))}\n")
    print("*"*40)
    print(f"Number of iterrations (bisection+linear optim) = {iterration}")
    print("="*40)

dx=L/(N-1)
x_grid=np.zeros(N)
for it in range(len(x_grid)):
    x_grid[it]=it*dx
solv_end_time=time.time()
save_and_plot(eigenvectors, eigenvalues, potential_type)

# ======================== ENGINE parabolic well ====================
start_time_PW=time.time()
#print(f"Using N given/from optimalization of infinite well\nCalculating solution for parabolic potential, given L={L}")
#potential_type="parabolic_well"
#eigenvalues, eigenvectors = shooting(N, L, potential_type, 0, V0)

#save_and_plot(eigenvectors, eigenvalues, potential_type)
#print("="*40)

# ======================== L TEST ================================

#L_min=1
#L_step=1
#L_max=5
#os.makedirs("./L_test", exist_ok=True)

#print(f"Calculating eigen values for parabolic pot. with L from range ({L_min}, {L_max}) with step {L_step}")
#for test_L in range(L_min, L_max, L_step):
#    eigenvalues, _ = shooting(N, test_L, potential_type, 0, V0)
#    filename_dir=f"./L_test/L={test_L}.txt"
#    np.savetxt(filename_dir, eigenvalues, fmt='%10.5f')
#
end_time_PW=time.time()
#============================ FINITE WELL ===========================

start_time_FW=time.time()
print(f"Using N given/from optimalization of infinite well\nCalculating solution for finite well, given L={L}...")
potential_type="finite_well"
eigenvalues, eigenvectors = shooting(N, L, potential_type, -1.5*V0, 0.5*V0)
end_time_FW=time.time()
save_and_plot(eigenvectors, eigenvalues, potential_type)
print("="*40)


end_time=time.time()

print(f"total time = {end_time-start_time} s")
#print(f"total time infinite well = {end_time-start_time} s")
print(f"solve time (infinite well) {solv_end_time-solv_start_time} s")
print(f"total time finite well = {end_time_FW-start_time_FW} s")