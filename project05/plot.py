import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt('E.dat')



E = data[:, 0]  
T = data[:, 1]

plt.xlim(E.min(), E.max())
plt.plot(E, T, color='red', linewidth=2, label="T(E)")
plt.xlabel("E")
plt.ylabel("T")
plt.legend()
plt.savefig('T(E).png', dpi=300)



L = 5.0
a_s = 1.0

a_d = 1.0
b = 0.5

V0_s = 10.0
sigma = 1.0

PSI_data=np.loadtxt("./outdir/PSI_E4.9605.dat")

x=PSI_data[:,0]
psi=PSI_data[:,1]

x_start = L/2 - a_s/2
x_end = L/2 + a_s/2

d=b/2+a_d/2
barrier2_start1=L/2-d-a_d/2
barrier2_end1 = L/2-d+a_d/2
barrier2_start2=L/2+d-a_d/2
barrier2_end2=L/2+d+a_d/2

plt.figure(figsize=(8,6))
plt.xlim(x.min(), x.max())

#plt.axvspan(x_start, x_end, color='lightgray', alpha=0.5, label=f'V={V0_s}')    #single
plt.plot(x, 2.0*np.exp(-((x-L/2)/sigma)**2), color='lightgrey', linewidth=2, label=f'V={V0_s}')    #gauss

#plt.axvspan(barrier2_start1, barrier2_end1, color='lightgray', alpha=0.5, label=f'V={V0_s}')
#plt.axvspan(barrier2_start2, barrier2_end2, color='lightgray', alpha=0.5)


plt.plot(x,psi*psi, color='red', linewidth=2, label=r"$|\psi|^2$")
plt.xlabel("x")
plt.ylabel(r"$|\psi(x)|^2$")
plt.legend()
plt.savefig('psi.png', dpi=300)