import numpy as np
import os
import matplotlib.pyplot as plt

dataV=np.loadtxt("./outdir/V.dat")
dataG=np.loadtxt("./outdir/G.dat")

x=dataV[:,0]
y=dataV[:,1]
V=dataV[:,2]

E=dataG[:,0]
G=dataG[:,1]

plt.figure(figsize=(16,8))
plt.title("V")
plt.xlabel("x")
plt.ylabel("y")
heatmap=plt.tripcolor(x,y,V, cmap='viridis', shading='gouraud')
plt.colorbar(heatmap, label='V(x,y)')
plt.savefig("./outdir/V.png")


plt.figure(figsize=(8,8))
plt.title("G")
plt.xlabel("E_f")
plt.ylabel("G")
plt.plot(E,G,color='red')
plt.legend("G(E_f)")
plt.grid()
plt.savefig("./outdir/G.png")