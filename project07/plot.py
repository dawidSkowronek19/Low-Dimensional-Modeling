import numpy as np
import matplotlib.pyplot as plt

data_RS=np.loadtxt("bands.dat")
data_TBA=np.loadtxt("./epsilon.dat")
rgb=["red", "blue", "black", "orange"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)

for i in range(1,5):
    ax1.plot(data_RS[:,0], data_RS[:,i], color=rgb[i-1], label=f"Pasmo {i}")

ax1.set_title("a")
ax1.set_xlabel("k")
ax1.set_ylabel("E")

ax1.legend()
ax1.grid(True)
for i in range(1,4):
    ax2.plot(data_TBA[:,0], data_TBA[:,i], color=rgb[i-1], label=f"Pasmo {i}")

ax2.set_title("b")
ax2.set_xlabel("k")
ax2.set_ylabel("E")
ax2.grid(True)
ax2.legend()

plt.savefig("./bands.png")


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

for i in range(1,5):
    ax1.plot(data_RS[:,0], data_RS[:,i], color=rgb[i-1], label=f"RS", linestyle='--')
    if (i==4): continue
    ax1.plot(data_TBA[:,0], data_TBA[:,i], color=rgb[i-1], label=f"TBA", linewidth=3, alpha=0.4)

ax1.set_title("a")
ax1.set_xlabel("k")
ax1.set_ylabel("E")

ax1.legend()
ax1.grid(True)


ax2.plot(data_RS[:,0], data_RS[:,1], color=rgb[0], label=f"RS", linestyle='--')
ax2.plot(data_TBA[:,0], data_TBA[:,1], color=rgb[0], label=f"TBA", linewidth=3, alpha=0.4)

ax2.set_title("b")
ax2.set_xlabel("k")
ax2.set_ylabel("E")
ax2.grid(True)
ax2.legend()

plt.savefig("./bands_porownanie.png")


change_val=1

min_bands=[]
max_bands=[]

for band in data_RS.T:
    minE=np.min(band)
    maxE=np.max(band)
    min_bands.append(minE)
    max_bands.append(maxE)
    print(maxE-minE)
print()
for i in range(3):
    print(min_bands[i+1]-max_bands[i])





