import numpy as np
import matplotlib.pyplot as plt

with open("1Dpotential/1/dump.xyz","r") as file:
    lines = file.readlines()
xs = []
ys = []
for i in range(2,len(lines),3):
    xs.append(float(lines[i].split()[1]))
    ys.append(float(lines[i].split()[2]))

H, xedges, yedges = np.histogram2d(xs, ys, bins=(20,200))
H = H.T
fes = -np.log(H)
plt.imshow(fes-fes.min(), interpolation='nearest',aspect="auto", origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
cbar = plt.colorbar()
plt.show()
plt.plot(xs)
plt.show()