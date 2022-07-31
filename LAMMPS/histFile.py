import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("filename",nargs=1,type=str)
args = parser.parse_args()
filename = args.filename[0]

fpts = np.loadtxt(filename,skiprows=1)
logbins = np.logspace(0,np.log10(max(fpts)),50)
xs = 0.5*(logbins[1:] + logbins[:-1])
c,b = np.histogram(fpts, bins=logbins, density=True)
plt.plot(xs/1e6,c,c="blue")
plt.xscale("log")
plt.ylabel("Probability")
plt.xlabel("FPT (ns)")
plt.show()