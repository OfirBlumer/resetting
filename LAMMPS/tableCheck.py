import argparse
import numpy as np
import os
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("dir",nargs=1,type=str)

args = parser.parse_args()
dir = args.dir[0]


datalist = []
for filename in os.listdir(dir):
    if "FPT" in filename:
        file = os.path.join(dir,filename)
        with open(file, "r") as newfile:
            lines = newfile.readlines()
        params = eval(lines[0])
        mean = int(np.loadtxt(file,skiprows=1).mean())
        std = int(np.loadtxt(file,skiprows=1).std())
        N = len(np.loadtxt(file,skiprows=1))
        datalist.append(pd.DataFrame(dict(file=filename,mean=[mean],CV=std/mean,N=N,**params)))

print(pd.concat(datalist).drop(columns=["passageCoordinate"]))
