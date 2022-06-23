import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import scipy as scp
from scipy import optimize
from dataManager import dataManager
from itertools import product
import os
import seaborn as sns
myData = dataManager("resetting","simplePotentials1")

def pade4_2(x,A,U1,D1,U2,D2,D3,D4):
    return (U2*x**2+x*U1+A)/(D4*x**4+D3*x**3+D2*x**2+x*D1+A)

def pade5_3(x,A,U1,D1,U2,D2,U3,D3,D4,D5):
    return (U3*x**3+U2*x**2+x*U1+A)/(D5*x**5+D4*x**4+D3*x**3+D2*x**2+x*D1+A)

def pade6_4(x,A,U1,D1,U2,D2,U3,D3,U4,D4,D5,D6):
    return (U4*x**4+U3*x**3+U2*x**2+x*U1+A)/(D6*x**6+D5*x**5+D4*x**4+D3*x**3+D2*x**2+x*D1+A)

def moonWalkPrediction(func,values,rate,stepSize,forwardStep=1e-6,forwardN=100):
    zs = np.linspace(0,forwardStep,forwardN)
    mus = []
    for z in zs:
        Tr = sum(np.exp(-z*values))/len(values)
        mus.append((1-Tr)/(z*Tr))
    mus[0] = np.mean(values)
    for i in range(int(rate/stepSize)):
        zs = zs + stepSize
        Trs = 1/(mus*zs+1)
        zs = np.insert(zs,0,0)
        Trs = np.insert(Trs,0,1)
        try:
            fit = scp.optimize.curve_fit(func, zs, Trs,bounds=(0,np.inf))
            mus.insert(0,(fit[0][2]-fit[0][1])/fit[0][0])
        except:
            mus.insert(0,mus[0])

    zs = zs + rate - stepSize*int(rate/stepSize)
    Trs = 1/(mus*zs+1)
    zs = np.insert(zs,0,0)
    Trs = np.insert(Trs,0,1)
    try:
        fit = scp.optimize.curve_fit(func, zs, Trs,bounds=(0,np.inf))
        mus.insert(0,(fit[0][2]-fit[0][1])/fit[0][0])
    except:
        mus.insert(0,mus[0])
    return zs, mus

for doc in list(myData.collection.find({"N":10000})):
    if doc["restartRate"]!=0 and doc["fileSeed"]!='46717817747070256':
        newDocParams = {}
        for key in doc.keys():
            if key not in ["_id","path","N","dataType"]:
                newDocParams[key] = doc[key]
        print(newDocParams)
        fpts = np.loadtxt(doc["path"],skiprows=1)
        datalist = []
        Nfptss = [1000,5000,10000]
        forwardSteps = [1e-5,1e-4,1e-3]
        forwardNs = [50,500,5000]
        stepSizes = [5e-8,5e-7,5e-6]
        count = 10
        for Nfpts, stepSize, forwardStep, forwardN in product(Nfptss,stepSizes,forwardSteps,forwardNs):
            count += 1
            print(Nfpts,stepSize, forwardStep, forwardN)
            for func, funcname in zip([pade4_2,pade5_3,pade6_4],["pade4_2","pade5_3","pade6_4"]):
                zs, mus = moonWalkPrediction(func=func,values=fpts[:Nfpts],rate=doc["restartRate"],stepSize=stepSize,forwardStep=forwardStep,forwardN=forwardN)
                newD = pd.DataFrame(dict(Nfpts=Nfpts,stepSize=[stepSize],forwardStep=forwardStep,forwardN=forwardN,func=funcname,prediction=mus[0]))
                datalist.append(newD)
        data = pd.concat(datalist,ignore_index=True)
        data.sort_values("prediction")
        filenumber = str(np.random.uniform()).split(".")[1]
        data.to_csv(f"predictions/pr{filenumber}")
        fileDict = dict(path=os.path.abspath(f"predictions/pr{filenumber}"),dataType="predictions",**newDocParams)
        myData.collection.insert_one(fileDict)