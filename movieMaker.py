import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import numpy as np

xs = np.linspace(-150,150,200)
ys = np.linspace(-50,50,200)
kbT = 0.0003
Cxy = 1
Cx = 4
Cy = 1
levels=[i for i in range(16)]
datalist = []

for xtag in xs:
    x = xtag/100
    newRow = []
    for ytag in ys:
        y = ytag/15
        newRow.append(kbT*(x**4+y**4-2*x**2-4*y**2+x*y+0.3*x+0.1*y))
    datalist.append(pd.DataFrame({"x":xtag,"y":ys,"V":newRow}))
surface = pd.concat(datalist).reset_index()
surface.V -= surface.V.min()
fig,ax = plt.subplots()
ax.tricontourf(surface.x,surface.y,surface.V/0.00025,levels=levels)#levels=[surface.V.min()+i*(surface.V.max()-surface.V.min())/100 for i in range(100)])
ax.set_xlim(-150,150)
ax.set_ylim(-50,50)

with open("positionsResults/ps3539959879060586","r") as file:
    lines = file.readlines()
xs = []
ys = []

for line in lines[1:409]:
    xs.append(float(line.split()[0]))
    ys.append(float(line.split()[1]))
for i in range(30):
    xs.append(xs[-1])
    ys.append(ys[-1])

for line in lines[409:]:
    xs.append(float(line.split()[0]))
    ys.append(float(line.split()[1]))
for i in range(30):
    xs.append(xs[-1])
    ys.append(ys[-1])

positionSeed = 0
scatter=ax.scatter(xs[0],ys[0],edgecolors='black',c="white")
count = 1
def update(frame_number):

   scatter.set_offsets([xs[0],ys[0]])
   xs.pop(0)
   ys.pop(0)
   return scatter,

anim = FuncAnimation(fig, update,save_count=1000)
anim.save('WolfeQuapp.gif', writer='PillowWriter',fps=30)
# plt.plot(xs)
# plt.show()