import numpy as np
import os
import pandas as pd

for rate, seed in zip([0.0],[10000]):
    params = dict(resetRate=rate,
                  mass=40,
                  A1=28,
                  A2=28,
                  x01=2.5,
                  x02=-2.5,
                  y01=0.,
                  y02=0.,
                  alpha1=1.3,
                  alpha2=1.3,
                  lambda1=200,
                  lambda2=1,
                  Cx=4,
                  Cy = 0.005,
                  passageCoordinate="x",
                  passageValue=-1.,
                  passageSign="<",
                  startx=1.35,
                  starty=0,
                  round=1000,
                  seed=seed)

    fileSeed = str(np.random.uniform()).split('.')[1]

    np.random.seed(params["seed"])
    rounding = params["round"]
    alpha21 = params["alpha1"]**2
    alpha22 = params["alpha2"]**2
    lambda21 = params["lambda1"]**2
    lambda22 = params["lambda2"]**2

    fpts = f"{params}\n"

    pointsy = np.linspace(-40,40,101)
    pointsx = np.linspace(0,3,100)

    xs = []
    ys = []
    Vs = []
    levels = np.linspace(0,15,11)
    for x in pointsx:
        for y in pointsy:
            W_1 = params["A1"]*np.exp(-(x-params["x01"])**2/(2*params["alpha1"]**2)-(y-params["y01"])**2/(2*params["lambda1"]**2))
            W_2 = params["A2"]*np.exp(-(x-params["x02"])**2/(2*params["alpha2"]**2)-(y-params["y02"])**2/(2*params["lambda2"]**2))
            Vs.append((-(W_1 + W_2) + 4*x**2 + 0.005*y**2))
            xs.append(x)
            ys.append(y)

    surface = pd.DataFrame({"x":xs,"y":ys,"V":Vs})
    surface.V -= surface.V.min()
    surface["H"] = np.exp(-surface.V/0.59616)
    flat = surface.H
    flat = flat/sum(flat)
    std = np.sqrt((8.31445e-7*300 / params["mass"]))

    for i in range(10000):
        resetTimes = []
        if params["resetRate"] == 0:
            resetTimes.append(int(100000000 / rounding))
        else:
            nsteps = 0
            while nsteps < 100000000:
                newResetTime = np.random.exponential(1 / params["resetRate"])
                nsteps += newResetTime
                resetTimes.append(int(newResetTime / rounding))
        VxString = ""
        VyString = ""
        PxString = ""
        PyString = ""
        resetTimesString = ""
        for time in resetTimes:
            resetTimesString += f"{time} "
            sample_index = np.random.choice(a=flat.size, p=flat)
            PxString += f"{surface.x[sample_index]} "
            PyString += f"{surface.y[sample_index]} "
            V = np.sqrt(np.random.normal(0, std) ** 2 + np.random.normal(0, std) ** 2 + np.random.normal(0, std) ** 2)
            phi = np.random.uniform() * 2 * np.pi
            VxString += f"{np.cos(phi) * V} "
            VyString += f"{np.sin(phi) * V} "
        sample_index = np.random.choice(a=flat.size, p=flat)
        V = np.sqrt(np.random.normal(0, std) ** 2 + np.random.normal(0, std) ** 2 + np.random.normal(0, std) ** 2)
        phi = np.random.uniform() * 2 * np.pi
        lammpsFile = "units        real\n"\
                     "atom_style    atomic\n"\
                     "dimension 2\n"\
                     "atom_modify map yes\n"\
                     "region        box block -1000 1000 -1000 1000 -0.1 0.1\n"\
                     "create_box    1 box\n"\
                     f"create_atoms    1 single {surface.x[sample_index]} {surface.y[sample_index]} 0.0\n"\
                     f"mass        1 {params['mass']}\n"\
                     f"velocity    all set {np.cos(phi) * V} {np.sin(phi) * V} 0 sum yes\n" \
                     "pair_style      none\n" \
                     "fix        1 all nve\n"\
                     f"fix        2  all langevin 300 300 100.0 {params['seed']+i}\n" \
                     f"variable W_1 atom {params['A1']}*exp(-((x-{params['x01']})^2)/(2*{alpha21})-((y-{params['y01']})^2)/(2*{lambda21}))\n" \
                     f"variable W_2 atom {params['A2']}*exp(-((x-{params['x02']})^2)/(2*{alpha22})-((y-{params['y02']})^2)/(2*{lambda22}))\n" \
                     f"variable harm1d atom (-(v_W_1+v_W_2)+{params['Cx']}*x^2+{params['Cy']}*y^2)\n"\
                     f"variable fx atom (-(x-{params['x01']})*v_W_1/{alpha21}-(x-{params['x02']})*v_W_2/{alpha22}-2*{params['Cx']}*x)\n"\
                     f"variable fy atom (-(y-{params['y01']})*v_W_1/{lambda21}-(y-{params['y02']})*v_W_2/{lambda22}-2*{params['Cy']}*y)\n"\
                     "fix harm all addforce v_fx v_fy 0.0 energy v_harm1d\n"\
                     "fix_modify harm energy yes\n"\
                     "fix        3 all enforce2d\n"\
                     f"variable resetTimes index {resetTimesString}\n" \
                     f"variable Pxs index {PxString}\n" \
                     f"variable Pys index {PyString}\n" \
                     f"variable Vxs index {VxString}\n" \
                     f"variable Vys index {VyString}\n" \
                     f'variable pos equal "{params["passageCoordinate"]}[1]"\n'\
                     "label loop\n"\
                     f"variable a loop {len(resetTimes)}\n"\
                     "label innerLoop\n"\
                     "variable b loop ${resetTimes}\n"\
                     f"run {rounding}\n"\
                     'if "(${pos} '\
                     f'{params["passageSign"]} {params["passageValue"]})" then &\n'\
                     '   "jump SELF break"\n'\
                     "next b\n"\
                     "jump SELF innerLoop\n"\
                     "set atom 1 x ${Pxs} y ${Pys} vx ${Vxs} vy ${Vys}\n" \
                     "next a\n"\
                     "next resetTimes\n" \
                     "next Pxs\n" \
                     "next Pys\n" \
                     "next Vxs\n" \
                     "next Vys\n" \
                     "jump SELF loop\n"\
                     "label break\n"\
                     'print "ALL DONE"\n' \
                     'write_dump all xyz dump.atom'

        with open("newRun.lmp", "w") as newFile:
            newFile.write(lammpsFile)

        os.system("mpirun -np 2 /home/ofirblumer/mylammps/src/lmp_g++_openmpi -in newRun.lmp")

        with open("dump.atom", "r") as file:
            lines = file.readlines()
        fpts += f"{lines[1].split()[-1]}\n"

        with open(f"../bothFPT{fileSeed}","w") as newFile:
            newFile.write(fpts)
