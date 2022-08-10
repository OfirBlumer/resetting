import numpy as np
import os
import pandas as pd

for rate, seed in zip([0.0],[15000]):
    params = dict(resetRate=rate,
              mass=40,
              height=0.59616*1.5,
              Cx = 2.,
              Cx4 = 1,
              xscale=15,
              yscale=1.,
              xslope=1.5,
              yslope=1.2,
              startx=-14.9,
              starty=-1.4,
              passageCoordinate="y",
              passageValue=1,
              passageSign=">",
              round=1000,
              seed=10000)

    fileSeed = str(np.random.uniform()).split('.')[1]

    np.random.seed(params["seed"])
    rounding = params["round"]

    fpts = f"{params}\n"

    pointsy = np.linspace(0,-2.5,100)
    pointsx = np.linspace(-2.5,2.5,100)

    xs = []
    ys = []
    Vs = []
    for xt in pointsx:
        x = xt/params["xscale"]
        for yt in pointsy:
            y = yt/params["yscale"]
            Vs.append(params["height"]*(params["Cx4"]*x**4+y**4-params["Cx"]*x**2-4*y**2+x*y+params["xslope"]*x+params["yslope"]*y))
            xs.append(xt)
            ys.append(yt)

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
            VxString += f"{np.random.normal(0, std)} "
            VyString += f"{np.random.normal(0, std)} "
        sample_index = np.random.choice(a=flat.size, p=flat)
        lammpsFile = "units        real\n"\
                     "atom_style    atomic\n"\
                     "dimension 2\n"\
                     "atom_modify map yes\n"\
                     "region        box block -1000 1000 -1000 1000 -0.1 0.1\n"\
                     "create_box    1 box\n"\
                     f"create_atoms    1 single {surface.x[sample_index]} {surface.y[sample_index]} 0.0\n"\
                     f"mass        1 {params['mass']}\n"\
                     f"velocity    all set {np.random.normal(0, std)} {np.random.normal(0, std)} 0 sum yes\n" \
                     "pair_style      none\n" \
                     "fix        1 all nve\n"\
                     f"fix        2  all langevin 300 300 100.0 {params['seed']+i}\n" \
                     f"variable harm1d atom {params['height']}*({params['Cx4']}*(x/{params['xscale']})^4+(y/{params['yscale']})^4-{params['Cx']}*(x/{params['xscale']})^2"\
                     f"-4*(y/{params['yscale']})^2+(x/{params['xscale']})*(y/{params['yscale']})+" \
                     f"{params['xslope']}*(x/{params['xscale']})+{params['yslope']}*(y/{params['yscale']}))\n"\
                     f"variable fx atom ({params['height']}/{params['xscale']})*(-4*{params['Cx4']}*(x/{params['xscale']})^3+2*{params['Cx']}*(x/{params['xscale']})-(y/{params['yscale']})-{params['xslope']})\n"\
                     f"variable fy atom ({params['height']}/{params['yscale']})*(-4*(y/{params['yscale']})^3+8*(y/{params['yscale']})-(x/{params['xscale']})-{params['yslope']})\n"\
                     f"fix harm all addforce v_fx v_fy 0.0 energy v_harm1d\n"\
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
