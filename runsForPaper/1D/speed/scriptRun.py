import numpy as np
import os
import pandas as pd

for rate, seed in zip([0.0],[200000]):
    params = dict(resetRate=rate,
                  mass=40,
                  A = 0.59616/10000,
                  B = 0.59616,
                  C = 1,
                  startx=3,
                  starty=0,
                  passageCoordinate="x",
                  passageValue=-3,
                  passageSign="<",
                  round=1000,
                  seed=seed)
    seed += 10000

    fileSeed = str(np.random.uniform()).split('.')[1]
    np.random.seed(params["seed"])

    rounding = params["round"]

    fpts = f"{params}\n"
    std = np.sqrt((8.31445e-7 * 300 / params["mass"]))

    for i in range(10000):
        resetTimes = []
        if params["resetRate"] == 0:
            resetTimes.append(int(50000000 / rounding))
        else:
            nsteps = 0
            while nsteps < 50000000:
                newResetTime = np.random.exponential(1 / params["resetRate"])
                nsteps += newResetTime
                resetTimes.append(int(newResetTime / rounding))
        resetTimesString = ""
        VxString = ""
        for time in resetTimes:
            resetTimesString += f"{time} "
            VxString += f"{np.random.normal(0, std)} "

        V = np.random.normal(0, std)
        lammpsFile = "units        real\n" \
                     "atom_style    atomic\n" \
                     "dimension 2\n" \
                     "atom_modify map yes\n" \
                     "region        box block -1000 1000 -1000 1000 -0.1 0.1\n" \
                     "create_box    1 box\n" \
                     f"create_atoms    1 single {params['startx']} 0. 0.0\n" \
                     f"mass        1 {params['mass']}\n" \
                     f"velocity    all set {V} 0 0 sum yes\n" \
                     "pair_style      none\n" \
                     "fix        1 all nve\n" \
                     f"fix        2  all langevin 300 300 100.0 {params['seed'] + i}\n" \
                     f"variable myexp atom exp(-{params['C']}*x^2)\n" \
                     f"variable harm1d atom {params['A']}*x^2+{params['B']}*v_myexp\n" \
                     f"variable fx atom 2*x*({params['B']}*{params['C']}*v_myexp-{params['A']})\n" \
                     f"variable fy atom -300*y\n" \
                     f"fix harm all addforce v_fx v_fy 0.0 energy v_harm1d\n" \
                     "fix_modify harm energy yes\n" \
                     "fix        3 all enforce2d\n" \
                     f"variable resetTimes index {resetTimesString}\n" \
                     f"variable Vxs index {VxString}\n" \
                     f'variable pos equal "{params["passageCoordinate"]}[1]"\n' \
                     "label loop\n" \
                     f"variable a loop {len(resetTimes)}\n" \
                     "label innerLoop\n" \
                     "variable b loop ${resetTimes}\n" \
                     f"run {rounding}\n" \
                     'if "(${pos} ' \
                     f'{params["passageSign"]} {params["passageValue"]})" then &\n' \
                     '   "jump SELF break"\n' \
                     "next b\n" \
                     "jump SELF innerLoop\n" \
                     f"set atom 1 x {params['startx']} y 0 " \
                     "vx ${Vxs} vy 0\n" \
                     "next a\n" \
                     "next Vxs\n" \
                     "next resetTimes\n" \
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

        with open(f"../speedFPT{fileSeed}","w") as newFile:
            newFile.write(fpts)



