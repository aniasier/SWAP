import os
from loader import *
from plotter import *

def main():
    name = "RUN_V0_0.06"
    dir = f"/home/alina/Documents/SWAP-results/Time_evolution/120125/tests/{name}/"
    xyz = 2 # 0 -> x, 1 -> y, 2 -> z

    e1 = LoadEnergy('Energies1', dir)
    e2 = LoadEnergy('Energies2', dir)
    exp1 = LoadExpectations1('Expectations_1',dir) 
    exp2 = LoadExpectations2('Expectations_2',dir) 

    PlotEnergies1_Orbital(e1, exp1, name)
    PlotEnergies2_Orbital(e2, exp2, name)
    PlotEnergies1_Spin(e1, exp1, name)
    PlotEnergies2_Spin(e2, exp2, name)

    for n in range(1,3):
        psi_LR = LoadLRPsi(n,dir)
        PlotLRPsi(psi_LR,n, name)

    for n in range(1,11):
        psi_1 = LoadSingleElectronPsi(n,dir)
        PlotSingleElectronPsi(psi_1,n, name)

    potential = LoadPotential(dir)
    PlotPotential(potential, name)
    spin = LoadSpinTime(name,dir)  
    PlotSpinTime(spin, name, xyz)
    

if __name__ == "__main__":
    main()
