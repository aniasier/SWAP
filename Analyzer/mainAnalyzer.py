import os
from loader import *
from plotter import *

def main():
    name = "RUN_Bz_0.01"
    dir = f"/home/alina/Documents/SWAP-results/Time_evolution/130126/wyniki/noSO/{name}/"
    x = 0 # 0 -> x, 1 -> y, 2 -> z
    y = 0
    z = 1
    xyz = 2

    e1 = LoadEnergy('Energies1', dir)
    e2 = LoadEnergy('Energies2', dir)
    exp1 = LoadExpectations1('Expectations_1',dir) 
    exp2 = LoadExpectations2('Expectations_2',dir) 

    PlotEnergies1_Orbital(e1, exp1, name, output_folder=f"{dir}/Plots")
    PlotEnergies2_Orbital(e2, exp2, name, output_folder=f"{dir}Plots")
    PlotEnergies1_Spin(e1, exp1, name,xyz, output_folder=f"{dir}.Plots")
    PlotEnergies2_Spin(e2, exp2, name,xyz, output_folder=f"{dir}/Plots")

    for n in range(1,3):
        psi_LR = LoadLRPsi(n,dir)
        PlotLRPsi(psi_LR,n, name, output_folder=f"{dir}/Plots")

    for n in range(1,11):
        psi_1 = LoadSingleElectronPsi(n,dir)
        PlotSingleElectronPsi(psi_1,n, name, output_folder=f"{dir}/Plots")

    potential = LoadPotential(dir)
    PlotPotential(potential, name, output_folder=f"{dir}/Plots")
    spin = LoadSpinTime(name,dir)  
    PlotSpinTime(spin, name, x, y,z, output_folder=f"{dir}/Plots")
    

if __name__ == "__main__":
    main()
