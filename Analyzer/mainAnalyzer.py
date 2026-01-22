import os
from loader import *
from plotter import *

def main():
    name = "RUN_Bz_0.01"
    dir = f"/home/alina/Documents/SWAP-results/220126/small/period/{name}/"
    # dir = f"/home/alina/Documents/SWAP/"
    x = 0 # 0 -> x, 1 -> y, 2 -> z
    y = 0
    z = 1
    xyz = 2

    # Bz = np.linspace(0.0,20.0,21)
    # v = np.linspace(10.0e-3,100.0e-3,10)
    # print(v)
    
    # # e1_all = []
    # # exp1_all =[]
    # # dso = np.array([0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1])
    # # namelist = []
    # # time = []
    # # for i in v:
    # #     i = np.floor(i * 100 + 1e-12) / 100

    # #     dir1 = f"/home/alina/Documents/SWAP-results/140126/V0_spectrum_SO/RUN_V0_{i}/"
    # #     e1_all.append(LoadEnergy('Energies2', dir1))
    # #     exp1_all.append(LoadExpectations2('Expectations_2',dir1))
    # #     # spin = spin = LoadSpinTime(name,dir1)
    # #     # time.append(FindTime(spin, name, x, y,z))
    # # print(time)
    # # PlotEnergyGap_S_minus2_to_0(e1_all, exp1_all, v, name, xyz=2, output_folder="../Plots")

    # # PlotEnergies1_Spectrum_Spin(e1_all, exp1_all, Bz, name, xyz,output_folder="../Plots")
    # # PlotSwitchingTime(time,Bz)
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
    # spin_density = LoadSpinDensity(dir)
    # PlotSpinDensity(spin_density, name, output_folder=f"{dir}/Plots")
    spin = LoadSpinTime(name,dir)  
    PlotSpinTime(spin, name, x, y,z, output_folder=f"{dir}/Plots")
    

if __name__ == "__main__":
    main()
