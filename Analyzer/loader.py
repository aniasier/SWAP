import os
import pandas as pd

def LoadEnergy(filename,dir):
    # dir = #"/home/ania/Documents/Doktorat/test_results/spektrum/RUN_Bz_1.0/"
    psiPath = os.path.join(
        dir + "OutputData", f"{filename}.dat"
    )
    #print(f'Loading Psi_1_n{n}.dat')
    if os.path.exists(psiPath):
        psi1 = pd.read_fwf(
            psiPath,
            skiprows=1,
            infer_nrows=100,
            names=[
                "n",
                "energy",
            ],
        )
    else:
        print("File does not exists, skipping: ", psiPath)
        
    return psi1

def LoadExpectations1(filename,dir):
    psiPath = os.path.join(
        dir + "OutputData", f"{filename}.dat"
    )
    if os.path.exists(psiPath):
        psi1 = pd.read_fwf(
            psiPath,
            skiprows=1,
            infer_nrows=100,
            names=[
                "n",
                "sx",
                "sy",
                "sz",
                "dxy_up",
                "dxy_down",
                "dxz_up",
                "dxz_down",
                "dyz_up",
                "dyz_down",
                "parity",
                "x",
                "y",
                "sz_L",
                "sz_R",
            ],
        )
    else:
        print("File does not exists, skipping: ", psiPath)
        
    return psi1
def LoadExpectations2(filename,dir):
    psiPath = os.path.join(
        dir + "OutputData", f"{filename}.dat"
    )
    if os.path.exists(psiPath):
        psi1 = pd.read_fwf(
            psiPath,
            skiprows=1,
            infer_nrows=100,
            names=[
                "n",
                "x",
                "sx",
                "sy",
                "sz",
                "dxy_up",
                "dxy_down",
                "dxz_up",
                "dxz_down",
                "dyz_up",
                "dyz_down",
                "parity",
                "sz_L",
                "sz_R",
            ],
        )
    else:
        print("File does not exists, skipping: ", psiPath)
        
    return psi1

def LoadLRPsi(n,dir):
    psiPath = os.path.join(
        dir + "OutputData", f"Psi_LR_n{n}.dat"
    )
    if os.path.exists(psiPath):
        psi1 = pd.read_fwf(
            psiPath,
            skiprows=1,
            infer_nrows=100,
            names=[
                "kx",
                "ky",
                "rePsi_xy_up",
                "imPsi_xy_up",
                "rePsi_xy_down",
                "imPsi_xy_down",
                "rePsi_xz_up",
                "imPsi_xz_up",
                "rePsi_xz_down",
                "imPsi_xz_down",
                "rePsi_yz_up",
                "imPsi_yz_up",
                "rePsi_yz_down",
                "imPsi_yz_down",
            ],
        )
    else:
        print("File does not exists, skipping: ", psiPath)
        
    return psi1

def LoadSingleElectronPsi(n,dir):
    psiPath = os.path.join(
        dir + "OutputData", f"Psi_1_n{n}.dat"
    )
    if os.path.exists(psiPath):
        psi1 = pd.read_fwf(
            psiPath,
            skiprows=1,
            infer_nrows=100,
            names=[
                "kx",
                "ky",
                "rePsi_xy_up",
                "imPsi_xy_up",
                "rePsi_xy_down",
                "imPsi_xy_down",
                "rePsi_xz_up",
                "imPsi_xz_up",
                "rePsi_xz_down",
                "imPsi_xz_down",
                "rePsi_yz_up",
                "imPsi_yz_up",
                "rePsi_yz_down",
                "imPsi_yz_down",
            ],
        )
    else:
        print("File does not exists, skipping: ", psiPath)
        
    return psi1

def LoadPotential(dir):
    psiPath = os.path.join(
        dir + "OutputData", f"potential.dat"
    )
    if os.path.exists(psiPath):
        psi1 = pd.read_fwf(
            psiPath,
            skiprows=1,
            infer_nrows=100,
            names=[
                "kx",
                "ky",
                "potential",
            ],
        )
    else:
        print("File does not exists, skipping: ", psiPath)
        
    return psi1


def LoadSpinDensity(dir):
    psiPath = os.path.join(
        dir + "OutputData", f"spin_densitytswitch.dat"
    )
    if os.path.exists(psiPath):
        psi1 = pd.read_fwf(
            psiPath,
            skiprows=1,
            infer_nrows=100,
            names=[
                "kx",
                "ky",
                "potential",
            ],
        )
    else:
        print("File does not exists, skipping: ", psiPath)
        
    return psi1
def LoadSpinTime(filename,dir):
    psiPath = os.path.join(
        dir + "OutputData/" + "Spin_time_evolution.dat"
    )
    if os.path.exists(psiPath):
        psi1 = pd.read_fwf(
            psiPath,
            skiprows=1,
            infer_nrows=10000,
            names=[
                "t",
                "sx_L",
                "sx_R",
                "sy_L",
                "sy_R",
                "sz_L",
                "sz_R",
            ],
        )
    else:
        print("File does not exists, skipping: ", psiPath)
        
    return psi1
