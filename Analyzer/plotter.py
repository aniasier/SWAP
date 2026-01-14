import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import os
import sys
import re

from matplotlib import rcParams
from matplotlib.cm import get_cmap
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.signal import find_peaks
from matplotlib.lines import Line2D

## Font
#rcParams['font.family'] = 'Times New Roman'
rcParams['font.family']= 'serif'
rcParams['font.serif']= ['Times New Roman', 'Times', 'DejaVu Serif']
rcParams['font.size'] = 14

## Lines
rcParams['lines.solid_joinstyle'] = 'miter'  # other options: 'round' or 'bevel'
rcParams['lines.antialiased'] = True  # turning on/off of antialiasing for sharper edges
rcParams['lines.linewidth'] = 1.25

## Legend
rcParams['legend.loc'] = 'upper left'
rcParams['legend.frameon'] = False

## Ticks
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = True
rcParams['ytick.right'] = True

rcParams['xtick.minor.visible'] = True
rcParams['ytick.minor.visible'] = True

## Resolution
rcParams['figure.dpi'] = 150

## Colors
### cmaps
cm_inferno = get_cmap("inferno")
cm_viridis = get_cmap("viridis")
cm_seismic = get_cmap("seismic")
cm_jet = get_cmap("jet")
cm_tab10 = get_cmap("tab10")
### Palettes from color-hex.com/
c_google = ['#008744', '#0057e7', '#d62d20', '#ffa700'] # G, B, R, Y # https://www.color-hex.com/color-palette/1872
c_twilight = ['#363b74', '#673888', '#ef4f91', '#c79dd7', '#4d1b7b'] # https://www.color-hex.com/color-palette/809

import numpy as np
import matplotlib.pyplot as plt
import os

def fill_maxwell(ax, corner_labels=('B', 'G', 'R')):
    Nlignes = 300
    Ncol = 300
    img = np.zeros((Nlignes, Ncol, 4))
    dx = 2.0 / (Ncol - 1)
    dy = 1.0 / (Nlignes - 1)
    for i in range(Ncol - 1):
        for j in range(Nlignes - 1):
            x = -1.0 + i * dx
            y = j * dy
            v = y
            r = (x + 1 - v) / 2.0
            b = 1.0 - v - r
            if 0 <= r <= 1 and 0 <= v <= 1 and 0 <= b <= 1:
                img[j][i] = np.array([r, v, b, 1.0])
            else:
                img[j][i] = np.array([1.0, 1.0, 1.0, 0.0])
    
    a = 1.0 / np.sqrt(3)
    ax.imshow(img, origin='lower', extent=[-a, a, 0.0, 1.0])
    ax.axis('off')

    pozycje = [(-a-0.3, 0), (a+0.3, 0), (0, 1+0.15)]
    for (x, y), label in zip(pozycje, corner_labels):
        ax.text(x, y, label, ha='center', va='center', fontsize=10, color='black')


def PlotEnergies1_Orbital(e1, exp1, name, output_folder="../Plots"):
    """
    Scatter plot energii z kolorami RGB zależnymi od orbitalów (dxy, dxz, dyz).
    Trójkąt Maxwella w dolnym prawym rogu, 2 razy mniejszy.
    """
    os.makedirs(output_folder, exist_ok=True)

    orbital_array = np.vstack((
        exp1['dxy_down'] + exp1['dxy_up'], 
        exp1['dxz_down'] + exp1['dxz_up'], 
        exp1['dyz_down'] + exp1['dyz_up']
    )).T

    orbital_array = np.clip(orbital_array, 0, 1)
    energies = e1.iloc[:, 1].values
    fig, ax_main = plt.subplots(figsize=(8, 6))

    # Scatter główny
    ax_main.scatter(np.arange(len(energies)), energies, color=orbital_array, s=20)
    ax_main.set_xlabel('Index')
    ax_main.set_ylabel('E [meV]')

    # Mały trójkąt Maxwella w dolnym prawym rogu (2x mniejszy)
    inset_ax = fig.add_axes([0.85, 0.15, 0.08, 0.08])  # left, bottom, width, height
    fill_maxwell(inset_ax, corner_labels=('d$_{yz}$', 'd$_{xy}$', 'd$_{zx}$'))

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"Energies1_dxdy_{name}.png"), dpi=300)
    plt.close()


def PlotEnergies1_Spectrum(e1_list, exp1_list, dso, name, output_folder="../Plots"):
    """
    Scatter energii vs DSO.
    e1_list – lista np.arrayów (po jednym na każdą wartość DSO)
    dso     – np.array z wartościami DSO (tej samej długości co e1_list)
    Kolory RGB zależne od orbitalów (dxy, dxz, dyz).
    """

    os.makedirs(output_folder, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    for i, dso_val in enumerate(dso):

        # --- energie ---
        energies = e1_list[i]

        # pandas DataFrame
        if hasattr(energies, "iloc"):
            energies = energies.iloc[:, 1].values

        # numpy array
        elif energies.ndim > 1:
            energies = energies[:, 1]

        # --- kolory orbitalne dla tego DSO ---
        exp1 = exp1_list[i]
        orbital_array = np.vstack((
            exp1['dxy_down'] + exp1['dxy_up'],
            exp1['dxz_down'] + exp1['dxz_up'],
            exp1['dyz_down'] + exp1['dyz_up']
        )).T

        orbital_array = np.clip(orbital_array, 0, 1)

        # --- oś X = stałe DSO ---
        x = np.full(len(energies), dso_val)

        ax.scatter(
            x,
            energies,
            color=orbital_array,
            s=20
        )

    ax.set_xlabel("DSO")
    ax.set_ylabel("E [meV]")

    # --- trójkąt Maxwella ---
    inset_ax = fig.add_axes([0.82, 0.15, 0.08, 0.08])
    fill_maxwell(inset_ax, corner_labels=('d$_{yz}$', 'd$_{xy}$', 'd$_{zx}$'))
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"Energies1_dxdy_DSO_{name}.png"), dpi=300)
    plt.close()

def PlotEnergies1_Spectrum_Spin(
    e1_list, exp1_list, dso, name, xyz, output_folder="../Plots"
):
    """
    Scatter energii vs DSO.
    Kolorowanie ciągłe po spinach:
    blue (-2) -> black (0) -> red (+2)
    """

    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    os.makedirs(output_folder, exist_ok=True)

    # --- colormap spinowy ---
    spin_cmap = LinearSegmentedColormap.from_list(
        "spin_bkr",
        ["blue", "black", "red"]
    )
    norm = Normalize(vmin=-2, vmax=2)

    fig, ax = plt.subplots(figsize=(8, 6))

    for i, dso_val in enumerate(dso):

        # --- energie ---
        energies = e1_list[i]
        if hasattr(energies, "iloc"):
            energies = energies.iloc[:, 1].values
        elif energies.ndim > 1:
            energies = energies[:, 1]

        # --- spiny ---
        exp1 = exp1_list[i]
        if xyz == 0:
            spin = exp1['sx']
        elif xyz == 1:
            spin = exp1['sy']
        elif xyz == 2:
            spin = exp1['sz']
        else:
            raise ValueError("xyz must be 0 (sx), 1 (sy) or 2 (sz)")

        # --- oś X = DSO ---
        x = np.full(len(energies), dso_val)

        sc = ax.scatter(
            x,
            energies,
            c=spin,
            cmap=spin_cmap,
            norm=norm,
            s=20
        )

    ax.set_xlabel("DSO")
    ax.set_ylabel("E [meV]")
    ax.grid(True)

    # --- colorbar ---
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(r"$\langle s_{%s} \rangle$" % ["x", "y", "z"][xyz])

    plt.tight_layout()
    plt.savefig(
        os.path.join(output_folder, f"Energies1_spin_cont_s{xyz}_DSO_{name}.png"),
        dpi=300
    )
    plt.close()



def PlotEnergies2_Orbital(e2, exp2, name, output_folder="../Plots"):

    os.makedirs(output_folder, exist_ok=True)

    orbital_array = np.vstack((
        exp2['dxy_down'] + exp2['dxy_up'], 
        exp2['dxz_down'] + exp2['dxz_up'], 
        exp2['dyz_down'] + exp2['dyz_up']
    )).T
    orbital_array = orbital_array / 2
    orbital_array = np.clip(orbital_array, 0, 1)
    energies = e2.iloc[:, 1].values

    fig, ax_main = plt.subplots(figsize=(8, 6))

    # Scatter główny
    ax_main.scatter(np.arange(len(energies)), energies, color=orbital_array, s=20)
    ax_main.set_xlabel('Index')
    ax_main.set_ylabel('E [meV]')

    # Mały trójkąt Maxwella w dolnym prawym rogu (2x mniejszy)
    inset_ax = fig.add_axes([0.85, 0.15, 0.08, 0.08])  # left, bottom, width, height
    fill_maxwell(inset_ax, corner_labels=('d$_{yz}$', 'd$_{xy}$', 'd$_{zx}$'))
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"Energies2_dxdy_{name}.png"), format='png', dpi=300)
    plt.close()

def PlotEnergies1_Spin(e1, exp1, name,xyz, output_folder="../Plots"):

    os.makedirs(output_folder, exist_ok=True)


    if(xyz==0):
        sz_down = np.where(exp1['sx']<=-0.99)
        sz_up = np.where(exp1['sx']>=.99)
    elif(xyz==1):
        sz_down = np.where(exp1['sy']<=-0.99)
        sz_up = np.where(exp1['sy']>=.99)
    elif(xyz==2):
        sz_down = np.where(exp1['sz']<=-0.99)
        sz_up = np.where(exp1['sz']>=.99)

    plt.scatter(sz_down[0], e1.iloc[sz_down[0], 1], color=c_google[1], label='down')
    plt.scatter(sz_up[0], e1.iloc[sz_up[0], 1], color=c_google[2],label='up')
    plt.xlabel('E [meV]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder,f"Energies1_s{xyz}_{name}.png"), format='png', dpi=300)
    plt.close()

def PlotEnergies2_Spin(e2, exp2, name,xyz, output_folder="../Plots"):

    os.makedirs(output_folder, exist_ok=True)
    eps = 1e-10
    if(xyz==0):
        sz_down = np.where(exp2['sx']<=-1.99)
        sz_up = np.where(exp2['sx']>=1.99)
        sz_zero = np.where(np.abs(exp2['sx']) < eps)
    elif(xyz==1):
        sz_down = np.where(exp2['sy']<=-1.99)
        sz_up = np.where(exp2['sy']>=1.99)
        sz_zero = np.where(np.abs(exp2['sy']) < eps)
    elif(xyz==2):
        sz_down = np.where(exp2['sz']<=-1.99)
        sz_up = np.where(exp2['sz']>=1.99)
        sz_zero = np.where(np.abs(exp2['sz']) < eps)


    plt.scatter(sz_down[0], e2.iloc[sz_down[0], 1], color=c_google[1], label='-2')
    plt.scatter(sz_up[0], e2.iloc[sz_up[0], 1], color=c_google[2],label='+2')
    plt.scatter(sz_zero[0], e2.iloc[sz_zero[0], 1], color='gray',label='0')
    plt.xlabel('E [meV]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder,f"Energies2_s{xyz}_{name}.png"), format='png', dpi=300)
    plt.close()

def PlotSingleElectronPsi(psi1,n, name, output_folder="../Plots"):
    os.makedirs(output_folder, exist_ok=True)
    # liczba orbitali (xy, xz, yz)
    norb = 3  
    
    density = np.zeros(len(psi1))
    for i in range(norb):
        re_up = psi1.iloc[:, 2 + i*4]
        im_up = psi1.iloc[:, 3 + i*4]
        re_down = psi1.iloc[:, 4 + i*4]
        im_down = psi1.iloc[:, 5 + i*4]
        density += re_up**2 + im_up**2 + re_down**2 + im_down**2

    plt.scatter(psi1["kx"], psi1["ky"], c=density, cmap="inferno", s=30)
    plt.colorbar(label=r"$|\psi(k_x, k_y)|^2$")
    plt.axis('equal')

    plt.xlabel(r"$k_x$")
    plt.ylabel(r"$k_y$")
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder,f"Psi_1_{n}_{name}.png"), format='png', dpi=300)
    plt.close()

def PlotLRPsi(psi1,n, name, output_folder="../Plots"):
    os.makedirs(output_folder, exist_ok=True)
    # liczba orbitali (xy, xz, yz)
    norb = 3  
    
    density = np.zeros(len(psi1))
    for i in range(norb):
        re_up = psi1.iloc[:, 2 + i*4]
        im_up = psi1.iloc[:, 3 + i*4]
        re_down = psi1.iloc[:, 4 + i*4]
        im_down = psi1.iloc[:, 5 + i*4]
        density += re_up**2 + im_up**2 + re_down**2 + im_down**2

    plt.scatter(psi1["kx"], psi1["ky"], c=density, cmap="inferno", s=30)
    plt.colorbar(label=r"$|\psi(k_x, k_y)|^2$")
    plt.axis('equal')

    plt.xlabel(r"$k_x$")
    plt.ylabel(r"$k_y$")
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder,f"Psi_LR_{n}_{name}.png"), format='png', dpi=300)
    plt.close()

def PlotPotential(potential, name, output_folder="../Plots"):
    os.makedirs(output_folder, exist_ok=True)
       
    plt.scatter(potential["kx"], potential["ky"], c=potential["potential"], cmap="inferno", s=30)
    plt.colorbar(label="V [eV]")
    plt.axis('equal')

    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.savefig(os.path.join(output_folder,f"potential_{name}.png"), format='png', dpi=300)
    plt.close()

def PlotSpinTime(spin, name,x, y,z, output_folder="../Plots"):
    os.makedirs(output_folder, exist_ok=True)

    if (x ==1):
        plt.plot(spin.iloc[:,0], spin.iloc[:,1], label='$S_x$ (L)')
        plt.plot(spin.iloc[:,0], spin.iloc[:,2], label='$S_x$ (R)')
        spin_time = spin.iloc[:,2].values
    if(y ==1):
        plt.plot(spin.iloc[:,0], spin.iloc[:,3], label='$S_y$ (L)')
        plt.plot(spin.iloc[:,0], spin.iloc[:,4], label='$S_y$ (R)')
        spin_time = spin.iloc[:,4].values
    if(z == 1):
        plt.plot(spin.iloc[:,0], spin.iloc[:,5], label='$S_z$ (L)')
        plt.plot(spin.iloc[:,0], spin.iloc[:,6], label='$S_z$ (R)')
        spin_time = spin.iloc[:,6].values
    xyz = x + y + z
    plt.ylim(-1,1)
    plt.xlim(spin.iloc[0,0],spin.iloc[-1,0])
    plt.xlabel('t [ns]')
    plt.ylabel('$<S>$ [$\hbar$/2]')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder,f"spin_swap{xyz}_{name}.png"), format='png', dpi=300)
    plt.close()

    peaks, _ = find_peaks(spin_time,prominence=0.1)
    if len(peaks) == 0:
        print("No peaks found")
        return
    print(peaks)
    first_peak = peaks[0]
    print('Value at maximum: ', spin_time[first_peak])
    print('Switching time: ', spin.iloc[:,0][first_peak])


def FindTime(spin, name,x, y,z,):
    if (x ==1):
        spin_time = spin.iloc[:,2].values
    if(y ==1):
        spin_time = spin.iloc[:,4].values
    if(z == 1):
        spin_time = spin.iloc[:,6].values
    peaks, _ = find_peaks(spin_time,prominence=0.1)
    if len(peaks) == 0:
        print("No peaks found")
        return np.nan
    # print(peaks)
    first_peak = peaks[0]
    # print(spin.iloc[:,0][first_peak])
    return spin.iloc[:,0][first_peak]


def PlotSwitchingTime(time, v0, output_folder="../Plots"):
    plt.plot(v0, time,'.')
    # plt.xlim(v0[0],v0[-1])
    plt.xlabel('V$_0$ [eV]')
    plt.ylabel('switching time')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder,f"switching_time.png"), format='png', dpi=300)
    plt.close()