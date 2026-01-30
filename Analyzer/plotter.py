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
from matplotlib.colors import ListedColormap

## Font
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
rcParams['font.size'] = 9

rcParams['axes.labelsize'] = 9
rcParams['axes.titlesize'] = 9
rcParams['legend.fontsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8

## Lines
rcParams['lines.linewidth'] = 1.1
rcParams['lines.solid_joinstyle'] = 'miter'
rcParams['lines.antialiased'] = True
rcParams['lines.markersize'] = 4

## Axes
rcParams['axes.linewidth'] = 0.8

## Legend
rcParams['legend.frameon'] = False
rcParams['legend.loc'] = 'best'

## Ticks
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = True
rcParams['ytick.right'] = True

rcParams['xtick.major.size'] = 4
rcParams['ytick.major.size'] = 4
rcParams['xtick.minor.size'] = 2
rcParams['ytick.minor.size'] = 2

rcParams['xtick.major.width'] = 0.8
rcParams['ytick.major.width'] = 0.8
rcParams['xtick.minor.width'] = 0.6
rcParams['ytick.minor.width'] = 0.6

rcParams['xtick.minor.visible'] = True
rcParams['ytick.minor.visible'] = True

## Figure
rcParams['figure.figsize'] = (3.35, 2.5)
rcParams['figure.dpi'] = 300
rcParams['savefig.dpi'] = 300

## Colormaps
cm_inferno = get_cmap("inferno")
cm_plasma = get_cmap("plasma")
cm_viridis = get_cmap("viridis")
cm_seismic = get_cmap("seismic")
cm_tab10 = get_cmap("tab10")

### Palettes from color-hex.com/ 
c_google = ['#008744', '#0057e7', '#d62d20', '#ffa700'] # G, B, R, Y # https://www.color-hex.com/color-palette/1872 
c_twilight = ['#363b74', '#673888', '#ef4f91', '#c79dd7', '#4d1b7b'] # https://www.color-hex.com/color-palette/809
c_palette = ["#780000","#c1121f","#fdf0d5","#003049","#669bbc"]
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

    pozycje = [(-a-0.4, 0), (a+0.4, 0), (0, 1+0.25)]
    for (x, y), label in zip(pozycje, corner_labels):
        ax.text(x, y, label, ha='center', va='center', fontsize=8, color='black')


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
    fig, ax_main = plt.subplots()

    # Scatter główny
    ax_main.scatter(np.arange(len(energies)), energies, color=orbital_array, s=20)
    ax_main.set_xlabel('Index')
    ax_main.set_ylabel('E [meV]')

    # Mały trójkąt Maxwella w dolnym prawym rogu (2x mniejszy)
    inset_ax = fig.add_axes([0.82, 0.25, 0.08, 0.08])  # left, bottom, width, height
    fill_maxwell(inset_ax, corner_labels=('d$_{yz}$', 'd$_{xy}$', 'd$_{zx}$'))

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"Energies1_dxdy_{name}.png"))
    plt.close()


def PlotEnergies1_Spectrum(e1_list, exp1_list, dso, name, output_folder="../Plots"):
    """
    Scatter energii vs DSO.
    e1_list – lista np.arrayów (po jednym na każdą wartość DSO)
    dso     – np.array z wartościami DSO (tej samej długości co e1_list)
    Kolory RGB zależne od orbitalów (dxy, dxz, dyz).
    """

    os.makedirs(output_folder, exist_ok=True)

    fig, ax = plt.subplots()

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
    plt.savefig(os.path.join(output_folder, f"Energies1_dxdy_DSO_{name}.png"))
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
        ["tab:blue", "black", "tab:red"]
    )
    norm = Normalize(vmin=-2, vmax=2)

    fig, ax = plt.subplots()

    # --- collect energies and spins into matrices (shape: n_dso x n_levels)
    n_dso = len(dso)
    # determine number of levels from first entry
    first_energies = e1_list[0]
    if hasattr(first_energies, "iloc"):
        first_arr = first_energies.iloc[:, 1].values
    elif getattr(first_energies, 'ndim', 1) > 1:
        first_arr = first_energies[:, 1]
    else:
        first_arr = np.asarray(first_energies).ravel()

    n_levels = len(first_arr)

    energies_mat = np.zeros((n_dso, n_levels))
    spins_mat = np.zeros((n_dso, n_levels))

    for i in range(n_dso):
        energies = e1_list[i]
        if hasattr(energies, "iloc"):
            energies = energies.iloc[:, 1].values
        elif getattr(energies, 'ndim', 1) > 1:
            energies = energies[:, 1]
        energies_mat[i, :] = energies

        exp1 = exp1_list[i]
        if xyz == 0:
            spin = np.asarray(exp1['sx'])
        elif xyz == 1:
            spin = np.asarray(exp1['sy'])
        elif xyz == 2:
            spin = np.asarray(exp1['sz'])
        else:
            raise ValueError("xyz must be 0 (sx), 1 (sy) or 2 (sz)")
        spins_mat[i, :] = spin

    from matplotlib.collections import LineCollection

    x = np.array(dso)

    if n_dso < 2:
        # nothing to connect: fallback to single scatter
        for lvl in range(n_levels):
            ax.scatter(x, energies_mat[:, lvl], c=spins_mat[:, lvl], cmap=spin_cmap, norm=norm, s=20)
        sc_for_cb = None
    else:
        sc_for_cb = None
        # for each energy level (index), draw a LineCollection whose segments
        # are colored according to the spin values along DSO
        for lvl in range(n_levels):
            y = energies_mat[:, lvl]
            spin_vals = spins_mat[:, lvl]
            points = np.column_stack([x, y])
            segments = np.stack([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap=spin_cmap, norm=norm, linewidths=0.9)
            # color segments by spin at the left endpoint of each segment
            lc.set_array(spin_vals[:-1])
            ax.add_collection(lc)
            sc_for_cb = lc

    ax.set_xlabel("B [T]")
    ax.set_ylabel("E [meV]")
    # ax.grid(True)

    # --- colorbar ---
    if sc_for_cb is not None:
        cbar = plt.colorbar(sc_for_cb, ax=ax)
        cbar.set_label(r"$\langle s_{%s} \rangle$" % ["x", "y", "z"][xyz])

    ax.set_xlim(np.min(x), np.max(x))
    ax.autoscale_view()
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_folder, f"Energies1_spin_cont_s{xyz}_DSO_{name}.png")
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

    fig, ax_main = plt.subplots()

    # Scatter główny
    ax_main.scatter(np.arange(len(energies)), energies, color=orbital_array, s=20)
    ax_main.set_xlabel('Index')
    ax_main.set_ylabel('E [meV]')

    # Mały trójkąt Maxwella w dolnym prawym rogu (2x mniejszy)
    inset_ax = fig.add_axes([0.82, 0.25, 0.08, 0.08])  # left, bottom, width, height
    fill_maxwell(inset_ax, corner_labels=('d$_{yz}$', 'd$_{xy}$', 'd$_{zx}$'))
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"Energies2_dxdy_{name}.png"), format='png')
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
    plt.savefig(os.path.join(output_folder,f"Energies1_s{xyz}_{name}.png"), format='png')
    plt.close()

def PlotEnergies2_Spin(e2, exp2, name,xyz, output_folder="../Plots"):

    os.makedirs(output_folder, exist_ok=True)
    eps = 1e-10
    if(xyz==0):
        sz_down = np.where(exp2['sx']<=-1.9)
        sz_up = np.where(exp2['sx']>=1.9)
        sz_zero = np.where(np.abs(exp2['sx']) < eps)
    elif(xyz==1):
        sz_down = np.where(exp2['sy']<=-1.9)
        sz_up = np.where(exp2['sy']>=1.9)
        sz_zero = np.where(np.abs(exp2['sy']) < eps)
    elif(xyz==2):
        sz_down = np.where(exp2['sz']<=-1.9)
        sz_up = np.where(exp2['sz']>=1.9)
        sz_zero = np.where(np.abs(exp2['sz']) < eps)


    plt.scatter(sz_down[0], e2.iloc[sz_down[0], 1], color=c_google[1], label='-2')
    plt.scatter(sz_up[0], e2.iloc[sz_up[0], 1], color=c_google[2],label='+2')
    plt.scatter(sz_zero[0], e2.iloc[sz_zero[0], 1], color='gray',label='0')
    plt.xlabel('E [meV]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder,f"Energies2_s{xyz}_{name}.png"), format='png')
    plt.close()

def PlotEnergyGap_S_minus2_to_0(e_list, exp_list, dso, name, xyz=2, output_folder="../Plots"):
    """
    For each value in `dso`, find the first two states with spin approximately 
    zero (|spin| < eps_zero), and plot their energy difference vs DSO.
    """

    os.makedirs(output_folder, exist_ok=True)

    gaps = []
    dso_vals = []
    eps_zero = 1e-4
    for i, d in enumerate(dso):
        # extract energies similar to other functions
        energies = e_list[i]
        if hasattr(energies, "iloc"):
            e_arr = energies.iloc[:, 1].values
        elif getattr(energies, 'ndim', 1) > 1:
            e_arr = energies[:, 1]
        else:
            e_arr = np.asarray(energies).ravel()

        exp = exp_list[i]
        if xyz == 0:
            spin = np.asarray(exp['sx'])
        elif xyz == 1:
            spin = np.asarray(exp['sy'])
        elif xyz == 2:
            spin = np.asarray(exp['sz'])
        else:
            raise ValueError("xyz must be 0 (sx), 1 (sy) or 2 (sz)")

        # find indices for spin ~ 0
        idx_zero = np.where(np.abs(spin) < eps_zero)[0]

        if idx_zero.size < 2:
            # need at least 2 zero-spin states
            gaps.append(np.nan)
            dso_vals.append(d)
            continue

        # first two zero-spin energies
        idx_1 = idx_zero[0]
        idx_2 = idx_zero[1]
        e_1 = e_arr[idx_1]
        e_2 = e_arr[idx_2]

        gap = np.abs(e_2 - e_1)
        gaps.append(gap)
        dso_vals.append(d)

    # convert to arrays
    dso_vals = np.array(dso_vals)
    gaps = np.array(gaps)

    # plot
    fig, ax = plt.subplots()
    ax.plot(dso_vals, gaps, linestyle='-')
    # ax.set_yscale('log')
    ax.set_xlabel('$V_0$ [eV]')
    ax.set_ylabel(r'J [meV]')
    ax.set_xlim(np.min(dso_vals), np.max(dso_vals))
    ax.set_ylim(np.min(gaps), np.max(gaps))
    # ax.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"gap_zero_spin_s{xyz}_DSO_{name}.png"))
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
    # rectangular figure with 2:1 width:height
    fig, ax = plt.subplots()

    sc = ax.scatter(psi1["kx"], psi1["ky"], c=density, cmap="inferno", s=30)
    # fig.colorbar(sc, ax=ax, label=r"$|\psi(k_x, k_y)|^2$")
    # set limits exactly to data bounds to avoid white margins
    kx = psi1["kx"].values
    ky = psi1["ky"].values
    ax.set_xlim(np.min(kx), np.max(kx))
    ax.set_ylim(np.min(ky), np.max(ky))

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(r"$x$ [nm]")
    ax.set_ylabel(r"$y$ [nm]")
    fig.tight_layout()
    fig.savefig(os.path.join(output_folder, f"Psi_1_{n}_{name}.png"), format='png')
    plt.close(fig)

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
    # rectangular figure with 2:1 width:height
    fig, ax = plt.subplots()

    sc = ax.scatter(psi1["kx"], psi1["ky"], c=density, cmap="inferno", s=30)
    # fig.colorbar(sc, ax=ax, label=r"$|\psi(k_x, k_y)|^2$")
    # set limits exactly to data bounds to avoid white margins
    kx = psi1["kx"].values
    ky = psi1["ky"].values
    ax.set_xlim(np.min(kx), np.max(kx))
    ax.set_ylim(np.min(ky), np.max(ky))

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(r"$x$ [nm]")
    ax.set_ylabel(r"$y$ [nm]")
    fig.tight_layout()
    fig.savefig(os.path.join(output_folder, f"Psi_LR_{n}_{name}.png"), format='png')
    plt.close(fig)

def PlotPotential(potential, name, output_folder="../Plots"):
    os.makedirs(output_folder, exist_ok=True)
    # Use a rectangular figure with 2:1 width:height ratio (keeps similar overall size)
    fig, ax = plt.subplots()

    sc = ax.scatter(potential["kx"], potential["ky"], c=potential["potential"], cmap="inferno", s=30)
    fig.colorbar(sc, ax=ax, label="V [eV]")

    # set limits exactly to data bounds to avoid white margins
    kx = np.asarray(potential["kx"])
    ky = np.asarray(potential["ky"])
    ax.set_xlim(np.min(kx), np.max(kx))
    ax.set_ylim(np.min(ky), np.max(ky))

    # keep data aspect equal but the figure itself is 2:1
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel(r"$x$ [nm]")
    ax.set_ylabel(r"$y$ [nm]")
    fig.tight_layout()
    fig.savefig(os.path.join(output_folder, f"potential_{name}.png"), format='png')
    plt.close(fig)

    # --- profile at ky == 0: plot V(kx) ---
    try:
        ky_arr = np.asarray(potential["ky"])
        mask = np.isclose(ky_arr, 0.0, atol=1e-8)
    except Exception:
        mask = None

    if mask is None or not np.any(mask):
        # nothing exactly at ky==0 -- skip profile (user can adjust tolerance if needed)
        print(f"PlotPotential: no entries with ky==0 for profile_potential_{name}; skipping profile plot")
    else:
        kx_profile = np.asarray(potential["kx"])[mask]
        V_profile = np.asarray(potential["potential"])[mask]
        # sort by kx
        order = np.argsort(kx_profile)
        kx_profile = kx_profile[order]
        V_profile = V_profile[order]

        fig2, ax2 = plt.subplots(figsize=(8, 4))
        ax2.plot(kx_profile, V_profile, linestyle='-')
        ax2.set_xlabel(r"$x$ [nm]")
        ax2.set_ylabel('V [eV]')
        ax2.set_xlim(np.min(kx_profile), np.max(kx_profile))
        ax2.set_ylim(np.min(V_profile), np.max(V_profile))
        ax2.grid(False)
        fig2.tight_layout()
        fig2.savefig(os.path.join(output_folder, f"profile_potential_{name}.png"))
        plt.close(fig2)
    

def PlotSpinTime(spin, name,x, y,z, output_folder="../Plots"):
    os.makedirs(output_folder, exist_ok=True)

    # Determine number of lines to plot
    num_lines = (x + y + z) * 2
    
    # Select visually appealing colors from inferno
    if num_lines == 2:
        color_indices = [0.2, 0.85]
    elif num_lines == 4:
        color_indices = [0.15, 0.4, 0.65, 0.9]
    else:  # num_lines == 6
        color_indices = [0.1, 0.3, 0.5, 0.7, 0.85, 0.95]
    
    colors = [cm_plasma(idx) for idx in color_indices]
    color_idx = 0

    if (x ==1):
        plt.plot(spin.iloc[:,0], spin.iloc[:,1], color=c_palette[1])
        color_idx += 1
        plt.plot(spin.iloc[:,0], spin.iloc[:,2], color=c_palette[4])
        color_idx += 1
        spin_time = spin.iloc[:,2].values
    if(y ==1):
        plt.plot(spin.iloc[:,0], spin.iloc[:,3], color=c_palette[1])
        color_idx += 1
        plt.plot(spin.iloc[:,0], spin.iloc[:,4], color=c_palette[4])
        color_idx += 1
        spin_time = spin.iloc[:,4].values
    if(z == 1):
        plt.plot(spin.iloc[:,0], spin.iloc[:,5], color=c_palette[1])
        color_idx += 1
        plt.plot(spin.iloc[:,0], spin.iloc[:,6], color=c_palette[4])
        color_idx += 1
        spin_time = spin.iloc[:,6].values
    
    xyz = x + y + z
    plt.ylim(-1,1)
    plt.xlim(spin.iloc[0,0], spin.iloc[-1,0])
    # plt.xlim(spin.iloc[0,0], 0.2)
    plt.xlabel('t [ns]')
    plt.ylabel('$S_z$ [$\hbar$/20]')
    
    # Add text labels instead of legend (can be positioned manually)
    ax = plt.gca()
    label_text = []

    label_text.append(f'(L)')
    label_text.append(f'(R)')
    
    # Position labels individually for easy customization
    ax.text(0.02, 0.95, label_text[0], transform=ax.transAxes, 
            fontsize=8, verticalalignment='top', color=c_palette[1])
    ax.text(0.02, 0.1, label_text[1], transform=ax.transAxes, 
            fontsize=8, verticalalignment='top', color=c_palette[4])
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder,f"spin_swap_small_period_{name}.png"), format='png')
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

def PlotSpinDensity(spin_up, spin_down, name, output_folder="../Plots"):
    os.makedirs(output_folder, exist_ok=True)

    red_cmap = ListedColormap([cm_seismic(i) for i in np.linspace(0.5, 1.0, 256)])
    blue_cmap = ListedColormap([cm_seismic(i) for i in np.linspace(0.5, 0.0, 256)])

    fig, ax = plt.subplots()
    
    # Plot spin_down with blue spectrum (plotted first so it's behind)
    scatter_down = ax.scatter(spin_down["kx"], spin_down["ky"], c=spin_down["potential"], 
                             cmap=blue_cmap, s=30, alpha=0.6, label='Spin Down')
    
    # Plot spin_up with red spectrum (plotted after so it appears on top)
    scatter_up = ax.scatter(spin_up["kx"], spin_up["ky"], c=spin_up["potential"], 
                           cmap=red_cmap, s=30, alpha=0.6, label='Spin Up')
    
    # plt.colorbar(scatter_up, ax=ax, label="Spin Density [$\hbar$/2]")
    plt.axis('equal')

    plt.xlabel(r"$x$ [nm]")
    plt.ylabel(r"$y$ [nm]")
    plt.legend()
    plt.savefig(os.path.join(output_folder,f"spin_density_t-half_{name}.png"), format='png')
    plt.close()

def PlotSpinDensityUp(spin_up, name, output_folder="../Plots"):
    os.makedirs(output_folder, exist_ok=True)

    red_cmap = ListedColormap([cm_seismic(i) for i in np.linspace(0.5, 1.0, 256)])

    fig, ax = plt.subplots()
    
    # Plot spin_up with red spectrum
    scatter_up = ax.scatter(spin_up["kx"], spin_up["ky"], c=spin_up["potential"], 
                           cmap=red_cmap, s=30)
    
    plt.colorbar(scatter_up, ax=ax, label="Spin Density [$\hbar$/2]")
    plt.axis('equal')

    plt.xlabel(r"$x$ [nm]")
    plt.ylabel(r"$y$ [nm]")
    plt.savefig(os.path.join(output_folder,f"spin_density_up_switch_{name}.png"), format='png')
    plt.close()

def PlotSpinDensityDown(spin_down, name, output_folder="../Plots"):
    os.makedirs(output_folder, exist_ok=True)

    blue_cmap = ListedColormap([cm_seismic(i) for i in np.linspace(0.5, 0.0, 256)])

    fig, ax = plt.subplots()
    
    # Plot spin_down with blue spectrum
    scatter_down = ax.scatter(spin_down["kx"], spin_down["ky"], c=spin_down["potential"], 
                             cmap=blue_cmap, s=30)
    
    plt.colorbar(scatter_down, ax=ax, label="Spin Density [$\hbar$/2]")
    plt.axis('equal')

    plt.xlabel(r"$x$ [nm]")
    plt.ylabel(r"$y$ [nm]")
    plt.savefig(os.path.join(output_folder,f"spin_density_down_switch_{name}.png"), format='png')
    plt.close()

def PlotSwitchingTime(time, v0, output_folder="../Plots"):
    plt.plot(v0, time,'.')
    # plt.xlim(v0[0],v0[-1])
    plt.xlabel('V$_0$ [eV]')
    plt.ylabel('switching time [ns]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder,f"switching_time.png"), format='png')
    plt.close()

def PlotExchangeEnergy(e2, v0, output_folder="../Plots"):
    plt.plot(v0, time,'.')
    # plt.xlim(v0[0],v0[-1])
    plt.xlabel('V$_0$ [eV]')
    plt.ylabel('switching time [ns]')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder,f"switching_time.png"), format='png')
    plt.close()