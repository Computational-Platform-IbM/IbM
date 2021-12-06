import collections
from typing import Dict, List
import os
from tqdm import tqdm
import argparse
import imageio
from matplotlib.lines import Line2D
import matplotlib.collections
import matplotlib.pyplot as plt
import h5py
import numpy as np


def HEX2RGBsplit(c):
    """Convert HEX code of colours to RGB, spliting in Red, Green and Blue values.

    Args:
        c (string): HEX value of colours (it can be an array)

    Returns:
        rC, gC, bC (float): Red, Gren and Blue values of RGB code, respectively
    """

    rC, gC, bC = [0]*len(c), [0]*len(c), [0]*len(c)
    for cSet in range(len(c)):
        C = c[cSet]
        RGB = tuple(int(C.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
        rC[cSet] = RGB[0] / 255
        gC[cSet] = RGB[1] / 255
        bC[cSet] = RGB[2] / 255

    return(rC, gC, bC)


def muRatio(mu, s, inc):
    """Calculate ratio between mu[i] and max(mu[i,s]) for each bacterium [i] of specific specie [s].

    Args:
        mu (float): actual mu values
        s (int): specie ID
        inc(float): value to modify the darkening range of bacteria colour

    Returns:
        muA (float): alpha values to represent the ratio of mu (with darkening)
    """

    if inc == 'NoAlpha':
        return np.ones(np.size(mu))

    mu_noneg = np.maximum(mu, 0.0)
    dct = {}
    for i, j in zip(s, mu_noneg):
        dct.setdefault(i, []).append(j)
    sort_dct = collections.OrderedDict(sorted(dct.items()))
    max_mu = [max(sort_dct[key]) for key in sort_dct]
    max_mu = [1.0 if i == 0.0 else i for i in max_mu]
    muR = [mu/max_mu[s-1] for mu, s in zip(mu_noneg, s)]

    if inc == 0:
        muA = muR
    else:
        muA = (np.array(muR) + inc) / (1 + inc)

    return muA


def save_plot(i: int, xlim: List[float], ylim: List[float], bac: Dict, bacNames: list):
    """Create and save a figure of the bacteria at a certain point in time.

    Args:
        i (int): [description]
        xlim (List[float]): [description]
        ylim (List[float]): [description]

    Returns:
        string: name of the file in which the figure is saved.
    """

    nBacs = bac['nBacs'][i]
    x = bac['x'][0:nBacs, i] * 1e6
    y = bac['y'][0:nBacs, i] * 1e6
    r = bac['radius'][0:nBacs, i] * 1e6
    s = bac['species'][0:nBacs, i]
    a = bac['active'][0:nBacs, i]
    mu = bac['mu'][0:nBacs, i]

    # Calculus of mu/max(mu)
    # Recommended: 1.0 - 2.0  (inc = 'NoAlpha' -> alpha = 1; inc = 0 -> alpha = mu/max_mu)
    inc = 'NoAlpha'
    muAlpha = muRatio(mu, s, inc)

    # Colours: HEX code
    # c = ['#CC66FF', '#00B04F', '#FFA200', '#FF1482']  # Old colours
    # Colourblind-friendly: ( https://www.color-hex.com/color-palette/49436 )
    # c = ['#0072B2', '#D55E00', '#F0E442', '#CC79A7']
    c = ['#D81B60', '#1E88E5', '#FFC107', '#004D40']  # colors Chiel
    c = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00']

    rC, gC, bC = HEX2RGBsplit(c)  # HEX to RGB

    patches = [plt.Circle((x, y), radius) for x, y, radius in zip(x, y, r)]

    fig, ax = plt.subplots()

    coll = matplotlib.collections.PatchCollection(patches)
    coll.set_facecolor(
        [(rC[species-1]*muA, gC[species-1]*muA, bC[species-1]*muA) if active else '#000000' for muA, species, active in zip(muAlpha, s, a)])
    coll.set_alpha([1.0 if active else 0.5 for active in a])
    coll.set_edgecolor('k')
    coll.set_linewidth(0.05)
    ax.add_collection(coll)

    plt.xlim(xlim * 1e6)
    plt.ylim(ylim * 1e6)
    plt.xlabel('Position along x-axis [μm]')
    plt.ylabel('Position along y-axis [μm]')
    ax.set_aspect(1/ax.get_data_ratio(), adjustable='box')
    # plt.margins(0.01)

    # Legend
    legend_list = []
    for ii, bacName in enumerate(bacNames):
        legend_list.append(Line2D([0], [0], linestyle="none", marker="o",
                                  markersize=10, markerfacecolor=c[ii], markeredgecolor=c[ii]))

    plt.legend(legend_list, bacNames,
               numpoints=1, loc="upper left", frameon=False)

    filename = f'{directory}/{i}.png'
    # TODO: set dynamic value of dT_save from constants!!!
    plt.title(f'Time = {i*96}')
    plt.savefig(filename)
    plt.close()

    return filename


# %%
def generate_gif(args: Dict):
    # load data from results file
    with h5py.File(f'{directory}/results1D.mat', 'r') as f:
        print('Loading results file...')
        bac = {}
        for k in f['bac_saved'].keys():
            bac[k] = np.array(f['bac_saved'][k]).squeeze()
            print(f'Loaded bac.{k}')

    with h5py.File(simulation_file, 'r') as f:
        print('Loading simulation file...')
        grid = {}
        for k in f['grid'].keys():
            grid[k] = np.array(f['grid'][k]).squeeze()
            print(f'Loaded grid.{k}')

        for columns in f['constants']['speciesNames']:
            bacNames = []
            for row in range(len(columns)):
                bacNames.append(''.join(map(chr, f[columns[row]][:])))

    # %%
    # calculate boundaries for the plot
    lastNonzero = np.max(np.nonzero(bac['nBacs']))
    print(lastNonzero)
    final_nBacs = bac['nBacs'][lastNonzero]
    print(f'final number of bacteria: {final_nBacs}')
    xlim = np.array([bac['x'][0:final_nBacs, lastNonzero].min() - 5*grid['dx'],
                     bac['x'][0:final_nBacs, lastNonzero].max() + 5*grid['dx']])
    ylim = np.array([bac['y'][0:final_nBacs, lastNonzero].min() - 5*grid['dy'],
                     bac['y'][0:final_nBacs, lastNonzero].max() + 5*grid['dy']])
    print(xlim, ylim)

    # %%
    # create figure per timepoint
    filenames = []
    for i in tqdm(range(lastNonzero), desc='Generation frames'):
        if bac['nBacs'][i]:
            filenames.append(save_plot(i, xlim, ylim, bac, bacNames))

    # build gif
    with imageio.get_writer(f'{directory}/bacteria.gif', mode='I', fps=4) as writer:
        for filename in tqdm(filenames, desc='Gif'):
            image = imageio.imread(filename)
            writer.append_data(image)

    # %%
    # remove files afterwards
    for filename in tqdm(set(filenames), desc='Removing images'):
        os.remove(filename)

    print('DONE!')


# parse command line input
parser = argparse.ArgumentParser(
    description="Create an animation of the bacteria over time.")
parser.add_argument("simulationNumber",
                    help="[int] Simulation number for which to create the animation (1-9999)", type=int)
parser.add_argument('-nf', '--notFinished', dest='finished', default=True, action='store_false',
                    help="[bool] Has the simulation finished?")
args = parser.parse_args()

# set directory and resultsFile
sim = f'{args.simulationNumber:04d}'
directory = f'Results/{sim}'
if args.finished:
    simulation_file = f'{directory}/sim_{sim}.mat'
else:
    simulation_file = f'sim_{sim}.mat'

generate_gif(args)
