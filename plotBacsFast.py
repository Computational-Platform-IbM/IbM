# %%
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.collections
import imageio
import argparse
from tqdm import tqdm
import os
from typing import Dict, List

# %%


def save_plot(i: int, xlim: List[float], ylim: List[float], bac:Dict):
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

    c = ['#D81B60', '#1E88E5', '#FFC107', '#004D40']
    patches = [plt.Circle((x, y), radius) for x, y, radius in zip(x, y, r)]

    fig, ax = plt.subplots()

    coll = matplotlib.collections.PatchCollection(patches)
    coll.set_facecolor(
        [c[species-1] if active else '#000000' for species, active in zip(s, a)])
    coll.set_alpha([1 if active else 0.5 for active in a])
    coll.set_edgecolor('k')
    coll.set_linewidth(0.1)
    ax.add_collection(coll)

    plt.xlim(xlim * 1e6)
    plt.ylim(ylim * 1e6)
    plt.xlabel('Position along x-axis [μm]')
    plt.ylabel('Position along y-axis [μm]')
    ax.set_aspect(1/ax.get_data_ratio(), adjustable='box')
    # plt.margins(0.01)

    filename = f'{directory}/{i}.png'
    plt.savefig(filename)
    plt.close()

    return filename


# %%
def generate_gif(args:Dict):
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
        dtBac = np.array(f['constants']['dT_bac']).squeeze()

    # %%
    # calculate boundaries for the plot
    lastNonzero = np.max(np.nonzero(bac['nBacs']))
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
        filenames.append(save_plot(i, xlim, ylim, bac))

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
