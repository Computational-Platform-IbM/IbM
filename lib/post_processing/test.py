# -*- coding: utf-8 -*-
import numpy as np
from scipy.ndimage import convolve
# import h5py
# import matplotlib.pyplot as plt
# import matplotlib.collections
# from matplotlib.lines import Line2D
# import imageio
# import argparse
# from tqdm import tqdm
# import os
# from typing import Dict, List
# import collections

# x = np.arange(3., 10.)
# y = np.arange(3., 10.)
# x = [7, 8, 9, 10, 11, 12, 13, 14, 10, 10, 10, 10] #, 10, 10, 10, 10, 9, 12, 9, 12]
# y = [10, 10, 10, 10, 10, 10, 10, 10, 7, 8, 9, 10] #, 11, 12, 13, 14, 9, 12, 12, 9]
x = [7, 8, 9, 10, 11, 12, 13, 14, 10, 10, 10, 10, 10, 10, 10, 10, 9, 12, 9, 12]
y = [10, 10, 10, 10, 10, 10, 10, 10, 7, 8, 9, 10, 11, 12, 13, 14, 9, 12, 12, 9]
dx = 2
dy = 2
nX = 11;
nY = 11;
blayer = 2;

ix = np.floor_divide(x, dx)
iy = np.floor_divide(y, dy)

grid2nBacs = np.zeros((nX, nY))
grid2Bac = np.zeros((4, nX, nY))
# print(grid2Bac)
diffRegion = np.zeros((nX, nY), dtype=bool)

# print(diffRegion)
# print(grid2nBacs)
for ii, xx, yy in zip(np.arange(1, len(x)+1), ix, iy):
    grid2Bac[int(grid2nBacs[int(xx), int(yy)]), int(xx), int(yy)] += ii
    grid2nBacs[int(xx), int(yy)] += 1 
# print('--')
# # print(grid2nBacs)
# print(grid2Bac)

##### def getDiffusionNodes(x, y, dx, dy, nX, nY, blayer):
####### ...
# getDifusionNodes()
x_min = np.min(x)
x_max = np.max(x)
y_min = np.min(y)
y_max = np.max(y)

# define an offset, such that the boundary layer always falls in the right node
offsetX = blayer + dx
offsetY = blayer + dy

# create an array (columnvector) with all end coordinates of the nodes
nodeEndCoordinatesX = np.arange(0, nX) * dx
nodeEndCoordinatesY = np.arange(0, nY) * dy

# check which noded have either the boundary layer of granule in them
diffNodesX = (nodeEndCoordinatesX > (x_min - offsetX)) & ((nodeEndCoordinatesX) < (x_max + offsetY))
diffNodesY = (nodeEndCoordinatesY > (y_min - offsetY)) & ((nodeEndCoordinatesY) < (y_max + offsetY))
# print('diffNodes')
# print(diffNodesX)
# print(diffNodesY)
# print(' ')
# print('mNodes')
sX = np.argwhere(diffNodesX)[-1] - np.argwhere(diffNodesX)[0] + 1
sY = np.argwhere(diffNodesY)[-1] - np.argwhere(diffNodesY)[0] + 1
# diffNodesX = diffNodesX[:,np.newaxis]
# mNodesX = np.tile(diffNodesX,(1,nY))
# mNodesY = np.tile(diffNodesY,(nX,1))
# mNodes = np.logical_and(mNodesX, mNodesY)
mNodes = np.logical_and(np.tile(diffNodesX[:,np.newaxis], (1,nY)), np.tile(diffNodesY, (nX,1)))
# print(mNodesX)
# print('--')
# print(mNodesY)
# print('--')
# print(mNodes)

####### return(diffNodesX, diffNodesY, sX, sY, mNodes)

# in the diffusion region, apply convolution to find boundary of grid cell
#   with bacteria
kernel = np.full([3,3], -1/8)
kernel[1,1] = 1
hasBac = np.reshape(grid2nBacs[mNodes], (int(sX),int(sY))) > 0
isBacBoundary = convolve(hasBac.astype(float), kernel) > 0
# print('--')
print(hasBac.astype(float))
print('--')
print(isBacBoundary)

# for the boundary of bacterial grid cells, compute for neighbouring grid
#   cells whetehr they are in diffusion region
maxOffsetX = np.ceil(blayer / dx)
dX = np.argwhere(diffNodesX)[0] - 1
maxOffsetY = np.ceil(blayer / dy)
dY =  np.argwhere(diffNodesY)[0] - 1

[iBB, jBB] = np.where(isBacBoundary == True)
isDiffRegion = hasBac # working matrix: 1 == diffregion, 0 == bulk

for i, j in zip(iBB, jBB):
    for di in np.arange(-maxOffsetX, maxOffsetX + 1):
        for dj in np.arange(-maxOffsetY, maxOffsetY + 1):
            if isDiffRegion[int(i + di), int(j + dj)] != 0:
                continue
            else:
                # perform actual check between gridcell and bacteria
                # bacs = np.argwhere(grid2Bac[:, (int(i) + int(di)), (int(j) + int(dj))])
                bacs = grid2Bac[np.argwhere(grid2Bac[:, (int(i) + int(di)), (int(j) + int(dj))]), (int(i) + int(di)), (int(j) + int(dj))]
                # print(bacs)
                gridcell_center = np.array([(dx * (dX + i + di) - dx / 2), (dy * (dY + j + di) - dy / 2)])

                for b in bacs:
                    ##### def isWithinBoundaryLayer(bac_x, bac_y, gridcell_center, dx, dy, di, dj, blayer):
                    print(b)
                    
                    isBLayer = True
                    ##### return(isBLayer)
                    if isBLayer:
                        print('Y')
                        isDiffRegion[int(i + di), int(j + dj)] = True
                        
print(isDiffRegion)