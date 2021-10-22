import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.collections



f = h5py.File('./Results/0003/results1D.mat', 'r')
bac = {}
bac_x = np.array(f.get('bac_saved/x'))
bac_y = np.array(f.get('bac_saved/y'))
bac_r = np.array(f.get('bac_saved/radius'))
bac_s = np.array(f.get('bac_saved/species'))
nbacs = np.array(f.get('bac_saved/nBacs')).ravel()
bac_a = np.array(f.get('bac_saved/active'))
i = np.max(np.nonzero(nbacs))

x = bac_x[0:nbacs[i], i]
y = bac_y[0:nbacs[i], i]
r = bac_r[0:nbacs[i], i]
s = bac_s[0:nbacs[i], i]
a = bac_a[0:nbacs[i], i]
xy = np.stack([x, y], axis=1)

c = ['#D81B60', '#1E88E5', '#FFC107', '#004D40']
patches = [plt.Circle(center, radius) for center, radius in zip(xy, r)]

fig, ax = plt.subplots()

coll = matplotlib.collections.PatchCollection(patches)
coll.set_facecolor([c[species-1] if active else '#000000' for species, active in zip(s, a)])
coll.set_alpha([1 if active else 0.5 for active in a])
coll.set_edgecolor('k')
coll.set_linewidth(0.1)
ax.add_collection(coll)

ax.margins(0.01)
plt.axis('square')
plt.show()



