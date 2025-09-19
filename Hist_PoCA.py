### To draw the 2D histogram of poca points

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import rcParams
# import warnings
# warnings.filterwarnings('error')

plt.style.use('classic')

SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 10, 14, 16
#plt.rcParams["font.family"] = "serif"
#plt.rcParams['font.serif'] = ['Century Schoolbook']
plt.rc('font', size=MEDIUM_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

file_i = open("POCA_output.txt", 'r')

# Minimum value of scattering angle for a POCA point to be considered in the analysis.
# Modify this parameter to visualize its effects on the choice of POCA points for further analysis.
thmin = 0.57 # degrees ~= 0.01 radian
thmax = 30.0 # degrees ~= 0.523599 radian

x_poca, y_poca, z_poca, dev_th = [], [], [], []
for i, line in enumerate(file_i):
    i = i + 1
   # if i > 100 : break
    if i%1000==0 : print ("Work in progress :", i)
    data = np.array([float(_) for _ in line.split()])
    # data = np.array([float(_) for _ in line.split()])
    # print('linedata: ', linedata)
    xp, yp, zp, th = data[0], data[1], data[2], data[3]
    x_poca.append(xp)
    y_poca.append(yp)
    z_poca.append(zp)
    dev_th.append(th)
    # if (thmax > data[3] > thmin):
       # x_poca.append(xp)
       # y_poca.append(yp)
       # z_poca.append(zp)
       # dev_th.append(th)

color, s, alpha, marker = 'red', 1, 0.1, 's'
#fig = plt.figure(figsize=(8,8))
plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.8, wspace=0.1, hspace=0.1)
# Use 3D scatter plot to visualize the distribution of POCA points.
# This visualization will be useful to choose parameters for the histogram plots.
# Modify the min, max and bins parameters to check the effects on final results.
xmin, xmax = -300, 300
ymin, ymax = -300, 300
zmin, zmax = -300, 300
bins = 120
"""
color, s, alpha, marker = 'red', 1, 0.1, 's'
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(projection='3d')
plt.subplots_adjust(left=0.2, bottom=0.2, right=0.8, top=0.8, wspace=0.1, hspace=0.1)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, zmax)
ax.scatter(x_poca, y_poca, z_poca, marker=marker, color=color, s=s, alpha = alpha)
ax.grid(True)
plt.show()
"""
#ax.scatter(x_poca, y_poca, z_poca, marker=marker, color=color, s=s, alpha = alpha)
#ax.grid(True)

plt.hist2d(x_poca, y_poca, bins=120, range=[[-50, 50], [-50, 50]])#, cmin = 80, cmax = 1000)
plt.colorbar()
plt.savefig('Histogram_poca_ideal.png', dpi=150, bbox_inches='tight')
plt.show()
