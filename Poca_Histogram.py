#2d Histogram Plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

file_i = open("POCA_output.txt", 'r')


thmin = 0.57 
thmax = 30.0 

x_poca, y_poca, z_poca, angle = [], [], [], []
for i, line in enumerate(file_i):
    i = i + 1
  
    print (i,"th Iteration is done")
    data = np.array([float(_) for _ in line.split()])
    xp, yp, zp, th = data[0], data[1], data[2], data[3]
    if (True):
       x_poca.append(xp)
       y_poca.append(yp)
       z_poca.append(zp)
       angle.append(th)

  
# Use 3D scatter plot to visualize the distribution of POCA points.
# This visualization will be useful to choose parameters for the histogram plots.
# Modify the min, max and bins parameters to check the effects on final results.
xmin, xmax = -300, 300
ymin, ymax = -300, 300
zmin, zmax = -300, 300
bins = 120



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
plt.hist2d(x_poca, y_poca, bins=120, range=[[-50, 50], [-50, 50]])#, cmin = 80, cmax = 1000)
plt.colorbar()
plt.savefig('Histogram_2d.png', dpi=150, bbox_inches='tight')
plt.show()
