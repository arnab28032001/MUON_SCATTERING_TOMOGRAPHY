#2d Histogram Plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

file_i = open("POCA_output_120cm_2days_100mum.txt", 'r')


thmin = 0.57
thmax = 30.0 

x_poca, y_poca, z_poca, angle = [], [], [], []
for i, line in enumerate(file_i):
    i = i + 1
  
    print (i,"th Iteration is done")
    data = np.array([float(_) for _ in line.split()])
    
    if(data[3]>=thmin and data[3]<=thmax):
        x_poca.append(data[0])
        y_poca.append(data[1])
        z_poca.append(data[2])
       #angle.append(data[3])

  
# Use 3D scatter plot to visualize the distribution of POCA points.
# This visualization will be useful to choose parameters for the histogram plots.
# Modify the min, max and bins parameters to check the effects on final results.
"""
xmin, xmax = -30, 30
ymin, ymax = -30, 30
zmin, zmax = -100, 100
bins = 120
"""

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

plt.hist2d(x_poca, y_poca, bins=300, range=[[-300, 300], [-300, 300]],cmap="viridis")
plt.title("Average Deviation Angle_120cm_2days_(resolution=0.001)")
plt.xlabel("X-->")
plt.ylabel("Y--->")
plt.colorbar()
plt.savefig('Histogram_POCA_120cm_2days_(resolution=0.001)_viridis.jpg', dpi=300, bbox_inches='tight')
plt.show()

