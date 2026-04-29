
# Lot of hard-coded numbers.
# Parameters need to be properly defined.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# np.seterr(invalid='ignore')
plt.rc('font', size=15)
plt.rcParams["font.family"] = "serif"

def write_file(den_data, dev_data, label, filename):
    for den, dev in zip(den_data.flat, dev_data.flat):
        out_file = open(filename, 'a+')
        out_file.write(str(den) + ' ' + str(round(dev, 3)) + ' ' + str(label) + '\n')
        out_file.close()
    return

def write_file_small(den_data, label, filename):
    for den in zip(den_data.flat):
        out_file = open(filename, 'a+')
        out_file.write(str(den) + ' ' + str(label) + '\n')
        out_file.close()
    return

file_name = "GMTE_output.txt"
file_i = open(file_name, 'r')
data_file_all = file_i.readlines()
print('Total data in file :', len(data_file_all))
file_i.close()

Cluster_Data_file =  'Cluster_sep_10cm_det_120_2days_mu_trec_100um.txt'
Con_data_file = 'Con_data_48hrs_rotation_poca_100um.txt'
U_data_file = 'U_data_48hrs_rotation_poca_100um.txt'
Pb_data_file = 'Pb_data_48hrs_rotation_poca_100um.txt'
SS_data_file = 'SS_data_48hrs_rotation_poca_100um.txt'
Air_data_file = 'Air_data_48hrs_rotation_poca_100um.txt'

sample_size = len(data_file_all) - 1
print('Chosen sample size:', sample_size)

# Minimum value of scattering angle for a POCA point to be considered in the analysis.
# Modify this parameter to visualize its effects on the choice of POCA points for further analysis.
#thmin = 0.57 # degrees ~= 0.01 radian as in Hist_POCAnew.py
for s in range(sample_size, len(data_file_all), sample_size):
    data_file = data_file_all[s-sample_size:s]
    print()
    print ('Calculated : ', s)
    XD = []
    for i, line in enumerate(data_file):
        # if i < 000000 : continue
        # if i > 100000 : break
        if i%100000 == 0 : print("Progress : ", i)
        data = [float(_) for _ in line.split()]
        # print(data)
        XD.append(data)
    XD = np.array(XD)
    x, y,z, do= XD[:, 0], XD[:, 1], XD[:, 2], XD[:, 3]
    print('Useful data :', len(XD))
    bins = 120
    # Remove extreme outliers (VERY IMPORTANT)
    x_min, x_max = np.percentile(x, [1, 99])
    y_min, y_max = np.percentile(y, [1, 99])

    # ------------------------------------------------------------
    # 1. COUNT MAP
    # ------------------------------------------------------------
    counts, xedges, yedges = np.histogram2d(
    x, y,
    bins=bins,
    range=[[x_min, x_max], [y_min, y_max]]
    )

    # ------------------------------------------------------------
    # 2. SCATTERING SUM MAP
    # ------------------------------------------------------------
    scatter_sum, _, _ = np.histogram2d(
    x, y,
    bins=bins,
    range=[[x_min, x_max], [y_min, y_max]],
    weights=do
    )

    # ------------------------------------------------------------
    # 3. SAFE AVERAGE (FIXES DIVISION ISSUE)
    # ------------------------------------------------------------
    average_scatter = np.zeros_like(scatter_sum)

    mask = counts > 0
    average_scatter[mask] = scatter_sum[mask] / counts[mask]

    # Optional: remove noise
    average_scatter[~mask] = 0

    # ------------------------------------------------------------
    # 4. LOG SCALING (CRUCIAL FOR VISIBILITY)
    # ------------------------------------------------------------
    avg_plot = np.log1p(average_scatter)

    # ------------------------------------------------------------
    # 5. PLOT (FINAL RESULT)
    # ------------------------------------------------------------
    plt.figure(figsize=(7.5, 6))

    plt.title('Scattering-based Hit Map (GMTE)')

    plt.imshow(
    avg_plot.T,
    origin='lower',
    extent=[x_min, x_max, y_min, y_max],
    cmap='jet',
    aspect='auto'
    )

    plt.colorbar(label="log(1 + scattering)")

    plt.xlabel('X (cm)', fontsize=17)
    plt.ylabel('Y (cm)', fontsize=17)

    plt.tight_layout()

    plt.savefig('GMTE_hit_map_fixed.png', dpi=150)
    plt.show()
    plt.close()


