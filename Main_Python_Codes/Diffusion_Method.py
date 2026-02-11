# Lot of hard-coded numbers.
# Parameters need to be properly defined.
import torch
import torch.nn as nn
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
    
class DiffusionUNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(2, 32, 3, padding=1),
            nn.ReLU(),
            nn.Conv2d(32, 64, 3, padding=1),
            nn.ReLU(),
            nn.Conv2d(64, 32, 3, padding=1),
            nn.ReLU(),
            nn.Conv2d(32, 2, 3, padding=1)
        )

    def forward(self, x):
        return self.net(x)

file_name = "POCA_output_120cm_2days_100mum.txt"
file_i = open(file_name, 'r')
data_file_all = file_i.readlines()
print('Total data in file :', len(data_file_all))
file_i.close()

device = "cpu"
diff_model = DiffusionUNet().to(device)
diff_model.load_state_dict(
    torch.load("diffusion_prm.pth", map_location=device)
)
diff_model.eval()


Cluster_Data_file =  'Cluster_sep_10cm_det_120_2days_poca_100um_DIFFUSION.txt'
Con_data_file = 'Con_data_48hrs_rotation_poca_100um_DIFFUSION.txt'
U_data_file = 'U_data_48hrs_rotation_poca_100um_DIFFUSION.txt'
Pb_data_file = 'Pb_data_48hrs_rotation_poca_100um_DIFFUSION.txt'
SS_data_file = 'SS_data_48hrs_rotation_poca_100um_DIFFUSION.txt'
Air_data_file = 'Air_data_48hrs_rotation_poca_100um_DIFFUSION.txt'

sample_size = len(data_file_all) - 1
print('Chosen sample size:', sample_size)





# Minimum value of scattering angle for a POCA point to be considered in the analysis.
# Modify this parameter to visualize its effects on the choice of POCA points for further analysis.
thmin = 0.57 # degrees ~= 0.01 radian as in Hist_POCAnew.py
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
    x, y,z, do = XD[:, 0], XD[:, 1], XD[:, 2], XD[:, 3]
    print('Useful data :', len(XD))
    XD = np.array(XD)
    x, y,z, do = XD[:, 0], XD[:, 1], XD[:, 2], XD[:, 3]
    print('Useful data :', len(XD))
    
    bins = 120
    xmin, xmax = -300, 300
    ymin, ymax = -300, 300
    
    cluster_density = plt.hist2d(x,y, bins=bins, range=[[xmin, xmax], [ymin, ymax]])
    cluster_density_s = cluster_density[0].reshape(bins*bins)
    plt.close()
    plt.figure(figsize=(7.5, 6))

    plt.title('Cluster density')
    plt.imshow(cluster_density[0].T, origin='lower', cmap='jet')
    plt.grid()
    plt.colorbar()
    plt.xlabel('X bin (5.0 mm)', fontsize=17) # change value in brackets as needed
    plt.ylabel('Y bin (5.0 mm)', fontsize=17)
    xticks = [0.5, 1.5, 2.5, 3.5]
    xticks_labels = ['Al', 'Fe', 'Pb', 'U']
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig('clusden_sep_10cm_det_120_2days_poca_100um_DIFFUSION.png', dpi=150, bbox_inches='tight')	
    plt.show()
    plt.close()

    cluster_density_dev = plt.hist2d(x,y, weights = do,  bins=bins, range=[[xmin, xmax], [ymin, ymax]])
    cluster_density_dev = cluster_density_dev[0].reshape(bins*bins)
    plt.close()
    average_deviation = cluster_density_dev/cluster_density_s
    average_deviation = average_deviation.reshape((bins,bins))
    # ================= DIFFUSION DENOISING (NEW) =================
    PRM_density = np.nan_to_num(cluster_density[0])
    PRM_deviation = np.nan_to_num(average_deviation)
    # Normalize for diffusion
    PRM_density = PRM_density / (PRM_density.max() + 1e-6)
    PRM_deviation = PRM_deviation / (PRM_deviation.max() + 1e-6)
    # Convert to torch tensor
    inp = torch.tensor(
    np.stack([PRM_density, PRM_deviation]),
    dtype=torch.float32
    ).unsqueeze(0).to(device)
    with torch.no_grad():
        denoised = diff_model(inp)

    # Replace original PRM maps
    cluster_density = denoised[0,0].cpu().numpy()
    average_deviation = denoised[0,1].cpu().numpy()

    # =============================================================
        # print('cluster_density : ', average_deviation)
    plt.title('Average deviation')
    plt.imshow(average_deviation.T, origin='lower', cmap='jet')# cmin=2500, cmax=3000)
    plt.grid()
    plt.colorbar()
    plt.savefig('avgdev_sep_10cm_det_60_3days_poca_1000um_DIFFUSION.png', dpi=150, bbox_inches='tight')	
    plt.show()
    plt.close()

# Create dataframes from cluster density and deviation data 
    df=pd.DataFrame(cluster_density[0])
    df_d=pd.DataFrame(average_deviation)

    # Following min and max values need to be cross-checked with relevant histos.
    # In fact, they can be considered as input parameters that the use provides after carefully
    # analyzing the histograms plotted just ahead.
    xminXmean, xmaxXmean = 50, 70
    yminYmean, ymaxYmean = 50, 70    
    xminCon, xmaxCon = 50, 70 # Note that  there is a rotation involved that modifies the geometry 
    yminCon, ymaxCon = 95, 105 # quite significantly (as is apparent from the histograms).
    xminU, xmaxU = 70, 90
    yminU, ymaxU = 50, 70
    xminAir, xmaxAir = 30, 50
    yminAir, ymaxAir = 50, 70
    xminPb, xmaxPb = 50, 70
    yminPb, ymaxPb = 70, 90
    xminSS, xmaxSS = 50, 70
    yminSS, ymaxSS = 30, 50
    
    # data = np.array(df.iloc[xminXmean:xmaxXmean, yminYmean:ymaxYmean])    # Mid Position
    data_All = np.array(df.iloc[35:45, 35:45])    # Mid Position
    # print('filter :', data)

    # data_Con = np.array(df.iloc[xminCon:xmaxCon, yminCon:ymaxCon])
    # data_U = np.array(df.iloc[xminU:xmaxU, yminU:ymaxU])
    # data_Pb = np.array(df.iloc[xminPb:xmaxPb, yminPb:ymaxPb])
    # data_SS = np.array(df.iloc[xminSS:xmaxSS, yminSS:ymaxSS])
    # data_Air = np.array(df.iloc[xminAir:xmaxAir, yminAir:ymaxAir])

    # data_Con_d = np.array(df.iloc[xminCon:xmaxCon, yminCon:ymaxCon])
    # data_U_d = np.array(df.iloc[xminU:xmaxU, yminU:ymaxU])
    # data_Pb_d = np.array(df.iloc[xminPb:xmaxPb, yminPb:ymaxPb])
    # data_SS_d = np.array(df.iloc[xminSS:xmaxSS, yminSS:ymaxSS])
    # data_Air_d = np.array(df.iloc[xminAir:xmaxAir, yminAir:ymaxAir])

    data_Con = np.array(df.iloc[74:91, 76:95])
    data_U = np.array(df.iloc[50:69, 20:39])
    data_Pb = np.array(df.iloc[49:69, 81:100])
    data_SS = np.array(df.iloc[18:39, 50:69])
    data_Air = np.array(df.iloc[79:89, 54:67])
    
    data_Con_d = np.array(df_d.iloc[74:91, 76:95])  # consistent with earlier numbers and works well
    data_U_d = np.array(df_d.iloc[50:70, 20:40])  # consistent with earlier numbers and works well
    data_Pb_d = np.array(df_d.iloc[50:70, 80:100])  # consistent with earlier numbers and works well
    data_SS_d = np.array(df_d.iloc[17:40, 49:70])
    data_Air_d = np.array(df_d.iloc[80:99, 49:70])

# Output cluster data for further processing using logistic regression / SVM
    write_file(data_Con, data_Con_d, 1, Cluster_Data_file)
    write_file(data_U, data_U_d, 2, Cluster_Data_file)
    write_file(data_Pb, data_Pb_d, 3, Cluster_Data_file)
    write_file(data_SS, data_SS_d, 4, Cluster_Data_file)
    write_file(data_Air, data_Air_d, 5, Cluster_Data_file)
    """
    plt.title('data_Con')
    plt.imshow(data_Con.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig("Count_concrete.png")
    plt.show()
    plt.title('data_Con_d')
    plt.imshow(data_Con_d.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig("Dev_concrete.png")
    plt.show()

    plt.title('data_U')
    plt.imshow(data_U.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig("Count_U.png")
    plt.show()
    plt.title('data_U_d')    
    plt.imshow(data_U_d.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig("Dev_U.png")
    plt.show()

    plt.title('data_SS')    
    plt.imshow(data_SS.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig("Count_SS.png")
    plt.show()
    plt.title('data_SS_d')    
    plt.imshow(data_SS_d.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig("Dev_SS.png")
    plt.show()

    plt.title('data_Pb')    
    plt.imshow(data_Pb.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig("Count_Pb.png")
    plt.show()
    plt.title('data_Pb_d')    
    plt.imshow(data_Pb_d.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig("Dev_Pb.png")
    plt.show()

    plt.title('data_Air')  
    plt.imshow(data_Air.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig("Count_Air.png")
    plt.show()
    plt.title('data_Air_d')  
    plt.imshow(data_Air_d.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.savefig("Dev_air.png")
    plt.show()
    """

    # break
    
    # break

    mean_All = data_All.mean()
    std_deviation_All = data_All.std()
    print('For All, Mean :', mean_All, ', STD: ', std_deviation_All)  

    mean_Con = data_Con.mean()
    std_deviation_Con = data_Con.std()
    print('For Con, Mean :', mean_Con, ', STD: ', std_deviation_Con)  

    mean_Pb = data_Pb.mean()
    std_deviation_Pb = data_Pb.std()
    print('For Pb, Mean :', mean_Pb, ', STD: ', std_deviation_Pb)  

    mean_U = data_U.mean()
    std_deviation_U = data_U.std()
    print('For U, Mean :', mean_U, ', STD: ', std_deviation_U)  
    
    mean_SS = data_SS.mean()
    std_deviation_SS = data_SS.std()
    print('For SS, Mean :', mean_SS, ', STD: ', std_deviation_SS)  

    mean_Air = data_Air.mean()
    std_deviation_Air = data_Air.std()
    print('For Air, Mean :', mean_Air, ', STD: ', std_deviation_Air)  

    


    
    
    

