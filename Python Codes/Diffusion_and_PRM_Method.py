# Lot of hard-coded numbers.
# Parameters need to be properly defined.
import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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
    
class DiffusionUNet(nn.Module):#defining a neural network
    def __init__(self):#constructor
        super().__init__()#calls constructor of parent class nn.Module
        self.net = nn.Sequential( #creating a container 
            nn.Conv2d(2, 32, 3, padding=1), #adding a Convolutional 2Dlayer
            nn.ReLU(), #activation function
            nn.Conv2d(32, 64, 3, padding=1),
            nn.ReLU(),
            nn.Conv2d(64, 32, 3, padding=1),#reduces feature channels
            nn.ReLU(),
            nn.Conv2d(32, 2, 3, padding=1)#Final Convolution layer
        )

    def forward(self, x):#defining forward pass
        return self.net(x)

file_name = "POCA_output_120cm_2days_100mum.txt"
file_i = open(file_name, 'r')
data_file_all = file_i.readlines()
print('Total data in file :', len(data_file_all))
file_i.close()




Cluster_Data_file =  'Cluster_DIFFUSION.txt'
Con_data_file = 'Con_data_48hrs_rotation_poca_100um_DIFFUSION.txt'
U_data_file = 'U_data_48hrs_rotation_poca_100um_DIFFUSION.txt'
Pb_data_file = 'Pb_data_48hrs_rotation_poca_100um_DIFFUSION.txt'
SS_data_file = 'SS_data_48hrs_rotation_poca_100um_DIFFUSION.txt'
Air_data_file = 'Air_data_48hrs_rotation_poca_100um_DIFFUSION.txt'

sample_size = len(data_file_all) - 1
print('Chosen sample size:', sample_size)


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

cluster_density_dev = plt.hist2d(x,y, weights = do,  bins=bins, range=[[xmin, xmax], [ymin, ymax]])
cluster_density_dev = cluster_density_dev[0].reshape(bins*bins)
plt.close()
average_deviation = cluster_density_dev/cluster_density_s
average_deviation = average_deviation.reshape((bins,bins))

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
plt.show()
plt.close()

df=pd.DataFrame(cluster_density[0])
df_d=pd.DataFrame(average_deviation)

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

# save diffusion cluster data


print("Cluster density file saved:", Cluster_Data_file)


#  DIFFUSION DENOISING (NEW) 

device = "cpu" #using only cpu
diff_model = DiffusionUNet().to(device)
diff_model.load_state_dict(#prepares to load trained parameters in the model
    torch.load("diffusion_prm.pth", map_location=device)
)
diff_model.eval()


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

   
cluster_density = denoised[0,0].cpu().numpy()
average_deviation = denoised[0,1].cpu().numpy()

    # Following min and max values need to be cross-checked with relevant histos.
    # In fact, they can be considered as input parameters that the use provides after carefully
    # analyzing the histograms plotted just ahead.

# save diffusion cluster data

"""    
plt.title('Diffusion Average deviation')
plt.imshow(average_deviation.T, origin='lower', cmap='jet')# cmin=2500, cmax=3000)
plt.grid()
plt.colorbar()
plt.savefig('Diffusion_Average_Deviation.png', dpi=150, bbox_inches='tight')	
plt.show()
plt.close()
plt.title('Diffusion Cluster Density')
plt.imshow(cluster_density.T, origin='lower', cmap='jet')# cmin=2500, cmax=3000)
plt.grid()
plt.colorbar()
plt.savefig('Diffusion_Cluster_Density.png', dpi=150, bbox_inches='tight')	
plt.show()
plt.close()  
"""

# Adjusting cutoffs to show each material
#PRM AFTER DIFFUSION 


from scipy.signal import convolve2d

#GLOBAL KERNEL
ker_All = -np.ones((3,3))

output_All = convolve2d(cluster_density, ker_All, mode='same')

All_cutoff = np.mean(convolve2d(ker_All, ker_All, mode='valid'))
print("All Cutoff:", All_cutoff)

# ROI EXTRACTION 
df = pd.DataFrame(cluster_density)
df_d = pd.DataFrame(average_deviation)

data_Con = np.array(df.iloc[74:91, 76:95])
data_U   = np.array(df.iloc[50:69, 20:39])
data_Pb  = np.array(df.iloc[49:69, 81:100])
data_SS  = np.array(df.iloc[18:39, 50:69])
data_Air = np.array(df.iloc[79:89, 54:67])

mean_Con = data_Con.mean()
std_Con  = data_Con.std()

mean_U = data_U.mean()
std_U  = data_U.std()

mean_Pb = data_Pb.mean()
std_Pb  = data_Pb.std()

mean_SS = data_SS.mean()
std_SS  = data_SS.std()

mean_Air = data_Air.mean()
std_Air  = data_Air.std()

# MATERIAL KERNELS

ker_Con = np.random.normal(mean_Con, std_Con, (3,3))
ker_U   = np.random.normal(mean_U, std_U, (3,3))
ker_Pb  = np.random.normal(mean_Pb, std_Pb, (3,3))
ker_SS  = np.random.normal(mean_SS, std_SS, (3,3))
ker_Air = np.random.normal(mean_Air, std_Air, (3,3))

output_Con = convolve2d(cluster_density, ker_Con, mode='same')
output_U   = convolve2d(cluster_density, ker_U, mode='same')
output_Pb  = convolve2d(cluster_density, ker_Pb, mode='same')
output_SS  = convolve2d(cluster_density, ker_SS, mode='same')
output_Air = convolve2d(cluster_density, ker_Air, mode='same')

#CUTOFFS 
Con_cutoff = np.mean(convolve2d(ker_Con, ker_Con, mode='valid'))
U_cutoff   = np.mean(convolve2d(ker_U, ker_U, mode='valid'))
Pb_cutoff  = np.mean(convolve2d(ker_Pb, ker_Pb, mode='valid'))
SS_cutoff  = np.mean(convolve2d(ker_SS, ker_SS, mode='valid'))
Air_cutoff = np.mean(convolve2d(ker_Air, ker_Air, mode='valid'))
Con_cutoff = Con_cutoff + 0.2 * np.std(output_Con)
U_cutoff   = U_cutoff   + 0.2 * np.std(output_U)
Pb_cutoff  = Pb_cutoff  + 0.2 * np.std(output_Pb)
SS_cutoff  = SS_cutoff  + 0.2 * np.std(output_SS)
Air_cutoff = Air_cutoff + 0.2 * np.std(output_Air)

print("Concrete cutoff:", Con_cutoff)
print("U cutoff:", U_cutoff)
print("Pb cutoff:", Pb_cutoff)
print("SS cutoff:", SS_cutoff)
print("Air cutoff:", Air_cutoff)

# THRESHOLDING 
PRM_Density = np.where(output_All < All_cutoff, cluster_density, np.nan)
PRM_Deviation = np.where(output_All < All_cutoff, average_deviation, np.nan)
Con_cluster = np.where(output_Con > Con_cutoff, output_Con, np.nan)
U_cluster   = np.where(output_U > U_cutoff, output_U, np.nan)
Pb_cluster  = np.where(output_Pb > Pb_cutoff, output_Pb, np.nan)
SS_cluster  = np.where(output_SS > SS_cutoff, output_SS, np.nan)
Air_cluster = np.where(output_Air > Air_cutoff, output_Air, np.nan)

#Calculating Accuracy
from sklearn.metrics import accuracy_score

true_map = np.zeros_like(cluster_density)

true_map[74:91, 76:95] = 1
true_map[50:69, 20:39] = 2
true_map[49:69, 81:100] = 3
true_map[18:39, 50:69] = 4
true_map[79:89, 54:67] = 5
pred_map = np.zeros_like(cluster_density)

pred_map[~np.isnan(Con_cluster)] = 1
pred_map[~np.isnan(U_cluster)] = 2
pred_map[~np.isnan(Pb_cluster)] = 3
pred_map[~np.isnan(SS_cluster)] = 4
pred_map[~np.isnan(Air_cluster)] = 5
accuracy = accuracy_score(true_map.flatten(), pred_map.flatten())
print("Accuracy:", accuracy)
"""
from sklearn.metrics import confusion_matrix, accuracy_score

y_true = true_map.flatten()
y_pred = pred_map.flatten()

cm = confusion_matrix(y_true, y_pred)
accuracy = accuracy_score(y_true, y_pred)
"""
#VISUALIZATION 

plt.title('Diff_PRM Average deviation')
plt.imshow(PRM_Deviation.T, origin='lower', cmap='jet')# cmin=2500, cmax=3000)
plt.grid()
plt.colorbar()
plt.text(5,110,
         f'Accuracy = {accuracy:.3f}',
         fontsize=16,
         color='white',
         bbox=dict(facecolor='black', alpha=0.6))
plt.savefig('Diffusion_PRM_Average_Deviation.png', dpi=150, bbox_inches='tight')	
plt.show()
plt.close()
plt.title('Diff_PRM Cluster Density')
plt.imshow(PRM_Density.T, origin='lower', cmap='jet')# cmin=2500, cmax=3000)
plt.grid()
plt.colorbar()
plt.text(5,110,
         f'Accuracy = {accuracy:.3f}',
         fontsize=16,
         color='white',
         bbox=dict(facecolor='black', alpha=0.6))
plt.savefig('Diffusion_PRM_Cluster_Density.png', dpi=150, bbox_inches='tight')	
plt.show()
plt.close()  
plt.figure()
plt.imshow(Con_cluster.T, origin='lower', cmap='jet')
#plt.title("Concrete (Diffusion + PRM)")
plt.colorbar()
plt.savefig("Concrete_Diffusion_PRM.png")
plt.show()

plt.figure()
plt.imshow(U_cluster.T, origin='lower', cmap='jet')
#plt.title("Uranium (Diffusion + PRM)")
plt.colorbar()
plt.savefig("Uranium_Diffusion_PRM.png")
plt.show()

plt.figure()
plt.imshow(Pb_cluster.T, origin='lower', cmap='jet')
#plt.title("Lead (Diffusion + PRM)")
plt.colorbar()
plt.savefig("Lead_Diffusion_PRM.png")
plt.show()

plt.figure()
plt.imshow(SS_cluster.T, origin='lower', cmap='jet')
#plt.title("Steel (Diffusion + PRM)")
plt.colorbar()
plt.savefig("Steel_Diffusion_PRM.png")
plt.show()

plt.figure()
plt.imshow(Air_cluster.T, origin='lower', cmap='jet')
#plt.title("Air (Diffusion + PRM)")
plt.colorbar()
plt.savefig("Air_Diffusion_PRM.png")
plt.show()
    
    

