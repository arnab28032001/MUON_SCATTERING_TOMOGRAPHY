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

file_name = "POCA_output_120cm_2days_100mum.txt"
file_i = open(file_name, 'r')
data_file_all = file_i.readlines()
print('Total data in file :', len(data_file_all))
file_i.close()

Cluster_Data_file =  'Cluster_sep_10cm_det_120_2days_poca_100um.txt'
Con_data_file = 'Con_data_48hrs_rotation_poca_50um.txt'
U_data_file = 'U_data_48hrs_rotation_poca_50um.txt'
Pb_data_file = 'Pb_data_48hrs_rotation_poca_50um.txt'
SS_data_file = 'SS_data_48hrs_rotation_poca_50um.txt'
Air_data_file = 'Air_data_48hrs_rotation_poca_50um.txt'

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
    plt.savefig('clusden_sep_10cm_det_120_2days_poca_100um.png', dpi=150, bbox_inches='tight')	
    plt.show()
    plt.close()

    cluster_density_dev = plt.hist2d(x,y, weights = do,  bins=bins, range=[[xmin, xmax], [ymin, ymax]])
    cluster_density_dev = cluster_density_dev[0].reshape(bins*bins)
    plt.close()
    average_deviation = cluster_density_dev/cluster_density_s
    average_deviation = average_deviation.reshape((bins,bins))
    # print('cluster_density : ', average_deviation)
    plt.title('Average deviation')
    plt.imshow(average_deviation.T, origin='lower', cmap='jet')# cmin=2500, cmax=3000)
    plt.grid()
    plt.colorbar()
    plt.savefig('avgdev_sep_10cm_det_60_3days_poca_1000um.png', dpi=150, bbox_inches='tight')	
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

    plt.title('data_Con')
    plt.imshow(data_Con.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.show()
    plt.title('data_Con_d')
    plt.imshow(data_Con_d.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.show()

    plt.title('data_U')
    plt.imshow(data_U.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.show()
    plt.title('data_U_d')    
    plt.imshow(data_U_d.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.show()

    plt.title('data_SS')    
    plt.imshow(data_SS.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.show()
    plt.title('data_SS_d')    
    plt.imshow(data_SS_d.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.show()

    plt.title('data_Pb')    
    plt.imshow(data_Pb.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.show()
    plt.title('data_Pb_d')    
    plt.imshow(data_Pb_d.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.show()

    plt.title('data_Air')  
    plt.imshow(data_Air.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.show()
    plt.title('data_Air_d')  
    plt.imshow(data_Air_d.T, origin='lower', cmap='jet')
    plt.colorbar()
    plt.show()

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


# Preparation for filtering operations
    from scipy.signal import convolve2d 
    
# Filter for the entire figure    
    # ker_All = np.random.normal(mean_All, std_deviation_All, (3,3))
    # ker_All = [[-1,-1,-1],[-1,-1,-1],[-1,-1,-1]]
    # print(ker_All)
    # print(ker_All[::-1, ::-1])
    ker_All = [[-1,-1,-1],[-1,-1,-1],[-1,-1,-1]]
    ker_All = np.array(ker_All)
    output_All=convolve2d(cluster_density[0], ker_All, 'same')   #kernel is moved over all the counts
    print('Dimension of cluster density[0]: ', cluster_density[0].shape)
    print('Dimension of output_All after convolution with ker_All: ', output_All.shape)

# Concrete filter
    ker_Con = np.random.normal(mean_Con, std_deviation_Con, (3,3))
    ker_Con = np.array(ker_Con)
    output_Con=convolve2d(cluster_density[0], ker_Con, 'same')
    
# Uranium filter
    ker_U = np.random.normal(mean_U, std_deviation_U, (3,3))
    ker_U = np.array(ker_U)
    output_U = convolve2d(cluster_density[0], ker_U, 'same')

# Lead filter
    ker_Pb = np.random.normal(mean_Pb, std_deviation_Pb, (3,3))
    ker_Pb = np.array(ker_Pb)
    output_Pb = convolve2d(cluster_density[0], ker_Pb, 'same')

# SS filter
    ker_SS = np.random.normal(mean_SS, std_deviation_SS, (3,3))
    ker_SS = np.array(ker_SS)
    output_SS = convolve2d(cluster_density[0], ker_SS, 'same')
    
# Air filter
    ker_Air = np.random.normal(mean_Air, std_deviation_Air, (3,3))
    ker_Air = np.array(ker_Air)
    output_Air = convolve2d(cluster_density[0], ker_Air, 'same')

 # Convolute each filter with itself to find different cut-offs
    self_All=convolve2d(ker_All, ker_All, 'valid')
    All_cutoff = self_All.mean()
    print('All cut off :', All_cutoff)
    self_Con=convolve2d(ker_Con, ker_Con, 'valid')
    Con_cutoff = self_Con.mean()
    print('Concrete cut off :', Con_cutoff)
    self_U=convolve2d(ker_U, ker_U, 'valid')
    U_cutoff = self_U.mean()
    print('U cut off :', U_cutoff)
    self_Pb=convolve2d(ker_Pb, ker_Pb, 'valid')
    Pb_cutoff = self_Pb.mean()
    print('Pb cut off :', Pb_cutoff)
    self_SS=convolve2d(ker_SS, ker_SS, 'valid')
    SS_cutoff = self_SS.mean()
    print('SS cut off :', SS_cutoff)
    self_Air=convolve2d(ker_Air, ker_Air, 'valid')
    Air_cutoff = self_Air.mean()
    print('Air cut off :', Air_cutoff)

    cluster_density = cluster_density[0].reshape(len(cluster_density[0])*len(cluster_density[0]))
    average_deviation = average_deviation.reshape(len(average_deviation)*len(average_deviation))
    cluster_output_All  = output_All.reshape(len(output_All)*len(output_All))
    plt.title('All cluster output')
    plt.plot(cluster_output_All)
    plt.show()
    cluster_output_Con  = output_Con.reshape(len(output_Con)*len(output_Con))
    plt.title('Con cluster output')
    plt.plot(cluster_output_Con)
    plt.show()
    cluster_output_U  = output_U.reshape(len(output_U)*len(output_U))
    plt.title('U cluster output')
    plt.plot(cluster_output_U)
    plt.show()
    cluster_output_Pb  = output_Pb.reshape(len(output_Pb)*len(output_Pb))
    plt.title('Pb cluster output')
    plt.plot(cluster_output_Pb)
    plt.show()
    cluster_output_SS  = output_SS.reshape(len(output_SS)*len(output_SS))
    plt.title('SS cluster output')
    plt.plot(cluster_output_SS)
    plt.show()
    cluster_output_Air  = output_Air.reshape(len(output_Air)*len(output_Air))
    plt.title('Air cluster output')
    plt.plot(cluster_output_Air)
    plt.show()

# The limiting values necessary for the following analysis can be guessed from the mean and deviation,
# the cluster output plots and the cutoff values computed by convoluting different kernels, as done
# above.
    PRM_Density = []
    PRM_Deviation = []
    Con_cluster = []
    U_cluster = []
    Pb_cluster = []
    SS_cluster = []
    Air_cluster = []

    cutoff = -300
    for i, x in enumerate(cluster_output_All):
        if x > cutoff:
            PRM_Density.append(0)
            PRM_Deviation.append(0)
        else:
            PRM_Density.append(cluster_density[i])
            PRM_Deviation.append(average_deviation[i])
    PRM_Density = np.array(PRM_Density).reshape((bins,bins))
    PRM_Density[PRM_Density <= 0] = None
    PRM_Deviation= np.array(PRM_Deviation).reshape((bins,bins))
    print('Dimension of PRM_Density and PRM_Deviation: ', PRM_Density.shape, PRM_Deviation.shape)

    for i, x in enumerate(cluster_output_Con):
        if (1.2e5) > x > (4.0e4):
            Con_cluster.append(cluster_output_Con[i])
        else:
            Con_cluster.append(0)
    Con_cluster = np.array(Con_cluster).reshape((bins,bins))
    Con_cluster[Con_cluster <= 0] = None
    # write_file_small(Con_cluster, 1, Con_data_file) # lot of nan data in these files!
    # print('Con data file completed')    
 
    for i, x in enumerate(cluster_output_U):
        if x >= 1.5e5:
            U_cluster.append(cluster_output_U[i])
        else:
            U_cluster.append(0)
    U_cluster = np.array(U_cluster).reshape((bins,bins))
    U_cluster[U_cluster <= 0] = None
    # write_file_small(U_cluster, 2, U_data_file) # lot of nan data in these files!
    # print('U data file completed')
    
    for i, x in enumerate(cluster_output_Pb):
        if  4.2e5 > x > 1.2e5 :
            Pb_cluster.append(cluster_output_Pb[i])
        else:
            Pb_cluster.append(0)
    Pb_cluster = np.array(Pb_cluster).reshape((bins,bins))
    Pb_cluster[Pb_cluster <= 0] = None
    # write_file_small(Pb_cluster, 3, Pb_data_file) # lot of nan data in these files!
    # print('Pb data file completed')
   
    for i, x in enumerate(cluster_output_SS):
        if 1.8e5 > x > 6.0e4:
            SS_cluster.append(cluster_output_SS[i])
        else:
            SS_cluster.append(0)
    SS_cluster = np.array(SS_cluster).reshape((bins,bins))
    SS_cluster[SS_cluster <= 0] = None
    # write_file_small(SS_cluster, 4, SS_data_file) # lot of nan data in these files!    
    # print('SS data file completed')
    
    for i, x in enumerate(cluster_output_Air):
        if 1.5e5 > x  > 2.0e4:
            Air_cluster.append(cluster_output_Air[i])
        else:
            Air_cluster.append(0)
    Air_cluster = np.array(Air_cluster).reshape((bins,bins))
    Air_cluster[Air_cluster <= 0] = None
    # write_file_small(Air_cluster, 5, Air_data_file) # lot of nan data in these files!
    # print('Air data file completed')    

    plt.imshow(PRM_Density.T, origin='lower', cmap='jet', extent=[0,100, 0, 100] )#vmax = -300, vmin = -1500)
    plt.title('PRM density')
    plt.colorbar()
    plt.rc('font' , size = 20)
    plt.savefig('PRMsep_10cm_det_120_2days_poca_100um.png', dpi=150, bbox_inches='tight')	
    plt.show()
    plt.close()

    plt.imshow(Con_cluster.T, origin='lower', cmap='jet', extent=[0,100, 0, 100])#, vmin = 130000, vmax = 165000)
    plt.title('Con cluster')
    plt.colorbar()
    plt.rc('font' , size = 20)
    plt.savefig('concrete_cluster_sep_10cm_det_120_2days_poca_100um.png', dpi=150, bbox_inches='tight')	
    plt.show()
    plt.close()

    plt.imshow(U_cluster.T, origin='lower', cmap='jet', extent=[0,100, 0, 100] )# vmin = 440000, vmax = 550000)#vmax = -300, vmin = -1500)
    plt.title('U cluster')
    plt.colorbar()
    plt.rc('font' , size = 20)
    plt.savefig('uranium_cluster_sep_10cm_det_120_2days_poca_100um.png', dpi=150, bbox_inches='tight')	
    plt.show()
    plt.close()

    plt.imshow(Pb_cluster.T, origin='lower', cmap='jet', extent=[0,100, 0, 100])# vmin = 330000, vmax = 380000 )#vmax = -300, vmin = -1500)
    plt.title('Pb cluster')
    plt.colorbar()
    plt.rc('font' , size = 20)
    plt.savefig('lead_cluster_sep_10cm_det_120_2days_poca_100um.png', dpi=150, bbox_inches='tight')	
    plt.show()
    plt.close()
    
    plt.imshow(SS_cluster.T, origin='lower', cmap='jet', extent=[0,100, 0, 100])
    plt.title('SS cluster')
    plt.colorbar()
    plt.rc('font' , size = 20)
    plt.savefig('SS_cluster_sep_10cm_det_120_2days_poca_100um.png', dpi=150, bbox_inches='tight')	
    plt.show()
    plt.close()

    plt.imshow(Air_cluster.T, origin='lower', cmap='jet', extent=[0,100, 0, 100])
    plt.title('Air cluster')
    plt.colorbar()
    plt.rc('font' , size = 20)
    plt.savefig('air_cluster_sep_10cm_det_120_2days_poca_100um.png', dpi=150, bbox_inches='tight')	
    plt.show()
    plt.close()

    break
