# ===================== IMPORTS =====================
import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import convolve2d 

plt.rc('font', size=15)
plt.rcParams["font.family"] = "serif"

# ===================== FILE WRITERS =====================
def write_file(den_data, dev_data, label, filename):
    with open(filename, 'a+') as f:
        for d, v in zip(den_data.flat, dev_data.flat):
            f.write(f"{d} {v:.3f} {label}\n")

# ===================== DIFFUSION MODEL =====================
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

# ===================== LOAD MODEL =====================
device = "cpu"
diff_model = DiffusionUNet().to(device)
diff_model.load_state_dict(torch.load("diffusion_prm.pth", map_location=device))
diff_model.eval()

# ===================== LOAD POCA DATA =====================
file_name = "POCA_output_120cm_2days_100mum.txt"
data = np.loadtxt(file_name)

x, y, z, theta = data[:,0], data[:,1], data[:,2], data[:,3]

# ===================== HISTOGRAM (PRM) =====================
bins = 120
rng = [[-300, 300], [-300, 300]]

density, _, _ = np.histogram2d(x, y, bins=bins, range=rng)
weighted, _, _ = np.histogram2d(x, y, bins=bins, range=rng, weights=theta)

average_deviation = np.divide(
    weighted, density,
    out=np.zeros_like(weighted),
    where=density > 0
)

# ===================== DIFFUSION DENOISING =====================
density_n = density / (density.max() + 1e-6)
dev_n = average_deviation / (average_deviation.max() + 1e-6)

inp = torch.tensor(
    np.stack([density_n, dev_n]),
    dtype=torch.float32
).unsqueeze(0)

with torch.no_grad():
    out = diff_model(inp)

density_denoised = out[0,0].numpy()
deviation_denoised = out[0,1].numpy()

# ===================== FEATURE EXTRACTION =====================
df = pd.DataFrame(density_denoised)
df_d = pd.DataFrame(deviation_denoised)

data_Con = np.array(df.iloc[74:91, 76:95])
data_U   = np.array(df.iloc[50:69, 20:39])
data_Pb  = np.array(df.iloc[49:69, 81:100])
data_SS  = np.array(df.iloc[18:39, 50:69])
data_Air = np.array(df.iloc[79:89, 54:67])

data_Con_d = np.array(df_d.iloc[74:91, 76:95])
data_U_d   = np.array(df_d.iloc[50:69, 20:39])
data_Pb_d  = np.array(df_d.iloc[49:69, 81:100])
data_SS_d  = np.array(df_d.iloc[18:39, 50:69])
data_Air_d = np.array(df_d.iloc[79:89, 54:67])

# ===================== WRITE OUTPUT =====================
write_file(data_Con, data_Con_d, 1, "Cluster_DIFF_PRM.txt")
write_file(data_U,   data_U_d,   2, "Cluster_DIFF_PRM.txt")
write_file(data_Pb,  data_Pb_d,  3, "Cluster_DIFF_PRM.txt")
write_file(data_SS,  data_SS_d,  4, "Cluster_DIFF_PRM.txt")
write_file(data_Air, data_Air_d, 5, "Cluster_DIFF_PRM.txt")

# ===================== PRM AFTER DIFFUSION =====================

# Use denoised maps directly
PRM_Density = density_denoised.copy()
PRM_Deviation = deviation_denoised.copy()

# Optional: mask empty regions
PRM_Density[PRM_Density <= 0] = np.nan
PRM_Deviation[PRM_Deviation <= 0] = np.nan

# ----------------- PRM FILTERING -----------------
# Global PRM kernel (same philosophy as original code)
ker_All = -np.ones((3, 3))

PRM_response = convolve2d(
    np.nan_to_num(PRM_Density),
    ker_All,
    mode='same'
)

# ----------------- Thresholding -----------------
cutoff = np.nanmean(PRM_response)

PRM_mask = PRM_response < cutoff

PRM_Density_filtered = np.where(PRM_mask, PRM_Density, np.nan)
PRM_Deviation_filtered = np.where(PRM_mask, PRM_Deviation, np.nan)

# ----------------- VISUALIZATION -----------------
plt.figure(figsize=(6,5))
plt.imshow(PRM_Density_filtered.T, origin='lower', cmap='jet')
plt.colorbar(label="PRM Density (Diffusion)")
plt.title("Diffusion and PRM Density")
plt.savefig("Cluster_Density_Diffusion_plus_PRM.png")
plt.show()

plt.figure(figsize=(6,5))
plt.imshow(PRM_Deviation_filtered.T, origin='lower', cmap='jet')
plt.colorbar(label="PRM Deviation (Diffusion)")
plt.title("Diffusion and PRM Deviation")
plt.savefig("Avg_Deviation_Diffusion_plus_PRM.png")
plt.show()


