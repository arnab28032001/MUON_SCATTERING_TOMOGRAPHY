import torch
import torch.nn as nn
import numpy as np

#Model 
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

#Fake PRM generator
def generate_clean_prm(size=120):
    x, y = np.meshgrid(
        np.linspace(-1, 1, size),
        np.linspace(-1, 1, size)
    )

    density = np.exp(-4 * (x**2 + y**2))
    deviation = density * 0.8

    return density, deviation

def add_noise(density, deviation):
    noisy_density = density + np.random.normal(0, 0.1, density.shape)
    noisy_deviation = deviation + np.random.normal(0, 0.1, deviation.shape)
    return noisy_density, noisy_deviation

#Training
device = "cpu"
model = DiffusionUNet().to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
loss_fn = nn.MSELoss()

for epoch in range(200):
    clean_d, clean_dev = generate_clean_prm()
    noisy_d, noisy_dev = add_noise(clean_d, clean_dev)

    inp = torch.tensor(
        np.stack([noisy_d, noisy_dev]),
        dtype=torch.float32
    ).unsqueeze(0)

    target = torch.tensor(
        np.stack([clean_d, clean_dev]),
        dtype=torch.float32
    ).unsqueeze(0)

    optimizer.zero_grad()
    out = model(inp)
    loss = loss_fn(out, target)
    loss.backward()
    optimizer.step()

    if epoch % 20 == 0:
        print(f"Epoch {epoch}, Loss = {loss.item():.6f}")

#Save model 
torch.save(model.state_dict(), "diffusion_prm.pth")
print("Saved diffusion_prm.pth")

