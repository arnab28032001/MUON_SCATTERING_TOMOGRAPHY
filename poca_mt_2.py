import numpy as np
from multiprocessing import Pool
import Lab_3419 as lb

def add_resolution(points, res=0.0):
    # Apply resolution only to x and y (z remains unchanged)
    noise = np.random.normal(0, res, (len(points), 2))
    points[:, :2] += noise
    return points

def fit_3D(points):
    # Fit line: x(z), y(z) using direct slope formula (faster than np.polyfit)
    z = points[:, 2]
    x = points[:, 0]
    y = points[:, 1]

    dz = z[-1] - z[0]
    if dz == 0:  # avoid division by zero
        dz = 1e-10

    m_zx = (x[-1] - x[0]) / dz
    m_zy = (y[-1] - y[0]) / dz

    fit_z = np.array([z[0], z[-1]])
    fit_x = x[0] + m_zx * (fit_z - z[0])
    fit_y = y[0] + m_zy * (fit_z - z[0])

    return np.column_stack((fit_x, fit_y, fit_z))

def POCA_Point(line1, line2):
    P0, P1 = line1
    Q0, Q1 = line2
    u, v = P1 - P0, Q1 - Q0
    w = P0 - Q0

    u_norm = np.linalg.norm(u)
    v_norm = np.linalg.norm(v)
    cos_theta = np.dot(u, v) / (u_norm * v_norm + 1e-12)
    theta_deg = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

    a, b, c = np.dot(u, u), np.dot(u, v), np.dot(v, v)
    d, e = np.dot(u, w), np.dot(v, w)
    denom = a * c - b * b
    if denom == 0.0:
        denom = 1e-10

    sc = (b * e - c * d) / denom
    tc = (a * e - b * d) / denom
    M1 = P0 + sc * u
    M2 = Q0 + tc * v
    M = 0.5 * (M1 + M2)

    return M, theta_deg

def calculate(data, sigma=0.1):
    data_c = np.fromstring(data, sep=' ', count=18)
    rpc_hits = data_c.reshape(6, 3)
    rpc_hits = add_resolution(rpc_hits, sigma)
    top_hits, bottom_hits = rpc_hits[:3], rpc_hits[3:]
    top_line = fit_3D(top_hits)
    bottom_line = fit_3D(bottom_hits)
    (xp, yp, zp), theta_deg = POCA_Point(top_line, bottom_line)
    return xp, yp, zp, theta_deg

data_file = "hits_output.txt"
poca_file = "POCA_output.txt"
print("Data file:", data_file)

with open(data_file, 'r') as file_d:
    data_d = file_d.readlines()

print('Process started...')

# Choose single-thread or multiprocessing based on data size
if len(data_d) > 5000:  # threshold for multiprocessing benefit
    with Pool() as process_pool:
        output_data = process_pool.map(calculate, data_d)
else:
    output_data = list(map(calculate, data_d))

print('Reading completed..')
output_data = np.asarray(output_data)
np.savetxt(poca_file, output_data, fmt='%.4f')
print("Poca file saved:", poca_file)
print('Analysis completed..')
