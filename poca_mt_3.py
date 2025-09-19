import numpy as np
from multiprocessing import Pool
import Lab_3419 as lb

def add_resolution(points, res=0.0):
    # Only modify x and y; z stays unchanged
    if res > 0:
        points[:, :2] += np.random.normal(0, res, points[:, :2].shape)
    return points

def fit_3D(points):
    # Direct two-point linear fit for x(z) and y(z)
    z0, z1 = points[0, 2], points[-1, 2]
    x0, x1 = points[0, 0], points[-1, 0]
    y0, y1 = points[0, 1], points[-1, 1]

    dz = z1 - z0 if z1 != z0 else 1e-10
    m_zx = (x1 - x0) / dz
    m_zy = (y1 - y0) / dz

    fit_z = np.array([z0, z1])
    fit_x = x0 + m_zx * (fit_z - z0)
    fit_y = y0 + m_zy * (fit_z - y0)
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

def calculate_line(data, sigma=0.1):
    # Convert line to 18 floats efficiently
    data_c = np.fromstring(data, sep=' ', count=18)
    if data_c.size != 18:
        return np.nan, np.nan, np.nan, np.nan  # skip malformed lines
    rpc_hits = data_c.reshape(6, 3)
    rpc_hits = add_resolution(rpc_hits, sigma)
    top_line = fit_3D(rpc_hits[:3])
    bottom_line = fit_3D(rpc_hits[3:])
    (xp, yp, zp), theta_deg = POCA_Point(top_line, bottom_line)
    return xp, yp, zp, theta_deg

def process_file(input_file, output_file, sigma=0.1, chunk_size=100000, use_multiprocessing=True):
    print("Processing:", input_file)
    results = []

    with open(input_file, 'r') as f:
        lines_buffer = f.readlines(chunk_size)  # read chunk by chunk
        while lines_buffer:
            if use_multiprocessing and len(lines_buffer) > 5000:
                with Pool() as pool:
                    chunk_results = pool.map(calculate_line, lines_buffer)
            else:
                chunk_results = map(calculate_line, lines_buffer)

            results.extend(chunk_results)
            lines_buffer = f.readlines(chunk_size)

    results = np.array(results, dtype=np.float32)
    np.savetxt(output_file, results, fmt='%.4f')
    print("Saved POCA file:", output_file)
    return results

# Run
data_file = "hits_output.txt"
poca_file = "POCA_output.txt"
process_file(data_file, poca_file, sigma=0.1)

