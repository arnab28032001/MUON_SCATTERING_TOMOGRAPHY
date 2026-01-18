import numpy as np
from multiprocessing import Pool
#import Lab_3419 as lb

def add_resolusion(points_, res_=0.0):
    data_x_, data_y_, data_z_ = points_.T
    data_x_e_ = data_x_ + np.random.normal(0, res_, len(data_x_))
    data_y_e_ = data_y_ + np.random.normal(0, res_, len(data_y_))
    data_z_e_ = data_z_
    data_xyz_e_ = np.array([data_x_e_, data_y_e_, data_z_e_]).T
    return data_xyz_e_

def fit_3D(points_):
    data_x_, data_y_, data_z_ = points_.T
    fit_z_ = np.array([data_z_[0], data_z_[-1]])
    m_zx_, c_zx_, *_ = np.polyfit(data_z_, data_x_, 1)
    m_zy_, c_zy_, *_ = np.polyfit(data_z_, data_y_, 1)
    fit_x_ = m_zx_ * fit_z_ + c_zx_
    fit_y_ = m_zy_ * fit_z_ + c_zy_
    fit_xyz_ = np.array([fit_x_, fit_y_, fit_z_]).T
    return fit_xyz_
    
def POCA_Point(line_1_, line_2_):
    P0_, P1_ = line_1_[0], line_1_[1]
    Q0_, Q1_ = line_2_[0], line_2_[1]
    u_, v_, w_ = (P1_ - P0_), (Q1_ - Q0_), (P0_ - Q0_)
    cos_theta_ = np.dot(u_, v_) / (np.linalg.norm(u_) * np.linalg.norm(v_))
    theta_rad_ = np.arccos(np.clip(cos_theta_, -1.0, 1.0))
    theta_deg_ = np.rad2deg(theta_rad_)
    a_, b_, c_ = np.dot(u_, u_), np.dot(u_, v_), np.dot(v_, v_)
    d_, e_, f_ = np.dot(u_, w_), np.dot(v_, w_), (a_ * c_ - b_ * b_)
    if f_ == 0.0: f_ = 1.0e-10  # to avoid zero division error
    sc_, tc_ = (b_ * e_ - c_ * d_) / f_, (a_ * e_ - b_ * d_) / f_
    M1_, M2_ = (P0_ + sc_ * u_), (Q0_ + tc_ * v_)
    M_ = (M1_ + M2_) / 2.0
    return M_, theta_deg_


def calculate(data_, sigma_=0.100):
    data_c_ = np.array([float(_) for _ in data_.split()])
    data_c_ = data_c_[0:18]
    rpc_hits_ = data_c_.reshape(6, 3)
#    sigma_ = 0.100
    rpc_hits_ = add_resolusion(rpc_hits_, sigma_)
    top_hits_, bottom_hits_ = rpc_hits_[0:3], rpc_hits_[3:6]
    top_line_ = fit_3D(top_hits_)
    bottom_line_ = fit_3D(bottom_hits_)
    (xp_, yp_, zp_), theta_deg_ = POCA_Point(top_line_, bottom_line_)
    return xp_, yp_, zp_, theta_deg_

data_file = "Arnab_Hits_120cm_2days.txt"
poca_file = "POCA_output_120cm_2days_100mum.txt"
print("Data file :", data_file)

file_d = open(data_file, 'r')
data_d = file_d.readlines()
file_d.close()

print('Process started...')
process_pool = Pool()
output_data = process_pool.map(calculate, data_d)
process_pool.close()
process_pool.join()
print('Reading completed ..')
output_data = np.array(output_data)
np.savetxt(poca_file, output_data, fmt='%.4f')
print("Poca file saved:", poca_file)
print('Analysis completed ..')
