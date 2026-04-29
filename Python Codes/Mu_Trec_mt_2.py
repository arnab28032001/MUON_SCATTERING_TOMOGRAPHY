from __future__ import annotations

import numpy as np
import jax
import jax.numpy as jnp
from typing import Tuple
from dataclasses import dataclass


# ============================================================
# PARAMETERS (No external import → fixes ReconParams issue)
# ============================================================

@dataclass(frozen=True)
class ReconParams:
    start_point: Tuple[float, float, float]
    cube_size: float
    voxel_size: float
    X0: float
    p_muon: float
    DE: float
    const1: float


# ============================================================
# NUMPY PART (GEOMETRY + DETECTOR HANDLING)
# ============================================================

def add_resolution(points_, res_=0.0):
    x, y, z = points_.T
    x_e = x + np.random.normal(0, res_, len(x))
    y_e = y + np.random.normal(0, res_, len(y))
    return np.array([x_e, y_e, z]).T


def fit_3D(points_):
    # Sort by z (important for consistency)
    points_ = points_[np.argsort(points_[:, 2])]

    x, y, z = points_.T

    # Linear fit: x(z), y(z)
    m_zx, c_zx = np.polyfit(z, x, 1)
    m_zy, c_zy = np.polyfit(z, y, 1)

    z_fit = np.array([z[0], z[-1]])

    x_fit = m_zx * z_fit + c_zx
    y_fit = m_zy * z_fit + c_zy

    return np.array([x_fit, y_fit, z_fit]).T  # (2,3)


def prepare_gmte_input(data_line, sigma_=0.001):
    data_c = np.array([float(_) for _ in data_line.split()])
    data_c = data_c[0:18]

    rpc_hits = data_c.reshape(6, 3)

    # Add detector resolution
    rpc_hits = add_resolution(rpc_hits, sigma_)

    # Split detectors
    top_hits = rpc_hits[0:3]
    bottom_hits = rpc_hits[3:6]

    # Fit tracks using all detectors
    top_line = fit_3D(top_hits)
    bottom_line = fit_3D(bottom_hits)

    p1, p2 = top_line[0], top_line[1]
    p3, p4 = bottom_line[0], bottom_line[1]

    return p1, p2, p3, p4


def prepare_batch(data_lines, sigma_=0.001):
    p1_list, p2_list, p3_list, p4_list = [], [], [], []

    for line in data_lines:
        p1, p2, p3, p4 = prepare_gmte_input(line, sigma_)
        p1_list.append(p1)
        p2_list.append(p2)
        p3_list.append(p3)
        p4_list.append(p4)

    return (
        jnp.array(p1_list),
        jnp.array(p2_list),
        jnp.array(p3_list),
        jnp.array(p4_list),
    )


# ============================================================
# JAX PART (GMTE CORE)
# ============================================================

def _angles_xy(p1, p2):
    dz = p2[..., 2] - p1[..., 2]
    dz = jnp.where(dz == 0, 1e-12, dz)

    theta_y = jnp.arctan2(p2[..., 1] - p1[..., 1], dz)
    theta_x = jnp.arctan2(p2[..., 0] - p1[..., 0], dz)

    return theta_x, theta_y


def scattering_angle_signal(p1, p2, p3, p4):
    entry_x, entry_y = _angles_xy(p1, p2)
    exit_x, exit_y = _angles_xy(p3, p4)

    dx = jnp.abs(exit_x - entry_x)
    dy = jnp.abs(exit_y - entry_y)

    return jnp.sqrt((dx**2 + dy**2) / 2.0)


def _variance_terms(xeta, X0, p_muon, DE, const1):
    ratio = jnp.maximum(xeta / X0, 1e-12)
    denom = jnp.maximum(p_muon**2 - DE * p_muon * xeta, 1e-12)

    variance_theta = const1 * (xeta / X0) * (1.0 / denom) * (1.0 + 0.038 * jnp.log(ratio))**2
    variance_y = (xeta**2) * variance_theta / 3.0
    variance_th_y = 0.5 * jnp.sqrt(3.0) * jnp.sqrt(variance_theta * variance_y)

    return variance_theta, variance_y, variance_th_y


def mutrec_reshma(p1, p2, p3, p4, params: ReconParams):

    start_z = params.start_point[2]
    T = int(round(params.cube_size / params.voxel_size))

    z0 = start_z - 1e-3
    z = z0 + params.voxel_size * jnp.arange(T)

    entry_x, entry_y = _angles_xy(p1, p2)
    exit_x, exit_y = _angles_xy(p3, p4)

    YA_y = jnp.stack([p2[1], entry_y])
    YB_y = jnp.stack([p3[1], exit_y])
    YA_x = jnp.stack([p2[0], entry_x])
    YB_x = jnp.stack([p3[0], exit_x])

    RA = jnp.stack([
        jnp.stack([jnp.ones_like(z), z - p2[2]], axis=-1),
        jnp.stack([jnp.zeros_like(z), jnp.ones_like(z)], axis=-1)
    ], axis=-2)

    RB = jnp.stack([
        jnp.stack([jnp.ones_like(z), p3[2] - z], axis=-1),
        jnp.stack([jnp.zeros_like(z), jnp.ones_like(z)], axis=-1)
    ], axis=-2)

    # Covariance matrices
    xeta1 = z - p1[2]
    vt1, vy1, vth1 = _variance_terms(xeta1, params.X0, params.p_muon, params.DE, params.const1)

    S1 = jnp.stack([
        jnp.stack([vy1, vth1], axis=-1),
        jnp.stack([vth1, vt1], axis=-1)
    ], axis=-2)

    xeta2 = p4[2] - z
    vt2, vy2, vth2 = _variance_terms(xeta2, params.X0, params.p_muon, params.DE, params.const1)

    S2 = jnp.stack([
        jnp.stack([vy2, vth2], axis=-1),
        jnp.stack([vth2, vt2], axis=-1)
    ], axis=-2)

    def solve_component(YA, YB):
        YA_b = jnp.broadcast_to(YA, (z.shape[0], 2))
        YB_b = jnp.broadcast_to(YB, (z.shape[0], 2))

        invS1 = jnp.linalg.inv(S1)
        S2_RB = jnp.linalg.solve(S2, RB)
        RB_T = jnp.swapaxes(RB, -1, -2)

        M = invS1 + RB_T @ S2_RB

        rhs1 = invS1 @ (RA @ YA_b[..., None])
        rhs2 = RB_T @ (jnp.linalg.solve(S2, YB_b[..., None]))

        rhs = (rhs1 + rhs2)[..., 0]

        out = jnp.linalg.solve(M, rhs[..., None]).squeeze(-1)

        return out[..., 0], out[..., 1]

    y_gmte, theta_y = solve_component(YA_y, YB_y)
    x_gmte, theta_x = solve_component(YA_x, YB_x)

    return x_gmte, y_gmte, z, theta_x, theta_y


# ============================================================
# BATCHED VERSION (FAST)
# ============================================================

_mutrec_vmap = jax.jit(
    jax.vmap(mutrec_reshma, in_axes=(0, 0, 0, 0, None)),
    static_argnums=(4,)
)


def mutrec_reshma_batch(p1, p2, p3, p4, params):
    x, y, z, thx, thy = _mutrec_vmap(p1, p2, p3, p4, params)
    return x, y, z[0], thx, thy


# ============================================================
# MAIN EXECUTION
# ============================================================

if __name__ == "__main__":

    data_file = "Arnab_Hits_120cm_2days.txt"

    print("Reading data...")
    with open(data_file, "r") as f:
        data_lines = f.readlines()

    print("Preparing GMTE inputs...")
    p1, p2, p3, p4 = prepare_batch(data_lines, sigma_=0.001)

    # IMPORTANT: keep units consistent (MeV here)
    params = ReconParams(
        start_point=(0.0, 0.0, 45.0),
        cube_size=90.0,
        voxel_size=2.0,
        X0=10.0,
        p_muon=3000.0,   # MeV
        DE=1.5,
        const1=13.6**2   # (MeV)^2
    )

    print("Running GMTE reconstruction...")
    x, y, z, thx, thy = mutrec_reshma_batch(p1, p2, p3, p4, params)

    print("Reconstruction complete.")

    # Extract midpoint
    mid_idx = len(z) // 2

    # Scattering signal
    s = scattering_angle_signal(p1, p2, p3, p4)

    output = np.column_stack([
        np.array(x[:, mid_idx]),
        np.array(y[:, mid_idx]),
        np.full(len(x), z[mid_idx]),
        np.array(s)
    ])

    np.savetxt("GMTE_output.txt", output, fmt="%.4f")

    print("Saved GMTE output.")