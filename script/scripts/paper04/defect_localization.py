"""
defect_localization.py
=============================

Numerical/variational verification for Appendix E.6:
  Reduced 1D AX-torsion defect benchmark eta = eta(z) on slice-wise
  homogeneous backgrounds.

  The benchmark operator is
      -d/dz [K_geo(z) df/dz] + M_geo(z) f = E f,
      K_geo(z) = M_geo(z) = c_geo rho(z),
  where current dppu gives the exact mass datum
      M_geo = d^2 V_eff / d eta^2 |_{eta=0},
  and the reduced kinetic datum matches it for the AX contortion gradient.

  Topology-dependent coefficients:
      c_T3 = c_Nil3 = c_Sol3 = 96 pi^4 L / kappa^2
      c_S3 = 12 pi^2 L / kappa^2

  For a Gaussian radial dip rho(z) = rho0 (1 - A exp(-z^2/w^2)),
  the Rayleigh quotient satisfies
      E_var(lambda) / M0
        = (1 + lambda^2) [1 - A lambda w sqrt(pi) exp(lambda^2 w^2) erfc(lambda w)],
      M0 = c_geo rho0,
  so every sampled geometry/profile pair should localize.

Sample profiles (DPPU canonical scan):
   (A, w/rho0) = (0.3, 1.0), (0.5, 2.0), (0.7, 1.0)
Expected outcome:
   4 topologies x 3 profiles = 12/12 LOCALIZED.

Run:
  python -X utf8 script/scripts/paper04/defect_localization.py

Author: Muacca
Date: 2026-04-24
"""

import math
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
from sympy import S, diff, simplify

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

_DATA_DIR = Path(__file__).resolve().parents[3] / "data"
_DATA_DIR.mkdir(exist_ok=True)
from dppu.utils.tee_logger import setup_log
setup_log(__file__, log_dir=str(_DATA_DIR))

from dppu.topology.unified import DOFConfig, TopologyType, UnifiedEngine
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant


TOPOLOGIES = [
    ("T3", TopologyType.T3),
    ("Nil3", TopologyType.NIL3),
    ("S3", TopologyType.S3),
    ("Sol3", TopologyType.SOL3),
]


def eta_geometry_coefficient(topology):
    """
    Return c_geo from the exact homogeneous datum
        M_geo = d^2 V_eff / d eta^2 |_{eta=0} = c_geo * rho.

    Numeric convention used in this script: L = kappa = 1.
    """
    cfg = DOFConfig(
        topology=topology,
        torsion_mode=Mode.AX,
        enable_squash=False,
        enable_shear=False,
        ny_variant=NyVariant.FULL,
    )
    eng = UnifiedEngine(cfg)
    eng.run()

    params = eng.data["params"]
    rho = params.get("R", params["r"])
    veff = eng.data["potential"].subs({
        params["V"]: S.Zero,
        params["theta_NY"]: S.Zero,
    })
    m_geo = simplify(diff(veff, params["eta"], 2).subs(params["eta"], 0))
    c_geo = simplify(m_geo / rho).subs({
        params["L"]: S.One,
        params["kappa"]: S.One,
    })
    return float(c_geo), str(m_geo)


def R_dip(z, R0, A, w):
    """R-dip profile: R(z) = R0 (1 - A exp(-z^2/w^2))."""
    return R0 * (1.0 - A * np.exp(-(z / w) ** 2))


def gaussian_variational_upper_bound(rho0, A, w, c_geo, lam):
    """
    Variational upper bound using the trial function f_lambda = exp(-lambda |z|).

    For the Gaussian dip profile R(z) = R0 (1 - A exp(-z^2 / w^2)) one finds

        E_var(lambda) / M0
          = (1 + lambda^2)
            [1 - A lambda w sqrt(pi) exp(lambda^2 w^2) erfc(lambda w)],

    where M0 = c_geo * rho0 is the continuum threshold.

    As lambda -> 0,

        E_var(lambda) / M0 = 1 - A sqrt(pi) w lambda + O(lambda^2),

    so any Gaussian dip with A > 0 admits E_var(lambda) < M0 for sufficiently
    small lambda, proving the existence of at least one bound state.
    """
    M0 = c_geo * rho0
    x = lam * w
    dip_overlap = w * math.sqrt(math.pi) * math.exp(x * x) * math.erfc(x)
    return M0 * (1.0 + lam * lam) * (1.0 - A * lam * dip_overlap)


def best_variational_upper_bound(rho0, A, w, c_geo):
    """Numerically minimize the explicit variational upper bound."""
    lam_grid = np.logspace(-4, 1, 400) / w
    vals = [gaussian_variational_upper_bound(rho0, A, w, c_geo, float(lam))
            for lam in lam_grid]
    idx = int(np.argmin(vals))
    return float(lam_grid[idx]), float(vals[idx])


def sl_min_eigenvalue(rho_arr, c_geo, dz):
    """
    Smallest eigenvalue of the SL operator
        -d/dz [K_geo(z) df/dz] + M_geo(z) f = E f,  K_geo = M_geo = c_geo * rho.

    Discretization: standard 2nd-order finite difference with K
    evaluated at half-integer points.
    """
    N = len(rho_arr)
    K_arr = c_geo * rho_arr
    M_arr = c_geo * rho_arr

    # K at half-grid points: K_{i+1/2} = (K_i + K_{i+1}) / 2
    Kh = 0.5 * (K_arr[:-1] + K_arr[1:])  # length N-1

    # Diagonal: (K_{i-1/2} + K_{i+1/2}) / dz^2 + M_i
    diag = np.zeros(N)
    diag[1:-1] = (Kh[:-1] + Kh[1:]) / dz**2 + M_arr[1:-1]
    # Dirichlet boundary on edges:
    diag[0]  = Kh[0]  / dz**2 + M_arr[0]
    diag[-1] = Kh[-1] / dz**2 + M_arr[-1]

    # Off-diagonal: -K_{i+1/2} / dz^2
    off = -Kh / dz**2  # length N-1

    H = diags([off, diag, off], offsets=[-1, 0, 1], format="csr")
    eigs = eigsh(H, k=1, which="SA", return_eigenvectors=False)
    return float(eigs[0])


def main():
    print("=" * 64)
    print("Verification: Reduced 1D eta-mode defect localization")
    print("Target: Appendix E.6")
    print("=" * 64)

    rho0 = 3.0

    print()
    print("  Reduced 1D AX benchmark:")
    print("    K_geo(z) = M_geo(z) = c_geo rho(z)")
    print("    c_T3 = c_Nil3 = c_Sol3 = 96 pi^4   (L = kappa = 1)")
    print("    c_S3 = 12 pi^2                     (L = kappa = 1)")
    print(f"  rho0 = {rho0}")
    print("  Variational benchmark for Gaussian R-dips:")
    print("    E_var(lambda) / M0 = 1 - A sqrt(pi) w lambda + O(lambda^2)")
    print("    so any Gaussian dip with A > 0 has at least one bound state.")

    # Grid: z in [-z_max, z_max], z_max = 6 R0, N = 4001
    z_max = 6.0 * rho0
    N = 4001
    z = np.linspace(-z_max, z_max, N)
    dz = z[1] - z[0]

    profiles = [
        ("dip-A03-w1.0R0", 0.3, 1.0),
        ("dip-A05-w2.0R0", 0.5, 2.0),
        ("dip-A07-w1.0R0", 0.7, 1.0),
    ]

    print()
    print(f"  Grid: N = {N}, z in [{-z_max:.1f}, {z_max:.1f}], dz = {dz:.5f}")
    print()
    print("  Exact eta-mass data from current dppu:")
    coeffs = []
    for label, topo in TOPOLOGIES:
        c_geo, m_expr = eta_geometry_coefficient(topo)
        coeffs.append((label, c_geo, m_expr))
        print(f"    {label:4}: c_geo = {c_geo:.4f}   M_geo = {m_expr}")
    print()
    print("-" * 64)
    print(f"  {'topology':8}  {'profile':18}  {'A':4}  {'w/r0':5}  "
          f"{'E_min':12}  {'E_var':12}  {'M_0':12}  {'localized?'}")
    print("-" * 64)

    n_loc = 0
    n_var = 0
    for topo_label, c_geo, _ in coeffs:
        M0 = c_geo * rho0
        for name, A_amp, w_ratio in profiles:
            w = w_ratio * rho0
            rho_arr = R_dip(z, rho0, A_amp, w)
            E_min = sl_min_eigenvalue(rho_arr, c_geo, dz)
            lam_var, E_var = best_variational_upper_bound(rho0, A_amp, w, c_geo)
            localized = E_min < M0
            localized_var = E_var < M0
            flag = "LOCALIZED" if localized else "NOT LOCALIZED"
            print(f"  {topo_label:8}  {name:18}  {A_amp:4}  {w_ratio:5}  "
                  f"{E_min:12.4f}  {E_var:12.4f}  {M0:12.4f}  {flag}")
            print(f"    variational minimizer: lambda*w = {lam_var * w:.4e}   "
                  f"E_var < M0 ? {'YES' if localized_var else 'NO'}")
            if localized:
                n_loc += 1
            if localized_var:
                n_var += 1

    print("-" * 64)
    print()
    total_cases = len(coeffs) * len(profiles)
    print(f"  Localization count: {n_loc} / {total_cases}")
    print(f"  Variational upper-bound count: {n_var} / {total_cases}")
    overall_ok = (n_loc == total_cases == n_var)

    print()
    print("=" * 64)
    if overall_ok:
        print("RESULT: PASS  (4 topologies x 3 Gaussian R-dip profiles: "
              "variational proof + numerical bound state)")
    else:
        print("RESULT: FAIL  "
              f"(num={n_loc}/{total_cases}, var={n_var}/{total_cases})")
    print("=" * 64)
    return overall_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
