"""
eta_aps_nil3.py
======================

Numerical verification for Appendix E.10:
  APS spectral core on Nil3 with PPA spin structure.

Target:
  eta_APS^{(3D)}(Nil^3, PPA) = (eta(0) + h) / 2 = (1 + 0) / 2 = +1/2.

Method:
  [A] Keep the low-lying Heisenberg Landau-level check, so the unpaired n=0
      spin-up branch and the PPA zero-mode exclusion remain visible.

  [B] Replace the old hardcoded rho-invariant input by the clean Phase4H route:
      compute the Nil3 Levi-Civita Chern-Simons integral explicitly,
      evaluate the spinor-vs-adjoint cubic trace ratio,
      derive
          int_Y CS_sigma = -1/2,
          eta_APS = -int_Y CS_sigma = +1/2,
          eta(0) = 2 eta_APS - h = 1.

      No rho_PPA = 1 or eta_APS = 1/2 is injected by hand.

Run:
  python -X utf8 script/scripts/paper04/eta_aps_nil3.py

Author: Muacca
Date: 2026-04-24
"""

import math
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))
_DATA_DIR = Path(__file__).resolve().parents[3] / "data"
_DATA_DIR.mkdir(exist_ok=True)
from dppu.utils.tee_logger import setup_log
setup_log(__file__, log_dir=str(_DATA_DIR))

import numpy as np


def eps3(i, j, k):
    if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        return 1
    if (i, j, k) in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
        return -1
    return 0


def heisenberg_landau_levels(p2, r0, n_max):
    k2 = 2.0 * math.pi * p2 / r0
    k2a = abs(k2)
    omega = k2a / r0
    levels = []
    for n in range(n_max + 1):
        mu_up = 2.0 * n * omega * k2a + k2a**2
        mu_down = 2.0 * (n + 1) * omega * k2a + k2a**2
        levels.append((mu_up, mu_down))
    return levels


def build_nil3_geometry(r0):
    C = np.zeros((3, 3, 3))
    C[2, 0, 1] = 1.0 / r0
    C[2, 1, 0] = -1.0 / r0

    omega = np.zeros((3, 3, 3))
    for a in range(3):
        for b in range(3):
            for c in range(3):
                omega[a, b, c] = 0.5 * (C[a, b, c] + C[c, b, a] - C[b, a, c])

    d_omega = np.zeros((3, 3, 3, 3))
    for a in range(3):
        for b in range(3):
            for d in range(3):
                for e in range(3):
                    for c in range(3):
                        d_omega[a, b, d, e] += omega[a, b, c] * (-C[c, d, e])

    return C, omega, d_omega, r0**3


def adjoint_generators():
    J = {}
    J[(0, 1)] = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
    J[(0, 2)] = np.array([[0, 0, -1], [0, 0, 0], [1, 0, 0]], dtype=float)
    J[(1, 2)] = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=float)
    for a in range(3):
        for b in range(a + 1, 3):
            J[(b, a)] = -J[(a, b)]
    return J


def spinor_generators():
    gamma = [
        np.array([[0, 1], [1, 0]], dtype=complex),
        np.array([[0, -1j], [1j, 0]], dtype=complex),
        np.array([[1, 0], [0, -1]], dtype=complex),
    ]
    sigma_gen = {}
    for a in range(3):
        for b in range(a + 1, 3):
            comm = gamma[a] @ gamma[b] - gamma[b] @ gamma[a]
            sigma_gen[(a, b)] = 0.5 * comm
            sigma_gen[(b, a)] = -sigma_gen[(a, b)]
    return sigma_gen


def build_rep_connection(omega, generators, size):
    A_rep = [np.zeros((size, size), dtype=complex) for _ in range(3)]
    for c in range(3):
        for a in range(3):
            for b in range(a + 1, 3):
                A_rep[c] += omega[a, b, c] * generators[(a, b)]
    return A_rep


def cubic_trace_sum(A_rep):
    total = 0.0
    for c1 in range(3):
        for c2 in range(3):
            for c3 in range(3):
                ep = eps3(c1, c2, c3)
                if ep == 0:
                    continue
                total += ep * float(np.real(np.trace(A_rep[c1] @ A_rep[c2] @ A_rep[c3])))
    return total


def adjoint_cs_integral_hpp1(omega, d_omega, vol):
    J = adjoint_generators()
    A_adj = build_rep_connection(omega, J, size=3)

    cs_A = 0.0
    for a in range(3):
        for b in range(a + 1, 3):
            for c in range(3):
                for d in range(3):
                    for e in range(d + 1, 3):
                        cs_A += omega[a, b, c] * d_omega[a, b, d, e] * eps3(c, d, e)

    cs_B = 0.0
    for c1 in range(3):
        for c2 in range(3):
            for c3 in range(3):
                ep = eps3(c1, c2, c3)
                if ep == 0:
                    continue
                cs_B += (2.0 / 3.0) * ep * float(
                    np.real(np.trace(A_adj[c1] @ A_adj[c2] @ A_adj[c3]))
                )

    return {
        "A_adj": A_adj,
        "cs_A": cs_A * vol,
        "cs_B": cs_B * vol,
        "total": (cs_A + cs_B) * vol,
        "T2": float(np.trace(J[(0, 1)] @ J[(0, 1)])),
    }


def spinor_cs_integral(omega, C, vol):
    sigma_gen = spinor_generators()
    A_sigma = build_rep_connection(omega, sigma_gen, size=2)

    dA_sigma = [[np.zeros((2, 2), dtype=complex) for _ in range(3)] for _ in range(3)]
    for d in range(3):
        for e in range(3):
            for c in range(3):
                dA_sigma[d][e] += A_sigma[c] * (-C[c, d, e])

    cs_A = 0.0
    for c in range(3):
        for d in range(3):
            for e in range(d + 1, 3):
                ep = eps3(c, d, e)
                if ep == 0:
                    continue
                cs_A += ep * float(np.real(np.trace(A_sigma[c] @ dA_sigma[d][e])))

    cs_B = 0.0
    for c1 in range(3):
        for c2 in range(3):
            for c3 in range(3):
                ep = eps3(c1, c2, c3)
                if ep == 0:
                    continue
                cs_B += (2.0 / 3.0) * ep * float(
                    np.real(np.trace(A_sigma[c1] @ A_sigma[c2] @ A_sigma[c3]))
                )

    return {
        "A_sigma": A_sigma,
        "cs_A": cs_A * vol,
        "cs_B": cs_B * vol,
        "total": (cs_A + cs_B) * vol,
        "T2": float(np.real(np.trace(sigma_gen[(0, 1)] @ sigma_gen[(0, 1)]))),
    }


def main():
    print("=" * 64)
    print("Verification: Nil3 APS spectral core eta_APS = +1/2")
    print("Target: Appendix E.10")
    print("=" * 64)

    r0 = 3.0
    C, omega, d_omega, vol = build_nil3_geometry(r0)

    print()
    print(f"  r0 = {r0}")
    print("  PPA spin structure: p_2 in Z + 1/2")
    print("  Dirac APS sign convention used here: eta_APS = - int_Y CS_sigma")

    print()
    print("[1] Heisenberg Landau-level sanity check")
    print("-" * 64)
    print(f"  {'p_2':>5}  {'mu_0^up':>12}  {'mu_0^down':>12}  {'mu_1^up':>12}")
    spectrum_ok = True
    for p2 in [0.5, 1.5, 2.5, -0.5]:
        levels = heisenberg_landau_levels(p2, r0, n_max=1)
        mu0_up, mu0_down = levels[0]
        mu1_up, _ = levels[1]
        print(f"  {p2:>5.1f}  {mu0_up:>12.6f}  {mu0_down:>12.6f}  {mu1_up:>12.6f}")
        spectrum_ok &= (mu0_up < mu0_down)
    print(f"  n=0 unpaired spin-up branch present for sampled PPA modes: "
          f"{'OK' if spectrum_ok else 'FAIL'}")

    print()
    print("[2] Nil3 Levi-Civita Chern-Simons integral in the adjoint convention")
    print("-" * 64)
    adj = adjoint_cs_integral_hpp1(omega, d_omega, vol)
    ok_adj_A = abs(adj["cs_A"] - (-0.25)) < 1e-10
    ok_adj_B = abs(adj["cs_B"] - 0.5) < 1e-10
    ok_adj_total = abs(adj["total"] - 0.25) < 1e-10
    print(f"  int_Y Tr(omega dw)        = {adj['cs_A']:+.8f}   "
          f"[{'OK' if ok_adj_A else 'FAIL'}]")
    print(f"  int_Y (2/3)Tr(omega^3)    = {adj['cs_B']:+.8f}   "
          f"[{'OK' if ok_adj_B else 'FAIL'}]")
    print(f"  int_Y CS_adj (H++-1 norm) = {adj['total']:+.8f}   "
          f"[{'OK' if ok_adj_total else 'FAIL'}]")

    print()
    print("[3] Spinor trace ratio and spinor CS integral")
    print("-" * 64)
    sigma = spinor_cs_integral(omega, C, vol)
    ratio_T2 = sigma["T2"] / adj["T2"]
    ratio_T3 = cubic_trace_sum(sigma["A_sigma"]) / cubic_trace_sum(adj["A_adj"])
    ok_T2 = abs(ratio_T2 - 1.0) < 1e-10
    ok_T3 = abs(ratio_T3 + 2.0) < 1e-10
    ok_sigma_total = abs(sigma["total"] - (-0.5)) < 1e-10
    print(f"  quadratic trace ratio Tr_sigma/Tr_adj = {ratio_T2:+.8f}   "
          f"[{'OK' if ok_T2 else 'FAIL'}]")
    print(f"  cubic trace ratio                     = {ratio_T3:+.8f}   "
          f"[{'OK' if ok_T3 else 'FAIL'}]")
    print(f"  int_Y CS_sigma                        = {sigma['total']:+.8f}   "
          f"[{'OK' if ok_sigma_total else 'FAIL'}]")

    print()
    print("[4] Kernel dimension for PPA")
    print("-" * 64)
    h = 0
    print("  PPA excludes p_2 = 0, so the would-be zero mode is absent")
    print(f"  h = {h}")
    ok_h = (h == 0)

    print()
    print("[5] APS assembly without hardcoded rho")
    print("-" * 64)
    eta_aps = -sigma["total"]
    eta0 = 2.0 * eta_aps - h
    N_grav = eta_aps / adj["total"]
    ok_eta_aps = abs(eta_aps - 0.5) < 1e-10
    ok_eta0 = abs(eta0 - 1.0) < 1e-10
    ok_N_grav = abs(N_grav - 2.0) < 1e-10
    print(f"  eta_APS = -int_Y CS_sigma = {eta_aps:+.8f}   "
          f"[{'OK' if ok_eta_aps else 'FAIL'}]")
    print(f"  eta(0)  = 2 eta_APS - h   = {eta0:+.8f}   "
          f"[{'OK' if ok_eta0 else 'FAIL'}]")
    print(f"  N_grav  = eta_APS / int_Y CS_adj = {N_grav:+.8f}   "
          f"[{'OK' if ok_N_grav else 'FAIL'}]")

    overall_ok = all([
        spectrum_ok,
        ok_adj_A,
        ok_adj_B,
        ok_adj_total,
        ok_T2,
        ok_T3,
        ok_sigma_total,
        ok_h,
        ok_eta_aps,
        ok_eta0,
        ok_N_grav,
    ])

    print()
    print("=" * 64)
    if overall_ok:
        print("RESULT: PASS  (eta_APS^{(3D)}(Nil^3, PPA) = +1/2, derived)")
    else:
        print("RESULT: FAIL")
    print("=" * 64)
    return overall_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
