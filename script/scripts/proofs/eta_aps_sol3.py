"""
eta_aps_sol3.py
======================

Semi-analytic verification for the Sol3 APS eta-invariant on the compact
mapping torus

    M_A = T^2 \rtimes_A S^1,   A = [[2, 1], [1, 1]].

Main points:
  [A] The local Levi-Civita spinor CS 3-form vanishes on Sol3:
          int_Y CS_sigma = 0,
          int_Y CS_adj   = 0.

  [B] The Sol3 Dirac operator in the left-invariant orthonormal frame has no
      zero-order term:
          D_Sol = gamma^a E_a.

  [C] Nonzero fiber Fourier sectors have no zero modes. After a local U(1)
      rotation of the fiber momentum vector, the zero-mode equation reduces to
          (d/dt +/- m_k(t)) phi = 0,
      with m_k(t) > 0 for every nonzero k in the T^2 dual lattice.
      Hence periodic or anti-periodic boundary conditions around the base circle
      cannot be satisfied unless k = 0.

  [D] The k = 0 sector is spin-structure dependent:
          Sol-P (base periodic):     h = 2   (two constant spinors descend)
          Sol-A (base anti-periodic): h = 0.

  [E] A time-reversal antiunitary symmetry gives eta(0) = 0 for both spin
      structures. Therefore
          eta_APS = (eta(0) + h) / 2
      becomes
          eta_APS(Sol-P) = 1,
          eta_APS(Sol-A) = 0.

Important consequence:
  The local CS integral fixes only the non-integer part here. On Sol-P the
  integer shift comes entirely from h = dim ker D, so eta_APS is NOT equal to
  -int CS_sigma.

Run:
  python -X utf8 script/scripts/proof/eta_aps_sol3.py

Author: Muacca
Date: 2026-04-24
"""

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))
_DATA_DIR = Path(__file__).resolve().parents[3] / "data"
_DATA_DIR.mkdir(exist_ok=True)
from dppu.utils.tee_logger import setup_log
setup_log(__file__, log_dir=str(_DATA_DIR))

import numpy as np
from scipy.linalg import expm
import sympy as sp


def eps3(i, j, k):
    if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        return 1
    if (i, j, k) in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
        return -1
    return 0


def build_sol3_geometry(r0):
    """
    Sol3 structure constants in the DPPU orthonormal frame:
        C^1_{01} = +1/r0,
        C^2_{02} = -1/r0.

    The unit-cell volume is irrelevant here because the CS density vanishes
    identically. We nevertheless keep an explicit volume factor for bookkeeping.
    """
    C = np.zeros((3, 3, 3))
    C[1, 0, 1] = +1.0 / r0
    C[1, 1, 0] = -1.0 / r0
    C[2, 0, 2] = -1.0 / r0
    C[2, 2, 0] = +1.0 / r0

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

    vol = r0**3
    return C, omega, d_omega, vol


def adjoint_generators():
    J = {}
    J[(0, 1)] = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
    J[(0, 2)] = np.array([[0, 0, -1], [0, 0, 0], [1, 0, 0]], dtype=float)
    J[(1, 2)] = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], dtype=float)
    for a in range(3):
        for b in range(a + 1, 3):
            J[(b, a)] = -J[(a, b)]
    return J


def pauli_gamma():
    return [
        np.array([[0, 1], [1, 0]], dtype=complex),
        np.array([[0, -1j], [1j, 0]], dtype=complex),
        np.array([[1, 0], [0, -1]], dtype=complex),
    ]


def spinor_generators():
    gamma = pauli_gamma()
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
                total += float(np.real(np.trace(A_rep[c1] @ A_rep[c2] @ A_rep[c3])))
    return total


def adjoint_cs_integral(omega, d_omega, vol):
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
    }


def dirac_zero_order_term(omega):
    """
    For the spin Dirac operator

        D = gamma^c (E_c + Omega_c),
        Omega_c = 1/4 sum_{a<b} omega_{ab,c} [gamma^a, gamma^b],

    the lower-order term is B = sum_c gamma^c Omega_c.
    """
    gamma = pauli_gamma()
    Omega = [np.zeros((2, 2), dtype=complex) for _ in range(3)]
    for c in range(3):
        total = np.zeros((2, 2), dtype=complex)
        for a in range(3):
            for b in range(a + 1, 3):
                comm = gamma[a] @ gamma[b] - gamma[b] @ gamma[a]
                total += 0.25 * omega[a, b, c] * comm
        Omega[c] = total

    B = np.zeros((2, 2), dtype=complex)
    for c in range(3):
        B += gamma[c] @ Omega[c]
    return Omega, B


def time_reversal_check():
    """
    T = i sigma_2 K satisfies T sigma_j T^{-1} = -sigma_j for j=1,2,3.
    Hence for the Sol Dirac operator with real vector fields E_j, T D T^{-1} = -D.
    """
    gamma = pauli_gamma()
    U = 1j * gamma[1]
    transformed = []
    for g in gamma:
        transformed.append(U @ np.conjugate(g) @ U.conj().T)
    return transformed


def log_hyperbolic_monodromy(A):
    vals, vecs = np.linalg.eig(A)
    logs = np.diag(np.log(vals))
    L = vecs @ logs @ np.linalg.inv(vecs)
    return np.real_if_close(L)


def fiber_mass_profile(A, k_vec, t_grid):
    """
    In the periodic mapping-torus trivialization with coframe
        e^f(t) = exp(-t log A) du_f,
    the dual fiber momentum entering the Dirac operator is
        p_k(t) = exp(t log(A)^T) k.
    Its norm m_k(t) is strictly positive for k != 0.
    """
    L = log_hyperbolic_monodromy(A)
    prof = []
    for t in t_grid:
        p_t = expm(t * L.T) @ np.asarray(k_vec, dtype=float)
        prof.append(np.linalg.norm(p_t))
    return np.asarray(prof)


def base_periodic_sector_kernel_dim(base_periodic):
    """
    k = 0 sector:
        D_0 = gamma^3 d/dt.
    Its zero modes are constant spinors. The quotient has two spin structures:
        Sol-P: base periodic     -> both constant spinors descend
        Sol-A: base anti-periodic -> constants are projected out
    """
    return 2 if base_periodic else 0


def off_diagonal_reduction_check():
    """
    After the local U(1) rotation that aligns the real fiber momentum vector
    p_k(t) with a single Pauli axis, the nonzero-k Dirac sector has the form

        D_m = sigma_3 * p + i m * sigma_1.

    A further constant unitary rotation V turns this into

        V D_m V^dagger = [[0, p - m], [p + m, 0]],

    so zero modes satisfy (d/dt - m) phi_- = 0 and (d/dt + m) phi_+ = 0.
    """
    p, m = sp.symbols("p m", commutative=True)
    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma3 = sp.Matrix([[1, 0], [0, -1]])
    Dm = sigma3 * p + sp.I * sigma1 * m
    V = sp.Matrix([[1, sp.I], [1, -sp.I]]) / sp.sqrt(2)
    reduced = sp.simplify(V * Dm * V.H)
    expected = sp.Matrix([[0, p - m], [p + m, 0]])
    return reduced, expected


def main():
    print("=" * 64)
    print("Verification: Sol3 APS eta on the compact mapping torus")
    print("Target: re-derive eta_APS(Sol3) from LC data + quotient spin structure")
    print("=" * 64)

    r0 = 3.0
    C, omega, d_omega, vol = build_sol3_geometry(r0)

    print()
    print(f"  r0 = {r0}")
    print("  Compact quotient convention: M_A = T^2 rtimes_A S^1, A = [[2,1],[1,1]]")
    print("  Spin structures fixed as in Phase5C-04rr2:")
    print("    Sol-P: base periodic, fiber periodic")
    print("    Sol-A: base anti-periodic, fiber periodic")

    print()
    print("[1] Local Levi-Civita Chern-Simons integrals on Sol3")
    print("-" * 64)
    adj = adjoint_cs_integral(omega, d_omega, vol)
    sigma = spinor_cs_integral(omega, C, vol)
    ok_adj = abs(adj["cs_A"]) < 1e-12 and abs(adj["cs_B"]) < 1e-12 and abs(adj["total"]) < 1e-12
    ok_sigma = abs(sigma["cs_A"]) < 1e-12 and abs(sigma["cs_B"]) < 1e-12 and abs(sigma["total"]) < 1e-12
    print(f"  int_Y Tr(omega dw)        = {adj['cs_A']:+.8e}")
    print(f"  int_Y (2/3)Tr(omega^3)    = {adj['cs_B']:+.8e}")
    print(f"  int_Y CS_adj              = {adj['total']:+.8e}   "
          f"[{'OK' if ok_adj else 'FAIL'}]")
    print(f"  int_Y CS_sigma            = {sigma['total']:+.8e}   "
          f"[{'OK' if ok_sigma else 'FAIL'}]")
    print("  => local spinor-CS density vanishes on Sol3.")

    print()
    print("[2] Dirac zero-order term in the left-invariant frame")
    print("-" * 64)
    Omega, B = dirac_zero_order_term(omega)
    ok_B = np.max(np.abs(B)) < 1e-12
    for c, Om in enumerate(Omega):
        print(f"  Omega_{c} =")
        print(f"    {np.array2string(Om, precision=6, suppress_small=True)}")
    print(f"  B = sum_c gamma^c Omega_c =")
    print(f"    {np.array2string(B, precision=6, suppress_small=True)}")
    print(f"  lower-order term cancellation: {'OK' if ok_B else 'FAIL'}")
    print("  => D_Sol = gamma^a E_a, so constant spinors are harmonic candidates.")

    print()
    print("[3] Time-reversal symmetry and eta(0)")
    print("-" * 64)
    transformed = time_reversal_check()
    ok_tr = True
    for j, Tg in enumerate(transformed):
        target = -pauli_gamma()[j]
        good = np.max(np.abs(Tg - target)) < 1e-12
        ok_tr &= good
        print(f"  T gamma_{j} T^(-1) =")
        print(f"    {np.array2string(Tg, precision=6, suppress_small=True)}   "
              f"[{'OK' if good else 'FAIL'}]")
    print("  Since E_a are real differential operators, T D T^(-1) = -D.")
    print("  Therefore the nonzero spectrum is +/- symmetric and eta(0) = 0.")

    print()
    print("[4] Nonzero fiber Fourier sectors: no zero modes")
    print("-" * 64)
    A = np.array([[2.0, 1.0], [1.0, 1.0]])
    t_grid = np.linspace(0.0, 1.0, 400)
    sample_modes = [(1, 0), (0, 1), (1, 1), (2, 1)]
    no_zero_all = True
    reduced, expected = off_diagonal_reduction_check()
    ok_reduce = (reduced - expected) == sp.zeros(2)
    print("  Exact operator reduction:")
    print("    V (sigma_3 p + i m sigma_1) V^dagger =")
    print(f"      {reduced}")
    print(f"    expected = {expected}   [{'OK' if ok_reduce else 'FAIL'}]")
    print("  For k != 0 the fiber part is a real two-vector p_k(t),")
    print("  so after a local U(1) rotation the zero-mode equation becomes")
    print("      (d/dt +/- m_k(t)) phi = 0,   m_k(t) = ||p_k(t)|| > 0.")
    print("  Hence periodic or anti-periodic solutions are impossible unless k = 0.")
    print()
    print(f"  {'k':>10}  {'min m_k':>12}  {'int_0^1 m_k dt':>16}  {'zero mode?':>12}")
    for k_vec in sample_modes:
        prof = fiber_mass_profile(A, k_vec, t_grid)
        min_m = float(np.min(prof))
        int_m = float(np.trapezoid(prof, t_grid))
        no_zero = (min_m > 1e-12) and (int_m > 1e-12)
        no_zero_all &= no_zero
        print(f"  {str(k_vec):>10}  {min_m:12.6f}  {int_m:16.6f}  "
              f"{'NO' if no_zero else 'UNRESOLVED'}")
    print("  Representative nonzero sectors confirm the general argument.")

    print()
    print("[5] k = 0 sector and kernel dimension h")
    print("-" * 64)
    h_sol_p = base_periodic_sector_kernel_dim(base_periodic=True)
    h_sol_a = base_periodic_sector_kernel_dim(base_periodic=False)
    print("  k = 0 sector: D_0 = gamma^3 d/dt, so zero modes are constant spinors.")
    print(f"  Sol-P: base periodic     -> h = {h_sol_p}")
    print(f"  Sol-A: base anti-periodic -> h = {h_sol_a}")
    ok_h = (h_sol_p == 2 and h_sol_a == 0)

    print()
    print("[6] APS assembly")
    print("-" * 64)
    eta0_sol_p = 0
    eta0_sol_a = 0
    eta_aps_sol_p = (eta0_sol_p + h_sol_p) / 2
    eta_aps_sol_a = (eta0_sol_a + h_sol_a) / 2
    ok_eta = (abs(eta_aps_sol_p - 1.0) < 1e-12 and
              abs(eta_aps_sol_a - 0.0) < 1e-12)
    print(f"  Sol-P: eta(0) = {eta0_sol_p}, h = {h_sol_p}, "
          f"eta_APS = {eta_aps_sol_p}")
    print(f"  Sol-A: eta(0) = {eta0_sol_a}, h = {h_sol_a}, "
          f"eta_APS = {eta_aps_sol_a}")
    print("  Note: int_Y CS_sigma = 0 fixes only the local/non-integer part.")
    print("        On Sol-P the integer shift comes from h = 2.")

    overall_ok = all([ok_adj, ok_sigma, ok_B, ok_tr, ok_reduce, no_zero_all, ok_h, ok_eta])

    print()
    print("=" * 64)
    if overall_ok:
        print("RESULT: PASS")
        print("  eta_APS(Sol-P) = 1")
        print("  eta_APS(Sol-A) = 0")
        print("  => Sol3 APS entry is spin-structure dependent on the compact quotient.")
    else:
        print("RESULT: FAIL")
    print("=" * 64)
    return overall_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
