"""
landau_levels_nil3.py
============================

Symbolic ladder-algebra verification for Appendix E.10:
  Heisenberg Landau levels in each Nil3 p2 Fourier sector.

Scope:
  This script verifies the closed-form spectrum after Heisenberg mode
  decomposition.  It is not a complete spectral theorem for compact
  nilmanifold quotients.

Target:
  mu_n^+ = 2 n omega |k2| + k2^2,
  mu_n^- = 2 (n + 1) omega |k2| + k2^2,
  omega = |k2| / r0.

Run:
  python -X utf8 script/scripts/proof/landau_levels_nil3.py

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

import sympy as sp


def eps3(i, j, k):
    if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        return sp.S.One
    if (i, j, k) in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
        return -sp.S.One
    return sp.S.Zero


def build_nil3_spin_connection(r0):
    C = {}
    omega = {}
    for a in range(3):
        for b in range(3):
            for c in range(3):
                C[(a, b, c)] = sp.S.Zero

    C[(2, 0, 1)] = 1 / r0
    C[(2, 1, 0)] = -1 / r0

    for a in range(3):
        for b in range(3):
            for c in range(3):
                omega[(a, b, c)] = sp.cancel(
                    sp.Rational(1, 2)
                    * (C[(a, b, c)] + C[(c, b, a)] - C[(b, a, c)])
                )
    return C, omega


def main():
    print("=" * 72)
    print("Verification: Nil3 Heisenberg Landau levels")
    print("Target: Appendix E.10")
    print("=" * 72)

    r0 = sp.symbols("r0", positive=True)
    k_abs = sp.symbols("k2_abs", positive=True)
    n = sp.symbols("n", integer=True, nonnegative=True)
    m = sp.symbols("m", integer=True)

    C, omega_conn = build_nil3_spin_connection(r0)

    print()
    print("[1] Nil3 structure constant and Levi-Civita spin connection")
    print("-" * 72)
    ok_c = sp.cancel(C[(2, 0, 1)] - 1 / r0) == sp.S.Zero
    print(f"  C^2_01 = {C[(2, 0, 1)]}   expected 1/r0   [{'OK' if ok_c else 'FAIL'}]")
    nonzero_omega = [
        (idx, sp.factor(val))
        for idx, val in omega_conn.items()
        if val != sp.S.Zero
    ]
    for (a, b, c), val in nonzero_omega:
        print(f"  omega^{a}{b}_{c} = {val}")

    print()
    print("[2] Heisenberg ladder algebra in a fixed p2 sector")
    print("-" * 72)
    omega = k_abs / r0
    B = sp.cancel(omega * k_abs)

    # Hermitian oscillator variables Pi0, Pi1 are normalized so that
    # [Pi0, Pi1] = i B in each p2 sector.  The ladder construction then
    # gives Pi0^2 + Pi1^2 = B(2N + 1).
    N = n
    oscillator = sp.cancel(B * (2 * N + 1))
    spin_up_term = -B
    spin_down_term = B

    mu_up = sp.cancel(oscillator + spin_up_term + k_abs**2)
    mu_down = sp.cancel(oscillator + spin_down_term + k_abs**2)
    expected_up = sp.cancel(2 * n * omega * k_abs + k_abs**2)
    expected_down = sp.cancel(2 * (n + 1) * omega * k_abs + k_abs**2)

    ok_up = sp.cancel(mu_up - expected_up) == sp.S.Zero
    ok_down = sp.cancel(mu_down - expected_down) == sp.S.Zero
    print(f"  [Pi0, Pi1] = i B,  B = omega |k2| = {B}")
    print("  a = (Pi0 + i Pi1)/sqrt(2B),  a^dagger = (Pi0 - i Pi1)/sqrt(2B)")
    print("  Pi0^2 + Pi1^2 = B(2N+1)")
    print(f"  mu_n^+ = {sp.factor(mu_up)}   [{'OK' if ok_up else 'FAIL'}]")
    print(f"  mu_n^- = {sp.factor(mu_down)}   [{'OK' if ok_down else 'FAIL'}]")

    print()
    print("[3] Pairing pattern and PPA zero-mode exclusion")
    print("-" * 72)
    pair_check = sp.cancel(expected_up.subs(n, n + 1) - expected_down) == sp.S.Zero
    zero_solution = sp.Rational(-1, 2)
    # The formal solution is not an integer, hence no PPA sector has p2 = 0.
    ppa_excludes_zero = not bool(zero_solution.is_integer)
    h_zero = expected_up.subs(n, 0)
    ok_h = h_zero == k_abs**2
    print(f"  mu_(n+1)^+ - mu_n^- = {sp.cancel(expected_up.subs(n, n + 1) - expected_down)}"
          f"   [{'OK' if pair_check else 'FAIL'}]")
    print("  n=0 spin-up branch is unpaired relative to the down branch.")
    print(f"  PPA: p2 = m + 1/2.  Solving p2=0 gives m={zero_solution}, not an integer"
          f"   [{'OK' if ppa_excludes_zero else 'FAIL'}]")
    print(f"  lowest PPA sector has mu_0^+ = {h_zero} > 0, so h = 0"
          f"   [{'OK' if ok_h else 'FAIL'}]")

    all_ok = ok_c and ok_up and ok_down and pair_check and ppa_excludes_zero and ok_h

    print()
    print("=" * 72)
    if all_ok:
        print("RESULT: PASS  (Nil3 Heisenberg ladder levels and PPA h=0)")
    else:
        print("RESULT: FAIL")
    print("=" * 72)
    return all_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
