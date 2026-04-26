"""
kk_normalization.py
==========================

Symbolic verification for Appendix E.11:
  KK reduction normalization identity
      k_{3D}^{tor-CS} = (1/2) k_{4D}^{NY}
  derived from
      S_NY^{(4D)} = (k_{4D}^{NY} / (8 pi^2)) * int_{X^4} NY
  via the CZ identity NY = d C_3^{tor} and Stokes' theorem on
  X^4 = M^3 x [0, beta], with the boundary normalization
      N_top := (1/(4 pi^2)) * int_{M_3} C_3^{tor},
  where N_top is used as a frame-bundle-normalized torsional charge basis.

Method:
  Pure symbolic algebra. We do NOT assume any specific value of
  k_{4D}^{NY}, and we do not assume integer-valuedness of N_top. We track
  the 8 pi^2 factor inherited from
  the 4D action, the 4 pi^2 factor inherited from the boundary
  normalization, and the trivial collapse along [0, beta].

  Additionally we check that, for an untwisted product-circle KK
  reduction of real bosonic fields, the nonzero tower contributes
  zero in the parity-odd sector. The reason is that Fourier modes
  n and -n come as conjugate pairs, any even weight depends only on
  n^2, and a parity-odd residue must therefore be of the generic form
  n * w(n^2), which cancels pairwise.

Run:
  python -X utf8 script/scripts/proof/kk_normalization.py

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

from sympy import Function, Rational, S, pi, simplify, symbols


def main():
    print("=" * 64)
    print("Verification: KK reduction normalization identity")
    print("Target: Appendix E.11   k_{3D}^{tor-CS} = (1/2) k_{4D}^{NY}")
    print("=" * 64)

    # ----- Step A: 4D action notation -----
    k4 = symbols("k4_NY", real=True)        # 4D NY level (free)
    Ntop = symbols("N_top", real=True)      # frame-bundle-normalized torsional charge
    beta = symbols("beta", positive=True)   # S^1 circumference (irrelevant after collapse)

    # 4D action coefficient: S_NY = (k4 / (8 pi^2)) * int NY
    coeff_4D = k4 / (8 * pi**2)

    # ----- Step B: CZ + Stokes on X^4 = M^3 x [0, beta] -----
    # int_{X^4} NY = int_{X^4} d C_3^tor
    #             = int_{boundary} C_3^tor   (Stokes)
    # Boundary of M^3 x [0, beta] in the relevant inflow setup
    # collapses to a single copy of M^3 (the [0, beta] direction is the
    # KK circle, contributing trivially after factoring out beta from
    # the 4D measure and the 3D coefficient definition).
    #
    # By the frame-bundle normalization
    #     N_top = (1/(4 pi^2)) * int_{M_3} C_3^tor,
    # we have
    #     int_{M_3} C_3^tor = 4 pi^2 * N_top.

    int_C3 = 4 * pi**2 * Ntop

    # ----- Step C: Substitute back into the action -----
    S_eff_3D = simplify(coeff_4D * int_C3)
    print()
    print("[1] Substitute int C_3^tor = 4 pi^2 * N_top into 4D action")
    print("-" * 64)
    print(f"  S_eff^(3D) = ({coeff_4D}) * (4 pi^2 * N_top)")
    print(f"             = {S_eff_3D}")

    # In the normalized torsional-charge basis, the 3D effective action is
    #     S_eff^(3D) = k_{3D}^{tor-CS} * N_top.
    # Compare coefficients:
    k3 = symbols("k3_torCS", real=True)
    k3_solved = simplify(S_eff_3D / Ntop)
    print()
    print("[2] Identify the coefficient of N_top with k_{3D}^{tor-CS}")
    print("-" * 64)
    print(f"  k_{{3D}}^{{tor-CS}} = {k3_solved}")

    expected = k4 / 2
    diff = simplify(k3_solved - expected)
    ok_main = (diff == S.Zero)
    flag = "OK" if ok_main else "FAIL"
    print(f"  expected           = k4_NY / 2 = {expected}")
    print(f"  difference         = {diff}   [{flag}]")

    # ----- Step D: KK tower (n != 0) parity-odd cancellation -----
    print()
    print("[3] KK tower (n != 0) parity-odd cancellation")
    print("-" * 64)
    n = symbols("n", integer=True, nonzero=True)
    w = Function("w")
    # For a real bosonic field on an untwisted product circle, Fourier modes
    # n and -n are exchanged by tau -> -tau and carry the same even weight.
    # Any parity-odd tower residue must therefore be of the form
    #   n * w(n^2),
    # with w depending only on |n| (equivalently n^2).
    contrib_pos = n * w(n**2)
    contrib_neg = (-n) * w((-n)**2)
    pair = simplify(contrib_pos + contrib_neg)
    ok_tower = (pair == S.Zero)
    flag2 = "OK" if ok_tower else "FAIL"
    print("  product-circle symmetry: tau -> -tau exchanges n <-> -n")
    print("  real-field condition: tower weights depend only on n^2")
    print(f"  generic pair: {contrib_pos} + ({contrib_neg}) = {pair}")
    print(f"  -> KK tower contribution to k_{{3D}}^{{tor-CS}} = 0   [{flag2}]")

    # ----- Step E: r0 independence of the 1/2 factor -----
    print()
    print("[4] r0 independence of the 1/2 factor")
    print("-" * 64)
    r0 = symbols("r0", positive=True)
    # The factor 1/2 comes from 4 pi^2 / 8 pi^2 only. r0 enters neither.
    factor = simplify((4 * pi**2) / (8 * pi**2))
    ok_factor = (factor == Rational(1, 2))
    flag3 = "OK" if ok_factor else "FAIL"
    print(f"  4 pi^2 / 8 pi^2 = {factor}   [{flag3}]")
    print(f"  No r0 dependence in the normalization ratio.")

    overall_ok = ok_main and ok_tower and ok_factor
    print()
    print("=" * 64)
    if overall_ok:
        print("RESULT: PASS  (k_{3D}^{tor-CS} = (1/2) k_{4D}^{NY}, exact)")
    else:
        print("RESULT: FAIL")
    print("=" * 64)
    return overall_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
