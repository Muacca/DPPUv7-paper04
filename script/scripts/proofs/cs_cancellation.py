"""
cs_cancellation.py
=========================

Symbolic verification for Appendix E.3 and E.4:
  Sol3 Chern-Simons activation vanishes algebraically in the profile-local
  KK ansatz,
  and the off-diagonal CS rule (E.4) classifies the 4 geometries
  exactly into (0, 1, 3, 0) CS direction count for (T3, Nil3, S3, Sol3).

Target statements:
  E.3:  CS = 0 for Sol3 for arbitrary scalar self-coupling coefficients
        c01(z), c02(z) multiplying A1 in F01 and A2 in F02.  This supports
        profile-local inhomogeneous openings whose only effect is to change
        those scalar coefficients; it is not a full derivative-dependent
        inhomogeneous field-equation proof.
  E.4:  CS direction count = (0, 1, 3, 0) for (T3, Nil3, S3, Sol3),
        determined by whether the F-tilde correction leg lies inside
        or outside the field-strength index pair {i, j}.

Method:
  Two complementary checks are performed:

  (a) Structural classification: count off-diagonal correction directions
      directly from the F-tilde correction dictionary {(i,j): A_k}, using
      the criterion "leg k not in {i, j}". This is the literal statement
      of the off-diagonal rule.

  (b) General Sol3 self-coupling check: run the dppu KK extractor with
      independent scalar coefficients c01 and c02 in
          Ftilde_01 = F_01 + c01 A_1,
          Ftilde_02 = F_02 + c02 A_2.
      The algebraic CS coefficient must vanish without using c01 = -c02.

  (c) Pipeline confirmation: run the current Sol3 correction dictionary with
      constant scale and with Rinv promoted to a free profile-local symbol.

Run:
  python -X utf8 script/scripts/proof/cs_cancellation.py

Author: Muacca
Date: 2026-04-24
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

_DATA_DIR = Path(__file__).resolve().parents[3] / "data"
_DATA_DIR.mkdir(exist_ok=True)
from dppu.utils.tee_logger import setup_log
setup_log(__file__, log_dir=str(_DATA_DIR))

from sympy import Symbol, S, symbols, cancel

from dppu.kk.field_strength import (
    make_F_tilde, make_omega1,
    s3_corrections, nil3_corrections, t3_corrections,
    sol3_corrections, sol3_inhomogeneous_corrections,
)
from dppu.kk.gamma_gamma import gamma_gamma_ricci
from dppu.kk.extractor import extract_all


r0 = symbols("r0", positive=True)
L  = symbols("L",  positive=True)
A  = [Symbol(f"A{k}") for k in range(3)]
dA = [[Symbol(f"dA{j}{k}") for k in range(3)] for j in range(3)]


def kk_for(corrections):
    F_fn = lambda j, k: make_F_tilde(j, k, A, dA, corrections)
    omega1_fn = make_omega1(F_fn, r0, L)
    R_GG = gamma_gamma_ricci(omega1_fn)
    return extract_all(R_GG, A, dA)


def index_of_A(expr):
    """Recover the leg index k of an A_k inside a coefficient expression."""
    legs = []
    for k, sym in enumerate(A):
        if expr.has(sym):
            legs.append(k)
    return legs


def off_diagonal_count(corrections):
    """
    Count off-diagonal directions in a correction dict.
    A correction (i,j) -> coeff*A_k is OFF-DIAGONAL iff k not in {i,j}.

    Pairs (i,j) and (i',j') that are self-referential (k in {i,j})
    with opposite signs cancel and contribute 0.
    """
    off = 0
    self_ref = []
    for (i, j), coeff in corrections.items():
        legs = index_of_A(coeff)
        for k in legs:
            if k in (i, j):
                self_ref.append((i, j, k, coeff))
            else:
                off += 1
    # Self-referential pairs do not generate a 3D CS direction by E.4.
    return off, self_ref


def main():
    print("=" * 64)
    print("Verification: Chern-Simons cancellation and direction count")
    print("Target: Appendix E.3 (Sol3) and E.4 (off-diagonal rule)")
    print("=" * 64)

    overall_ok = True

    # ---- Part (a): structural off-diagonal count ----
    print()
    print("[1] Structural off-diagonal CS count (Appendix E.4)")
    print("-" * 64)
    print(f"  {'geo':5}  {'corrections':45}  {'off':3}  {'expected':8}  flag")
    print("  " + "-" * 70)
    table = [
        ("T3",   t3_corrections(),               0),
        ("Nil3", nil3_corrections(A),            1),
        ("S3",   s3_corrections(A),              3),
        ("Sol3", sol3_corrections(A),            0),
    ]
    for name, corr, expected in table:
        off, self_ref = off_diagonal_count(corr)
        ok = (off == expected)
        flag = "OK" if ok else "FAIL"
        corr_str = str(corr)[:43]
        print(f"  {name:5}  {corr_str:45}  {off:3}  {expected:8}  [{flag}]")
        overall_ok &= ok

    # ---- Part (b): arbitrary scalar self-couplings for Sol3 (E.3) ----
    print()
    print("[2] General Sol3 self-coupling check")
    print("-" * 64)
    c01, c02 = symbols("c01 c02", real=True)
    res_general = kk_for({
        (0, 1): c01 * A[1],
        (0, 2): c02 * A[2],
    })
    cs_general = res_general["cs"]
    ok_general = (cs_general is None or cs_general == S.Zero
                  or cancel(cs_general) == S.Zero)
    print("  Ftilde_01 = F_01 + c01 A_1")
    print("  Ftilde_02 = F_02 + c02 A_2")
    print(f"  CS = {cs_general}   expected 0 for arbitrary c01,c02   "
          f"[{('OK' if ok_general else 'FAIL')}]")
    overall_ok &= ok_general

    # ---- Part (c): pipeline confirmation for current Sol3 dictionaries ----
    print()
    print("[3] Pipeline confirmation: Sol3 CS = 0 algebraically (constant scale)")
    print("-" * 64)
    res_const = kk_for(sol3_corrections(A))
    cs_const = res_const["cs"]
    ok_const = (cs_const is None or cs_const == S.Zero
                or cancel(cs_const) == S.Zero)
    print(f"  CS = {cs_const}   expected 0   [{('OK' if ok_const else 'FAIL')}]")
    overall_ok &= ok_const

    print()
    print("[4] Pipeline confirmation: Sol3 CS = 0 for profile-local 1/R(z)")
    print("-" * 64)
    R_inv = symbols("Rinv", positive=True)
    res_inh = kk_for(sol3_inhomogeneous_corrections(A, R_inv))
    cs_inh = res_inh["cs"]
    ok_inh = (cs_inh is None or cs_inh == S.Zero
              or cancel(cs_inh) == S.Zero)
    print(f"  CS = {cs_inh}   expected 0 within the profile-local ansatz   "
          f"[{('OK' if ok_inh else 'FAIL')}]")
    overall_ok &= ok_inh

    print()
    print("=" * 64)
    if overall_ok:
        print("RESULT: PASS  (off-diagonal count (0,1,3,0); Sol3 profile-local CS=0)")
    else:
        print("RESULT: FAIL")
    print("=" * 64)
    return overall_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
