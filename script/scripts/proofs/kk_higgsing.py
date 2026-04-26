"""
kk_higgsing.py
=====================

Symbolic verification for Appendix E.2 and E.5:
  KK reduction Maxwell baseline universality, biaxial Higgsing on Sol3,
  and the four-type Higgsing classification (none / uniaxial / triaxial /
  biaxial) realized exactly once on (T3, Nil3, S3, Sol3).

Target statements:
  E.2:  m^2(A_0) = 0,   m^2(A_1) = m^2(A_2) = K_Mxw = -L^2/(2 r0^4)   (Sol3)
  E.5:  Higgsing pattern = (none, uniaxial, triaxial, biaxial)
        on (T3, Nil3, S3, Sol3) respectively.
        K_Mxw = -L^2/(2 r0^4) is universal across all four geometries.

Method:
  - Use the dppu KK pipeline (Gamma-Gamma shortcut).
  - For each geometry, build the F-tilde correction dictionary, compute
    R_GG = gamma_gamma_ricci(omega1_fn), then call extract_all to obtain
    {Maxwell, mass dict, CS}.
  - Compare against the closed-form expectations.

Run:
  python -X utf8 script/scripts/proof/kk_higgsing.py

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

from sympy import Rational, Symbol, S, symbols, cancel

from dppu.kk.field_strength import (
    make_F_tilde, make_omega1,
    s3_corrections, nil3_corrections, t3_corrections, sol3_corrections,
)
from dppu.kk.gamma_gamma import gamma_gamma_ricci
from dppu.kk.extractor import extract_all


# Shared symbols
r0 = symbols("r0", positive=True)
L  = symbols("L",  positive=True)
A  = [Symbol(f"A{k}") for k in range(3)]
dA = [[Symbol(f"dA{j}{k}") for k in range(3)] for j in range(3)]


def kk_for(corrections):
    """Run the KK extraction pipeline for one geometry."""
    F_fn = lambda j, k: make_F_tilde(j, k, A, dA, corrections)
    omega1_fn = make_omega1(F_fn, r0, L)
    R_GG = gamma_gamma_ricci(omega1_fn)
    return extract_all(R_GG, A, dA)


def main():
    print("=" * 64)
    print("Verification: KK Higgsing four-type classification")
    print("Target: Appendix E.2 and E.5")
    print("=" * 64)

    K_Mxw_expected = -L**2 / (2 * r0**4)

    rows = [
        ("T3",   t3_corrections(),       "none",     []),
        ("Nil3", nil3_corrections(A),    "uniaxial", [2]),
        ("S3",   s3_corrections(A),      "triaxial", [0, 1, 2]),
        ("Sol3", sol3_corrections(A),    "biaxial",  [1, 2]),
    ]
    expected_akk = {
        "T3": S.Zero,
        "Nil3": Rational(2, 3),
        "S3": S.Zero,
        "Sol3": Rational(2, 3),
    }

    overall_ok = True
    print()
    print(f"  Universal Maxwell baseline expected: K_Mxw = {K_Mxw_expected}")
    print()
    print("-" * 64)
    print(f"  {'geo':5}  {'pattern':10}  {'massive':12}  {'CS':4}  {'Maxwell':12}")
    print("-" * 64)

    for name, corr, expected_pattern, expected_massive in rows:
        res = kk_for(corr)
        mass_dict = res["mass"] or {}
        massive_dirs = sorted(mass_dict.keys())
        n_massive = len(massive_dirs)
        if n_massive == 0:
            pattern = "none"
        elif n_massive == 1:
            pattern = "uniaxial"
        elif n_massive == 2:
            pattern = "biaxial"
        elif n_massive == 3:
            pattern = "triaxial"
        else:
            pattern = f"({n_massive})"

        # Maxwell value (any one entry; universality is checked below)
        mxw_dict = res["maxwell"] or {}
        if mxw_dict:
            mxw_val = list(mxw_dict.values())[0]
            mxw_match = (cancel(mxw_val - K_Mxw_expected) == S.Zero)
        else:
            mxw_val, mxw_match = None, False

        cs = res["cs"]
        cs_zero = (cs is None or cs == S.Zero)

        ok_pattern = (pattern == expected_pattern)
        ok_massive = (massive_dirs == expected_massive)
        row_ok = ok_pattern and ok_massive and mxw_match
        flag = "OK" if row_ok else "FAIL"

        mass_values = [mass_dict.get(k, S.Zero) for k in [0, 1, 2]]
        mean_mass = sum(mass_values) / 3
        a_kk = cancel(sum((mk - mean_mass) ** 2 for mk in mass_values))
        a_kk_ratio = cancel(a_kk / K_Mxw_expected**2)
        ok_akk = cancel(a_kk_ratio - expected_akk[name]) == S.Zero

        print(f"  {name:5}  {pattern:10}  {str(massive_dirs):12}  "
              f"{('!=0' if not cs_zero else '0'):4}  "
              f"{str(mxw_val):12}  [{flag}]")
        print(f"         A_KK/K^2 = {a_kk_ratio}   expected {expected_akk[name]}"
              f"   [{'OK' if ok_akk else 'FAIL'}]")

        # Geometry-specific Sol3 mass values (E.2)
        if name == "Sol3":
            print()
            print("  [Sol3 mass dict detail (E.2)]")
            for k in [0, 1, 2]:
                mk = mass_dict.get(k, S.Zero)
                if k == 0:
                    expected = S.Zero
                else:
                    expected = K_Mxw_expected
                ok = cancel(mk - expected) == S.Zero
                fl = "OK" if ok else "FAIL"
                print(f"    m^2(A_{k}) = {mk}   expected {expected}   [{fl}]")
                row_ok &= ok
            print()

        overall_ok &= row_ok and ok_akk

    print("-" * 64)
    print()
    print("=" * 64)
    print("RESULT: " + ("PASS  (Maxwell universal; 4-type Higgsing and A_KK exact)"
                        if overall_ok else "FAIL"))
    print("=" * 64)
    return overall_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
