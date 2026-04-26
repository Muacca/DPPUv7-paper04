"""
torsional_charge.py
============================

Symbolic verification for Appendix E.9:
  Frame-bundle-normalized torsional-charge benchmark on S^3.

Target:
  N_top(S^3) := (1/(4 pi^2)) * int_{S^3} C_3^{tor} = 6 r0^2,
  where N_top denotes the frame-bundle-normalized torsional charge.

Method:
  [A] Derive the geometric factor int_{S^3} C_3^{tor} from the S^3 structure
      constants produced by the current dppu topology engine:
          C^i_{jk} = (4/r0) epsilon_{ijk}
          T^i = (1/2) C^i_{jk} e^j ^ e^k
          C_3^{tor} = e_i ^ T^i
      This yields the density coefficient 12/r0 without hardcoding it.

  [B] Apply the frame-bundle normalization
          (1 / (4 pi^2)) int_{S^3} C_3^{tor}
      associated with pi_3(SO(3)) = Z and vector-trace normalization T_R = 1.
      No integer-valuedness of N_top is assumed for continuous r0; the script
      verifies the normalized geometric value only.

Run:
  python -X utf8 script/scripts/paper04/torsional_charge.py

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

from sympy import Rational, S, pi, simplify

from dppu.topology.unified import DOFConfig, FiberMode, TopologyType, UnifiedEngine
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.utils.levi_civita import epsilon_symbol


def build_s3_structure_constants():
    cfg = DOFConfig(
        topology=TopologyType.S3,
        fiber_mode=FiberMode.NONE,
        torsion_mode=Mode.AX,
        ny_variant=NyVariant.FULL,
    )
    eng = UnifiedEngine(cfg)
    eng.run()
    return eng.data["structure_constants"], eng.data["params"]["r"]


def derive_cs_density_from_structure_constants(C):
    terms = []
    density = S.Zero
    for i in range(3):
        for j in range(3):
            for k in range(3):
                eps = epsilon_symbol(i, j, k)
                if eps == 0:
                    continue
                term = simplify(Rational(1, 2) * C[i, j, k] * eps)
                if term != S.Zero:
                    terms.append(((i, j, k), term))
                    density += term
    return terms, simplify(density)


def main():
    print("=" * 64)
    print("Verification: frame-bundle-normalized torsional charge N_top(S^3) = 6 r0^2")
    print("Target: Appendix E.9")
    print("=" * 64)

    C, r0 = build_s3_structure_constants()

    print()
    print("[1] S^3 structure constants from dppu topology engine")
    print("-" * 64)
    expected_checks = [
        ((0, 1, 2),  4 / r0),
        ((1, 2, 0),  4 / r0),
        ((2, 0, 1),  4 / r0),
        ((0, 2, 1), -4 / r0),
    ]
    ok_struct = True
    for (i, j, k), expected in expected_checks:
        actual = simplify(C[i, j, k])
        ok = simplify(actual - expected) == S.Zero
        print(f"  C[{i},{j},{k}] = {actual}   expected {expected}   "
              f"[{'OK' if ok else 'FAIL'}]")
        ok_struct &= ok

    print()
    print("[2] Derive C_3^{tor} = e_i ^ T^i from T^i = (1/2) C^i_{jk} e^j ^ e^k")
    print("-" * 64)
    terms, cs_density = derive_cs_density_from_structure_constants(C)
    for (i, j, k), term in terms:
        print(f"  (i,j,k)=({i},{j},{k}):  (1/2) C[{i},{j},{k}] epsilon_{{{i}{j}{k}}}"
              f" = {term}")
    expected_density = 12 / r0
    ok_density = simplify(cs_density - expected_density) == S.Zero
    print()
    print(f"  C_3^tor density coefficient = {cs_density}")
    print(f"  expected                    = {expected_density}   "
          f"[{'OK' if ok_density else 'FAIL'}]")

    print()
    print("[3] Integrate over S^3")
    print("-" * 64)
    vol_s3 = 2 * pi**2 * r0**3
    integral = simplify(cs_density * vol_s3)
    expected_integral = simplify(24 * pi**2 * r0**2)
    ok_integral = simplify(integral - expected_integral) == S.Zero
    print(f"  Vol(S^3)              = {vol_s3}")
    print(f"  int_{{S^3}} C_3^tor     = {integral}")
    print(f"  expected              = {expected_integral}   "
          f"[{'OK' if ok_integral else 'FAIL'}]")

    print()
    print("[4] Frame-bundle normalization")
    print("-" * 64)
    print("  Convention: N_top := (1/(4 pi^2)) int_{S^3} C_3^{tor}")
    print("  Origin: pi_3(SO(3)) = Z with vector-representation normalization T_R = 1")
    print("  No integrality claim is made for continuous r0 without extra quantization")
    norm = 1 / (4 * pi**2)
    n_top = simplify(norm * integral)
    expected_n_top = 6 * r0**2
    ok_n_top = simplify(n_top - expected_n_top) == S.Zero
    print(f"  N_top(S^3)            = {n_top}")
    print(f"  expected              = {expected_n_top}   "
          f"[{'OK' if ok_n_top else 'FAIL'}]")

    print()
    print("[5] Numeric spot checks")
    print("-" * 64)
    ok_numeric = True
    for r0_val in [3, 5]:
        n_val = simplify(n_top.subs(r0, r0_val))
        expected_val = 6 * r0_val**2
        ok = simplify(n_val - expected_val) == S.Zero
        print(f"  r0 = {r0_val}: N_top = {n_val}   expected {expected_val}   "
              f"[{'OK' if ok else 'FAIL'}]")
        ok_numeric &= ok

    overall_ok = ok_struct and ok_density and ok_integral and ok_n_top and ok_numeric
    print()
    print("=" * 64)
    if overall_ok:
        print("RESULT: PASS  (S^3 geometry derives frame-bundle-normalized torsional charge N_top = 6 r0^2)")
    else:
        print("RESULT: FAIL")
    print("=" * 64)
    return overall_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
