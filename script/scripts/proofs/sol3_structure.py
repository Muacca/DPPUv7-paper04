"""
sol3_structure.py
========================

Symbolic verification for Appendix E.1:
  Sol3 structure constants and frame rigidity.

Target statements (Appendix E.1):
  C^1_{01} = +1/R,  C^2_{02} = -1/R   (DPPU convention [Eb,Ec] = -C^a_{bc} Ea)
  d/d(eps) C = 0,  d/ds C = 0          (frame rigidity)
  d2V/deps2 = d2V/ds2 = d2V/depsds = 0 at eps=s=0

Method:
  - Build the Sol3 engine via dppu.topology.unified.UnifiedEngine.
  - Read the spatial structure constants and verify the closed-form values.
  - Differentiate the structure constants w.r.t. squash/shear deformation
    parameters and check exact zero (SymPy cancel).

Run:
  python -X utf8 script/scripts/proof/sol3_structure.py

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

from sympy import S, Symbol, cancel, diff

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant


def main():
    print("=" * 64)
    print("Verification: Sol3 structure constants and frame rigidity")
    print("Target: Appendix E.1")
    print("=" * 64)

    # Build Sol3 engine with deformation parameters enabled (squash + shear)
    cfg = DOFConfig(
        topology=TopologyType.SOL3,
        torsion_mode=Mode.AX,
        enable_squash=True,
        enable_shear=True,
        ny_variant=NyVariant.FULL,
    )
    eng = UnifiedEngine(cfg)
    eng.run()

    C_sym = eng.data["structure_constants"]
    R = eng.data["params"]["r"]
    eps = eng.data["params"].get("epsilon", Symbol("epsilon"))
    s   = eng.data["params"].get("s", Symbol("s"))

    print()
    print("[1] Closed-form values of nonzero spatial structure constants")
    print("-" * 64)
    expected = {
        (1, 0, 1): +1 / R,
        (1, 1, 0): -1 / R,
        (2, 0, 2): -1 / R,
        (2, 2, 0): +1 / R,
    }
    all_ok = True
    for (a, b, c), val in expected.items():
        actual = cancel(C_sym[a, b, c])
        diff_val = cancel(actual - val)
        ok = diff_val == S.Zero
        flag = "OK" if ok else "FAIL"
        print(f"  C[{a},{b},{c}] = {actual}   expected {val}   [{flag}]")
        all_ok &= ok

    print()
    print("[2] Frame rigidity: derivatives w.r.t. (epsilon, s)")
    print("-" * 64)
    rigid_ok = True
    for (a, b, c) in [(1, 0, 1), (2, 0, 2)]:
        d_eps = cancel(diff(C_sym[a, b, c], eps))
        d_s   = cancel(diff(C_sym[a, b, c], s))
        ok_eps = d_eps == S.Zero
        ok_s   = d_s == S.Zero
        flag_e = "OK" if ok_eps else "FAIL"
        flag_s = "OK" if ok_s else "FAIL"
        print(f"  d/d(eps) C[{a},{b},{c}] = {d_eps}  [{flag_e}]")
        print(f"  d/d(s)   C[{a},{b},{c}] = {d_s}  [{flag_s}]")
        rigid_ok &= ok_eps and ok_s

    print()
    print("[3] Effective-potential Hessian in (epsilon, s)")
    print("-" * 64)
    veff = eng.data["potential"]
    subs0 = {eps: S.Zero, s: S.Zero}
    d2_eps = cancel(diff(veff, eps, 2).subs(subs0))
    d2_s = cancel(diff(veff, s, 2).subs(subs0))
    d2_cross = cancel(diff(veff, eps, s).subs(subs0))
    ok_d2_eps = d2_eps == S.Zero
    ok_d2_s = d2_s == S.Zero
    ok_d2_cross = d2_cross == S.Zero
    potential_ok = ok_d2_eps and ok_d2_s and ok_d2_cross
    print(f"  d2V/d(eps)^2 |0 = {d2_eps}  [{'OK' if ok_d2_eps else 'FAIL'}]")
    print(f"  d2V/ds^2     |0 = {d2_s}  [{'OK' if ok_d2_s else 'FAIL'}]")
    print(f"  d2V/deps ds  |0 = {d2_cross}  [{'OK' if ok_d2_cross else 'FAIL'}]")

    print()
    print("=" * 64)
    if all_ok and rigid_ok and potential_ok:
        print("RESULT: PASS  (Sol3 structure constants exact; frame and Hessian rigid)")
    else:
        print("RESULT: FAIL")
    print("=" * 64)
    return all_ok and rigid_ok and potential_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
