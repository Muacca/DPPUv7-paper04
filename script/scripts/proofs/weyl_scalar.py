"""
weyl_scalar.py
=====================

Symbolic verification for Appendix E.7:
  Levi-Civita Weyl scalar C^2_LC for the four geometries.

Target statements:
  C^2_LC(T3)   = 0
  C^2_LC(Nil3) = 4 / (3 R^4)
  C^2_LC(S3)   = 0
  C^2_LC(Sol3) = 16 / (3 R^4)

Method:
  - Build each geometry via dppu.topology.unified.UnifiedEngine.
  - Read the precomputed Weyl scalar (data['weyl_scalar']).
  - Compare against the closed-form expectations symbolically.

Run:
  python -X utf8 script/scripts/proof/weyl_scalar.py

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

from sympy import S, Rational, cancel, simplify

from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant


def run_geom(topo_enum):
    cfg = DOFConfig(
        topology=topo_enum,
        torsion_mode=Mode.AX,
        enable_squash=False,
        enable_shear=False,
        ny_variant=NyVariant.FULL,
    )
    eng = UnifiedEngine(cfg)
    eng.run()
    return eng


def main():
    print("=" * 64)
    print("Verification: Levi-Civita Weyl scalar C^2_LC for 4 geometries")
    print("Target: Appendix E.7")
    print("=" * 64)

    cases = [
        (TopologyType.T3,   "T3",   lambda R: S.Zero),
        (TopologyType.NIL3, "Nil3", lambda R: Rational(4, 3) / R**4),
        (TopologyType.S3,   "S3",   lambda R: S.Zero),
        (TopologyType.SOL3, "Sol3", lambda R: Rational(16, 3) / R**4),
    ]

    overall_ok = True
    print()
    print(f"  {'geo':5}  {'computed':30}  {'expected':18}  flag")
    print("  " + "-" * 64)
    for topo_enum, name, expected_fn in cases:
        eng = run_geom(topo_enum)
        R = eng.data["params"]["r"]
        c2_lc = cancel(eng.data.get("weyl_scalar", S.Zero))
        expected = expected_fn(R)
        diff = cancel(c2_lc - expected)
        ok = (diff == S.Zero)
        flag = "OK" if ok else "FAIL"
        print(f"  {name:5}  {str(c2_lc):30}  {str(expected):18}  [{flag}]")
        overall_ok &= ok

    print()
    print("=" * 64)
    if overall_ok:
        print("RESULT: PASS  (4-geometry C^2_LC values exact)")
    else:
        print("RESULT: FAIL")
    print("=" * 64)
    return overall_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
