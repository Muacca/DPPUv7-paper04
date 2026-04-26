"""
eta_defect_coefficients.py
=================================

Symbolic verification for Appendix E.6:
  reduced 1D AX-torsion defect benchmark eta = eta(z).

Checks:
  1. Exact homogeneous mass datum from current dppu:
         M_geo = d^2 V_eff / d eta^2 |_{eta=0}
  2. Reduced kinetic datum from the AX contortion gradient:
         K_geo = 6 * Vol_geo(rho) / (kappa^2 rho^2)
  3. Equality K_geo = M_geo for all four geometries.

Expected coefficients:
  T3, Nil3, Sol3 : 96 pi^4 L rho / kappa^2
  S3             : 12 pi^2 L rho / kappa^2

Run:
  python -X utf8 script/scripts/paper04/eta_defect_coefficients.py

Author: Muacca
Date: 2026-04-24
"""

from pathlib import Path
import sys

import sympy as sp

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


def expected_mass_expr(label, params, rho):
    if label == "S3":
        return 12 * sp.pi**2 * params["L"] * rho / params["kappa"]**2
    return 96 * sp.pi**4 * params["L"] * rho / params["kappa"]**2


def reduced_k_expr(params, volume, rho):
    return sp.simplify(6 * volume / (params["kappa"]**2 * rho**2))


def main():
    print("=" * 72)
    print("Verification: eta-mode reduced defect coefficients")
    print("Target: Appendix E.6")
    print("=" * 72)
    print()
    print("  Benchmark assumptions:")
    print("    - AX torsion mode")
    print("    - slice-wise homogeneous background")
    print("    - reduced 1D ansatz eta = eta(z)")
    print("    - kinetic datum from contortion gradient")
    print()

    all_ok = True

    for label, topology in TOPOLOGIES:
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
            params["V"]: sp.S.Zero,
            params["theta_NY"]: sp.S.Zero,
        })
        mass_expr = sp.simplify(sp.diff(veff, params["eta"], 2).subs(params["eta"], 0))
        expected_expr = expected_mass_expr(label, params, rho)
        kinetic_expr = reduced_k_expr(params, eng.data["total_volume"], rho)

        ok_expected = sp.simplify(mass_expr - expected_expr) == 0
        ok_match = sp.simplify(mass_expr - kinetic_expr) == 0
        all_ok = all_ok and ok_expected and ok_match

        print("-" * 72)
        print(f"  Topology: {label}")
        print(f"    volume      = {eng.data['total_volume']}")
        print(f"    M_geo       = {sp.factor(mass_expr)}")
        print(f"    K_geo(red.) = {sp.factor(kinetic_expr)}")
        print(f"    expected    = {sp.factor(expected_expr)}")
        print(f"    M_geo == expected ? {'YES' if ok_expected else 'NO'}")
        print(f"    K_geo == M_geo  ? {'YES' if ok_match else 'NO'}")

    print("-" * 72)
    print()
    print("=" * 72)
    if all_ok:
        print("RESULT: PASS  (4 geometries: exact eta-mass datum matches reduced K_geo)")
    else:
        print("RESULT: FAIL")
    print("=" * 72)
    return all_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
