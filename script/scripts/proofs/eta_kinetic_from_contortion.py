"""
eta_kinetic_from_contortion.py
======================================

Independent symbolic derivation for Appendix E.6:
  reduced 1D AX-torsion eta-mode coefficient from the contortion gradient.

Checks:
  1. AX torsion T_ijk = (2 eta / rho) epsilon_ijk gives
         K_abc = (eta / rho) epsilon_abc
     under the DPPU contortion convention.
  2. The mass-side contraction K_abc K^abc gives the eta^2 coefficient
         M_geo = 6 Vol_geo(rho) / (kappa^2 rho^2).
  3. The kinetic-side contraction (d_z K_abc)(d_z K^abc) gives the
     (d_z eta)^2 coefficient
         K_geo = 6 Vol_geo(rho) / (kappa^2 rho^2).
  4. For all four geometries this equals c_geo rho and matches the
     exact homogeneous mass datum used by eta_defect_coefficients.py.

Run:
  python -X utf8 script/scripts/proof/eta_kinetic_from_contortion.py

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


def eps3(i, j, k):
    if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        return sp.S.One
    if (i, j, k) in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
        return -sp.S.One
    return sp.S.Zero


def build_ax_torsion_and_contortion(eta, rho, eta_z):
    """Return symbolic AX torsion, contortion, and z-gradient tensors."""
    T = {}
    K = {}
    dK = {}
    for a in range(3):
        for b in range(3):
            for c in range(3):
                T[(a, b, c)] = 2 * eta / rho * eps3(a, b, c)

    for a in range(3):
        for b in range(3):
            for c in range(3):
                K[(a, b, c)] = sp.cancel(
                    sp.Rational(1, 2)
                    * (T[(a, b, c)] + T[(b, c, a)] - T[(c, a, b)])
                )
                # Reduced 1D ansatz: the kinetic datum keeps the eta(z)
                # gradient and treats rho as slice-wise constant.
                dK[(a, b, c)] = sp.cancel(sp.diff(K[(a, b, c)], eta) * eta_z)

    return T, K, dK


def contraction_squared(tensor):
    return sp.cancel(sum(tensor[idx] ** 2 for idx in tensor))


def expected_mass_expr(label, params, rho):
    if label == "S3":
        return 12 * sp.pi**2 * params["L"] * rho / params["kappa"]**2
    return 96 * sp.pi**4 * params["L"] * rho / params["kappa"]**2


def engine_mass_expr(topology, params):
    cfg = DOFConfig(
        topology=topology,
        torsion_mode=Mode.AX,
        enable_squash=False,
        enable_shear=False,
        ny_variant=NyVariant.FULL,
    )
    eng = UnifiedEngine(cfg)
    eng.run()

    eta = eng.data["params"]["eta"]
    veff = eng.data["potential"].subs({
        eng.data["params"]["V"]: sp.S.Zero,
        eng.data["params"]["theta_NY"]: sp.S.Zero,
    })
    return sp.simplify(sp.diff(veff, eta, 2).subs(eta, 0)), eng


def main():
    print("=" * 72)
    print("Verification: K_geo from AX contortion gradient")
    print("Target: Appendix E.6")
    print("=" * 72)

    eta, eta_z, rho = sp.symbols("eta eta_z rho", positive=True)
    kappa = sp.symbols("kappa", positive=True)

    T, K, dK = build_ax_torsion_and_contortion(eta, rho, eta_z)

    torsion_contract = contraction_squared(T)
    contortion_contract = contraction_squared(K)
    gradient_contract = contraction_squared(dK)

    expected_torsion = 24 * eta**2 / rho**2
    expected_contortion = 6 * eta**2 / rho**2
    expected_gradient = 6 * eta_z**2 / rho**2

    ok_torsion = sp.cancel(torsion_contract - expected_torsion) == sp.S.Zero
    ok_contortion = sp.cancel(contortion_contract - expected_contortion) == sp.S.Zero
    ok_gradient = sp.cancel(gradient_contract - expected_gradient) == sp.S.Zero

    print()
    print("[1] AX torsion -> contortion contractions")
    print("-" * 72)
    print(f"  T_abc T^abc                  = {sp.factor(torsion_contract)}"
          f"   [{'OK' if ok_torsion else 'FAIL'}]")
    print(f"  K_abc K^abc                  = {sp.factor(contortion_contract)}"
          f"   [{'OK' if ok_contortion else 'FAIL'}]")
    print(f"  (d_z K_abc)(d_z K^abc)       = {sp.factor(gradient_contract)}"
          f"   [{'OK' if ok_gradient else 'FAIL'}]")

    all_ok = ok_torsion and ok_contortion and ok_gradient

    print()
    print("[2] Four-geometry volume factors and coefficient match")
    print("-" * 72)

    for label, topology in TOPOLOGIES:
        mass_from_engine, eng = engine_mass_expr(topology, None)
        params = eng.data["params"]
        rho_eng = params.get("R", params["r"])
        vol = eng.data["total_volume"]

        target = sp.cancel(6 * vol / (params["kappa"]**2 * rho_eng**2))
        mass_side = sp.cancel(
            vol / params["kappa"]**2
            * contortion_contract.coeff(eta, 2)
        ).subs({rho: rho_eng, kappa: params["kappa"]})
        kinetic_side = sp.cancel(
            vol / params["kappa"]**2
            * gradient_contract.coeff(eta_z, 2)
        ).subs({rho: rho_eng, kappa: params["kappa"]})
        expected = expected_mass_expr(label, params, rho_eng)

        ok_mass_self = sp.cancel(mass_side - target) == sp.S.Zero
        ok_kin_self = sp.cancel(kinetic_side - target) == sp.S.Zero
        ok_engine = sp.cancel(mass_from_engine - target) == sp.S.Zero
        ok_expected = sp.cancel(target - expected) == sp.S.Zero
        row_ok = ok_mass_self and ok_kin_self and ok_engine and ok_expected
        all_ok &= row_ok

        print(f"  {label}:")
        print(f"    Vol_geo(rho)                    = {vol}")
        print(f"    mass from K_abc K^abc            = {sp.factor(mass_side)}")
        print(f"    kinetic from d_z K contraction   = {sp.factor(kinetic_side)}")
        print(f"    6 Vol/(kappa^2 rho^2)            = {sp.factor(target)}")
        print(f"    homogeneous d2V/deta2            = {sp.factor(mass_from_engine)}")
        print("    mass/kinetic/engine/expected     = "
              f"{'PASS' if row_ok else 'FAIL'}")

    print()
    print("=" * 72)
    if all_ok:
        print("RESULT: PASS  (AX contortion gradient independently gives K_geo = M_geo)")
    else:
        print("RESULT: FAIL")
    print("=" * 72)
    return all_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
