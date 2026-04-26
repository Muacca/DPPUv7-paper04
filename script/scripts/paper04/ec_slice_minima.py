"""
ec_slice_minima.py
=========================

Symbolic verification for the EC slice-minimum branch used in paper04.

Target:
  On the homogeneous eta=V=0 slice, the current DPPU EC+NY+Weyl
  potential has an isolated EC-induced radial branch for Nil3 and Sol3,
  but not for T3 or round S3.

  For Nil3 and Sol3:
      R0 = 4 kappa sqrt(-alpha) / sqrt(3)        (alpha < 0)

  The full homogeneous Hessian in (R, eta, V) has the same spin-0
  determinant in both geometries:
      det H_(eta,V)
        = (262144 pi^8 L^2 alpha^2 / 9) * (1 - kappa^4 theta_NY^2)
  The off-diagonal radial/spin-0 Hessian entries vanish at the branch:
      H_Reta = H_RV = 0

  Hence the branch is a full local minimum iff |kappa^2 theta_NY| < 1,
  marginal at equality, and a saddle otherwise.  Sol3 differs from Nil3
  only by a factor 4 in the radial slice potential and H_RR.

Run:
  python -X utf8 script/scripts/paper04/ec_slice_minima.py

Author: Muacca
Date: 2026-04-24
"""

from pathlib import Path
import sys

from sympy import Matrix, Rational, S, cancel, diff, factor, pi, simplify, solve, sqrt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

_DATA_DIR = Path(__file__).resolve().parents[3] / "data"
_DATA_DIR.mkdir(exist_ok=True)
from dppu.utils.tee_logger import setup_log
setup_log(__file__, log_dir=str(_DATA_DIR))

from dppu.action.ec_action import build_veff_ec
from dppu.topology.unified import TopologyType
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant


def build_case(topology):
    veff_ec, _, params = build_veff_ec(
        topology,
        torsion_mode=Mode.MX,
        ny_variant=NyVariant.FULL,
    )
    rho = params.get("R", params["r"])
    eta = params["eta"]
    V = params["V"]
    theta = params["theta_NY"]
    alpha = params["alpha"]
    kappa = params["kappa"]
    L = params["L"]
    slice_veff = cancel(veff_ec.subs({eta: 0, V: 0, theta: 0}))
    d_slice = cancel(diff(slice_veff, rho))
    return {
        "veff": veff_ec,
        "params": params,
        "rho": rho,
        "eta": eta,
        "V": V,
        "theta": theta,
        "alpha": alpha,
        "kappa": kappa,
        "L": L,
        "slice": slice_veff,
        "d_slice": d_slice,
        "solutions": solve(d_slice, rho),
    }


def hessian_data(case, r0):
    variables = (case["rho"], case["eta"], case["V"])
    H = Matrix(
        [
            [cancel(diff(case["veff"], v1, v2)) for v2 in variables]
            for v1 in variables
        ]
    )
    H_at = simplify(H.subs({case["rho"]: r0, case["eta"]: 0, case["V"]: 0}))
    spin0 = H_at[1:3, 1:3]
    return {
        "H": H_at,
        "H_rr": simplify(H_at[0, 0]),
        "H_r_eta": cancel(H_at[0, 1]),
        "H_r_V": cancel(H_at[0, 2]),
        "trace_spin0": simplify(spin0.trace()),
        "H_eta_eta": simplify(spin0[0, 0]),
        "det_spin0": factor(simplify(spin0.det())),
    }


def classify_numeric(H_at, case, theta_value):
    numeric = H_at.subs(
        {
            case["alpha"]: -1,
            case["kappa"]: 1,
            case["L"]: 1,
            case["theta"]: theta_value,
        }
    ).evalf()
    eigs = sorted(float(ev) for ev in numeric.eigenvals().keys())
    if eigs[0] > 1e-8:
        label = "local minimum"
    elif abs(eigs[0]) <= 1e-8:
        label = "marginal"
    else:
        label = "saddle"
    return eigs, label


def main():
    print("=" * 72)
    print("Verification: EC slice minima across T3 / Nil3 / S3 / Sol3")
    print("Target: paper04 Appendix E.5b")
    print("=" * 72)

    cases = {
        "T3": build_case(TopologyType.T3),
        "Nil3": build_case(TopologyType.NIL3),
        "S3": build_case(TopologyType.S3),
        "Sol3": build_case(TopologyType.SOL3),
    }

    overall_ok = True

    print()
    print("[1] eta=V=0 slice potentials")
    print("-" * 72)
    for name, case in cases.items():
        print(f"  {name:5} V_slice = {case['slice']}")
        print(f"        dV/drho = {case['d_slice']}")
        print(f"        stationary solutions = {case['solutions']}")

    t3_ok = cases["T3"]["slice"] == S.Zero and cases["T3"]["d_slice"] == S.Zero
    s3_ok = cases["S3"]["solutions"] == []
    absent_ok = t3_ok and s3_ok

    print()
    print("[2] Nil3 and Sol3 branch formulae")
    print("-" * 72)
    for name, scale in [("Nil3", 1), ("Sol3", 4)]:
        case = cases[name]
        rho = case["rho"]
        alpha = case["alpha"]
        kappa = case["kappa"]
        L = case["L"]
        theta = case["theta"]

        expected_slice = cancel(
            scale * (4 * pi**4 * L * rho / kappa**2 - 64 * pi**4 * L * alpha / (3 * rho))
        )
        expected_r0 = 4 * sqrt(3) * kappa * sqrt(-alpha) / 3
        expected_v0 = scale * 32 * sqrt(3) * pi**4 * L * sqrt(-alpha) / (3 * kappa)
        expected_det = factor(
            Rational(262144, 9) * pi**8 * L**2 * alpha**2
            * (1 - kappa**4 * theta**2)
        )
        expected_hrr = scale * 2 * sqrt(3) * pi**4 * L / (kappa**3 * sqrt(-alpha))
        expected_trace = (
            Rational(128, 27) * sqrt(3) * pi**4 * L * sqrt(-alpha)
            * (27 - 16 * alpha * kappa**2) / kappa
        )

        r0 = case["solutions"][0] if case["solutions"] else None
        v0 = simplify(case["slice"].subs(rho, r0)) if r0 is not None else None
        H = hessian_data(case, r0)

        ok_slice = simplify(case["slice"] - expected_slice) == S.Zero
        ok_r0 = simplify(r0 - expected_r0) == S.Zero
        ok_v0 = simplify(v0 - expected_v0) == S.Zero
        ok_hrr = simplify(H["H_rr"] - expected_hrr) == S.Zero
        ok_det = simplify(H["det_spin0"] - expected_det) == S.Zero
        ok_block = H["H_r_eta"] == S.Zero and H["H_r_V"] == S.Zero
        ok_trace = simplify(H["trace_spin0"] - expected_trace) == S.Zero

        print(f"  {name}:")
        print(f"    V_slice                    = {case['slice']}")
        print(f"    R0                         = {r0}")
        print(f"    V0                         = {v0}")
        print(f"    H_RR                       = {H['H_rr']}")
        print(f"    H_Reta, H_RV               = {H['H_r_eta']}, {H['H_r_V']}"
              f"   [{'PASS' if ok_block else 'FAIL'}]")
        print(f"    tr H_(eta,V)               = {H['trace_spin0']}"
              f"   [{'PASS' if ok_trace else 'FAIL'}]")
        print(f"    det H_(eta,V)              = {H['det_spin0']}")
        print(f"    slice/R0/V0/H_RR/det check = "
              f"{'PASS' if all([ok_slice, ok_r0, ok_v0, ok_hrr, ok_det]) else 'FAIL'}")
        print(f"    block-diagonal/trace check = "
              f"{'PASS' if ok_block and ok_trace else 'FAIL'}")
        overall_ok &= (
            ok_slice and ok_r0 and ok_v0 and ok_hrr and ok_det
            and ok_block and ok_trace
        )

        print("    numeric Hessian labels:")
        for theta_value, expected_label in [(0.5, "local minimum"), (1.0, "marginal"), (1.2, "saddle")]:
            eigs, label = classify_numeric(H["H"], case, theta_value)
            ok_label = label == expected_label
            eig_text = ", ".join(f"{ev:+.4g}" for ev in eigs)
            print(f"      theta_NY={theta_value:>3}: [{eig_text}] -> {label}"
                  f"   [{'OK' if ok_label else 'FAIL'}]")
            overall_ok &= ok_label

    print()
    print("[3] Absent cases")
    print("-" * 72)
    print(f"  T3: flat zero slice potential, no isolated EC radial branch   [{'OK' if t3_ok else 'FAIL'}]")
    print(f"  S3: nonzero constant slope on the round slice, no stationary point   [{'OK' if s3_ok else 'FAIL'}]")
    overall_ok &= absent_ok

    print()
    print("=" * 72)
    if overall_ok:
        print("RESULT: PASS  (EC slice branch present for Nil3 and Sol3; absent for T3 and S3)")
    else:
        print("RESULT: FAIL")
    print("=" * 72)
    return overall_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
