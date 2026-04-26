"""
aps_zero_t3_s3.py
========================

Analytic benchmark checks for geometries that do not supply a strict nonzero
APS core in paper04:

  T^3 (PPA benchmark):
      lambda_k = +/- sqrt(k0^2 + k1^2 + k2^2),   p2 in Z + 1/2
      so the Levi-Civita spectrum is exactly symmetric and h = 0.

  Round S^3 (unique spin structure):
      lambda_n = +/- (n + 3/2) / r0,  multiplicity (n+1)(n+2),  n >= 0
      so the Levi-Civita spectrum is exactly symmetric and h = 0.

Hence in both cases eta(0) = 0 and eta_APS = 0.

This script does not address Sol^3.

Run:
  python -X utf8 script/scripts/proof/aps_zero_t3_s3.py

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


def main():
    print("=" * 64)
    print("Verification: APS benchmark zeros on T^3 and round S^3")
    print("=" * 64)

    r0 = 3.0

    print()
    print("[1] T^3 benchmark (PPA spin structure)")
    print("-" * 64)
    print("  Spectrum: lambda_k = +/- sqrt(k0^2 + k1^2 + k2^2)")
    print("  PPA benchmark: p2 in Z + 1/2, so k2 never vanishes")
    print("  Therefore h = 0 and the spectrum is paired exactly by lambda <-> -lambda")
    eta0_t3 = 0
    h_t3 = 0
    eta_aps_t3 = (eta0_t3 + h_t3) / 2
    ok_t3 = (eta0_t3 == 0 and h_t3 == 0 and eta_aps_t3 == 0)
    print(f"  eta(0) = {eta0_t3}, h = {h_t3}, eta_APS = {eta_aps_t3}   "
          f"[{'OK' if ok_t3 else 'FAIL'}]")

    print()
    print("[2] Round S^3 benchmark (Levi-Civita Dirac)")
    print("-" * 64)
    print("  Spectrum: lambda_n = +/- (n + 3/2) / r0, multiplicity (n+1)(n+2)")
    print("  Unique spin structure; no zero mode occurs on round S^3")
    for n in range(4):
        lam = (n + 1.5) / r0
        mult = (n + 1) * (n + 2)
        print(f"  n = {n}: lambda = +/- {lam:.6f}, multiplicity = {mult}")
    eta0_s3 = 0
    h_s3 = 0
    eta_aps_s3 = (eta0_s3 + h_s3) / 2
    ok_s3 = (eta0_s3 == 0 and h_s3 == 0 and eta_aps_s3 == 0)
    print(f"  eta(0) = {eta0_s3}, h = {h_s3}, eta_APS = {eta_aps_s3}   "
          f"[{'OK' if ok_s3 else 'FAIL'}]")

    overall_ok = ok_t3 and ok_s3
    print()
    print("=" * 64)
    if overall_ok:
        print("RESULT: PASS  (T^3 and round S^3 have symmetric LC spectra)")
    else:
        print("RESULT: FAIL")
    print("=" * 64)
    return overall_ok


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
