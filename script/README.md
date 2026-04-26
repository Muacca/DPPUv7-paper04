# DPPUv7 — Scripts Directory

⇒ [日本語](README_ja.md)

**Paper**: "Euclidean Geometric Response Dictionary and Selection Rules for Four Thurston Geometries" (paper04)

Python packages and execution scripts for numerical and symbolic computation in the Einstein-Cartan + Nieh-Yan framework across four Thurston geometries (S³×S¹, T³×S¹, Nil³×S¹, Sol³×S¹).

---

## Directory Structure

```
script/
├── docs/                      # Technical documentation and conventions
├── dppu/                      # Main Python package (DPPUv7)
│   ├── geometry/              # Metric, volume form, structure constants
│   ├── connection/            # Levi-Civita connection, Contortion, EC connection
│   ├── curvature/             # Riemann, Ricci, Hodge dual, Pontryagin, Weyl
│   ├── torsion/               # Torsion modes, Ansatz, Nieh-Yan density
│   ├── action/                # Lagrangian, effective potential, stability classification
│   ├── topology/              # Unified engine (S³×S¹, T³×S¹, Nil³×S¹, Sol³×S¹)
│   ├── engine/                # Computation pipeline, logging, checkpointing
│   ├── kk/                    # Kaluza-Klein photon effective theory (two-route pipeline)
│   └── utils/                 # Common utilities (Levi-Civita symbol, symbolic computation, visualization)
│
└── scripts/                   # Execution scripts
    ├── paper04/               # paper04 analysis scripts
    ├── proofs/                # Analytic and symbolic proof scripts
    └── visualize/             # Figures notebook and build script
```

### `docs/` — Documentation

Technical documentation and conventions:
- [DPPUv7 Engine CONVENTIONS](docs/CONVENTIONS.md) — Engine core conventions and specifications
- [DPPUv7 SymPy guideline](docs/SymPy_guideline.md) — SymPy usage guidelines and best practices

---

## Package Overview (dppu/)

| Module | Role | Key Classes / Functions |
|--------|------|------------------------|
| [`geometry`](dppu/geometry/README.md) | Metric and frame field definitions | `build_metric`, `frame_field` |
| [`connection`](dppu/connection/README.md) | EC connection construction | `levi_civita`, `contortion`, `ec_connection` |
| [`curvature`](dppu/curvature/README.md) | Curvature tensors, Pontryagin, Weyl | `RiemannTensor`, `compute_pontryagin_inner_product`, `WeylTensor` |
| [`torsion`](dppu/torsion/README.md) | Torsion structure | `Mode`, `NyVariant`, `build_torsion_tensor` |
| [`action`](dppu/action/README.md) | Action and stability analysis | `build_lagrangian`, `classify_stability` |
| [`topology`](dppu/topology/README.md) | Unified engine for all four Thurston geometries | `UnifiedEngine`, `DOFConfig`, `TopologyType`, `FiberMode` |
| [`engine`](dppu/engine/README.md) | 15-step computation pipeline | `BaseFrameEngine`, `ComputationLogger`, `CheckpointManager` |
| [`kk`](dppu/kk/README.md) | KK photon effective theory (Γ×Γ shortcut + full Riemann) | `extract_maxwell`, `extract_mass`, `extract_cs` |
| [`utils`](dppu/utils/README.md) | Common utilities | `epsilon_symbol`, `prove_zero`, `set_style` |

---

## Script Overview (scripts/)

### paper04/ — paper04 Analysis

#### EC Vacuum and Spectral Analysis

| Script | Description |
|--------|-------------|
| `ec_slice_minima.py` | Symbolic verification of EC slice-minimum branch on Nil³ and Sol³ (η=V=0 slice, isolated radial minimum, Hessian determinant) |
| `eta_aps_nil3.py` | Numerical verification of APS spectral core on Nil³ with PPA spin structure: η_APS = +1/2 via Heisenberg Landau-level route |

#### Defect and Topological Charge

| Script | Description |
|--------|-------------|
| `defect_localization.py` | Numerical/variational benchmark for AX-torsion defect localization across all four geometries (Rayleigh quotient, 12 profile cases) |
| `eta_defect_coefficients.py` | Symbolic verification of reduced 1D AX-torsion defect kinetic and mass coefficients: K_geo = M_geo across all four geometries |
| `torsional_charge.py` | Symbolic verification of frame-bundle-normalized torsional charge on S³: N_top = 6 r0² |

---

### proofs/ — Analytic and Symbolic Proofs

| Script | Theorem / Content |
|--------|------------------|
| `sol3_structure.py` | **Appendix E.1**: Sol³ structure constants and frame rigidity (C¹₀₁ = +1/R, C²₀₂ = −1/R; rigid under ε, s deformation) |
| `cs_cancellation.py` | **Appendix E.3–E.4**: Sol³ CS activation vanishes in profile-local KK ansatz; off-diagonal CS direction count = (0, 1, 3, 0) for (T³, Nil³, S³, Sol³) |
| `kk_higgsing.py` | **Appendix E.2, E.5**: KK Maxwell universality, biaxial Higgsing on Sol³; Higgsing pattern = (none, uniaxial, triaxial, biaxial) on (T³, Nil³, S³, Sol³) |
| `eta_kinetic_from_contortion.py` | **Appendix E.6**: Independent derivation of AX-torsion η-mode kinetic coefficient from contortion gradient |
| `weyl_scalar.py` | **Appendix E.7**: Levi-Civita Weyl scalar C²_LC = (0, 4/3R⁴, 0, 16/3R⁴) for (T³, Nil³, S³, Sol³) |
| `aps_zero_t3_s3.py` | **Appendix E.10**: APS benchmark zeros on T³ (PPA) and round S³: symmetric spectrum, η(0) = 0 |
| `landau_levels_nil3.py` | **Appendix E.10**: Heisenberg Landau-level spectrum on Nil³ (ladder-algebra verification) |
| `eta_aps_sol3.py` | **Appendix E.10**: Sol³ APS η-invariant on compact mapping torus M_A = T²⋊_A S¹; spin-structure dependence |
| `kk_normalization.py` | **Appendix E.11**: KK normalization identity k³D^{tor-CS} = (1/2) k⁴D^{NY} from CZ identity and Stokes' theorem |

---

### visualize/ — Figures

| File | Description |
|------|-------------|
| `DPPUv7_Paper04_Figures.ipynb` | Jupyter notebook generating paper04 figures across all four Thurston geometries |
| `_build_paper04_figures_notebook.py` | Script to build / regenerate the figures notebook |

---

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run all scripts (from the script/ directory)
bash run_paper04.sh

# Run proofs only
bash run_paper04.sh proofs

# Run paper04 analysis scripts only
bash run_paper04.sh paper04

# Custom output directory
bash run_paper04.sh --output-dir /path/to/logs
```
