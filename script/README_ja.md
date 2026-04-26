# DPPUv7 — スクリプトディレクトリ

⇒ [English](README.md)

**論文**: "Euclidean Geometric Response Dictionary and Selection Rules for Four Thurston Geometries"（paper04）

Einstein-Cartan + Nieh-Yan フレームワークにおける4つの Thurston 幾何（S³×S¹・T³×S¹・Nil³×S¹・Sol³×S¹）を対象とした数値・記号計算のための Python パッケージ群と実行スクリプト。

---

## ディレクトリ構成

```
script/
├── docs/                      # 技術ドキュメントと規約
├── dppu/                      # メイン Python パッケージ（DPPUv7）
│   ├── geometry/              # 計量・体積形式・構造定数
│   ├── connection/            # Levi-Civita 接続・Contortion・EC 接続
│   ├── curvature/             # Riemann・Ricci・Hodge 双対・Pontryagin・Weyl
│   ├── torsion/               # トーションモード・Ansatz・Nieh-Yan 密度
│   ├── action/                # ラグランジアン・有効ポテンシャル・安定性分類
│   ├── topology/              # 統一エンジン（S³×S¹・T³×S¹・Nil³×S¹・Sol³×S¹）
│   ├── engine/                # 計算パイプライン・ロギング・チェックポイント
│   ├── kk/                    # Kaluza-Klein 光子有効理論（二経路パイプライン）
│   └── utils/                 # 共通ユーティリティ（Levi-Civita 記号・記号計算・可視化）
│
└── scripts/                   # 実行スクリプト
    ├── paper04/               # paper04 固有の解析スクリプト
    ├── proofs/                # 解析的・記号的証明スクリプト
    └── visualize/             # 図表ノートブックとビルドスクリプト
```

### `docs/` — ドキュメント

技術ドキュメントと規約：
- [DPPUv7 Engine CONVENTIONS](docs/CONVENTIONS_ja.md) — エンジンコアの規約と仕様
- [DPPUv7 SymPy guideline](docs/SymPy_guideline_ja.md) — SymPy 使用ガイドラインとベストプラクティス

---

## パッケージ概要（dppu/）

| モジュール | 役割 | 主要クラス・関数 |
|-----------|------|----------------|
| [`geometry`](dppu/geometry/README_ja.md) | 計量・フレーム場定義 | `build_metric`, `frame_field` |
| [`connection`](dppu/connection/README_ja.md) | EC 接続の構築 | `levi_civita`, `contortion`, `ec_connection` |
| [`curvature`](dppu/curvature/README_ja.md) | 曲率テンソル群・Pontryagin・Weyl | `RiemannTensor`, `compute_pontryagin_inner_product`, `WeylTensor` |
| [`torsion`](dppu/torsion/README_ja.md) | トーション構造 | `Mode`, `NyVariant`, `build_torsion_tensor` |
| [`action`](dppu/action/README_ja.md) | 作用・安定性解析 | `build_lagrangian`, `classify_stability` |
| [`topology`](dppu/topology/README_ja.md) | 4 つの Thurston 幾何対応統一エンジン | `UnifiedEngine`, `DOFConfig`, `TopologyType`, `FiberMode` |
| [`engine`](dppu/engine/README_ja.md) | 15 ステップ計算パイプライン | `BaseFrameEngine`, `ComputationLogger`, `CheckpointManager` |
| [`kk`](dppu/kk/README_ja.md) | KK 光子有効理論（Γ×Γ ショートカット＋全 Riemann 検証） | `extract_maxwell`, `extract_mass`, `extract_cs` |
| [`utils`](dppu/utils/README_ja.md) | 共通ユーティリティ | `epsilon_symbol`, `prove_zero`, `set_style` |

---

## 実行スクリプト概要（scripts/）

### paper04/ — paper04 固有解析

#### EC 真空とスペクトル解析

| スクリプト | 説明 |
|-----------|------|
| `ec_slice_minima.py` | Nil³・Sol³ における EC slice minimum ブランチの記号的検証（η=V=0 スライス・孤立極小・Hessian 行列式） |
| `eta_aps_nil3.py` | PPA スピン構造をもつ Nil³ 上の APS スペクトルコアの数値検証：Heisenberg Landau 準位経由で η_APS = +1/2 |

#### デフェクトとトポロジカル電荷

| スクリプト | 説明 |
|-----------|------|
| `defect_localization.py` | 4 幾何全体にわたる AX トーション デフェクト局在の数値・変分ベンチマーク（Rayleigh 商・12 プロファイル） |
| `eta_defect_coefficients.py` | 縮約 1D AX トーション デフェクトの運動エネルギー係数と質量係数の記号的検証：全 4 幾何で K_geo = M_geo |
| `torsional_charge.py` | S³ 上のフレームバンドル正規化トーション電荷の記号的検証：N_top = 6 r0² |

---

### proofs/ — 解析的・記号的証明

| スクリプト | 定理・内容 |
|-----------|----------|
| `sol3_structure.py` | **付録 E.1**：Sol³ 構造定数とフレーム剛性（C¹₀₁ = +1/R, C²₀₂ = −1/R；ε, s 変形に対して剛性） |
| `cs_cancellation.py` | **付録 E.3–E.4**：プロファイルローカル KK Ansatz での Sol³ CS 活性化ゼロ；off-diagonal CS 方向数 = (0, 1, 3, 0) for (T³, Nil³, S³, Sol³) |
| `kk_higgsing.py` | **付録 E.2, E.5**：KK Maxwell 普遍性・Sol³ 上の biaxial Higgsing；Higgsing パターン = (none, uniaxial, triaxial, biaxial) for (T³, Nil³, S³, Sol³) |
| `eta_kinetic_from_contortion.py` | **付録 E.6**：Contortion 勾配からの AX トーション η モード運動係数の独立導出 |
| `weyl_scalar.py` | **付録 E.7**：Levi-Civita Weyl スカラー C²_LC = (0, 4/3R⁴, 0, 16/3R⁴) for (T³, Nil³, S³, Sol³) |
| `aps_zero_t3_s3.py` | **付録 E.10**：T³（PPA）と丸い S³ における APS ベンチマークゼロ：対称スペクトル, η(0) = 0 |
| `landau_levels_nil3.py` | **付録 E.10**：Nil³ における Heisenberg Landau 準位スペクトル（梯子代数による検証） |
| `eta_aps_sol3.py` | **付録 E.10**：コンパクト写像トーラス M_A = T²⋊_A S¹ 上の Sol³ APS η 不変量；スピン構造依存性 |
| `kk_normalization.py` | **付録 E.11**：KK 正規化恒等式 k³D^{tor-CS} = (1/2) k⁴D^{NY}（CZ 恒等式と Stokes の定理から導出） |

---

### visualize/ — 図表

| ファイル | 説明 |
|---------|------|
| `DPPUv7_Paper04_Figures.ipynb` | 4 つの Thurston 幾何にわたる paper04 図表生成 Jupyter ノートブック |
| `_build_paper04_figures_notebook.py` | 図表ノートブックのビルド・再生成スクリプト |

---

## クイックスタート

```bash
# 依存パッケージのインストール
pip install -r requirements.txt

# 全スクリプトの実行（script/ ディレクトリから）
bash run_paper04.sh

# 証明スクリプトのみ実行
bash run_paper04.sh proofs

# paper04 解析スクリプトのみ実行
bash run_paper04.sh paper04

# 出力ディレクトリを指定して実行
bash run_paper04.sh --output-dir /path/to/logs
```
