## Appendix A. 記号・規約

本付録では、本稿で用いる記号規約、4 geometry の略記、two-layer structure の用語、ならびに section 間で混同しやすい 3-form の区別をまとめる。本文の数式は compact に書いているため、比較の前提をこの付録で明示しておく。

### A.1 Geometry と共通 setup

- $T^3$: flat three-torus
- $Nil^3$: Heisenberg-type nil geometry
- $S^3$: round three-sphere
- $Sol^3$: solv geometry

本稿では `4 geometries` とは常にこの 4 つを指す。すべての比較は Euclidean product structure

$$
M^3 \times S^1
$$

の上で行い、 $S^1$ の周長を $L$ と書く。 $r_0$ は homogeneous background を特徴づける基準長、 $R$ は 3 次元空間部の構造定数に現れる幾何学的スケールである。

### A.2 二層構造と応答語彙

- **first layer**: boundary を導入せずに bulk 側で定義できる response
- **second layer**: boundary / interface / cobordism setting で定義される response
- **secondary reading**: 既に定まった Euclidean response dictionary を別言語で読む解釈層

本文ではさらに次の語を用いる。

- **allowed**: 幾何学的にその channel の定義が意味を持つ
- **protected**: 非零 perturbation を許しても当該 channel が自動的には開かない
- **activated**: geometry の構造によって非自明な response が立ち上がる
- **minimal**: active family を持たず benchmark 側に属する

### A.3 主要記号と代表値

- $K_{\rm Mxw}$: Maxwell baseline coefficient
- $S_{\rm CZ}$: Callias-Zanelli inflow action
- $\eta_{\rm APS}^{(3D)}$: APS spectral invariant
- $N_{\rm top}$: frame-bundle-normalized torsional-charge entry
- $C^2_{\rm LC}$: Levi-Civita Weyl scalar
- spin-2 rigidity tag: scaffold-layer label (`trivial / non-rigid / rigid`)
- $A_{\rm KK}$: KK anisotropy marker
- AX torsion ansatz: $T_{ijk}=(2\eta/\rho)\,\epsilon_{ijk}$（DPPU AX mode の paper03ec §3 規約。係数 $2/\rho$ はこの規約に固定）

代表的な benchmark value は次の通りである。導出スケッチおよび検証スクリプトは [Appendix E](DPPUv7-paper04_appE.md) を参照。

| quantity | value | role | derivation |
|---|---|---|---|
| $K_{\rm Mxw}$ | $-L^2/(2r_0^4)$ | universal bulk baseline | [E.2](DPPUv7-paper04_appE.md), [E.5](DPPUv7-paper04_appE.md) |
| $\eta_{\rm APS}^{(3D)}(Nil^3)$ | $+1/2$ | strict spectral core | [E.10](DPPUv7-paper04_appE.md) |
| $N_{\rm top}(S^3)$ | $6r_0^2$ | frame-bundle-normalized torsional-charge core | [E.9](DPPUv7-paper04_appE.md) |
| EC slice branch | present for $Nil^3,Sol^3$ | spin-0 bulk branch | [E.5b](DPPUv7-paper04_appE.md) |
| $m^2(A_1)=m^2(A_2)$ in $Sol^3$ | $-L^2/(2r_0^4)$ | biaxial KK splitting | [E.2](DPPUv7-paper04_appE.md), [E.5](DPPUv7-paper04_appE.md) |
| $C^2_{\rm LC}(Nil^3)$ | $4/(3R^4)$ | scaffold datum | [E.7](DPPUv7-paper04_appE.md) |
| $C^2_{\rm LC}(Sol^3)$ | $16/(3R^4)$ | scaffold datum | [E.7](DPPUv7-paper04_appE.md) |
| $k_{3D}^{\rm tor\text{-}CS}/k_{4D}^{\rm NY}$ | $1/2$ | NY-CS normalization identity | [E.11](DPPUv7-paper04_appE.md) |
| $Q_{\rm best}=\eta(0)/\Delta k_{\rm matter}$ | $1$ | normalization-free invariant | [E.12](DPPUv7-paper04_appE.md) |
| $k_q$ (formal) | $\in \mathbb{Z}$ | reduced CS level integrality | [E.13](DPPUv7-paper04_appE.md) |

### A.4 3-form の区別と connection

本稿では次の 3-form を厳密に分ける。

- $\mathcal{C}_3^{\rm tor}$: torsion-CS form
- $\mathcal{C}_3^{\sigma}$: spinor CS form
- $\mathcal{C}_3^{\rm adj}$: adjoint CS form

CZ inflow は $\mathcal{C}_3^{\rm tor}$ に、APS spectral entry は $\mathcal{C}_3^{\sigma}$ に対応する。本文で両者を混同しないことが重要である。概念的な関係は

$$
{\rm NY}=d\mathcal{C}_3^{\rm tor}
$$

で与えられ、 $\mathcal{C}_3^{\rm adj}$ は APS 変換の中間段階に現れる。

### A.5 Geometric scaffold quantities

observable family そのものではないが、background geometry の distinctness を記述するために次の量を用いる。

- $C^2_{\rm LC}$: Levi-Civita Weyl curvature の大きさ
- spin-2 rigidity tag: squash / shear Hessian の零パターンを要約する scaffold label
- $A_{\rm KK}/K^2$: KK anisotropy の大きさ。ここで

$$
A_{\rm KK}:=\sum_{i=0}^2(m_i^2-\bar m^2)^2,\qquad \bar m^2:=\frac13\sum_{i=0}^2 m_i^2,\qquad K^2:=K_{\rm Mxw}^2
$$

とする

この scaffold layer は additional observable family を定義しないが、 $T^3$ と $Sol^3$ のように local CS core を共有しない geometry の違いを区別する上で不可欠である。とくに $Sol^3$ は CS-inert であっても EC slice minimum と frame rigidity を持つ。

benchmark 値 $A_{\rm KK}/K^2$ は 4 geometry の mass dict から直接計算できる:

- $T^3$: $m_i^2=0$ なので $A_{\rm KK}/K^2=0$
- $Nil^3$ (uniaxial $m_2^2=K_{\rm Mxw}$): $\bar m^2=K/3$, $A_{\rm KK}=2(K/3)^2+(2K/3)^2=(2/3)K^2$
- $S^3$ (triplet $m_i^2=-2L^2/r_0^4$): 全 $i$ で同値なので $A_{\rm KK}/K^2=0$
- $Sol^3$ (biaxial $m_1^2=m_2^2=K_{\rm Mxw}$): $\bar m^2=2K/3$, $A_{\rm KK}=(2K/3)^2+2(K/3)^2=(2/3)K^2$

この symbolic 確認は [`proofs/kk_higgsing.py`](script/scripts/proofs/kk_higgsing.py) で行う。
