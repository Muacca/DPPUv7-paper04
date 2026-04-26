## Appendix E. Derivations and proof sketches

本付録では、本文および他付録で結果のみ提示した statement について、最小限の導出スケッチを集約する。本文の selection rule や observable benchmark がどの計算に支えられているかを、再現可能な形で固定することが目的である。

各小節は次の構成をとる:

- **Statement**: 本文での主張
- **Derivation**: 自己完結した導出スケッチ
- **Verification**: 同等の検証を行うスタンドアローン Python スクリプトへの参照（再現可能）

すべての検証スクリプトは `script/dppu/` に同梱した DPPU ライブラリ、`sympy`, `numpy`, `scipy` のみに依存する。スクリプトはそれぞれ `python -X utf8 <path>` で実行でき、終了コードと `RESULT: PASS/FAIL` 行で判定される。

---

### E.1 Sol³ structure constants and frame rigidity

**Statement** (§2.3, §3.3, App. A):

$$
C^1{}_{01} = +\frac{1}{R}, \qquad C^2{}_{02} = -\frac{1}{R}
$$

であり、これらは squash / shear deformation に対して非依存（frame rigidity）。

**Derivation**:
Sol Lie 代数（DPPU 規約 $[E_b, E_c] = -C^a{}_{bc} E_a$）から left-invariant coframe $\sigma^0 = dz, \sigma^1 = e^z dx, \sigma^2 = e^{-z} dy$ を取り、外微分

$$
d\sigma^0 = 0, \quad d\sigma^1 = \sigma^0 \wedge \sigma^1, \quad d\sigma^2 = -\sigma^0 \wedge \sigma^2
$$

を得る。volume-preserving 変形 $f_1 = (1+\varepsilon)^{2/3}(1+s),\ f_2 = (1+\varepsilon)^{-2/3}/(1+s)$ を入れた deformed coframe $e^1 = R f_1 \sigma^1$ について

$$
de^1 = R f_1 \, d\sigma^1 = R f_1 \, \sigma^0 \wedge \sigma^1 = \frac{e^0 \wedge e^1}{R}
$$

となり $f_1$ が完全にキャンセルする。 $de^2$ も同様。したがって構造定数は $\varepsilon, s$ に非依存。

**Verification**: [`proofs/sol3_structure.py`](script/scripts/proofs/sol3_structure.py)

---

### E.2 Sol³ biaxial Higgsing — `m²(A_i)` from KK reduction

**Statement** (§2.4, §3.5):

$$
m^2(A_0) = 0, \qquad m^2(A_1) = m^2(A_2) = -\frac{L^2}{2 r_0^4} = K_{\rm Mxw}
$$

**Derivation**:
KK 光子場 $A_i$ に対する modified field strength

$$
\tilde{F}_{01} = F_{01} + C^a{}_{01} A_a = F_{01} + \frac{A_1}{R}, \qquad
\tilde{F}_{02} = F_{02} - \frac{A_2}{R}
$$

を Maxwell action $-\tfrac{1}{4} \tilde{F}\_{ij}\tilde{F}^{ij}$ に代入し $S^1$ 方向を積分すると、 $A_1, A_2$ に対する mass 項が $K_{\rm Mxw}$ と同じ係数で出る（同脚自己結合）。一方 $A_0$ は $\tilde{F}$ に現れないため massless。

**Verification**: [`proofs/kk_higgsing.py`](script/scripts/proofs/kk_higgsing.py)

---

### E.2b Maxwell universality — common principal coefficient

**Statement** (§2.3, §5.1 R1):

$$
K_{{\rm Mxw},T^3}=
K_{{\rm Mxw},Nil^3}=
K_{{\rm Mxw},S^3}=
K_{{\rm Mxw},Sol^3}=-\frac{L^2}{2r_0^4}
$$

**Derivation**:
4 geometry の差は KK の modified field strength

$$
\tilde{F}_{ij}=F_{ij}+C^a{}_{ij}A_a
$$

に入る structure-constant 補正として現れるが、quadratic Maxwell principal part そのものは常に abelian 項

$$
{-}{\frac{1}{4}}{F_{ij}}{F^{ij}}
$$

から来る。したがって geometry 依存の構造定数は lower-order の mass / CS / parity-odd channel を変えうる一方、 $F_{ij}F^{ij}$ の係数自体は変えない。共通の product-circle 規約と $r_0$ scaling を用いて $S^1$ を積分すると、この principal coefficient は 4 geometry すべてで

$$
K_{\rm Mxw}=-\frac{L^2}{2r_0^4}
$$

に固定される。ゆえに Maxwell universality は `extractor の偶然の一致` ではなく、principal symbol が topology-blind であることの帰結である。

**Verification**: [`proofs/kk_higgsing.py`](script/scripts/proofs/kk_higgsing.py)

---

### E.3 Sol³ Chern-Simons cancellation — profile-local self-coupling rule

**Statement** (§2.3, §3.3, R2):

$$
{\rm CS}_{\rm Sol^3}=0
\quad \text{for profile-local self-couplings}\quad
\tilde F_{01}=F_{01}+c_{01}(z)A_1,\quad
\tilde F_{02}=F_{02}+c_{02}(z)A_2 .
$$

**Derivation**:
$Sol^3$ の non-abelian cubic correction は構造定数 $C^1{}\_{01} = +1/R, C^2{}\_{02} = -1/R$ に比例する 2 つの寄与を持つ。より一般に、profile-local な opening によってこれらが任意の scalar coefficients $c_{01}(z),c_{02}(z)$ に置き換わっても、CS は $\tilde{F}_{ij}$ が $\{i,j\}$ の **外側** の脚を持つ補正でのみ活性化する（off-diagonal rule, §E.4）。Sol³ の補正は脚が $\{i,j\}$ に含まれる **同脚 (self-referential)** 形をとる:

$$
A_1 \to \tilde{F}_{01}: \text{leg 1} \in \{0,1\}, \quad A_2 \to \tilde{F}_{02}: \text{leg 2} \in \{0,2\}
$$

したがって $c_{01},c_{02}$ の値や符号関係を用いなくても off-diagonal CS coefficient は 0 である。検証スクリプトでは $c_{01},c_{02}$ を独立記号として扱い、さらに current $Sol^3$ dictionary の constant scale と $R_{\rm inv}=1/R(z)$ を profile-local symbol として代入した場合にも CS coefficient が消えることを確認する。ここで主張しているのは reduced KK extractor の profile-local algebraic cancellation であり、 $R'(z)$ などの derivative terms を含む full inhomogeneous field equation の導出ではない。

**Verification**: [`proofs/cs_cancellation.py`](script/scripts/proofs/cs_cancellation.py)

---

### E.4 Off-diagonal CS rule

**Statement** (§2.3, §5.1 R2):
CS channel は補正項が $\tilde{F}_{ij}$ の脚 $\{i,j\}$ の外側にあるときのみ活性化する。

**Derivation**:
KK 補正 $\tilde{F}\_{ij} = F\_{ij} + (\text{coupling})\, A_k$ について、3 点結合の構造から CS 寄与は $\epsilon^{ijk}$ 型に展開される。 $k \in \{i, j\}$ なら $\epsilon^{ijk}$ は反対称性により消え、 $k \notin \{i,j\}$ のときのみ非自明な off-diagonal CS direction を寄与する。4 幾何の補正型を分類すると:

| 型 | 例 | 補正 | CS |
|---|---|---|---|
| 構造定数なし | $T^3$ | — | 0 |
| off-diagonal | $Nil^3$ | $A_2 \to \tilde{F}_{01}$（脚 2 ∉ {0,1}） | 1 方向 |
| off-diagonal (等方) | $S^3$ | $\varepsilon_{ijk}$ 項（全脚 off-diagonal） | 3 方向 |
| self-referential | $Sol^3$ | $A_1 \to \tilde{F}_{01}$（脚 1 ∈ {0,1}） | 0 |

これは「構造定数の存在」ではなく「補正型の幾何学的位置」が CS 活性化を決めることを示す。

**Verification**: [`proofs/cs_cancellation.py`](script/scripts/proofs/cs_cancellation.py)（off-diagonal count をスクリプトで構造的に確認）

---

### E.5 KK Higgsing four-type classification

**Statement** (§2.4, §5.1 R3):
KK Higgsing pattern は `none / uniaxial / triaxial / biaxial` の 4 type に分類され、 $T^3, Nil^3, S^3, Sol^3$ の順にちょうど一つずつ実現される。

**Derivation**:
KK reduction で得られる mass 行列 $m^2(A_i)$ を 4 幾何で並列に計算すると:

| 幾何 | mass dict | massive 数 | 型 |
|---|---|---|---|
| $T^3$  | $\{\}$ | 0 | none |
| $Nil^3$ | $\{2: K_{\rm Mxw}\}$ | 1 | uniaxial |
| $S^3$  | $\{0,1,2: -2L^2/r_0^4\}$ | 3 | triaxial |
| $Sol^3$ | $\{1, 2: K_{\rm Mxw}\}$ | 2 | biaxial |

4 type すべてが過不足なく実現されることが直接確認される。

**Verification**: [`proofs/kk_higgsing.py`](script/scripts/proofs/kk_higgsing.py)

---

### E.5b Nil³ / Sol³ EC slice minima

**Statement** (§2.4, §3.2, §5.1 R4):
current DPPU の homogeneous EC+NY+Weyl potential では、 $\eta=V=0$ 断面の EC-induced slice minimum が $Nil^3\times S^1$ と $Sol^3\times S^1$ に存在し、 $T^3$ と round $S^3$ には存在しない。

$Nil^3$ と $Sol^3$ の slice potentials は

$$
V_{\rm eff}^{Nil^3}(R;\alpha)=\frac{4\pi^4 LR}{\kappa^2}-\frac{64\pi^4L\alpha}{3R},
\qquad
V_{\rm eff}^{Sol^3}(R;\alpha)=4V_{\rm eff}^{Nil^3}(R;\alpha).
$$

この "4 倍" は $Sol^3$ の 2 個の対称な構造定数 $C^1{}\_{01}=+1/R$, $C^2{}\_{02}=-1/R$ が $Nil^3$ の 1 個の構造定数 $C^2{}\_{01}=+1/R$ に比べて Weyl curvature を 4 倍にすること（ $C^2\_{\rm LC}(Sol^3)=4C^2\_{\rm LC}(Nil^3)$, [E.7](DPPUv7-paper04_appE.md)）と整合する幾何学的比である。

したがって $\alpha<0$ では共通の stationary radius

$$
R_0=\frac{4\kappa}{\sqrt{3}}\sqrt{-\alpha}
$$

を持ち、

$$
V_0^{Nil^3}=\frac{32\sqrt{3}\,\pi^4L\sqrt{-\alpha}}{3\kappa},
\qquad
V_0^{Sol^3}=4V_0^{Nil^3}
$$

となる。さらに full homogeneous Hessian の spin-0 block は両 geometry で同じ determinant

$$
\det H_{(\eta,V)}=\frac{262144\pi^8L^2\alpha^2}{9}\left(1-\kappa^4\theta_{\rm NY}^2\right)
$$

を持つ。ここで停留点 $R=R_0,\eta=V=0$ では off-diagonal 成分が

$$
\partial_R\partial_\eta V|_{*}=0,\qquad
\partial_R\partial_V V|_{*}=0
$$

となり、3×3 Hessian は $(R)\oplus(\eta,V)$ に block-diagonal となる（[`paper04/ec_slice_minima.py`](script/scripts/paper04/ec_slice_minima.py) で symbolic に確認）。さらに spin-0 block は

$$
H_{\eta\eta}+H_{VV}>0 \quad\text{for}\quad |\kappa^2\theta_{\rm NY}|\le 1
$$

を満たす。したがって $\det H_{(\eta,V)}>0$ と合わせると spin-0 block は正定値、radial 成分 $H_{RR}^{Sol^3}=4H_{RR}^{Nil^3}>0$ と合わせて 3×3 Hessian 全体が正定値となる。この logical chain から、

$$
|\kappa^2\theta_{\rm NY}|<1 \Rightarrow \text{local minimum}, \qquad
|\kappa^2\theta_{\rm NY}|=1 \Rightarrow \text{marginal}, \qquad
|\kappa^2\theta_{\rm NY}|>1 \Rightarrow \text{saddle}
$$

が Nil³ と Sol³ に共通する EC slice-minimum criterion である。

**Derivation**:
`paper04/ec_slice_minima.py` は current `dppu` の `Mode.MX`, `NyVariant.FULL` homogeneous potential から上の slice potentials, stationary radius, $V_0$, full Hessian を直接計算する。 $T^3$ は $V_{\rm eff}|_{\eta=V=0}=0$ の flat zero slice で isolated radial branch を持たず、round $S^3$ は

$$
V_{\rm eff}^{S^3}(r)=-\frac{24\pi^2Lr}{\kappa^2}
$$

の slope が非零なので stationary point を持たない。 $Nil^3$ の criterion は revised paper03ec [4] §5.2 / Theorem 6 と一致するが、paper04 では同じ検証を $Sol^3$ へ拡張して 4 geometry inventory の row を固定する。

**Verification / source**: [`paper04/ec_slice_minima.py`](script/scripts/paper04/ec_slice_minima.py)

---

### E.6 Reduced 1D $\eta$-mode defect-localization benchmark

**Statement** (§3.4):
slice-wise homogeneous background 上の reduced 1D AX-torsion ansatz $\eta=\eta(z)$ に対し、defect 局在化の benchmark は

$$
\rho(z)=r(z)\quad (T^3,S^3), \qquad \rho(z)=R(z)\quad (Nil^3,Sol^3)
$$

を用いて Sturm-Liouville 形

$$
{-}{\frac{d}{dz}}\!\left[K_{\rm geo}(z)\frac{df}{dz}\right] + M_{\rm geo}(z)f = E f,
\qquad
K_{\rm geo}(z)=M_{\rm geo}(z)=c_{\rm geo}\rho(z)
$$

で与えられる。係数は

$$
c_{T^3}=c_{Nil^3}=c_{Sol^3}=\frac{96\pi^4 L}{\kappa^2},
\qquad
c_{S^3}=\frac{12\pi^2 L}{\kappa^2}
$$

であり、少なくとも Gaussian $\rho$-dip family では 4 geometry すべてに対して解析的に束縛状態が存在する。

**Derivation**:
ここで主張しているのは full inhomogeneous field equation の厳密導出ではなく、AX torsion mode に限定した reduced 1D defect benchmark である。mass 側は current `dppu` の homogeneous engine から exact に得られ、kinetic 側は AX contortion gradient を $z$ 方向にだけ開いた reduced ansatz から得られる。

まず exact homogeneous datum として

$$
M_{\rm geo}=
\left.\frac{\partial^2 V_{\rm eff}}{\partial\eta^2}\right|_{\eta=0}=
c_{\rm geo}\rho
$$

が current `dppu` から直接得られる。具体的には

$$
M_{T^3}=M_{Nil^3}=M_{Sol^3}=\frac{96\pi^4 L}{\kappa^2}\rho,
\qquad
M_{S^3}=\frac{12\pi^2 L}{\kappa^2}\rho
$$

である。

次に AX torsion ansatz

$$
T_{ijk}=\frac{2\eta}{\rho}\,\epsilon_{ijk}
$$

（係数 $2/\rho$ は paper03ec §3 の DPPU AX 規約; [App. A.3](DPPUv7-paper04_appA.md)）と contortion

$$
K_{abc}=\frac{1}{2}(T_{abc}+T_{bca}-T_{cab})
$$

を用いる。AX torsion は totally antisymmetric なので $K_{abc}=\tfrac{1}{2}T_{abc}\propto(\eta/\rho)\epsilon_{abc}$。EC 作用に代入したとき、 $\eta=\eta(z)$, $\rho=\rho(z)$ の slice-wise homogeneous ansatz では、 $K_{abc}K^{abc}$ contraction が mass term $M_{\rm geo}\eta^2$ を与え、 $(\partial_zK_{abc})(\partial_zK^{abc})$ contraction が kinetic term $K_{\rm geo}(\partial_z\eta)^2$ を与える。

後者の kinetic 部分について、 $\eta$ にだけ $z$ 依存を持たせ $\rho$ を slice-wise constant とみなす reduced ansatz の下では

$$
\partial_z K_{abc}\propto\frac{\partial_z\eta}{\rho}\epsilon_{abc}
$$

となる。その自己縮約から、index 縮約の組み合わせ係数 $\epsilon_{abc}\epsilon^{abc}=6$、 $\rho$ scaling 因子 $1/\rho^2$、体積積分 ${\rm Vol}_{\rm geo}(\rho)$ の 3 因子が出る。したがって gradient part は

$$
S^{(2)}_{\rm grad}=
\frac{1}{2}\int dz\,K_{\rm geo}(z)\,(\partial_z\eta)^2,
\qquad
K_{\rm geo}(z)=\frac{6\,{\rm Vol}_{\rm geo}(\rho)}{\kappa^2\rho^2}
$$

となる（[`proofs/eta_kinetic_from_contortion.py`](script/scripts/proofs/eta_kinetic_from_contortion.py) で sympy 検証）。4 geometry の volume factor を代入すると、ちょうど

$$
K_{\rm geo}(z)=c_{\rm geo}\rho(z)=M_{\rm geo}(z)
$$

を得る。 $K_{\rm geo}=M_{\rm geo}$ の一致は kinetic 側 (contortion gradient) と mass 側 (homogeneous engine の $\partial_\eta^2V_{\rm eff}$) が独立に計算され、それでも係数が一致するという非自明な事実である。したがって current codebase が支えるのは「general defect channel の exact theorem」ではなく、「reduced 1D $\eta$ benchmark における exact mass datum と reduced kinetic datum の一致」である。

Gaussian dip

$$
\rho(z) = \rho_0(1 - A e^{-z^2/w^2}), \qquad 0 < A < 1
$$

に対して試験関数 $f_\lambda(z) = e^{-\lambda |z|}$ を用いると、Rayleigh 商は

$$
\frac{E_{\rm var}(\lambda)}{M_0}=
(1+\lambda^2)\left[1 - A\lambda w\sqrt{\pi}\,e^{\lambda^2 w^2}\mathrm{erfc}(\lambda w)\right],
\qquad
M_0 = c_{\rm geo}\rho_0
$$

となる。したがって $\lambda \to 0$ で

$$
\frac{E_{\rm var}(\lambda)}{M_0}=
1 - A\sqrt{\pi}\,w\lambda + O(\lambda^2) < 1
$$

が従い、十分小さい $\lambda > 0$ に対して少なくとも 1 つの bound state が存在する。係数 $c_{\rm geo}$ は比 $E_{\rm var}/M_0$ の中で完全に相殺されるため、この結論は 4 geometry に共通である。

代表的な数値スキャンとして $(A, w/\rho_0) = (0.3, 1.0), (0.5, 2.0), (0.7, 1.0)$ の 3 profile を 4 geometry すべてに対して調べると、全 12 例で

$$
E_{\min}<M_0=c_{\rm geo}\rho_0
$$

が成り立ち **LOCALIZED** を得る。

ここで E.6 が与えるのはあくまで `kinematic localization` である点に注意する。すなわち、上の reduced benchmark から従うのは radial dip が bound state を作りうるという共通事実までであり、これだけでは torsional-charge / parity-odd activation は決まらない。実際に local parity-odd activated response になるためには CS-active structure が必要である。Nil³ は R2 の CS 活性と R4 の EC slice minimum を同時に持ち、S³ は R2 の CS 活性と triplet/nonlinear opening を持つ。Sol³ は R4 の EC slice minimum を持つが R2 の CS 活性を持たないため、EC-supported だが CS-inert な entry になる。T³ はいずれも持たない。したがって本文 table の `trivial / activated` は局在化の有無ではなく、この共通 $\eta$ benchmark が後段の active structure とどう結びつくかを表している。

**Verification**: [`paper04/eta_defect_coefficients.py`](script/scripts/paper04/eta_defect_coefficients.py), [`paper04/defect_localization.py`](script/scripts/paper04/defect_localization.py)

---

### E.7 Background Weyl scalar `C²_LC`

**Statement** (§4.5, App. A):

$$
C^2_{\rm LC}(T^3) = 0, \quad C^2_{\rm LC}(Nil^3) = \frac{4}{3R^4}, \quad C^2_{\rm LC}(S^3) = 0, \quad C^2_{\rm LC}(Sol^3) = \frac{16}{3R^4}
$$

**Derivation**:
各 geometry の left-invariant coframe から Levi-Civita 接続 $\omega^a{}_b$ を Koszul 公式で計算し、Riemann 曲率 $R^a{}\_b = d\omega^a{}\_b + \omega^a{}\_c \wedge \omega^c{}\_b$ から Weyl tensor $C\_{abcd}$ を抜き出して $C^2\_{\rm LC} = C\_{abcd}C^{abcd}$ を評価する。 $T^3$ は flat、 $S^3$ は等角平坦のため $C^2 = 0$。 $Nil^3$ は Bianchi II 型、 $Sol^3$ は可解幾何の Weyl がそれぞれ非零値を持つ。 $Sol^3$ が $Nil^3$ の **4 倍** という比は構造定数の符号対称性 (`±1/R`) から幾何学的に従う。

**Verification**: [`proofs/weyl_scalar.py`](script/scripts/proofs/weyl_scalar.py)

---

### E.8 Three CS 3-forms の区別

**Statement** (§4.1, App. A):
本稿で扱う 3-form は

- $\mathcal{C}_3^{\rm tor} = e_a \wedge T^a$ (torsion-CS, EC connection)
- $\mathcal{C}\_3^{\sigma} = {\rm Tr}\_\sigma(\omega \wedge d\omega + \tfrac{2}{3}\omega^3)$ (spinor CS, LC, APS final)
- $\mathcal{C}\_3^{\rm adj} = {\rm Tr}\_{\rm adj}(\omega \wedge d\omega + \tfrac{2}{3}\omega^3)$ (adjoint CS, LC, APS intermediate)

の 3 種であり、CZ inflow は $\mathcal{C}_3^{\rm tor}$、APS spectral entry は $\mathcal{C}_3^{\sigma}$ を用いる。

**Derivation**:
$\mathcal{C}\_3^{\rm tor}$ は Nieh-Yan の primitive: $d\mathcal{C}\_3^{\rm tor} = T^a \wedge T_a + R_{ab} \wedge e^a \wedge e^b = {\rm NY}$。Stokes により $\int\_{X^4} {\rm NY} = \int_{Y} \mathcal{C}_3^{\rm tor}$ が境界 inflow を与える。 $\mathcal{C}_3^\sigma, \mathcal{C}\_3^{\rm adj}$ は Levi-Civita 接続のみの関数で、APS index theorem の幾何学的部分（ $\eta(0)$ 計算）に現れる。ただし両者の関係は単純な「トレース正規化で常に 2 倍」ではない。 $Nil^3$ benchmark では quadratic trace ratio は 1 だが cubic trace ratio は $-2$ であり、さらに APS final sign には向き規約 $\eta\_{\rm APS} = -\int_Y \mathcal{C}_3^\sigma$ が入る。したがって $\mathcal{C}_3^\sigma$ と $\mathcal{C}_3^{\rm adj}$ の対応は representation と convention に依存し、本稿では §E.10 の明示計算に基づいて区別する。

したがって $\eta_{\rm APS}$ (spectral) と $N_{\rm top}$ (frame-bundle-normalized torsional charge) は同じ "parity-odd boundary quantity" でも異なる 3-form sector に属し、単純な差し引きは意味を持たない。

**Verification**: 本節は分類的議論（代数的トレース構造）のみで完結する。3-form の構造定義は [§4.1](DPPUv7-paper04_sec04.md) と [App. A.4](DPPUv7-paper04_appA.md) を参照。

---

### E.9 `N_top = 6r₀²` for S³ (frame-bundle-normalized torsional charge)

**Statement** (§4.2, §5.2 B3):

$$
N_{\rm top}(S^3) = \frac{1}{4\pi^2} \int_{S^3} \mathcal{C}_3^{\rm tor} = 6 r_0^2
$$

**Derivation**:
導出は 2 段階に分ける。まず幾何学的部分として、round $S^3$ の left-invariant structure constants

$$
C^i{}_{jk} = \frac{4}{r_0}\,\epsilon^i{}_{jk}
$$

から torsion 2-form

$$
T^i = \frac{1}{2} C^i{}_{jk} e^j \wedge e^k
$$

を作ると、

$$
\mathcal{C}_3^{\rm tor} = e_i \wedge T^i= \frac{1}{2} C^i{}_{jk}\, e_i \wedge e^j \wedge e^k= \frac{12}{r_0}\,{\rm vol}_{S^3}
$$

を得る。したがって $\mathrm{Vol}(S^3) = 2\pi^2 r_0^3$ を用いて

$$
\int_{S^3} \mathcal{C}_3^{\rm tor} = \frac{12}{r_0} \cdot 2\pi^2 r_0^3 = 24\pi^2 r_0^2
$$

ここまでは幾何計算である。次に位相正規化として、frame bundle $\pi_3(SO(3)) = \mathbb{Z}$ の生成元に対する Chern class 単位を $T_R = 1$ の vector/frame 規約で取ると、正規化は $1/(4\pi^2)$ に固定される（ $8\pi^2$ は adjoint $T_R = 2$ に対応）。したがって

$$
N_{\rm top} = \frac{24\pi^2 r_0^2}{4\pi^2} = 6 r_0^2
$$

本稿では vector/frame trace 規約 $T_R=1$ を採用する。これは frame bundle 上の Pontryagin class を 4D Nieh-Yan 作用 (§E.11) の $1/(8\pi^2)$ 正規化と整合させる選択であり、adjoint 規約 $T_R=2$ を採れば 3D 側正規化も $1/(8\pi^2)$ となる。本稿は $T_R=1$ 一系で全 normalization chain を通す。

本稿では $r_0$ に追加の量子化条件は課していないので、ここで得られるのは frame-bundle-normalized torsional charge の値 $N_{\rm top}=6r_0^2$ である。benchmark choice $r_0=3$ では $N_{\rm top}=54$ となるが、一般の連続 $r_0$ に対する整数性は別途の geometric quantization condition を要する。

**Verification**: [`paper04/torsional_charge.py`](script/scripts/paper04/torsional_charge.py)

---

### E.10 `η_APS^(3D)(Nil³) = +1/2` (local APS spectral core)

**Statement** (§4.2, §5.2 B2):

$$
\eta_{\rm APS}^{(3D)}(Nil^3) = +\frac{1}{2}
$$

**Derivation**:
導出は「低い Landau-level 構造の確認」と「CS 積分からの値決定」に分ける。まず $Nil^3$ 上の Levi-Civita Dirac 演算子を PPA spin structure $(\epsilon_0, \epsilon_1, \epsilon_2) = (+,+,-)$ で Heisenberg モード分解すると、 $p_2 \in \mathbb{Z} + 1/2$ により $p_2 = 0$ セクターが排除され、各 $p_2$ で

$$
\mu_n^+ = 2 n \omega |k_2| + k_2^2, \quad \mu_n^- = 2(n+1) \omega |k_2| + k_2^2, \quad \omega = |k_2|/r_0
$$

を得る。この固有値式は、 $Nil^3$ の Heisenberg group structure $[E_0,E_1]=(1/r_0)E_2$ から直接得られる。各 $p_2$ Fourier mode を fix すると、 $(E_0,E_1)$ は Heisenberg algebra に従う ladder 演算子 $a,a^\dagger$ を生成し、Levi-Civita Dirac 演算子は

$$
D=\gamma^0E_0+\gamma^1E_1+ip_2\gamma^2+(\text{spin connection})
$$

の形で書ける。 $D^\dagger D$ を ladder 基底で対角化すると、各 $p_2$ sector で Landau-level structure $\mu_n=2n\omega|k_2|+k_2^2$, $\omega=|k_2|/r_0$ が closed form で出る。spin up/down で zero-point shift が異なるため $\mu_n^\pm$ の対が現れる（[`proofs/landau_levels_nil3.py`](script/scripts/proofs/landau_levels_nil3.py) で sympy ladder 構成と固有値の closed-form 確認）。本稿が用いるのは Heisenberg モード分解後の各 sector における ladder algebra closed-form であり、compact quotient $\Gamma\backslash{\rm Nil}^3$ 全体の spectral theorem は別問題である。したがって $n = 0$ spin-up branch が spectral asymmetry の唯一の源であり、PPA では核は

$$
h = \dim \ker D = 0
$$

である。

値そのものは、 $Nil^3$ の structure constants $C^2{}\_{01}=+1/r_0$ から Levi-Civita connection $\omega^{ab}$ を作り、これを adjoint generators $J\_{ab}$ と spinor generators $\sigma_{ab}=[\gamma_a,\gamma_b]/2$ に結合して Chern-Simons 3-form を評価することで決める。具体的には、 $Nil^3$ benchmark $r_0=3$ で

$$
\int_Y {\rm Tr}_{\rm adj}(\omega\wedge d\omega)=-\frac{1}{4},\qquad
\int_Y \frac{2}{3}{\rm Tr}_{\rm adj}(\omega^3)=+\frac{1}{2},
$$

であり、合計

$$
\int_Y \mathcal{C}_3^{\rm adj} = +\frac{1}{4}
$$

を得る。この中間値は [`paper04/eta_aps_nil3.py`](script/scripts/paper04/eta_aps_nil3.py) で再現される。一方、spinor 側では quadratic trace ratio は $1$ だが cubic trace ratio が

$$
\frac{{\rm Tr}_\sigma(A^3)}{{\rm Tr}_{\rm adj}(A^3)} = -2
$$

となるため、

$$
\int_Y \mathcal{C}_3^{\sigma} = -\frac{1}{2}
$$

が従う。DPPU の向き規約では

$$
\eta_{\rm APS}^{(3D)} = -\int_Y \mathcal{C}_3^\sigma
$$

である。なお、この cubic ratio は生成子を明示すると直接確認できる。3D Euclidean gamma 代数から

$$
\sigma_{ab}=\frac{1}{2}[\gamma_a,\gamma_b]=i\epsilon_{ab}{}^c\gamma_c
$$

であり、spinor 側の基底 $T_i=i\gamma_i$ について

$$
{\rm Tr}_\sigma(T_iT_jT_k)=2\,\epsilon_{ijk}
$$

を得る。一方、adjoint 基底 $J_i$ では

$$
{\rm Tr}_{\rm adj}(J_iJ_jJ_k)=-\,\epsilon_{ijk}
$$

なので、同じ接続係数 $A=a^i T_i \leftrightarrow a^i J_i$ を代入すると

$$
\frac{{\rm Tr}_\sigma(A^3)}{{\rm Tr}_{\rm adj}(A^3)}=-2
$$

となる。

したがって

$$
\eta_{\rm APS}^{(3D)}(Nil^3) = +\frac{1}{2}
$$

である。したがって APS 公式を逆に用いれば

$$
\eta(0) = 2\eta_{\rm APS} - h = 1
$$

を得る。Landau-level 解析は「どの branch が asymmetry を担うか」と「PPA で $h=0$ になること」を与え、値の決定は Levi-Civita CS 積分と spinor/adjoint の cubic trace ratio によって固定される。

EC 補正 $\delta = -\kappa_{\rm tor}/(4r_0)$ は一様シフトのみで、通常パラメータ範囲では零交叉を生まないため $\eta^{\rm EC}(0) = \eta^{\rm LC}(0)$。

比較 benchmark として、 $T^3$ の PPA spin structure と round $S^3$ の Levi-Civita Dirac はともに固有値が厳密に $\lambda \leftrightarrow -\lambda$ で対になり $h=0$ なので $\eta_{\rm APS}=0$ である。これに対し $Sol^3$ では local CS 積分自体は零だが、compact quotient / spin structure に依存する global spectral branch が現れうる。これを次節で benchmark として明示する。

**Verification**: Nil³ strict value は [`paper04/eta_aps_nil3.py`](script/scripts/paper04/eta_aps_nil3.py)、 $T^3/S^3$ benchmark zero は [`proofs/aps_zero_t3_s3.py`](script/scripts/proofs/aps_zero_t3_s3.py)

---

### E.10b `η_APS^(3D)(Sol³)` on the compact mapping-torus benchmark

**Statement** (§4.2, §5.2 B4):

compact Sol benchmark

$$
M_A = T^2 \rtimes_A S^1,\qquad
A=\begin{pmatrix}2&1\\
1&1\end{pmatrix}
$$

と、追跡している 2 つの benchmark spin structure `Sol-P / Sol-A` に対して

$$
\eta_{\rm APS}^{(3D)}(Sol\text{-}P)=1,\qquad
\eta_{\rm APS}^{(3D)}(Sol\text{-}A)=0
$$

が成り立つ。ここで `Sol-P` は base periodic / fiber periodic、`Sol-A` は base anti-periodic / fiber periodic の benchmark を表す。これは universal な $Sol^3$ の値 ではなく、この compact quotient と spin structure choice に依存する global spectral branch である。

**Derivation**:
$Sol^3$ の structure constants

$$
C^1{}_{01}=+\frac{1}{r_0},\qquad
C^2{}_{02}=-\frac{1}{r_0}
$$

から Levi-Civita connection を作ると、adjoint / spinor のいずれの Chern-Simons 3-form も

$$
\int_Y \mathcal{C}_3^{\rm adj}=0,\qquad
\int_Y \mathcal{C}_3^\sigma=0
$$

となる。したがって $Sol^3$ では local CS 積分だけでは APS 値は決まらない。

次に left-invariant orthonormal frame で Dirac 演算子を調べると、spin connection の zero-order 項が厳密に相殺し

$$
D_{Sol}=\gamma^a E_a
$$

となる。mapping torus 上で fiber Fourier mode $k\in\Lambda^\ast$ ごとに分解すると、 $k\neq 0$ sector の zero-mode 方程式は局所的な $U(1)$ 回転で

$$
(\partial_t \pm m_k(t))\phi = 0,\qquad m_k(t)>0
$$

へ帰着される。したがって periodic / anti-periodic いずれの base boundary condition でも $k\neq 0$ から零モードは生じない。

残る $k=0$ sector では constant spinor だけが候補となる。base 方向 periodic の `Sol-P` では 2 つの constant spinor が quotient に descend するため

$$
h(Sol\text{-}P)=2,
$$

一方 base 方向 anti-periodic の `Sol-A` では descend せず

$$
h(Sol\text{-}A)=0
$$

である。ここで `Sol-P` における $h=2$ は、monodromy 

$$
A=\left(\begin{smallmatrix}2&1\\
1&1\end{smallmatrix}\right)
\in SL(2,\mathbb{Z})
$$

の standard Spin(2) lift（paper03ec から踏襲する benchmark choice）を採用した値である。 $A$ の Spin lift は符号 2 通りあり、別の lift では $h(Sol\text{-}P)=0$ となる場合がある。本稿の `Sol-P / Sol-A` 値は、この benchmark spin lift に固定したうえでの compact mapping-torus 値であり、 $Sol^3$ の universal な spectral entry ではないことに注意する。さらに

$$
T=i\sigma_2 K
$$

と取ると $T\sigma_jT^{-1}=-\sigma_j$ が成り立ち、 $E_a$ は実ベクトル場なので $TDT^{-1}=-D$ となる。 $T$ は antiunitary であり、 $D\psi=\lambda\psi$ ($\lambda$ 実) から $D(T\psi)=-\lambda(T\psi)$ が従う。すなわち nonzero spectrum は $\lambda\leftrightarrow-\lambda$ で pairing される。したがって反線形時間反転対称性により両 spin structure で

$$
\eta(0)=0
$$

が従う。したがって APS の定義

$$
\eta_{\rm APS}^{(3D)}=\frac{\eta(0)+h}{2}
$$

から

$$
\eta_{\rm APS}^{(3D)}(Sol\text{-}P)=1,\qquad
\eta_{\rm APS}^{(3D)}(Sol\text{-}A)=0
$$

を得る。

ゆえに $Sol^3$ の spectral entry は、現 benchmark では $Nil^3$ のような local Levi-Civita CS core ではなく、compact quotient / spin structure に支えられた global kernel branch である。より一般の hyperbolic monodromy quotient ではこの global branch がさらに変化しうるが、本稿では上の $M_A$ benchmark を基準として用いる。

**Verification**: [`proofs/eta_aps_sol3.py`](script/scripts/proofs/eta_aps_sol3.py)

---

### E.11 KK reduction normalization identity

**Statement** (§6.1, App. D.3):

$$
k_{3D}^{\rm tor\text{-}CS} = \frac{1}{2} k_{4D}^{\rm NY}
$$

**Derivation**:
DPPU 4D NY 作用を独立記号 $k_{4D}^{\rm NY}$ で書く:

$$
S_{\rm NY}^{(4D)} = \frac{k_{4D}^{\rm NY}}{8\pi^2} \int_{X^4} {\rm NY}
$$

CZ 恒等式 ${\rm NY} = d\mathcal{C}_3^{\rm tor}$ と Stokes 定理を $X^4 = M^3 \times [0,\beta]$ に適用すると

$$
\int_{X^4} {\rm NY} = \int_{Y} \mathcal{C}_3^{\rm tor} = 4\pi^2 \cdot N_{\rm top}
$$

（最後の等式は §E.9 の正規化）。代入して

$$
S_{\rm eff}^{(3D)} \supset \frac{k_{4D}^{\rm NY}}{8\pi^2} \cdot 4\pi^2 \cdot N_{\rm top} = \frac{k_{4D}^{\rm NY}}{2} \cdot N_{\rm top}
$$

正規化基底

$$
N_{\rm top}:=\frac{1}{4\pi^2}\int_Y \mathcal{C}_3^{\rm tor}
$$

にかかる係数が 3D torsion-CS coefficient と同定され、したがって

$$
k_{3D}^{\rm tor\text{-}CS} = \frac{k_{4D}^{\rm NY}}{2}
$$

を得る。 $1/2$ 因子は $4\pi^2/8\pi^2$ から来て $r_0$ 非依存である。ここで用いているのは frame-bundle-normalized torsional charge basis としての $N_{\rm top}$ であり、本稿は一般の連続 $r_0$ に対する整数性そのものは仮定しない。さらに untwisted product circle 上の実 bosonic KK 分解では、Fourier モード $n$ と $-n$ は複素共役な対であり、質量や even weight は $n^2$ のみで決まる。円周反転 $\tau \mapsto -\tau$ はこの 2 モードを交換するだけで 4D kinetic operator を保つので、parity-odd residue は一般に $n\,w(n^2)$ の形をとり、KK tower $(n \neq 0)$ の寄与は pairwise に厳密消滅する。

**Verification**: [`proofs/kk_normalization.py`](script/scripts/proofs/kk_normalization.py)

---

### E.12 Normalization-free quotient `Q_best = η(0)/Δk_matter = 1`

**Statement** (§6.1, App. D.3):

$$
Q_{\rm best} := \frac{\eta(0)}{\Delta k_{\rm matter}} = 1
$$

**Derivation**:
分子 $\eta(0)^{(3D)} = 1$ は §E.10 から。分母 $\Delta k_{\rm matter} = h_{\rm rep}/2 = 1$ は Redlich parity shift [9] で、表現次元 $h_{\rm rep} = 2$（DPPU minimal matter assumption）から theory-derived。両者は **doubly normalization-free**:

- 分子は spectral 不変量で 4D APS bridge を経由しない
- 分母は $h_{\rm rep}$ のみに依存し、KK normalization を経由しない

したがって商 $Q_{\rm best}$ は normalization 規約に依存しない strict invariant となる。

**注意**: ここでの $h_{\rm rep}=2$ は paper03ec 以来用いている minimal two-component matter convention に対応し、lower-dimensional language では標準的な Redlich parity shift の単位を与える。 $Q_{\rm best} = 1$ から $k_q = 1$ の同定は別の識別ステップ (DPPU APS 規約の spectral ↔ CS level 辞書) を要し、本付録の範囲を超える。

**Verification**: 本節は §E.10 の $\eta(0) = 1$ と $h_{\rm rep} = 2$ から $\Delta k_{\rm matter} = 1$ を代入する代数的恒等式のみ。

---

### E.13 CS level integrality `k_q ∈ ℤ` (single-valuedness)

**Statement** (§6.1, App. D.3):

$$
k_q \in \mathbb{Z} \quad [\text{from single-valuedness of } e^{i S_{\rm odd}}]
$$

**Derivation**:
3D 境界 $Y = M^3$ 上の parity-odd 有効作用を $S_{\rm odd}^{(3D)} = k_q \cdot W[A/\omega]$ と書く。ここで $W$ は parity-odd functional であり、large frame rotation (frame bundle $\pi_3(SO(3)) = \mathbb{Z}$) の unit winding に対して $\Delta W = 1$ となるよう正規化する。 $2\pi$-normalized formulation では

$$
S_{\rm odd}^{(3D)} = k_q \cdot (2\pi \cdot W), \quad \Delta W|_{n=1} = 1
$$

unit winding shift は $\Delta S_{\rm odd}|\_{n=1} = 2\pi k_q$。path integral measure の単価性 $e^{i\Delta S\_{\rm odd}} = 1$ から

$$
e^{2\pi i k_q} = 1 \iff k_q \in \mathbb{Z}
$$

**注意**: ここで用いているホモトピー入力 $\pi_3(SO(3))=\mathbb{Z}$ は §E.9 の frame-bundle normalization と同じ frame-bundle data であり、standard 3D Chern-Simons quantization の large-rotation 版に対応する。 native $\theta_{\rm NY}$ 正規化との対応は未固定で、 $k_q \in \mathbb{Z}$ 自体は 形式的な statement。

**Verification**: 本節は path integral measure の単価性に基づく formal な議論のみで、追加の数値検証を要しない。

---

### E.14 Summary table

| Statement | 本文位置 | Verification | Status |
|---|---|---|---|
| Sol³ 構造定数 + frame rigidity | §2.3, App. A | `proofs/sol3_structure.py` | exact (SymPy) |
| Sol³ rigidity: $\partial^2_{\varepsilon,s}V=0$ at $\varepsilon=s=0$ | §4.5, §E.1 | `proofs/sol3_structure.py` | direct potential Hessian check |
| Sol³ biaxial Higgsing $m^2(A_i)$ | §2.4, §3.5 | `proofs/kk_higgsing.py` | exact (KK extractor) |
| Sol³ CS = 0 cancellation | §2.3, §3.3 | `proofs/cs_cancellation.py` | exact within profile-local KK ansatz |
| Maxwell universality $K_{\rm Mxw}$ | §2.3, §5.1 R1 | `proofs/kk_higgsing.py` | common principal coefficient |
| Off-diagonal CS rule (R2) | §5.1 | `proofs/cs_cancellation.py` | classification |
| 4-type Higgsing 1対1 | §5.1 R3 | `proofs/kk_higgsing.py` | exhaustive |
| EC slice minima: Nil³/Sol³ present, T³/S³ absent | §2.4, §3.2, §5.1 R4 | `paper04/ec_slice_minima.py` | exact homogeneous Hessian criterion |
| K_geo from AX contortion gradient | §E.6 | `proofs/eta_kinetic_from_contortion.py` | independent mass/kinetic derivation |
| reduced 1D $\eta$-benchmark: $K_{\rm geo}=M_{\rm geo}=c_{\rm geo}\rho(z)$ | §3.4 | `paper04/eta_defect_coefficients.py`, `paper04/defect_localization.py` | exact mass datum + reduced kinetic match (4 geometries), variational + numerical (12/12 LOCALIZED) |
| `C²_LC` (4 geometries) | §4.5, App. A | `proofs/weyl_scalar.py` | exact (SymPy) |
| 3-form distinction | §4.1, App. A | (代数的, スクリプト不要) | classification |
| `N_top(S³) = 6r₀²` | §4.2 | `paper04/torsional_charge.py` | geometry-derived + frame-bundle normalization |
| `η_APS(T³)=η_APS(S³)=0` | §4.2, §5.2 | `proofs/aps_zero_t3_s3.py` | exact benchmark (symmetric LC spectra) |
| `η_APS^(3D)(Nil³) = +1/2` | §4.2 | `paper04/eta_aps_nil3.py` | local LC CS + spinor trace derivation |
| Heisenberg Landau levels for Nil³ Dirac | §E.10 | `proofs/landau_levels_nil3.py` | symbolic ladder construction |
| `η_APS^(3D)(Sol-P/Sol-A)=1/0` | §4.2, §5.2 | `proofs/eta_aps_sol3.py` | global/kernel-sourced compact benchmark |
| $k_{3D}^{\rm tor\text{-}CS} = k_{4D}^{\rm NY}/2$ | §6.1 | `proofs/kk_normalization.py` | strict (CZ+Stokes) |
| $Q_{\rm best} = 1$ | §6.1 | (代数的, §E.10 + Redlich) | strict (normalization-free) |
| $k_q \in \mathbb{Z}$ | §6.1 | (formal, single-valuedness) | formal |

すべての検証スクリプトは `script/scripts/proofs/` または `script/scripts/paper04/` に配置されており、論文に同梱した `script/dppu/` ライブラリに対する依存のみで自己完結する。
