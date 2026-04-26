## 5. Selection rules

本節では、§2-§4 で得られた結果を selection rule として圧縮する。ここで重要なのは、個別計算の詳細そのものよりも、「どの応答がどの幾何で allowed, protected, activated になるか」を一望できる rule set にある。selection rule は inventory の要約ではなく、4 geometry の比較を再利用可能な判定規則へ変換したものである。

### 5.1 Bulk rules

bulk layer に関する主要 rule は次のように要約できる。各 rule がどの幾何学的性質から生じるかは、まず次の対応表で整理できる。

| rule | geometric source |
|---|---|
| R1 | Maxwell action の全幾何共通性 |
| R2 | structure constant の off-diagonal 成分の有無 |
| R3 | KK vector sector の質量分裂パターン |
| R4 | scalar potential における EC slice-minimum branch |
| R5 | reduced 1D $\eta$-benchmark が開く局在化井戸と、その活性化条件 |

- **R1: Maxwell universality**  
  Maxwell baseline は 4 geometry で共通であり、

$$
{K}_{\rm Mxw}=-\frac{{L}^{2}}{{2r}_{0}^{4}}
$$
  
  が全系に対する master coefficient となる。したがって geometry 依存性は kinetic baseline ではなく、活性化される channel の側に現れる。  
  *Justify*: Maxwell の principal part は abelian quadratic term $F_{ij}F^{ij}$ に由来し、構造定数は modified field strength を通じて lower-order mass / CS 項にのみ入る。したがって共通 product-circle 規約の下では principal coefficient は 4 geometry すべてで $-L^2/(2r_0^4)$ に固定される（[Appendix E.2b](DPPUv7-paper04_appE.md); 検証 `proofs/kk_higgsing.py`）。

- **R2: CS activation rule**  
  Chern-Simons channel は、 ${\tilde{F}}\_{ij}$ に対する補正が $\{i,j\}$ の外側の脚を持つときにのみ活性化する。 $Nil^3$ の 1 方向、 $S^3$ の 3 方向はこの規則に従い、 $Sol^3$ の self-referential correction は厳密に相殺される。  
  *Justify*: 補正型を「構造定数なし / off-diagonal / self-referential」に分類することで、CS 方向数 0/1/3/0 が代数的に決まる（[Appendix E.4](DPPUv7-paper04_appE.md)）。 $Sol^3$ では $A_1\to{\tilde{F}}\_{01}$, $A_2\to{\tilde{F}}\_{02}$ の脚が field-strength index pair の内側にあるため、profile-local な任意係数 $c_{01}(z),c_{02}(z)$ に対して CS coefficient は 0 のまま保たれる（検証 `proofs/cs_cancellation.py`）。これは profile-local reduced KK ansatz の主張であり、full derivative-dependent inhomogeneous field equation の導出ではない。

- **R3: Higgsing classification**  
  KK Higgsing pattern は none / uniaxial / triaxial / biaxial の 4 type に分類される。これは $T^3$, $Nil^3$, $S^3$, $Sol^3$ の順にちょうど一つずつ実現され、bulk sector の mode dictionary を最も直接に圧縮する。  
  *Justify*: KK mass dict を 4 幾何で計算すると massive 数が $0, 1, 3, 2$ となり、4 type が過不足なく実現されることが確認できる（[Appendix E.5](DPPUv7-paper04_appE.md)）。

- **R4: EC slice-minimum rule**  
  current DPPU の homogeneous EC+NY+Weyl potential を 4 geometry で横断評価すると、 $\eta=V=0$ 断面の EC slice-minimum branch は $Nil^3$ と $Sol^3$ に存在し、 $T^3$ と round $S^3$ には存在しない。したがって $Nil^3$ は CS+EC が同時に立ち上がる spin-0 entry を担い、 $Sol^3$ は CS-inert だが EC-active な rigid entry を担う。
  *Justify*: revised paper03ec [4] では

$$
r_0=\frac{4\kappa}{\sqrt{3}}\sqrt{|\alpha|}, \qquad
|\kappa^2\theta_{\rm NY}|<1 \Rightarrow \text{full homogeneous local minimum}
$$
  
  が解析的に導かれている。 $Sol^3$ では同じ $r_0$ と同じ Hessian criterion が現れ、slice potential と $H_{RR}$ だけが $Nil^3$ の 4 倍になる。 $T^3$ は flat zero slice のため isolated radial branch を持たず、round $S^3$ は slice potential の slope が非零のため stationary point を持たない（[Appendix E.5b](DPPUv7-paper04_appE.md); 検証 `paper04/ec_slice_minima.py`）。

- **R5: reduced $\eta$-benchmark localization rule**  
  slice-wise homogeneous background 上の reduced 1D AX-torsion benchmark では、4 geometry すべてで Gaussian radial dip が局在化井戸を作る。だが local parity-odd activated response へ進むには CS-active structure が必要である。したがって $Nil^3$ は CS+EC により activated、 $S^3$ は CS+triplet/nonlinear opening により activated、 $T^3$ は kinematic / trivial、 $Sol^3$ は EC-supported だが CS-inert という分類になる。
  *Justify*: current `dppu` は exact homogeneous datum
  
$$
M_{\rm geo}=\left.\frac{\partial^2 V_{\rm eff}}{\partial\eta^2}\right|_{\eta=0}=c_{\rm geo}\rho
$$
  
  を与え、AX contortion gradient からの reduced 1D kinetic datum も
  
$$
K_{\rm geo}=c_{\rm geo}\rho
$$
  
  に一致する。ここで $c_{T^3}=c_{Nil^3}=c_{Sol^3}=96\pi^4L/\kappa^2$, $c_{S^3}=12\pi^2L/\kappa^2$ である。したがって Gaussian $\rho$-dip に対する変分評価で少なくとも 1 つの bound state が 4 geometry 全てに対して解析的に保証され、代表的 benchmark scan でも 4 geometry $\times$ 3 profile の全 12 例で LOCALIZED が確認される（[Appendix E.6](DPPUv7-paper04_appE.md); 検証 `paper04/eta_defect_coefficients.py`, `paper04/defect_localization.py`）。response の差は局在化そのものではなく、その後段に CS 活性 (R2), EC slice minimum (R4), triplet/nonlinear opening があるかどうかで決まる。 $Nil^3$ は CS+EC を、 $S^3$ は CS+triplet を、 $Sol^3$ は EC のみを持ち、 $T^3$ はいずれも持たない。

以上の bulk rule は、4 geometry が `universal baseline + geometry-dependent activation pattern` の形で比較できることを示している。とくに $Sol^3$ は CS-active family を導入しないが、EC slice minimum を持つため T³ と同じ minimal bulk class には置かれない。Sol³ は CS-inert / EC-active な rigid geometry として残る。

### 5.2 Boundary rules

boundary / cobordism layer に関する主要 rule は次の通りである。

- **B1: CZ inflow activation rule**  
  CZ inflow は CS-active geometry においてのみ activated である。実際、 $T^3$ と $Sol^3$ では torsion-CS sector が消え、 $Nil^3$ と $S^3$ でのみ nontrivial inflow が残る。  
  *Justify*: $S_{\rm CZ} \propto \int \mathcal{C}_3^{\rm tor}$ は bulk CS が活性な場合のみ非零。R2 の CS direction count 0/1/3/0 と一対一対応する。

- **B2: local spectral-core rule**  
  strict local APS spectral core は $Nil^3$ において代表的に現れ、その canonical value は

$$
\eta_{\rm APS}^{(3D)}(Nil^3)=+\frac{1}{2}
$$
  
  である。この entry は spinor spectral family に属し、torsion-CS family と混同しない。  
  *Justify*: Heisenberg モード分解は PPA で $h=0$ かつ $n=0$ spin-up branch のみが asymmetry を担うことを示す。値そのものは $Nil^3$ の Levi-Civita 接続から計算した spinor CS 積分 $\int_Y \mathcal{C}\_3^\sigma=-1/2$ と向き規約 $\eta\_{\rm APS}=-\int_Y\mathcal{C}\_3^\sigma$ から得られ、したがって $\eta\_{\rm APS}=+1/2$, $\eta(0)=1$ となる（[Appendix E.10](DPPUv7-paper04_appE.md); 検証 `paper04/eta_aps_nil3.py`）。比較として、 $T^3$ の PPA benchmark と round $S^3$ の Levi-Civita Dirac は厳密な $\lambda \leftrightarrow -\lambda$ 対称スペクトルを持つため $\eta_{\rm APS}=0$ である（検証 `proofs/aps_zero_t3_s3.py`）。一方、compact $Sol^3$ benchmark では local CS 積分は零のままでも $Sol\text{-}P / Sol\text{-}A$ に応じて $\eta_{\rm APS}=1/0$ が分岐するため、 $Sol^3$ は `local core` ではなく `global spectral branch` として扱う（[Appendix E.10b](DPPUv7-paper04_appE.md); 検証 `proofs/eta_aps_sol3.py`）。

- **B3: torsional-charge rule**
  torsional-charge entry は $S^3$ において代表的に現れ、
  
$$
N_{\rm top}=6r_0^2
$$
  
  がその canonical representative となる。pairwise cobordism dictionary で $S^3$ を含む pair が torsional-charge 主導になるのはこのためである。
  *Justify*: S³ の structure constants $C^i{}\_{jk}=(4/r_0)\epsilon^i{}\_{jk}$ から $T^i=(1/2)C^i{}\_{jk}e^j\wedge e^k$ を経て $\mathcal{C}\_3^{\rm tor}=(12/r_0)\,{\rm vol}\_{S^3}$、したがって $\int\_{S^3}\mathcal{C}\_3^{\rm tor}=24\pi^2 r_0^2$ が幾何学的に得られる。これを frame bundle 正規化 $1/(4\pi^2)$（ $\pi_{3}(SO(3))=\mathbb{Z}$, $T_R=1$）で割ることで $N\_{\rm top}=6r_0^2$ が従う（[Appendix E.9](DPPUv7-paper04_appE.md); 検証 `paper04/torsional_charge.py`）。

- **B4: spectral-source rule**  
  spectral family の内部では、observable の有無だけでなく source mechanism も区別する必要がある。 $Nil^3$ は nonzero Levi-Civita spinor CS 積分に支えられた local spectral core を持ち、 $Sol^3$ は compact quotient / spin structure に応じて kernel data から global spectral branch を持ちうる。  
  *Justify*: $Nil^3$ では $\int_Y \mathcal{C}\_3^\sigma=-1/2$ から直接 $\eta\_{\rm APS}=+1/2$ が固定される（B2）。これに対し compact mapping-torus benchmark $M_A=T^2\rtimes_A S^1$ では $Sol^3$ の local CS 積分は $\int_Y \mathcal{C}\_3^\sigma=0$ のままだが、 $k=0$ sector の constant spinor が $Sol\text{-}P$ でのみ descend するため $h=2/0$ と分岐し、反線形対称性で $\eta(0)=0$ だから $\eta\_{\rm APS}(Sol\text{-}P/Sol\text{-}A)=1/0$ を得る（[Appendix E.10b](DPPUv7-paper04_appE.md); 検証 `proofs/eta_aps_sol3.py`）。したがって Nil³ と Sol³ は同じ spectral family に属しても source は同じでない。

### 5.3 Pairwise rules

pairwise cobordism dictionary からは次のような rule が得られる。

- $S^3$ を含む pair は torsional-charge 主導になりやすい。これは $S^3$ が second layer で canonical frame-bundle-normalized torsional-charge core を担い（B3）、bulk 側でも等方的な活性化構造（CS 3 方向 + triplet degeneracy, R2/R3）を持つためである。ただし $S^3 \leftrightarrow Sol^3$ は Sol-A では torsional-charge 主導、Sol-P では torsional/global-spectral mixed pair になる。
- $Nil^3$ を含む pair は spectral 主導になりやすい。とくに $T^3 \leftrightarrow Nil^3$ では local APS core が比較の主語となり、 $Nil^3 \leftrightarrow Sol^3$ では同じ spectral family の内部で `local vs global` の source 差が比較の主語となる。これは $Nil^3$ が EC slice minimum (R4) と CS-active structure の両方を持ち、Levi-Civita Dirac の PPA spin structure で唯一 nonzero local CS core を持つからである（[Appendix E.10](DPPUv7-paper04_appE.md), [E.10b](DPPUv7-paper04_appE.md)）。なお $Sol^3$ も EC slice minimum を持つが、local CS core は持たないため source mechanism は $Nil^3$ と同一ではない。
- $S^3 \leftrightarrow Nil^3$ は torsional-charge family ( $\mathcal{C}_3^{\rm tor}$ ) と local spectral family ( $\mathcal{C}_3^{\sigma}$ ) の交差を与える代表的 mixed pair である。両者は別の 3-form sector に属するため単一の scalar observable に還元しない（3-form 区別: [Appendix E.8](DPPUv7-paper04_appE.md)）。
- $T^3 \leftrightarrow Sol^3$ は inert な最小 pair ではなく、compact quotient / spin structure に敏感な benchmark pair である。Sol-A では $T^3$ と同じ trivial spectral profile が現れる一方、Sol-P では $\eta_{\rm APS}=1$ の global spectral branch が立ち上がる。geometric scaffold の差（ $C^2_{\rm LC}$, spin-2 rigidity tag, $A_{\rm KK}/K^2$ ）はその両方で保たれる。

### 5.4 Spectral source と geometric distinctness

本稿の selection rule で特に注意すべきは、`同じ family に属すること` と `同じ source で値が決まること` を混同しないことである。 $Nil^3$ と $Sol^3$ はともに spectral family に現れうるが、 $Nil^3$ は local Levi-Civita CS 積分により、Sol³ は compact quotient / spin structure に依存する global kernel data により値が決まる。したがって same spectral family は same mechanism を意味しない。

$T^3$ は flat and inert であり、 $Sol^3$ は non-vanishing Weyl curvature, biaxial KK splitting, frame rigidity を持つ。Sol-A では boundary observable が $T^3$ の trivial benchmark に重なっても、背景幾何は一致しない。逆に Sol-P では global spectral branch が追加され、 $T^3$ と second layer の段階で分離する。この差は third observable family や fourth observable family の存在を意味するのではなく、observable を支える geometric scaffold と global spectral sourcing の差を意味する。したがって、本稿の selection rule は「observable class の圧縮」と「背景幾何の特徴づけ」を意図的に分離している。

### 5.5 本節のまとめ

selection rule を前面に出すことにより、inventory は geometry-dependent rule と結びついた比較構造として読まれる。次節では、この主構造を保ったまま Euclidean circle の二つの読みを述べる。
