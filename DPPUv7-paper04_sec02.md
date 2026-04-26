## 2. Homogeneous baseline

本節では、4 つの Thurston 幾何 $T^3$, $Nil^3$, $S^3$, $Sol^3$ の homogeneous baseline を整理し、比較軸を固定する。ここでの目的は、bulk / boundary の response dictionary を議論する前に、各 geometry の mode language, baseline invariant, minimal / activated distinction を同一規約で並べることである。後続節の selection rule は、すべて本節で固定する baseline の比較から読み取られる。

### 2.1 Four-geometry setup

本稿で比較する 4 geometry はすべて $M^3 \times S^1$ という共通 Euclidean product structure を持つが、3 次元空間部の構造定数、Weyl curvature, KK splitting pattern, frame rigidity は geometry ごとに異なる。したがって、同一の minisuperspace 変数を用いても、effective potential の安定性、ベクトル sector の質量分裂、Chern-Simons channel の活性化条件は幾何依存に変わる。

本稿では、まず homogeneous sector を comparison baseline として固定し、その上に first-layer bulk response と second-layer boundary response を重ねる。これにより、後続節で現れる `activated`, `protected`, `minimal` という語が、単なる記述上のラベルではなく、homogeneous background に対する明確な比較判定として読めるようになる。

### 2.2 基本比較表

現時点の baseline は次の表で要約される。

| 指標 | $T^3$ | $Nil^3$ | $S^3$ | $Sol^3$ |
|---|---|---|---|---|
| CS 方向数 | 0 | 1 | 3 | 0 |
| KK Higgsing | none | uniaxial | triaxial | biaxial |
| EC slice minimum | absent | present | absent | present |
| Maxwell baseline | universal | universal | universal | universal |
| spin-2 rigidity | non-rigid | non-rigid | non-rigid | rigid |
| $C^2_{\rm LC}$ | 0 | $4/(3R^4)$ | 0 | $16/(3R^4)$ |

この表は二つの事実を同時に示している。第一に、 $Nil^3$ と $S^3$ は CS direction count と Higgsing pattern の両方で非自明な activated geometry である。第二に、 $Sol^3$ は CS direction count では minimal 側に留まる一方、non-vanishing Weyl curvature によって $Nil^3$ と同型の EC slice-minimum branch を持つ。したがって $Sol^3$ は $T^3$ と同じ inert class ではなく、CS-inert だが spin-0 sector では active な rigid geometry として扱う。boundary layer では、この bulk distinction に加えて compact quotient / spin structure に依存した spectral branch の有無がさらに効く。

### 2.3 Maxwell universality と CS direction count law

bulk layer の最初の universal statement は Maxwell coefficient の geometry-independence である。本稿では

$$
K_{\rm Mxw} = -\frac{L^2}{2r_0^4}
$$

を 4 geometry に共通する baseline coefficient として採用する。この普遍性は、後続の inhomogeneous opening や KK / thermal reading においても主辞書の骨格を保つ。したがって、geometry 依存性は Maxwell sector そのものではなく、CS-active direction, Higgsing pattern, EC slice-minimum structure, defect response の側に現れる。

これに対して Chern-Simons activation は geometry-dependent であり、off-diagonal structure により方向数が決まる。 $Nil^3$ では 1 方向、 $S^3$ では 3 方向が活性化し、 $Sol^3$ では self-referential cancellation によりゼロになる。この違いは、単に「構造定数があるかどうか」ではなく、「補正が $\tilde{F}_{ij}$ に対して $\{i,j\}$ の外側の脚を持つかどうか」で判定される。

$Sol^3$ で cancellation が起こる理由は、基礎となる構造定数が

$$
C^1{}_{01}=+\frac{1}{R}, \qquad C^2{}_{02}=-\frac{1}{R}
$$

という対称な組をなしており、同脚自己結合型の補正だけが残るためである。具体的には $A_1 \to \tilde{F}_{01}$ の補正の脚 1 が $\{0,1\}$ に含まれ self-referential となるため、 $+1/R$ と $-1/R$ の寄与が代数的に厳密に相殺する（off-diagonal CS rule, [§5.1 R2](DPPUv7-paper04_sec05.md); 詳細は [Appendix E.3](DPPUv7-paper04_appE.md)）。このとき CS channel は新しい parity-odd direction を生まず、Maxwell baseline の上に biaxial Higgsing だけが残る。

### 2.4 Mode dictionary と Higgsing pattern

homogeneous sector を mode language で読むと、4 geometry の差は三つの sector に集約される。spin-0 sector では、 $\eta=V=0$ 断面の EC slice minimum が $Nil^3$ と $Sol^3$ に現れる。両者は同じ stationary radius と同じ full homogeneous Hessian criterion を持つが、 $Sol^3$ の slice potential と radial Hessian は $Nil^3$ の 4 倍になる。 $T^3$ は flat zero slice で isolated branch を持たず、 $S^3$ は round slice の slope が非零で stationary point を持たない（検証は [Appendix E.5b](DPPUv7-paper04_appE.md)）。spin-1 sector では Higgsing pattern が none / uniaxial / triaxial / biaxial に分かれ、これが bulk dictionary の主要比較軸になる。spin-2 sector では $Sol^3$ のみが frame rigidity を持ち、他の 3 geometry は non-rigid 側に属する。

とくに $Sol^3$ の biaxial Higgsing は

$$
m^2(A_1)=m^2(A_2)=-\frac{L^2}{2r_0^4}, \qquad m^2(A_0)=0
$$

として書ける。これは modified field strength $\tilde{F}\_{01} = F\_{01} + A_1/R$ から $A_1, A_2$ に対する mass 項が $K_{\rm Mxw}$ と同じ係数で誘起され（同脚自己結合）、 $A_0$ は $\tilde{F}$ に現れず massless であることから従う（KK 導出は [Appendix E.2](DPPUv7-paper04_appE.md)）。ここで二つの massive direction が Maxwell baseline と同じ係数を共有することは、 $Sol^3$ が CS-active direction を導入しない一方で、EC slice branch と rigid scaffold を持つ distinct な geometry であることを示している。

なお Higgsing の 4 type (`none / uniaxial / triaxial / biaxial`) は KK mass 行列の直接計算から $T^3, Nil^3, S^3, Sol^3$ にちょうど一つずつ過不足なく実現されることが確認できる（[Appendix E.5](DPPUv7-paper04_appE.md)）。

以上により、本節の役割は response dictionary の比較軸を固定する baseline を与えることにある。technical detail の大半、特に $Sol^3$ の構造定数、CS cancellation, frame rigidity の詳細は [付録 A](DPPUv7-paper04_appA.md) と [付録 E](DPPUv7-paper04_appE.md) にまとめる。
