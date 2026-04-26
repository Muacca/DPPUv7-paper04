## 6. Secondary readings of the Euclidean circle

本節では、Euclidean circle $S^1$ に対する二つの読みを整理する。これらは、既に得られた Euclidean response dictionary に対して reduced description と thermal description を与える。重要なのは、ここで述べる読みが主辞書を置き換えるのではなく、同じ Euclidean data に対する secondary interpretation を与える点である。

### 6.1 KK reduction reading

KK reading では、geometric response そのものと、それを lower-dimensional language で記述した reduced descriptor とを区別する。すなわち、lower-dimensional sector に現れる有効係数や parity-odd channel は、本稿で先に構成した dictionary の像であって、primary object そのものではない。

KK reduction normalization identity で用いる幾何設定を明示しておく。係数抽出では、Euclidean circle を長さ $\beta$ の product/collar direction として扱い、相対 cobordism

$$
X^4=M^3\times I_\beta,\qquad I_\beta=[0,\beta]
$$

上で CZ 恒等式 ${\rm NY}=d\mathcal{C}\_3^{\rm tor}$ と Stokes 定理を適用する。3D effective coefficient は $I\_\beta$ 方向の体積因子を吸収した単一の physical boundary slice $Y=M^3$ 上の functional として定義するため、ここでの $\int\_{\partial X^4}\mathcal{C}\_3^{\rm tor}$ はこの single-copy boundary contribution を表す。したがって、以下の $k\_{3D}^{\rm tor\text{-}CS}$ は closed thermal trace 上の二端点差ではなく、inflow/collar 正規化で得られる 3D boundary coefficient である。また nonzero KK tower の相殺を述べる箇所では、untwisted product circle 上の実 bosonic KK 分解を仮定する。

この reduced-side reading では、次の 3 つの consistency statement が有用である。

$$
k_{3D}^{\rm tor\text{-}CS}=\frac{1}{2}k_{4D}^{\rm NY},
\qquad
Q_{\rm best}:=\frac{\eta(0)}{\Delta k_{\rm matter}}=1,
\qquad
k_q \in \mathbb{Z}.
$$

前二者は torsional / matter-coupled parity-odd descriptor に対する normalization-free な statement を与え、最後の式は reduced Chern-Simons level に対する integrality condition を与える。第一式 $k_{3D}^{\rm tor\text{-}CS} = (1/2) k_{4D}^{\rm NY}$ は、4D Nieh-Yan action $S_{\rm NY} = (k_{4D}^{\rm NY}/8\pi^2)\int d\mathcal{C}\_3^{\rm tor}$ を $S^1$ に沿って KK 還元したときの 3D torsion-CS level として代数的に得られる（KK reduction normalization identity; [Appendix E.11](DPPUv7-paper04_appE.md); 検証 `proofs/kk_normalization.py`）。このとき untwisted product circle 上の実 bosonic KK 分解では Fourier モード $n$ と $-n$ は複素共役な対で、even weight は $n^2$ のみで決まるため、parity-odd residue は一般に $n\,w(n^2)$ の形になる。したがって nonzero tower は pairwise に相殺し、追加の parity-odd shift を残さない。第二式 $Q\_{\rm best}=1$ は spinor-CS と APS spectral entry の比較から、 $Nil^3$ で $\eta(0)=1$, $\Delta k_{\rm matter}=1$ となることに由来する（matter-minimal inflow audit; [Appendix E.12](DPPUv7-paper04_appE.md)）。第三式 $k_q \in \mathbb{Z}$ は frame bundle 上の large gauge transformation $\pi_3(SO(3))=\mathbb{Z}$ から要請される formal quantization condition である（formal quantization on frame bundle; [Appendix E.13](DPPUv7-paper04_appE.md)）。重要なのは、これらが second-layer の新しい observable family を定義するのではなく、Euclidean response を lower-dimensional language で読むときの consistency statement だという点である。

bulk 側では、universal Maxwell coefficient

$$
K_{\rm Mxw}=-\frac{L^2}{2r_0^4}
$$

が topology-independent な kinetic baseline として読み直され、activated Chern-Simons directions, $Nil^3$ の uniaxial Higgsing, $S^3$ の triplet degeneracy は reduced vector sector の topology-dependent channel として解釈される。したがって、reduced language は bulk inventory を lower-dimensional に読み替えた記述であり、bulk inventory そのものに先立つわけではない。

boundary 側でも同様に、torsional inflow action, APS spectral entry, frame-bundle-normalized torsional charge は lower-dimensional parity-odd language に写像できる。ただし、 $S_{\rm CZ}$ は $\mathcal{C}\_3^{\rm tor}$ を用いて書かれ、 $\eta\_{\rm APS}^{(3D)}$ は Levi-Civita spectral data に属するため、reduced language に落とした後も observable family の区別は保たれる。とくに $Nil^3$ の $\eta_{\rm APS}^{(3D)}=1/2$ と $S^3$ の $N_{\rm top}=6r_0^2$ は、どちらも reduced parity-odd descriptor を持つが、同一の quantity ではない。また、本稿では reduced-side statement と native geometric quantity の完全な同一視までは主張しない。

以上の意味で、KK language は response dictionary そのものを置き換えるのではなく、その reduced reading を与える。係数対応や normalization chain の補足は付録 D にまとめる。

### 6.2 Matsubara-style thermal reading

同じ Euclidean circle に対して thermal descriptor reading を重ねるため、本稿では

$$
\beta = L
$$

を採用し、Matsubara-style thermal reading を thermal interpretation として記述する。ここで $L$ は、すでに geometric formulation に現れている同じ Euclidean circle の周長である。したがって thermal reading は新しい circle を導入するのではなく、既存の Euclidean data に thermal meaning を与える。

この解釈のもとでも、Maxwell universality, CS direction-count law, $Nil^3$ の uniaxial Higgsing, $S^3$ の triplet degeneracy は Euclidean structural entry のままであり、それらに thermal descriptor language が追加されるだけである。また、Matsubara reading は元の product structure $M^3 \times S^1$ を壊さないため、それ自体が新しい parity-odd source を生成するわけではない。この意味で本稿が得るのは `P=0`-protected thermal sector であり、Lorentzian continuation そのものではない。

KK reading と Matsubara-style thermal reading は、同じ Euclidean circle に対する二つの secondary reading として共存する。ただし、両者が同じ circle に基づくからといって理論的同一性を主張するものではない。前者は lower-dimensional descriptor language、後者は thermal descriptor language であり、reduced gauge coefficient, boundary inflow action, phase label などが両方の読みを許すにすぎない。

本節の thermal reading は意図的に Euclidean scope に限定されている。したがって、full finite-temperature field theory の導出、partition function の厳密構成、Matsubara 和の詳細解析、real-time transport までは主張しない。KK と Matsubara の coexistence table は付録 C に、descriptor separation と normalization chain は付録 D にまとめる。
