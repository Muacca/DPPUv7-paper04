## Appendix D. Computational notes and normalization chain

本付録では、本文で用いる normalization chain と descriptor separation をまとめる。主辞書と secondary reading を混同しないための最小限の補助事項だけを記す。

### D.1 Geometric response と descriptor

本稿では、geometric response そのものと、それを別言語で表した descriptor とを区別する。本文の primary object は常に geometric response dictionary であり、KK reading と Matsubara-style thermal reading はそれに付加される secondary interpretation である。

- **geometric response**: 本文の primary object
- **reduced descriptor**: lower-dimensional language による記述
- **thermal descriptor**: thermal language による記述

この区別により、response dictionary 本体と読み替え言語を混同しない。たとえば $K_{\rm Mxw}$ は native geometric quantity であると同時に reduced kinetic baseline と thermal stiffness scale の両方を許すが、そのことは三つの quantity が相互に同一であることを意味しない。同様に、 $\eta_{\rm APS}^{(3D)}$ と $N_{\rm top}$ はどちらも lower-dimensional / thermal descriptor を持ちうるが、observable family としては別物である。

### D.2 Normalization chain

CZ boundary term の normalization patch は概念的には

$$
\mathrm{NY} = d\mathcal{C}_3^{\rm tor}
\quad \Longrightarrow \quad
\int_{X_4}\mathrm{NY} = \int_{M_3}\mathcal{C}_3^{\rm tor}
\quad \Longrightarrow \quad
N_{\rm top}:=\frac{1}{4\pi^2}\int_{M_3}\mathcal{C}_3^{\rm tor}
\quad \Longrightarrow \quad
S_{\rm CZ}[M_3] = \frac{k_{\rm cl}}{4\pi^2}\int_{M_3}\mathcal{C}_3^{\rm tor}
$$

という chain で読む。これにより、torsional boundary term, frame-bundle-normalized torsional charge, inflow action が同じ $\mathcal{C}\_3^{\rm tor}$ sector の異なる表現であることが明確になる。 $1/(4\pi^2)$ 正規化が $T_R = 1$ (frame bundle Chern class 単位) に対応すること、および $S^3$ で $N\_{\rm top} = 6r_0^2$ となることの導出は [Appendix E.9](DPPUv7-paper04_appE.md) を参照。

一方、APS spectral entry はこの chain には含まれない。 $\eta_{\rm APS}^{(3D)}$ は $\mathcal{C}\_3^\sigma$ を介する spectral family に属し、 $\mathcal{C}\_3^{\rm tor}$ から直接定義される $N\_{\rm top}$ や $S\_{\rm CZ}$ とは別の observable として扱う必要がある。本稿で 3 family taxonomy を採用した理由は、この normalization chain を通じてもなお family distinction が消えないからである。

### D.3 Reduced-side consistency statements

KK reduction reading を用いるとき、lower-dimensional parity-odd language には次の 3 つの consistency statement が現れる。

$$
k_{3D}^{\rm tor\text{-}CS}=\frac{1}{2}k_{4D}^{\rm NY},
\qquad
Q_{\rm best}:=\frac{\eta(0)}{\Delta k_{\rm matter}}=1,
\qquad
k_q \in \mathbb{Z}.
$$

第 1 の式は 4D Nieh-Yan coefficient と 3D torsion-Chern-Simons coefficient の対応を与える (KK reduction normalization identity; [Appendix E.11](DPPUv7-paper04_appE.md); 検証 `proofs/kk_normalization.py`)。第 2 の式は spectral quantity と matter shift を用いた normalization-free quotient を与える (matter-minimal inflow audit; [Appendix E.12](DPPUv7-paper04_appE.md))。第 3 の式は reduced Chern-Simons level に対する integrality condition である (formal quantization on the frame bundle; [Appendix E.13](DPPUv7-paper04_appE.md))。

本稿では、これらを boundary observable family に追加される新しい entry としては扱わない。むしろ、すでに本文で与えた Euclidean geometric response dictionary を lower-dimensional language で読む際の consistency statement と位置づける。また、これら reduced-side statement と native geometric quantity の完全な同一視は、本稿の主張範囲には含めない。
