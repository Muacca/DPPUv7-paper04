## Appendix C. KK vs Matsubara coexistence tables

本付録では、same Euclidean circle に対する二つの読みを表形式で並べる。目的は、KK reading と Matsubara-style thermal reading が同じ Euclidean data に対する coexisting secondary reading であることを一望できるようにすることである。

### C.1 Basic coexistence table

| Same Euclidean quantity | KK reading | Matsubara reading | Can coexist? | Caution |
|---|---|---|---|---|
| $S^1$ length $L$ | compactification scale | inverse-temperature circle with $\beta=L$ | yes | same circle, two readings |
| Maxwell baseline | reduced kinetic baseline | thermal stiffness scale | yes | coefficient sharing does not imply identity |
| CS-active structure | reduced active channel count | thermal parity-odd channel count | yes | interpretation differs |
| $S_{\rm CZ}$ | reduced inflow action | interface free-energy descriptor | yes | not a real-time statement |
| branch label | reduced branch language | thermal branch language | yes | not a full finite-temperature derivation |

この table の要点は、`shared coefficient` と `identical theory` を区別することである。同じ Euclidean quantity が二つの descriptor language を持つことは、その二つの理論的 reading が完全一致することを意味しない。

### C.2 Representative quantities and two readings

| native quantity | KK / reduced reading | Matsubara / thermal reading | remark |
|---|---|---|---|
| $K_{\rm Mxw}$ | reduced kinetic baseline | thermal stiffness scale | same coefficient, different interpretation |
| CS direction count | reduced parity-odd channel count | thermal parity-odd channel count | counting survives, meaning differs |
| $\eta_{\rm APS}^{(3D)}$ | lower-dimensional spectral asymmetry | thermal spectral descriptor | not a transport coefficient by itself |
| $N_{\rm top}$ | reduced interface CS jump | thermal torsional-sector label | remains frame-bundle-normalized torsional charge |
| $S_{\rm CZ}$ | reduced inflow functional | thermal interface free-energy descriptor | Euclidean, not real-time |

### C.3 Interpretative distinction and scope

KK reading は lower-dimensional descriptor language を与え、Matsubara reading は thermal descriptor language を与える。両者は同じ Euclidean circle に基づくが、互いに理論的同一性を主張するものではない。

また、Matsubara reading は $M^3 \times S^1$ という元の product structure を保つため、本稿の範囲では `P=0`-protected thermal sector を与える。したがって、ここでの thermal language は full finite-temperature field theory, Matsubara 和の完全解析, real-time transport の導出までは含まない。
