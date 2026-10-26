# Dual futile cycle MATLAB verification

MATLAB scripts accompanying the paper:  

“Sequential and distributive dual futile cycle: Hopf bifurcation can occur under parameter-rich kinetics but cannot occur under mass-action kinetics.”  

See the paper’s Implementation section for details.

All scripts were tested with MATLAB R2026a on macOS using an Apple M2 processor. 

## Optional exact Python verification

[Eight-parameter symmetry slice](verification/symmetry_slice/README.md): an independently written exact spectral classification and optional Python checkers. The unrestricted sixteen-parameter Conjecture 3.1 remains unresolved; priority of this symmetry corollary is unconfirmed.

[Full-parameter fourth-Hurwitz partial lemma](verification/README.md): an optional exact certificate proves `det(G)>=0 => Delta4>0` for sixteen independent positive parameters. This determinant-conditioned extension does not globally replace the source's mass-action Proposition5.6. Conjecture3.1 remains unresolved; no positivity proof for `Delta5` is supplied.
