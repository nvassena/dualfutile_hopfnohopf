# Optional exact Python companions

These independently written Python companions supplement the original MATLAB checks. The original scripts and the source paper retain their authorship and licenses.

- [Eight-parameter symmetry slice](symmetry_slice/README.md): exact spectral classification under eight specified parameter pairings.
- [Full-parameter fourth Hurwitz determinant](full-parameter-delta4/README.md): for sixteen independent strictly positive parameters, `det(G)>=0` implies `Delta4>0`. The full Conjecture3.1 remains unresolved; `Delta5>0` is not established.

The second result is determinant-conditioned. It enlarges the kinetic class on that region and does not globally supersede the source's mass-action Proposition5.6, which has no determinant-sign premise. Both companions retain the source's mass-action theorem and negative-core observation as prior work.

Run the new portable bundle from this repository root:

```sh
python -m pip install -r verification/full-parameter-delta4/requirements.txt
python -B verification/full-parameter-delta4/verify_certificate.py
python -B -O verification/full-parameter-delta4/verify_certificate.py
python -B verification/full-parameter-delta4/independent_reconstruction.py
python -B -O verification/full-parameter-delta4/independent_reconstruction.py
```

The portable bundle is byte-identical to the independently checked lab package; its embedded README gives location-independent commands to run from inside that directory. All checks print JSON to stdout and leave package inputs and archived evidence unchanged. A separate substantive AI review checked the partial lemma; external human review and priority confirmation remain pending. The new directory's MIT license applies only to those new files.
