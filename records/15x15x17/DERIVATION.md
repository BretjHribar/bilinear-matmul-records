# Provenance and verification of the `<15,15,17>` rank-2,256 certificate

The complete certificate was retained from the frozen search run and imported
without altering its 2,256 rank-one terms. Its embedded header records a
support-aware rank-47 outer construction with allocations

- `N = (3,4,4,4)`
- `M = (4,4,4,3)`
- `P = (4,5,4,4)`.

The search family combined exact GF(2) outer/leaf isotropies with support
clipping: basis changes preserve each source tensor but alter which local
coordinates fall outside smaller boundary blocks. Candidates were materialized,
zero terms removed, identical triples parity-cancelled, and only then passed to
an exact Brent-equation gate. This is the same general mechanism used for the
19 by 19 construction, applied to a rectangular allocation.

The imported file is not accepted merely on its header or claimed rank. Both
repository verifiers independently establish:

- 2,256 parsed nonzero, in-range, distinct triples;
- 6,564,657 sparse coefficient toggles;
- 3,825 target coefficients;
- zero residual coefficients.

This proves the stated GF(2) decomposition independently of the unavailable
original search program.
