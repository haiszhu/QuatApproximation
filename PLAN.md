# QuatApproximation -- Roadmap

Tracker for open items.  Closed items move to git history (commit log /
PR descriptions).

## Open

### Local-to-local translation operator
**Not yet implemented.**  A translation operator for solid spherical
harmonic local expansions: given coefficients `{a_{l,m}}` centered at
`c0`, produce coefficients `{b_{l,m}}` centered at `c1` representing the
same harmonic function in the new local frame.

Needed for:
- Refinement-based close-evaluation strategies where each sub-patch
  carries its own local harmonic expansion and the parent's expansion
  has to be transferred to children (and vice versa).
- A future r128 pipeline that mirrors the r64 close-eval chain at every
  step including local-expansion transfers.

Likely landing site: a new `localtransl_mod.f90` (prefix `qal_`) alongside
`harmonic_mod`, with paired r64 / r128 routines and the usual
`<pp><m>_mex.f90` wrappers + `@function` blocks in
`matlab/QuatApproximation.mw`.

### Unify triangle and quadrilateral patch routines
The quaternion approximation on **triangle** patches (Koornwinder /
Vioreanu-Rokhlin reference element) and on **quadrilateral** patches
(tensor-product Gauss-Legendre reference element) currently live in two
parallel module families that were written at different times and have
drifted in naming, API, and feature coverage:

| Aspect | Triangle | Quadrilateral |
|---|---|---|
| reference-element module | `koorn_geom_mod` (prefix `qakg_`) | `tensor_geom_mod` (prefix `qatg_`) |
| vals -> coeffs primitive  | `koorn_vals2coefs_coefs2vals` (Koornwinder) | inline `Legmat = (2k-1)/2 * w * P_{k-1}(t)` |
| boundary GL quadrature    | `line3quadr_3dline`                  | `line3quadr_3dline_T`                      |
| boundary chart            | private `tparam_to_uv` (triangle 0..2π → simplex) | private `tparam_to_uv_square` (square 0..2π → box) |
| r128 sibling              | partial (`circumcircle_transform_3d_r128`, `lqkg_setup_target_r128`) | `line3quadr_3dline_T_r128` |
| associated mex prefix     | `qakg_*_mex`                         | `qatg_*_mex`                                |

This split is a historical artifact and needs to be unified so:
- routine names follow one consistent convention across patch types
- API shape (argument order, output-buffer ownership, slot conventions)
  matches between the triangle and quadrilateral families
- feature coverage is symmetric (every routine that exists on one side
  has a counterpart on the other, at both r64 and r128)
- downstream callers can switch patch type by swapping a single
  family-prefix rather than rewriting call sites

Sketch of the unification: pick one canonical naming scheme (e.g.
`patch_<routine>_<tri|quad>_<r64|r128>` or a shared module with a
`patch_type_id` dispatch arg), then either rename the two families to
match it or fold both into a single `patch_geom_mod` whose internals
dispatch on the chart type.  No timeline; do this before the package
grows a third patch type or before downstream BIE-solver code is built
on top of the current fragmented API.

### Other follow-ups
- `qak_q*kl` and `qak_qnm_i`: `lptype = 'T'` (DLPn) branch is currently
  `error stop`-stubbed in both r64 and r128.  Port `q[0-3]kl_DLPn` from
  `qotential/private/utils/q[0-3]kl_DLPn.m` when a downstream caller
  needs it.
- `evaltensorproductharmonicgrad_r128` (harmonic_mod) — needed once a
  fully-r128 close-eval orchestration goes live.
