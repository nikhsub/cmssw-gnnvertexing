# Optional two-track vertex-fit features

## Design

`TrackGVFeaturesProducer` already builds one `reco::TransientTrack` per input
track, enumerates the saved `trackPair` rows, and computes the existing
selection/features. The fit is therefore integrated into that same selected
pair loop. This reuses the transient tracks and pair selection and avoids a
second O(Ntracks^2) enumeration. A separate producer would have to rebuild
transient tracks and consume/reconstruct the selected-pair identity.

The fit logic is isolated in the reusable
`PhysicsTools/SimpleVertexTable/interface/TwoTrackVertexFitter.h` helper. It
uses one stream-local `KalmanVertexFitter(false)` (unconstrained, no smoothing)
per producer instance and one fit per saved pair. The two-element input vector
is reserved once per event and reused. Fits are never gated on `edgeLabel` or
other generator truth. `KalmanVertexFitter::vertex` returns a
`CachingVertex<5>`; the helper reads its position, error, chi-squared, and
degrees of freedom directly and converts only the covariance matrix. It avoids
converting each result through `TransientVertex`, which would rebuild per-fit
track containers.

The pair table remains named `trackPair`. When enabled, fit columns are added
to that same table in precisely its row order. The row count, endpoints,
selection, legacy labels, and other existing pair columns are unchanged. Pairs
for which the existing closest-approach calculation fails are still retained
by the original selection (`dcaVal` stays -1 and therefore does not fail the
`dcaVal > 0.1` rejection); they receive a fit attempt and an explicit status.

## Fit and status semantics

The fit uses the two existing transient tracks with an unconstrained Kalman
vertex fit. No mass constraint, PV constraint, or track refit is requested.

`vtxFitStatus`: 1 = valid fit; 2 = invalid transient-track state; 3 = invalid
fitter result; 4 = exception or invalid fit numerics. `vtxFitValid` is 1 only
for status 1. Floating-point values unavailable for a row are `-9999`. Since
the current pair loop itself skips invalid transient tracks, status 2 is a
defensive code path and is expected to be rare/unreachable for selected pairs.
Fit failures do not remove pair rows.

Position is in cm; covariance in cm^2; uncertainties and displacement values
in cm. Fit outputs are `vtxChi2`, `vtxNdof`, `vtxChi2Ndof`, `vtxProb`, `vtxX`,
`vtxY`, `vtxZ`, `vtxCovXX`, `vtxCovXY`, `vtxCovXZ`, `vtxCovYY`, `vtxCovYZ`,
`vtxCovZZ`, `vtxSigmaX`, `vtxSigmaY`, and `vtxSigmaZ`.
The columns are stored as full `float` values without NanoAOD mantissa
reduction; covariance values are calculated in double precision before their
conversion to the flat-table float representation.

`vtxPVValid` indicates a nonempty, nonfake leading reconstructed PV. Existing
workflow semantics use the first vertex (`offlineSlimmedPrimaryVertices[0]`)
as PV. No PV refit is done. PV-dependent values are invalid when that PV is not
usable; raw fit values remain available. PV/SV covariance is added assuming
independent estimates, and their correlations are neglected. The PV may
include tracks from the candidate.

The displacement features are:

* `vtxLxy`, `vtxLxyErr`, `vtxLxySig`: transverse PV-to-SV length and uncertainty
  propagated along the radial transverse direction.
* `vtxSignedLxy`, `vtxSignedLxySig`: displacement projected along the sum of
  the two reconstructed track transverse momenta. Its direction uncertainty is
  treated as fixed.
* `vtxL3D`, `vtxL3DErr`, `vtxL3DSig`: 3D PV-to-SV length and radial uncertainty.
* `vtxCosPointing2D`, `vtxCosPointing3D`: cosine between displacement and
  reconstructed pair momentum. Values are clamped to [-1, 1] for roundoff.
* `vtxPerpDisplacement`, `vtxPerpDisplacementSig`: transverse displacement
  perpendicular to pair pT and its significance. The significance propagates
  position covariance and treats pair-momentum direction as fixed.

Each family has a corresponding `vtx...Valid` integer flag. A zero-length
flight vector, zero pair momentum, or nonpositive/nonfinite projected variance
invalidates only the affected derived quantity, not the fit.

## Timing and counters

Enabling the fit also emits one singleton `vertexFitSummary` table row per
event with track count, selected pair count, attempted/successful/failed fit
counts, and timing fields. `enableVertexFitTiming=True` enables steady-clock
measurements for transient-track construction, cumulative fitter calls, and
the producer's event-local wall time. Timing fields are `-1` when timing is
disabled. Counters have negligible work compared with the fits; no per-pair
logging is added.

The NanoAOD keep rule already used by the QFrag configs,
`keep nanoaodFlatTable_trackGVFeatures_*_*`, retains both the extended
`trackPair` columns and the optional summary table.

## Enable and run

Both existing configs expose the same options. Set the same `inputFiles` and
`maxEvents` for baseline and enabled comparisons:

```bash
cd CMSSW_15_1_1/src
cmsenv
scram b -j 8

# Small enabled smoke test (10 events)
cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_training_cfg.py \
  inputFiles=file:input.root maxEvents=10 enableTwoTrackVertexFit=True \
  enableVertexFitTiming=True outfile=vertex_fit_smoke.root

# Baseline run (same input, event count, and config)
/usr/bin/time -v -o baseline.time \
  cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_training_cfg.py \
  inputFiles=file:input.root maxEvents=100 enableTwoTrackVertexFit=False \
  enableVertexFitTiming=False outfile=baseline.root

# Fit-enabled comparison run
/usr/bin/time -v -o enabled.time \
  cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_training_cfg.py \
  inputFiles=file:input.root maxEvents=100 enableTwoTrackVertexFit=True \
  enableVertexFitTiming=True outfile=enabled.root
```

For a complete baseline timing record, run the disabled configuration under
`/usr/bin/time -v -o baseline.time` with the same event count and input file.
Use the matching full-table config instead of the compact training config if
that is the intended production workflow. The fit arguments are identical.

## ROOT validation and benchmark summary

From a CMSSW environment with PyROOT:

```bash
python3 PhysicsTools/SimpleVertexTable/test/validate_two_track_vertex_fit.py \
  --baseline baseline.root --enabled enabled.root \
  --baseline-time baseline.time --enabled-time enabled.time \
  --csv vertex_fit_check.csv
```

The script checks presence of fit branches, pair row counts and exact endpoint
identity against the baseline, fit/status consistency, finite valid values,
probability range, covariance principal minors, and per-event summary counts.
It writes a compact CSV with pair/fit counts and fit timing totals even when
consistency checks fail, and prints the failing check categories. The
external `/usr/bin/time -v` records provide total CPU/wall time and peak RSS;
the script includes these values in the CSV and computes incremental CPU per
attempted fit as `(enabled user+system seconds - baseline user+system seconds) /
enabled fit attempts`.

## Unverified items

This cloud environment cannot compile or execute CMSSW. The Kalman fitter,
table schema, thread behavior, and calculations have not been built or run
here. In particular, validate the exact CMSSW_15_1_1 API/build dependency,
the ROOT branch names, fit success rate, and numerical distributions locally.
