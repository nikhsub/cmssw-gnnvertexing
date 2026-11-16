# Optional two-track vertex fits on `simplelabel`

`TrackPairFeaturesProducer` can optionally fit every row that passes its
existing pair selection with the unconstrained two-track Kalman vertex fitter
from `TwoTrackVertexFitter.h`. The fit columns are added to the existing
`trackPair` table in row order. Track-pair selection, pair endpoint ordering,
per-track `hadidx`/`label`, and `edgeLabel` production are unchanged. Fits are
not selected or filtered using truth labels.

This branch stores directed `(i,j)` and `(j,i)` rows for the same physical
track pair because some existing features depend on endpoint order. The fit is
symmetric, so it is computed once per unordered pair and copied to both rows.
Accordingly, `vertexFitSummary.nFitAttempts` counts unique physical pairs,
while `nSelectedPairs` counts the existing directed rows.

The switch is off by default. This keeps existing reconstruction and training
output unchanged. Enable `enableTwoTrackVertexFit=True` on the `track`
collection in `main_cfg.py` to write fit results. `enableVertexFitTiming=True`
also enables event-level counters and timing in `vertexFitSummary`.

When enabled, each selected pair receives:

- `vtxFitStatus`: 1 success, 2 invalid track state, 3 invalid fitter result,
  4 exception or invalid fit numerics; `vtxFitValid` is 1 only for status 1.
- `vtxCovarianceValid` and the fit quality/position/covariance fields:
  `vtxChi2`, `vtxNdof`, `vtxChi2Ndof`, `vtxProb`, `vtxX`, `vtxY`, `vtxZ`,
  `vtxCovXX`, `vtxCovXY`, `vtxCovXZ`, `vtxCovYY`, `vtxCovYZ`, `vtxCovZZ`,
  `vtxSigmaX`, `vtxSigmaY`, `vtxSigmaZ`.
- PV-relative fields and their individual masks:
  `vtxPVValid`, `vtxLxyValid`, `vtxLxy`, `vtxLxyErr`, `vtxLxySig`,
  `vtxSignedLxyValid`, `vtxSignedLxy`, `vtxSignedLxySig`, `vtxL3DValid`,
  `vtxL3D`, `vtxL3DErr`, `vtxL3DSig`, `vtxPointing2DValid`,
  `vtxCosPointing2D`, `vtxPointing3DValid`, `vtxCosPointing3D`,
  `vtxPerpendicularValid`, `vtxPerpDisplacement`,
  `vtxPerpDisplacementSig`.

Position and displacement units are cm; covariance units are cm². Invalid
floating-point outputs use `-9999`. `vtxFitValid` and `vtxCovarianceValid`
are distinct: a fit can have a valid position and fit quality while its
covariance fails the positive-semidefinite check. Covariance-dependent
uncertainties and significances are masked independently. When timing is
enabled, `vertexFitSummary` records track/pair counts, fit attempts, successes,
failures, covariance-invalid fits, transient-track construction time, fit
time, and producer time. Timing values are `-1` when timing is disabled.

Example from the CMSSW source area:

```bash
scram b -j 8
cmsRun main_cfg.py collection=track isMC=True \
  enableTwoTrackVertexFit=True enableVertexFitTiming=True \
  outputFile=track_pairs_with_fits.root
```

`main_cfg.py` currently takes its input file and event limit from the values
in the configuration itself; replace those values there when using a different
sample or limiting the run.

The existing `main_cfg.py` keep rule for `trackPairFeatures` saves the added
fit columns and `vertexFitSummary`; the MC `trackPairLabels` extension remains
saved under its existing keep rule. The implementation has not been compiled
or run in this environment.
