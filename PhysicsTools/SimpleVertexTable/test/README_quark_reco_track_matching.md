# Stage 2: QGenCharged to reconstructed-track matching

`QuarkRecoTrackMatchProducer` writes one `QRecoTrack` row per row in
`unpackedTracksAndVertices`, the same `reco::TrackCollection` used by the
package's graph-vertex workflow. The configuration copies the unpacker inputs
and `minTrackPt=0.4` from `python/sv_reco_cff.py`. The producer consumes the
Stage 1 `QGenChargedTable` directly, so its `genChargedIdx` is the event-local
row index in that table. It does not repeat or modify Stage 1 ancestry logic.

The Stage 2 match is a kinematic reco–gen association, not a simulated-hit
truth match and not a measured tracking efficiency. It uses only charge,
kinematics, and the configured thresholds; ancestry and charm labels do not
enter the score. Defaults are `maxDeltaR=0.03` and
`maxRelPt=abs(recoPt-genPt)/genPt <= 0.5`. Charge must be nonzero and match
exactly, and both sides must have finite kinematics with positive pT. No
species filter is applied to `QGenCharged`.

Eligible pairs have cost
`(deltaR/maxDeltaR)^2 + (relativePt/maxRelPt)^2`, quantized to a 1e-12 integer
scale. A sparse min-cost max-flow solver augments until no path remains: this
maximizes the number of one-to-one matches, then minimizes total cost. Tracks
and generator rows are inserted in ascending event-local index; shortest-path
ties use stable node/edge ordering. This makes exact ties reproducible.

## QRecoTrack columns

| Column | Type | Meaning |
| --- | --- | --- |
| `trackIdx` | int | Position in the input `reco::TrackCollection`. |
| `pt`, `eta`, `phi` | float | Reconstructed-track kinematics; pT in GeV and phi in radians. |
| `charge` | int | Reconstructed charge in units of e. |
| `isHighPurity` | int | `reco::Track` high-purity quality bit, 0 or 1; not a matching requirement. |
| `genChargedIdx` | int | Event-local QGenCharged row, or `-1` when unmatched. |
| `matchDR` | float | Matched deltaR, or `-1` when unmatched. |
| `matchRelPt` | float | Matched absolute relative pT difference, or `-1` when unmatched. |
| `matchStatus` | int | Status enum below. |
| `nGenCandidates` | int | Number of admissible QGenCharged rows before assignment. |

`matchStatus` is 0 = no admissible candidate; 1 = matched with no competing
eligible pair; 2 = matched while the track had another eligible generator or
the selected generator was eligible for another track; 3 = had at least one
eligible candidate but remained unmatched after the global one-to-one
assignment. The unmatched sentinels are `genChargedIdx=-1`, `matchDR=-1`, and
`matchRelPt=-1`.

The validator writes six CSVs using a requested output prefix: overall and
physics counts, prompt-charm per-system counts, per-D-hadron counts,
representative ambiguity/failure cases, and deltaR/relative-pT histograms.
The prompt-D CSV separates accepted generator particles from matched tracks
that also pass reconstructed pT/eta acceptance. Its generator-defined
two-track baseline counts unique additional QGenCharged ancestry assignments
without promoting source certainty. The source-certainty column remains the
Stage 1 value; Stage 2 never changes it.

## Build, run, validate

From the CMSSW_15_1_1 `src` directory, after checking out this package:

```bash
scram b -j 8
```

The standalone assignment regression test can also be run without CMSSW:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -I. \
  PhysicsTools/SimpleVertexTable/test/test_quark_reco_track_assignment.cc \
  -o /tmp/test_quark_reco_track_assignment
/tmp/test_quark_reco_track_assignment
```

Run one file per sample (replace the input paths with the local MiniAODSIM
files):

```bash
cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_matching_cfg.py \
  inputFiles=file:/path/to/HWWc_MINIAODSIM.root maxEvents=1000 outfile=HWWc_stage2.root

cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_matching_cfg.py \
  inputFiles=file:/path/to/ttbarhad_MINIAODSIM.root maxEvents=1000 outfile=ttbarhad_stage2.root
```

With Python packages `uproot` and `awkward` available, run:

```bash
python3 PhysicsTools/SimpleVertexTable/test/validate_quark_reco_track_matching.py \
  --input HWWc_stage2.root --sample HWWc --output-prefix HWWc_stage2

python3 PhysicsTools/SimpleVertexTable/test/validate_quark_reco_track_matching.py \
  --input ttbarhad_stage2.root --sample ttbarhad --output-prefix ttbarhad_stage2
```

For a faster diagnostic pass, limit validation to the first 100 events in the
ROOT file:

```bash
python3 PhysicsTools/SimpleVertexTable/test/validate_quark_reco_track_matching.py \
  --input HWWc_stage2.root --sample HWWc --output-prefix HWWc_stage2_first100 \
  --max-events 100
```

`--max-events` limits ROOT reading as well as processing. It selects a
deterministic prefix, not a random subset. The summary CSV records the file's
total entries, the requested limit, and the number actually processed; the
historical 1,000-event baseline counts should only be compared with a full,
matching sample.

The validator checks event-local references, one-to-one uniqueness, status and
sentinel consistency, candidate counts, and saved residuals. Its reference
Stage 1 generator-only counts are 66/239 for HWWc and 99/227 for ttbarhad;
the observed counts are reported separately and may differ if the input event
sample or event selection differs from the prior 1,000-event validation.
