# GNN seed downstream context schema

This schema supplements the existing AssociationGNN graph and top-128 `Seed`
tables with event-level reconstructed PF, PV, and MiniAOD SV context. It does
not alter the model, its features, graph selection, seed ranking, or selected
seed fits. Objects are saved once per event; seed-relative neighborhoods and
interactions are left to downstream dataset construction.

The table prefix in ROOT is the NanoAOD flat-table name. Integer indices are
event-local. Four-momenta and positions are stored in GeV and cm respectively;
angles are radians. Unavailable scalar or derived float values use `-9999`
unless otherwise stated. Finite covariance elements are preserved even when
the covariance matrix validity flag is zero. Integer references that cannot be
resolved use `-1`.

## Existing frozen seeder tables

The config retains `track`, `trackPair`, `SeederNode`, `SeederEmbedding`,
`Seed`, and `SeedSummary` as produced by the existing modules.

- `track` rows are `filteredTracks` rows and are the graph node index space.
- `trackPair` rows retain the existing directed graph pairs and pair features;
  `trackPair_trk1Idx` and `trackPair_trk2Idx` reference `track` rows.
- `SeederNode_trackIdx` references `track`. Node scores and validity remain
  unchanged.
- `SeederEmbedding_trackIdx` references `track`; `SeederEmbedding_valid` tells
  whether `emb0` through `emb127` contain the inference embedding.
- `Seed_trk1Idx` and `Seed_trk2Idx` reference `track`. `Seed` rows are ranked
  unique physical pairs and include the existing two-track fit results.

To obtain the two seed embeddings, for a seed row `s`, look up the `SeederEmbedding`
rows whose `trackIdx` equals `Seed_trk1Idx[s]` and `Seed_trk2Idx[s]`, then require
`valid == 1`. An additional selected graph track uses the same direct lookup.
No embedding is defined for a charged PF row whose `filteredTrackIdx` is `-1`,
and embeddings are not copied into PF tables. With `saveEmbeddings=False`, the
embedding table is absent and every PF row's `embeddingValid` is zero.

## `CtxChargedPF`

Rows are selected independently from `packedPFCandidates` and `lostTracks`
using `pt >= minChargedPt`, `abs(eta) <= maxPFAbsEta`, and nonzero charge.
Defaults are 0.1 GeV and 3.0. No seeder-track, displacement, SV, or jet cut is
applied. Source ID 0 means packed PF; 1 means lost track. `sourceIdx` is the
original index in that input collection.

| Branch suffix | Type | Meaning |
|---|---|---|
| `sourceCollectionId`, `sourceIdx` | int | Original candidate identity; source IDs 0=packed, 1=lost. |
| `unpackedTrackIdx`, `filteredTrackIdx` | int | Exact source ValueMap -> unpacked -> filtered track row chain; `-1` means absent at that stage. |
| `embeddingValid` | int | `SeederEmbedding_valid` for the mapped filtered row; zero if no mapping/table/valid vector. |
| `pt`, `eta`, `phi`, `mass`, `energy`, `px`, `py`, `pz` | float | Raw packed-candidate four-momentum. |
| `charge`, `pdgId` | int | Candidate charge and PDG ID. |
| `puppiWeight`, `puppiWeightNoLep` | float | Raw stored PUPPI weights, not applied to the four-momentum. |
| `hasTrackDetails` | int | Whether packed track details are available. |
| `trackHighPurity`, `lostInnerHits`, `numberOfHits`, `numberOfPixelHits`, `numberOfStripHits` | int | Packed track details; hit fields are `-1` when unavailable. Strip hits are `numberOfHits-numberOfPixelHits`. |
| `fromPV`, `pvAssociationQuality` | int | PackedCandidate PV classification/quality; `fromPV=-1` when no usable source PV association is available. |
| `associatedPVIdx`, `associatedPVValid` | int | Original `offlineSlimmedPrimaryVertices` index from the candidate vertex reference, or `-1`/0 if it does not resolve into that collection. |
| `dxy`, `dz`, `dxyError`, `dzError` | float | PackedCandidate impact parameter values/uncertainties. `dxy` and `dz` are explicitly evaluated relative to PV0; errors are the stored packed-candidate errors. |
| `ipValid`, `ipErrorValid` | int | Validity flags for the corresponding value and error pairs. |
| `ipReference` | int | 0 means PV0 (`offlineSlimmedPrimaryVertices[0]`), -1 means unavailable. |

The `filteredTrackIdx` mapping uses `packedTrackIndex` or `lostTrackIndex`
followed by `filteredTracks:filteredIndex`. It is exact and does not use
kinematic matching. A mapped row is the input key for `SeederEmbedding`; require
both a nonnegative row and `embeddingValid == 1` before consuming the vector.

## `CtxNeutralPF`

Rows come only from neutral `packedPFCandidates` with `pt >= minNeutralPt` and
`abs(eta) <= maxPFAbsEta` (defaults 0.1 GeV and 3.0). `sourceIdx` is the original
packed-candidate index. Charge is zero by selection; `pdgId` distinguishes
photons, neutral hadrons, and other neutral categories.

| Branch suffix | Type | Meaning |
|---|---|---|
| `sourceIdx`, `pdgId`, `charge` | int | Original packed index and particle identity. |
| `pt`, `eta`, `phi`, `mass`, `energy`, `px`, `py`, `pz` | float | Raw four-momentum. |
| `puppiWeight`, `puppiWeightNoLep` | float | Raw stored PUPPI weights. |
| `caloFraction`, `hcalFraction`, `rawCaloFraction`, `rawHcalFraction` | float | PackedCandidate calorimeter fractions; `-9999` unless finite and meaningful. |
| `caloFractionsValid` | int | 1 when all four fraction values are finite and energy is positive. |

Neutral candidates have no GNN track index or embedding.

## `CtxPFSummary`

This is a singleton table. Its integer branches are `nPackedInput`, `nLostInput`,
`nUnpackedTracks`, `nFilteredTracks`, `nChargedPF`, `nNeutralPF`,
`nChargedMappedToFilteredTracks`, `nChargedUnmapped`, and
`nInvalidReferences`. Input sizes are before PF cuts. The charged counts refer
to rows in `CtxChargedPF`; mapped/unmapped counts partition that table.
`nInvalidReferences` counts missing/malformed ValueMap mappings or out-of-range
indices; normal negative selection results are not errors.

## `CtxPV`

Contains all rows of `offlineSlimmedPrimaryVertices`, in their original order;
`sourceIdx` is the original PV row. PV0 is row 0 and the same source PV used for
the seed-fit path. PF `associatedPVIdx` retains the original source index. PF
impact parameters and SV geometry in this context explicitly use PV0, not the
candidate-associated PV. `PATTrackFromVertexFilter` builds the PV collection
used by the seed fitter from this source: it preserves PV ordering, positions,
and covariances while remapping track references. The `track` and `Seed`
outputs retain their existing PV conventions.

| Branch suffix | Type | Meaning |
|---|---|---|
| `sourceIdx`, `isFake`, `nTracks` | int | Original PV index, fake flag, and referenced track count. |
| `x`, `y`, `z` | float | Position [cm]. |
| `chi2`, `ndof`, `chi2Ndof` | float | Fit quality; ratio unavailable if `ndof <= 0`. |
| `covXX`, `covXY`, `covXZ`, `covYY`, `covYZ`, `covZZ` | float | Position covariance [cm^2]. |
| `positionValid`, `covarianceValid`, `chi2NdofValid` | int | Corresponding validity flags. |

## `CtxSV`

Contains every `slimmedSecondaryVertices` row without the `mySV` displacement
significance selection. `sourceIdx` is the original collection index and, since
no rows are removed, equals the event-local `CtxSV` row. These are reconstructed
MiniAOD SVs, not the fitted two-track `Seed` vertices.

| Branch suffix | Type | Meaning |
|---|---|---|
| `sourceIdx`, `nDaughters`, `nMappedTracks` | int | Original SV index, full source daughter count, and daughters with a valid `filteredTrackIdx`. |
| `pt`, `eta`, `phi`, `mass`, `energy` | float | Source SV four-momentum. |
| `x`, `y`, `z` | float | Source SV position [cm]. |
| `chi2`, `ndof`, `chi2Ndof` | float | Source vertex fit quality; ratio requires positive ndof. |
| `covXX`, `covXY`, `covXZ`, `covYY`, `covYZ`, `covZZ` | float | Source SV position covariance [cm^2]. |
| `covarianceValid`, `chi2NdofValid` | int | Covariance must be finite and positive semidefinite; ratio validity requires positive ndof. |
| `pvGeometryValid` | int | Nonfake PV0, finite PV/SV positions, and valid PV/SV covariances were available for geometry calculations. |
| `lxy`, `lxyErr`, `lxySig`, `l3d`, `l3dErr`, `l3dSig` | float | `VertexDistanceXY`/`VertexDistance3D` PV0-to-SV measurement, uncertainty [cm], and significance. |
| `lxyValid`, `lxyErrValid`, `lxySigValid`, `l3dValid`, `l3dErrValid`, `l3dSigValid` | int | Independent validity flags for each distance, uncertainty, and significance. |
| `cosPointing2D`, `cosPointing3D` | float | PV0-to-SV displacement direction dotted with SV transverse/3D momentum direction. |
| `cosPointing2DValid`, `cosPointing3DValid` | int | Direction vectors and displacement are nonzero and finite. |

Invalid float values use `-9999`; rows are retained even when a validity flag
is zero. The distance measurements use CMSSW `VertexDistanceXY` and
`VertexDistance3D`; pointing cosines use the normalized displacement and SV
momentum and are bounded to [-1,1].

## `CtxSVConstituent`

One row is written for every daughter pointer of every `slimmedSecondaryVertices`
entry, including daughters without a GNN track. `svIdx` references `CtxSV`;
`sourceSVIdx` and `daughterIdx` preserve the original SV row and daughter order.

| Branch suffix | Type | Meaning |
|---|---|---|
| `svIdx`, `sourceSVIdx`, `daughterIdx` | int | Event-local CtxSV row, original source SV row, and source daughter ordinal. |
| `sourceCollectionId`, `sourceIdx` | int | 0=packed PF, 1=lost tracks, -1=unrecognized/null source; source key is preserved when pointer exists. |
| `unpackedTrackIdx` | int | Exact source-to-unpacked ValueMap result or -1. |
| `chargedPFIdx` | int | Exact `CtxChargedPF` row for charged selected source candidates, otherwise -1. |
| `filteredTrackIdx` | int | Exact filtered graph-track row, otherwise -1. |
| `mappingValid`, `chargedPFMappingValid`, `filteredTrackMappingValid` | int | Validity for source identity, charged-PF reference, and graph-track reference respectively. |
| `pt`, `eta`, `phi`, `charge` | float/ int | Original daughter kinematics and charge. |

For a valid charged-PF reference, both `(sourceCollectionId, sourceIdx)` pairs
must match exactly. `filteredTrackIdx` is obtained independently from the
source ValueMaps. A daughter omitted by PF thresholds may still have a valid
track mapping; a daughter omitted from the GNN graph remains present with
`filteredTrackIdx=-1`. No approximate matching is used.

## Pairwise information and offline work

The existing directed `trackPair` rows, their seeder logits/scores, and all
existing `Seed` fit observables remain unchanged. Pair rows missing from the
selected graph are absent; downstream code must use an explicit graph-edge
existence mask rather than treating absence as a negative score. Only the
selected top-K `Seed` pairs receive Kalman fits. `CtxSV` refers to original
MiniAOD SVs and must not be conflated with these fits.

Construct seed-relative angles, momentum fractions, local neighborhoods,
PF-to-seed or PF-to-SV compatibility, and other pairwise quantities offline.
There is no jet axis, seed-centered truncation, seed-to-extra-track vertex fit,
full-graph vertex fit, or duplicated per-seed PF collection in this output.

## Commands and validation

From a configured CMSSW release area:

```bash
cd "$CMSSW_BASE/src"
scram b -j 8
```

MC:

```bash
cmsRun PhysicsTools/SimpleVertexTable/test/trackGNN_downstream_context_cfg.py \
  inputFiles=file:/path/to/miniaod.root maxEvents=100 \
  outfile=downstream_context_mc.root isMC=True maxSeeds=128 saveEmbeddings=True
```

Collision data:

```bash
cmsRun PhysicsTools/SimpleVertexTable/test/trackGNN_downstream_context_cfg.py \
  inputFiles=file:/path/to/miniaod_data.root maxEvents=100 \
  outfile=downstream_context_data.root isMC=False maxSeeds=128 saveEmbeddings=True
```

Inspect table branches with ROOT:

```bash
python - <<'PY'
import uproot
f = uproot.open("downstream_context_mc.root")
for name in ("Events",):
    t = f[name]
    for key in t.keys():
        if any(x in key for x in ("Ctx", "Seed", "Seeder", "trackPair", "mySV", "myPV")):
            print(key)
PY
```

Run integrity checks; optionally compare frozen outputs against the prior
seeder file using event identities:

```bash
python PhysicsTools/SimpleVertexTable/test/validate_trackGNN_downstream_context.py \
  downstream_context_mc.root --baseline previous_seeder.root
```

The validator requires `uproot`, `awkward`, and NumPy. Runtime and CMSSW API
compatibility must be confirmed in the target CMSSW_10_6_30 area.
