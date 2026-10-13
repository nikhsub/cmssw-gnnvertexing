# AssociationGNN and QFrag ROOT output

Two independent configs share the same feature and truth producer chain:
`quark_reco_track_training_cfg.py` writes the compact training schema and
`quark_reco_track_matching_cfg.py` writes the full validation schema. The
compact output omits full intermediate generator and legacy GV tables after
the required targets have been computed. Neither config trains a model or
converts ROOT to HDF5.

## Build and run

From a CMSSW_15_1_1 runtime area with this checkout available:

```bash
scram b -j 8
# Compact training output
cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_training_cfg.py \
  inputFiles=file:/path/to/input_MINIAOD.root maxEvents=100 \
  outfile=qfrag_training.root compressionAlgorithm=ZSTD compressionLevel=5
# Full validation output
cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_matching_cfg.py \
  inputFiles=file:/path/to/input_MINIAOD.root maxEvents=100 \
  outfile=qfrag_full.root compressionAlgorithm=ZSTD compressionLevel=5
python3 -m pip install 'uproot>=5.2,<6' 'awkward>=2.5,<3'
python3 PhysicsTools/SimpleVertexTable/test/validate_association_gnn_output.py \
  --input qfrag_full.root
python3 PhysicsTools/SimpleVertexTable/test/compare_qfrag_training_outputs.py \
  --full qfrag_full.root --compact qfrag_training.root
```

The input must contain the MiniAOD `prunedGenParticles`, `packedGenParticles`,
`packedPFCandidates`, `lostTracks`, and vertex collections expected by the
existing producers. The production config retains the existing Stage 2 cuts
(`maxDeltaR=0.03`, `maxRelPt=0.5`) and uses the existing unpacker with
`minTrackPt=0.4 GeV`.

## Tables and indices

All row indices are event-local. `trackIdx` is the position in the unchanged
`unpackedTracksAndVertices` collection. No track is filtered or reordered in
the feature or truth join tables.

| Table | Contents and links |
| --- | --- |
| `QFrag`, `QHadron`, `QGenCharged` | Original Stage 1 outputs, unchanged. `QHadron.quarkIdx` links a hadron to `QFrag`; charged rows retain hadron/quark indices, quality, membership, and source-certainty fields. `QHadron.preHadPt/Eta/Phi` and `QFrag` kinematics preserve the inputs for auxiliary targets. |
| `QRecoTrack` | Original Stage 2 output, unchanged. One row per unpacked track; `trackIdx` is the input track row and `genChargedIdx` links to `QGenCharged`, or `-1`. Match diagnostics and one-to-one assignment status remain available. |
| `track` | Existing NanoAOD base track table plus the `TrackGVFeaturesProducer` extension (`hadidx`, `label`). The base table is row-aligned with `QRecoTrack`. Its `genPartIdx` branch is truth-only and must be excluded from model inputs. |
| `trackPair` | Existing `TrackGVFeaturesProducer` output: one row per physical pair, with `trk1Idx < trk2Idx`, the historical candidate selection, all pair features, and legacy `edgeLabel`. |
| `GVTable` | Original generator hadron rows used by the historical `track.hadidx` matching; persisted so legacy hadron indices remain resolvable. Daughter rows are used during production but are not saved because no output index refers to them. |
| `QGNNTrack` | QFrag truth join, one row per `QRecoTrack` row. `trackIdx`, `genChargedIdx`, `hadronIdx`, `quarkIdx`, `associationQuality`, `membershipType`, and `fragmentationSourceCertain` preserve the ancestry link and its certainty. It also provides `hadronClass`, `hadronTruthValid`, `isB`, `isPromptD`, and `isDFromB`. |
| `QGNNPairTruth` | Same row order as `trackPair`; `pairIdx` indexes that table. Its row endpoints are `trackPair.trk1Idx/trk2Idx`. It contains QFrag-derived `sameHadron`, per-class hadron pair labels, and source-certain `sameQuark` targets with validity masks. |
| `QGNNHadron` | Full-validation table in the same row order as `QHadron`; `hadronIdx`, accepted reconstructed decay-track count, N category, and hadron class. |
| `QTrainTrack` | Compact truth-only track targets and references. The compact config saves no duplicated reconstructed node features. `hadronIdx` references `QTrainHadron`; `quarkIdx` references `QTrainQuark`. It also stores Stage-2 matching diagnostics, ancestry quality, membership, source certainty, strict/weak validity, and B/prompt-D/D-from-B/D-unknown labels. |
| `QTrainPair` | Compact same-hadron, class-specific hadron, strict same-quark, and weak same-quark targets in exactly the `trackPair` row order. Invalid or ambiguous targets use -1 with a separate validity column. |
| `QTrainHadron` | One compact object row per original `QHadron`, including N0 objects, hadron class inputs, parent/hadron/quark row links, pre-hadronization proxy kinematics, N category, auxiliary targets, and validity masks. |
| `QTrainQuark` | One row per original `QFrag`; `pdgId` preserves quark flavor. Row number is the compact quark-system index. |
| `QTrainLegacyHadron` | Minimal legacy `GVTable` row map used to resolve `track.hadidx` in compact files. |

The full config saves all columns from those tables plus these intermediate
schemas (column names are prefixed in the ROOT `Events` tree by their table
name):

| Full table | Columns |
| --- | --- |
| `QFrag` | `pdgId`, `prunedGenIdx`, `pt`, `eta`, `phi`, `mass`, `hardPrunedIdx`, `hardPt`, `hardEta`, `hardPhi`. |
| `QHadron` | `prunedGenIdx`, `pdgId`, `pt`, `eta`, `phi`, `mass`, `charmOrigin`, `parentBHadronIdx`, `quarkIdx`, `associationQuality`, `fragmentationSourceCertain`, `preHadPrunedIdx`, `preHadPt`, `preHadEta`, `preHadPhi`, `preHadKinematicsValid`. |
| `QGenCharged` | `packedGenIdx`, `pdgId`, `charge`, `pt`, `eta`, `phi`, `mass`, `lastPrunedIdx`, `lastPrunedIsSelf`, `hadronIdx`, `quarkIdx`, `associationQuality`, `membershipType`, `fragmentationSourceCertain`, `isFromB`, `isFromPromptD`, `isFromBtoD`. |
| `QRecoTrack` | `trackIdx`, `pt`, `eta`, `phi`, `charge`, `isHighPurity`, `genChargedIdx`, `matchDR`, `matchRelPt`, `matchStatus`, `nGenCandidates`. |
| `QGNNTrack` | `trackIdx`, `genChargedIdx`, `hadronIdx`, `quarkIdx`, `associationQuality`, `membershipType`, `fragmentationSourceCertain`, `hadronClass`, `hadronTruthValid`, `isB`, `isPromptD`, `isDFromB`. |
| `QGNNPairTruth` | `pairIdx`, `hadronTruthValid`, `sameHadron`, `sameB`, `samePromptD`, `sameDFromB`, `quarkTruthValid`, `sameQuark`. |
| `QGNNHadron` | `hadronIdx`, `nAcceptedDecayTracks`, `nCategory`, `hadronClass`. |
| `GVTable` | `pt`, `eta`, `phi`, `x`, `y`, `z`, `x_i`, `y_i`, `z_i`, `Hadron_pdgId`, `isB`, `isD`, `fromHF`, `toHF`, `pdgClass`, `nDaughters`, `maxDaughterPairDeltaR`. |

The full base `track` additionally retains the unused `genPartIdx` and
`timeSig` columns from the original setup; the compact base table contains
only the 22 node features and the two legacy truth extensions.

### Compact versus full output

The compact config saves `track` (the 22 historical reconstructed node
features plus legacy `hadidx` and `label`), `trackPair` (the 14 edge features,
endpoints, and legacy `edgeLabel`), and `QTrainTrack`, `QTrainPair`,
`QTrainHadron`, and `QTrainQuark`. The compact base `track` table omits the
unused truth-only `genPartIdx` and unused `timeSig` columns; `trackGenMatch` is
not run there. The legacy `TrackGVFeaturesProducer` labels and edge features
are still produced with their unchanged matching and pair definitions.

The compact products and columns are:

| Product | Saved columns |
| --- | --- |
| `track` | The 22 feature columns listed above, then `hadidx`, `label`. |
| `trackPair` | `trk1Idx`, `trk2Idx`, `deltaR`, `deltaPhi`, `deltaEta`, `dca`, `dcaSig`, `cptopv`, `pvtoPCA1`, `pvtoPCA2`, `dotprod1`, `dotprod2`, `pairMom`, `pairInvMass`, `pt`, `eta`, `edgeLabel`. |
| `QTrainTrack` | `trackIdx`, `hadronIdx`, `quarkIdx`, `matchStatus`, `matchDR`, `matchRelPt`, `nGenCandidates`, `associationQuality`, `membershipType`, `fragmentationSourceCertain`, `quarkWeakValid`, `quarkStrictValid`, `hadronClass`, `hadronTruthValid`, `isB`, `isPromptD`, `isDFromB`, `isDUnknown`. |
| `QTrainPair` | `hadronTruthValid`, `sameHadron`, `sameB`, `samePromptD`, `sameDFromB`, `quarkTruthValid`, `sameQuark`, `quarkWeakTruthValid`, `sameQuarkWeak`. |
| `QTrainHadron` | `hadronIdx`, `pdgId`, `charmOrigin`, `pt`, `eta`, `phi`, `parentBHadronIdx`, `quarkIdx`, `associationQuality`, `fragmentationSourceCertain`, `preHadPt`, `preHadEta`, `preHadPhi`, `preHadKinematicsValid`, `nAcceptedDecayTracks`, `nCategory`, `zProj`, `zT`, `deltaR_DQuark`, `kTrel`, `observableInputsValid`, `auxiliaryTargetValid`. |
| `QTrainQuark` | `pdgId`; the row position is the quark-system index. |
| `QTrainLegacyHadron` | `hadronIdx`, `pdgId`, `isB`, `isD`, `fromHF`; row position matches legacy `GVTable`, which is the reference used by `track.hadidx`. |

The full config keeps the above training information and the full debugging
tables: `GVTable`, `QFrag`, `QHadron`, `QGenCharged`, `QRecoTrack`, `QGNNTrack`,
`QGNNPairTruth`, and `QGNNHadron`. It also keeps `QTrain*` products so the
comparison utility can test the compact targets directly against the same
events and producer calculation. `GVDaughtersTable` remains an in-memory input
to legacy label production and is not persisted because no saved reference
requires its rows.

`QTrainLegacyHadron` is the small exception to omitting the legacy generator
table: it retains the minimal flags needed to resolve the saved legacy
`track.hadidx` and validate `track.label` without writing the full legacy
`GVTable` or daughter collections in compact files.

`QTrainHadron` preserves one row for every QHadron, so N0 hadrons are retained.
Its `nAcceptedDecayTracks` and `nCategory` use the existing definition:
membership type 1 with matched reconstructed pT>=0.4 GeV and |eta|<=2.5;
categories 0, 1, 2, 3 mean N0, N1, N2, N3+. The auxiliary values use the
Stage-3 definitions in double precision and are stored as ROOT float columns;
invalid values are NaN and controlled by `observableInputsValid`. The stricter
`auxiliaryTargetValid` mask selects confirmed prompt D hadrons with a unique
charm-quark association and valid finite proxy momentum, matching the previous
Stage-3 training-target eligibility rule.

`QTrainPair` omits a redundant `pairIdx`: its row number is the corresponding
`trackPair` row, and `trackPair.trk1Idx/trk2Idx` provide the graph endpoints.
Strict same-quark targets require unique ancestry and source certainty. Weak
same-quark candidates require unique ancestry and membership type 1 or 2;
uncertain additional candidates can contribute only to this explicitly weak
target, never to `sameQuark`. The compact track table preserves each candidate's
membership, quality, and certainty separately.

The label vocabularies remain separate:

| Label family | Code | Meaning |
| --- | ---: | --- |
| Legacy `track.label` | 0 | Primary vertex |
|  | 1 | Pileup |
|  | 2 | From B |
|  | 3 | From B-to-charm/heavy flavor |
|  | 4 | Prompt charm |
|  | 5 | Other matched secondary |
|  | 6 | Unassociated |
| QFrag `charmOrigin` | 0 | Not a D hadron |
|  | 1 | Confirmed prompt D |
|  | 2 | D from B |
|  | 3 | D origin unknown |
| QFrag `hadronClass` | -1 | Unresolved/unassigned track target |
|  | 0 | Other hadron |
|  | 1 | B hadron |
|  | 2 | Prompt D |
|  | 3 | D from B |
|  | 4 | D with unknown origin |

`track.hadidx` indexes `GVTable`; `QTrainTrack.hadronIdx` indexes
`QTrainHadron`; `QTrainTrack.quarkIdx` indexes `QTrainQuark`. These are
event-local row numbers, not interchangeable generator-particle indices.

### Storage and throughput comparison

Both configs accept `compressionAlgorithm` and `compressionLevel`; ROOT's
lossless ZSTD at level 5 is the default. Compare algorithms on the same input
sample and event count, recording elapsed time, CPU time, and compressed bytes:

```bash
for mode in training matching; do
  if [ "$mode" = training ]; then
    cfg=PhysicsTools/SimpleVertexTable/test/quark_reco_track_training_cfg.py
  else
    cfg=PhysicsTools/SimpleVertexTable/test/quark_reco_track_matching_cfg.py
  fi
  for algo in ZLIB LZMA LZ4 ZSTD; do
    /usr/bin/time -v cmsRun "$cfg" inputFiles=file:/path/to/input_MINIAOD.root \
      maxEvents=1000 outfile="${mode}_${algo}.root" \
      compressionAlgorithm="$algo" compressionLevel=5 2>"${mode}_${algo}.time"
    stat -c '%n %s bytes' "${mode}_${algo}.root"
  done
done
/usr/bin/time -v cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_training_cfg.py \
  inputFiles=file:/path/to/input_MINIAOD.root maxEvents=1000 outfile=training.root \
  2>training.time
/usr/bin/time -v cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_matching_cfg.py \
  inputFiles=file:/path/to/input_MINIAOD.root maxEvents=1000 outfile=full.root \
  2>full.time
stat -c '%n %s bytes' training.root full.root
```

`trackPair` is still exactly the historical feature table: endpoints and all
float features are retained at their existing precision. No lossy quantization
or pair thinning is applied. The compact saving comes from omitting full
generator and legacy GV tables, excluding unused legacy truth columns in the
base track table, and omitting a redundant pair-row index. Pair-building CPU
cost stays quadratic and identical between configs because the existing
producer is reused unchanged.

No representative HWWc, ttbar, or HZZc files or CMSSW runtime are mounted in
this workspace, so actual ROOT size, runtime, compilation, and sample-by-sample
identity measurements have not been produced here. Run both configs against
the same event lists, then run `compare_qfrag_training_outputs.py`; it checks
event IDs, all 22 node inputs, all 14 edge inputs/endpoints, both legacy label
sets, QTrain tables, QGNN label equivalence, and independently recomputed
auxiliary targets.

## Feature order used by AssociationGNN

Use the following 22 existing `track` branches in this exact order for the
node tensor:

```text
track_isHighPurity, track_charge, track_nHits,
track_numberOfPixelHits, track_numberOfStripHits, track_passSeed,
track_chi2, track_energy, track_eta, track_phi, track_pt,
track_px, track_py, track_pz, track_qdotp, track_qoverp,
track_dz, track_dzSignificance, track_ip2dSignificance, track_ip2dValue,
track_ip3dSignificance, track_ip3dValue
```

Use these 14 `trackPair` branches in this exact order for the physical-edge
feature tensor:

```text
trackPair_deltaR, trackPair_deltaPhi, trackPair_deltaEta,
trackPair_dca, trackPair_dcaSig, trackPair_cptopv,
trackPair_pvtoPCA1, trackPair_pvtoPCA2,
trackPair_dotprod1, trackPair_dotprod2, trackPair_pairMom,
trackPair_pairInvMass, trackPair_pt, trackPair_eta
```

These are reconstructed features. Do not include `track_genPartIdx`, the
`track` truth extension, `QGNNTrack`, `QGNNPairTruth`, or generator tables in
the inference inputs. The pasted AssociationGNN implementation stores each
physical pair once; it creates both directed messages internally and predicts
one edge score per stored pair, so `trackPair` and `QGNNPairTruth` remain
undirected one-row-per-pair tables.

## Selection and feature definitions

The node collection is the full unpacked collection produced with
`minTrackPt=0.4 GeV`; no high-purity or eta filter is applied to node rows.
`TrackGVFeaturesProducer` is reused unchanged for the historical truth
extension and pair table. It applies its existing `pT>=0.4 GeV`,
`deltaR<0.02`, and relative-pT `<0.2` requirements only to its legacy
track-to-GV-daughter truth match. It does not filter node rows.

The historical pair loop considers `i<j`, requires valid transient tracks,
`deltaR<=1.5`, and pion-hypothesis invariant mass `<=5 GeV`; it drops a pair
only when a calculated DCA exceeds `0.1 cm`. If the closest-approach
calculation fails, the historical `-1` DCA sentinel is retained. It stores
signed `deltaPhi`, absolute `deltaEta`, DCA and DCA
significance, crossing-point/PCA distances to the leading PV, two direction
dot products, pair momentum, pion-hypothesis mass, pair pT, and pair eta. The
historical `edgeLabel` is preserved exactly: 1 when both legacy `hadidx`
values are nonnegative and equal, otherwise 0 (including unmatched legacy
tracks). Use `QGNNPairTruth.hadronTruthValid` when a masked label is needed.

## Truth labels and categories

The legacy `track.label` codes from `TrackGVFeaturesProducer` remain unchanged:
0 primary, 1 pileup, 2 from B, 3 from B-to-charm/heavy flavor, 4 prompt charm,
5 other matched secondary, and 6 unassociated. The legacy `track.hadidx` is a
row in the legacy GV table, matching the original training convention.
That index is not interchangeable with `QGNNTrack.hadronIdx`, which indexes
`QHadron`; the corresponding `GVTable` and `GVDaughtersTable` are also saved.

`QGNNTrack.hadronIdx` is instead a `QHadron` row. Its `hadronClass` codes are
`-1` unresolved, `0` other hadron, `1` B hadron, `2` confirmed prompt D,
`3` D from B, and `4` D with unknown origin. `hadronTruthValid=1` requires a
Stage-2 matched `QGenCharged` row with `membershipType=1` and a valid
`QHadron` row. `associationQuality` is preserved separately because it also
describes quark-lineage certainty; it does not erase an otherwise valid
heavy-hadron decay label.

For `QGNNPairTruth`, `sameHadron` is 1 for the same known QHadron, 0 for two
known but different QHadrons, and -1 when either endpoint lacks a valid decay
hadron assignment. `sameB`, `samePromptD`, and `sameDFromB` use the same mask.
The legacy `trackPair.edgeLabel` is retained alongside this Stage-2-derived
target because the two assignments intentionally come from their respective
validated producers.

`sameQuark` is valid only when both matched charged rows have
`associationQuality=1`, a valid quark row, and
`fragmentationSourceCertain=1`. Uncertain additional candidates remain in
`QGNNTrack` with their membership, indices, quality, and certainty fields, but
produce `sameQuark=-1` when paired with a track that cannot establish a
source-certain relation. This prevents them from becoming confirmed
same-quark positives.

`QGNNHadron.nAcceptedDecayTracks` counts distinct matched rows with
`membershipType=1`, valid `hadronIdx`, `pT>=0.4 GeV`, and `|eta|<=2.5`.
`nCategory` is 0 for N0, 1 for N1, 2 for N2, and 3 for N3+.
Additional candidates do not contribute to this decay multiplicity.

The legacy producer scans all unordered track pairs, so pair construction time
is quadratic in the number of unpacked tracks. It stores only pairs passing
the historical geometric, mass, and DCA requirements. The ROOT output adds
one compact truth row per stored pair; actual event size and compressed file
size depend on the input sample and must be measured on a representative
processing job. The edge feature, endpoint, and label payload is about 100
uncompressed bytes per retained pair, before ROOT branch overhead.

## Auxiliary fragmentation targets

The ROOT output already contains QHadron four-momentum (`pt`, `eta`, `phi`),
the pre-hadronization proxy (`preHadPt`, `preHadEta`, `preHadPhi`), and the
QHadron-to-QFrag index and quality fields needed for target construction.
Using the same definitions as the Stage 3 study:

```text
zProj = dot(pD, pQ) / |pQ|^2
zT = pT_D / pT_Q
deltaR_DQuark = deltaR(D, pre-hadronization quark proxy)
kTrel = |pD x pQ| / |pQ|
```

These generator-derived values are training targets only. They and all
ancestry-derived candidate counts must stay out of inference features.
