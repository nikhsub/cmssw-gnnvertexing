# AssociationGNN and QFrag ROOT output

`quark_reco_track_matching_cfg.py` now writes the historical hadron-level
AssociationGNN inputs together with the validated QFrag Stage 1/2 tables and
additional QFrag truth joins. It does not train a model or convert ROOT to
HDF5.

## Build and run

From a CMSSW_15_1_1 runtime area with this checkout available:

```bash
scram b -j 8
cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_matching_cfg.py \
  inputFiles=file:/path/to/input_MINIAOD.root maxEvents=100 \
  outfile=qfrag_association_gnn.root
python3 -m pip install 'uproot>=5.2,<6' 'awkward>=2.5,<3'
python3 PhysicsTools/SimpleVertexTable/test/validate_association_gnn_output.py \
  --input qfrag_association_gnn.root
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
| `QGNNHadron` | Same row order as `QHadron`; `hadronIdx`, accepted reconstructed decay-track count, N category, and hadron class. |

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
