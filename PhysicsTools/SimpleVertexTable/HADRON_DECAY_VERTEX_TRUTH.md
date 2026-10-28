# Generator decay-vertex truth and matched-track pair fits

This extension adds MC truth diagnostics alongside the existing QFrag and AssociationGNN products. The historical `track` and `trackPair` tables remain owned by `TrackGVFeaturesProducer`; truth-only rows are written to dedicated tables and never appended to `trackPair`.

## Existing relationships and persisted data

`QHadron` has one row per selected final-copy B or D hadron. Its `prunedGenIdx` is an index into the event's `prunedGenParticles`, not a row in `QTrainHadron` or the legacy `GVTable`. When `saveDecayVertexTruth=True`, it also stores an explicit `hadronIdx`, generator status and prompt/hard-process/last-copy flags, `p`, `px/py/pz`, energy, and the existing `pt`, `eta`, `phi`, and mass in GeV where appropriate. `QTrainHadron` deliberately retains one row for every `QHadron`, in the same event-local order; its `hadronIdx` is that row number. The new `QTruthMatchedHadron` makes this mapping explicit with both `hadronIdx` and `qTrainHadronIdx`.

`QGenCharged` has one row per selected stable packed charged generator particle. Its `lastPrunedIdx` is a pruned-generator index (or -1), `hadronIdx` is a `QHadron` row, and `membershipType=1` means the existing Stage-1 logic assigned it to a selected heavy-hadron decay chain. `QRecoTrack.genChargedIdx` is a `QGenCharged` row. The shared reconstructed track index is `QRecoTrack.trackIdx`; it is also the row index in the base track feature table and the index used by `trackPair.trk1Idx/trk2Idx`.

In the compact production config, the existing base track table, the legacy track-label extension, all original `trackPair` features/endpoints, `QFrag`, `QHadron`, `QGenCharged`, `QRecoTrack`, `QTrainTrack`, `QTrainPair`, `QTrainHadron`, and `QTrainQuark` are retained. `QTrainPair` has the same row order as `trackPair`; its `pairIdx` references the ordinary pair row. `QTrainLegacyHadron.hadronIdx` refers to the legacy `GVTable` row, while QTrainHadron rows refer to QHadron. These namespaces are never joined by assuming integer equality.

The old `track.hadidx` still indexes `GVTable`. It is a different namespace from all QHadron/QTrainHadron indices. This implementation does not compare those raw integers. Existing `QTrainTrack.hadronIdx` links to QTrainHadron and `QTrainTrack.matchStatus` and `nGenCandidates` preserve Stage-2 match ambiguity.

## Generator hadron geometry

Enable with `saveDecayVertexTruth=True` through the production configuration. When enabled, the Stage-1 `QHadron` table gets production and decay coordinates in cm. `prodX/Y/Z` come from the selected hadron record's own production vertex. The decay position is obtained from the production vertices of its immediate generator daughters. All finite immediate-daughter vertices must agree within `1e-4 cm`; otherwise `decayVertexValid=0`, and no production-coordinate fallback is used. `decayVertexSource=1` identifies this immediate-daughter method; zero means unavailable or inconsistent. `decayVertexResolved` records the same unambiguous resolution condition explicitly.

`trueLxy` and `trueL3D` are the transverse and 3D distances between the hadron production and resolved decay positions in cm. `trueCt` is `m (L dot p) / |p|^2` in cm. Each has a separate validity flag. No relation to a reconstructed PV is assumed. The same columns are copied to `QTrainHadron` when present, allowing compact ROOT production to retain them. If a hadron has no usable daughter vertex, it remains in the table with invalid flags and -999 coordinate/observable sentinels.

The use of pruned-record daughter links and daughter production coordinates is a MiniAOD truth approximation. Missing or inconsistent daughters remain explicitly invalid. The GenParticle interface exposes coordinates but no separate “coordinate was populated” bit; a finite origin-valued coordinate is therefore treated as present and cannot be distinguished from an unset zero solely from this record. This output does not infer a vertex from packed-particle ancestry, does not merge daughter branches, and does not claim that a generator record has a uniquely resolved decay point when the daughter records disagree.

## Matched tracks and daughter topology

With `enableTruthMatchedPairFit=True`, `QTruthMatchedTrack` contains one row for each `QRecoTrack` row. It preserves `trackIdx`, assigned `genChargedIdx`, selected `hadronIdx`, generator PDG ID/charge, generator production coordinates, `matchStatus`, and `nGenCandidates`. `daughterRelation` is 0 for no selected-hadron relationship, 1 when the selected QHadron record is an immediate mother of the matched pruned particle, 2 when the hadron is reached through an intermediate retained descendant, and 3 when a selected hadron exists but the retained mother chain cannot establish direct versus indirect ancestry. `immediateParentPrunedIdx` is stored only when exactly one retained immediate mother is available.

Only Stage-2 assigned generator rows with existing `QGenCharged.membershipType=1` and a valid selected `QHadron` row enter a hadron's reconstructed-track group. Both unopposed (`matchStatus=1`) and competing (`matchStatus=2`) assignments are retained with their status; unmatched tracks and unassigned/ambiguous QGenCharged rows remain visible in `QTruthMatchedTrack` but cannot be attributed to one hadron group. Stage-2 rows with no assigned generator particle (including status 3 conflicts) cannot be assigned to a hadron here. This is not a detector-hit truth match. Raw `lastPrunedIdx` and `immediateParentPrunedIdx` remain indices into the original MiniAOD pruned collection, which is not duplicated wholesale into NanoAOD; the selected particle PDG/charge, production point, parent PDG, hadron row, and computed direct/indirect status are saved so the relevant joins do not require that collection.

`QTruthMatchedHadron` contains every QHadron, including zero- and one-track groups, and stores the number of assigned reconstructed tracks plus `nMatchedRecoTracks choose 2`. Multiplicity codes elsewhere retain the existing N0/N1/N2/N3+ convention. Categories are 0 other, 1 B, 2 confirmed prompt D, 3 D-from-B, and 4 D-origin-unknown. Category 4 preserves uncertainty; missing pruned B ancestry is not used to claim prompt D.

## Supplementary pair fits

For every unordered pair in each hadron's assigned-track group, `QTruthMatchedPairFit` stores a hadron-specific row. Its `trackIdx1/2` endpoints refer to the shared unpacked reconstructed tracks, `genChargedIdx1/2` refer to QGenCharged, and `hadronIdx` refers to QHadron. `qTrainHadronIdx` is the explicit corresponding QTrainHadron row (currently equal because both tables retain one row per QHadron). The table preserves endpoint match statuses, direct/indirect classifications, `commonImmediateDecayVertex`, and N2/N3+ group multiplicity. That flag is set only when both retained generator particles have the selected hadron as immediate mother, their finite recorded production coordinates agree within `1e-4 cm`, and—when the saved hadron decay point is enabled—both coordinates agree with that resolved point. A pair of two descendants of one hadron is not assumed to originate at the same physical vertex.

`passedOrdinaryPairSelection` and `trackPairIdx` tell whether the pair appears in the unchanged ordinary edge table. A valid reference has exactly the ordinary pair's endpoint identities. The historical producer already builds one transient-track vector per event. When truth-pair production is enabled, it transfers that in-memory vector to the supplementary producer, so the truth-pair path does not rebuild transient tracks. If the existing selected-pair table contains fit branches, those fitted observables are copied and `fitReusedFromTrackPair=1`; no repeated fit is done. For any other physical pair, the existing `TwoTrackVertexFitter` is called and its result is stored. A small event-local cache prevents a second fit if the same physical track pair is referenced through another hadron ancestry row. `fitReusedFromEarlierTruthPair` marks that case. Invalid fits retain explicit status codes and sentinel-valued observables. The transient-track event product is only an in-memory handoff and is not included in NanoAOD output.

`ordinarySelectionFailure` uses 0 selected, 1 invalid transient track, 2 deltaR above 1.5, 3 pion-hypothesis pair mass above 5 GeV, 4 DCA above 0.1 cm, and 5 failed DCA/PCA calculation. The historical producer represents failed DCA calculation as -1 and does not reject it, so code 5 identifies a diagnostic edge case rather than asserting that the old producer rejected it. `pairDeltaR`, `pairInvMass`, and `pairDca` are included for offline diagnosis.

All `vtx*` outputs use the same names, status and validity semantics as ordinary `trackPair` fits: vertex positions and lengths are cm, covariance terms cm^2, and fit probability/dimensionless observables have no units. Covariance validity uses the helper's float-serialized positive-semidefinite test. Position/fit validity is separate from covariance validity; covariance-dependent uncertainties and significances remain governed by their individual validity flags.

`QTruthPairFitSummary` reports the CMS run/luminosity-block/event identifiers, the leading non-fake reconstructed PV position and covariance (with `recoPVValid`), event-local hadron/group/pair counts, selected-pair reuse availability, actual extra fitter calls, supplementary fit outcomes and producer elapsed time. This keeps reconstructed-PV-relative observables distinct from the saved generator hadron production position. `nPairsRequiringFit` counts physical pair fit calls after selected-fit reuse and event-local physical-pair caching.

## Configuration and execution

### Dedicated prompt-D train/validate configuration

`test/twotrackfit_train_validate_cfg.py` is the minimized configuration for
prompt-D decay-geometry and momentum studies. It is MC-only because it consumes
MiniAOD pruned and packed generator particles. The `train` mode uses the same
compact reconstructed-track inputs and edge selection as the package training
path. The `validate` mode adds the legacy `GVTable` and the full `QGNN*`
cross-check tables. Both modes retain `trackPair` with its original endpoints,
14 pair features, truth label and enabled two-track fit columns, and both retain
all supplementary truth-associated pair-fit rows. Neither mode writes unrelated
NanoAOD physics tables.

This path omits fragmentation-regression fields (`preHad*`, `zProj`, `zT`,
`deltaR_DQuark`, `kTrel`, and their masks), detailed QFrag copy-history
branches, and derived `trueLxy`, `trueL3D`, and `trueCt` columns. It retains
minimal `QFrag` signed PDG ID, four-momentum and origin category. The existing
`QHadron` and `QTrainHadron` association-quality and source-certainty fields
remain available even when the quark association is ambiguous or missing.
The new options are opt-in through this dedicated config; standard production
customizations keep their prior defaults.

Relevant event-local joins are: `QRecoTrack.trackIdx` to the shared
`trackTable` row; `QRecoTrack.genChargedIdx` to `QGenCharged`; `QGenCharged.hadronIdx`
to `QHadron`; `QTrainTrack.trackIdx` to that same reconstructed-track index;
`QTrainTrack.hadronIdx` and `quarkIdx` to `QTrainHadron` and `QTrainQuark`;
`QTrainPair.pairIdx` to the ordinary `trackPair` row; and
`QTruthMatchedPairFit.trackIdx1/2` to reconstructed tracks, with its `hadronIdx`
and `qTrainHadronIdx` naming their separate hadron namespaces. Legacy
`track.hadidx` remains a `GVTable` index and is not interchangeable with either
QHadron index.

Compile in a configured CMSSW release area:

```bash
scram b -j 8
```

Run a compact training ROOT production:

```bash
cmsRun PhysicsTools/SimpleVertexTable/test/twotrackfit_train_validate_cfg.py \
  mode=train inputFiles=file:/path/to/HZZc_MINIAODSIM.root maxEvents=10 \
  outfile=promptd_train_10.root
```

Run the larger validation-table variant:

```bash
cmsRun PhysicsTools/SimpleVertexTable/test/twotrackfit_train_validate_cfg.py \
  mode=validate inputFiles=file:/path/to/HZZc_MINIAODSIM.root maxEvents=10 \
  outfile=promptd_validate_10.root
```

Inspect either output with the existing lightweight validator:

```bash
python PhysicsTools/SimpleVertexTable/test/validate_hadron_decay_vertex_truth.py \
  promptd_train_10.root
python PhysicsTools/SimpleVertexTable/test/sanity_check_hadron_decay_vertex_truth.py \
  promptd_train_10.root --prefix promptd10
```

The validator checks row/index consistency, truth-pair completeness, vertex
validity, and fit fields. Its checks should be run on the local CMSSW output;
configuration syntax alone cannot verify the produced ROOT schema.

Compile in a configured CMSSW release area:

```bash
scram b -j 8
```

Run a 10-event smoke sample, replacing the input with a local or XRootD MiniAODSIM file:

```bash
cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_training_cfg.py \
  inputFiles=file:/path/to/HZZc_MINIAODSIM.root maxEvents=10 \
  enableTwoTrackVertexFit=True enableVertexFitTiming=True \
  enableTruthMatchedPairFit=True saveDecayVertexTruth=True \
  outfile=hzzc_decay_vertex_10.root
```

Representative HZZc production (start with 100 events, then use 1000 after checking resources):

```bash
/usr/bin/time -v -o hzzc_decay_vertex_100.time cmsRun \
  PhysicsTools/SimpleVertexTable/test/quark_reco_track_training_cfg.py \
  inputFiles=root://cms-xrd-global.cern.ch//store/mc/REPLACE_WITH_HZZC_MINIAODSIM.root \
  maxEvents=100 enableTwoTrackVertexFit=True enableVertexFitTiming=True \
  enableTruthMatchedPairFit=True saveDecayVertexTruth=True \
  outfile=hzzc_decay_vertex_100.root
```

Set `maxEvents=1000` and update output/time names for the larger sample. The configuration accepts the same unpacked tracks, original historical pair features, QFrag matching, and NanoAOD output rules as the existing compact `quark_reco_track_training_cfg.py`.

Install `uproot`, `awkward`, `numpy`, and optionally `matplotlib` in the analysis Python environment. Then run:

```bash
python PhysicsTools/SimpleVertexTable/test/validate_hadron_decay_vertex_truth.py hzzc_decay_vertex_10.root
python PhysicsTools/SimpleVertexTable/test/sanity_check_hadron_decay_vertex_truth.py hzzc_decay_vertex_10.root --prefix hzzc10
```

For a regression check, give the validator a baseline file created with the same input, event limit, track selection and original fitting switch before enabling truth output:

```bash
python PhysicsTools/SimpleVertexTable/test/validate_hadron_decay_vertex_truth.py \
  hzzc_decay_vertex_100.root --baseline hzzc_baseline_100.root
```

For a runtime comparison, run baseline and extended commands over the same input and event limit using `/usr/bin/time -v -o ...`. The baseline should use `enableTruthMatchedPairFit=False`; keep `enableTwoTrackVertexFit` identical in both runs. Compare original `trackPair` row counts/identities and `vertexFitSummary` timing, then inspect `QTruthPairFitSummary` for incremental pair counts and `supplementaryProducerMs`. Runtime depends on input access, CPU contention and selected truth multiplicities; no estimate is assumed here.

## Limitations and interpretation

The hadron groups are ancestry-based and inherit the available pruned/packed MiniAOD graph and Stage-2 kinematic matching. Competing Stage-2 matches are kept but are not promoted to unambiguous truth. The immediate-daughter flag is based on retained direct mother references, not a detector simulation record. Cascade descendants can share a QHadron ancestor while originating at distinct vertices. A generator decay vertex can be unavailable or inconsistent. These cases remain represented through validity and relation codes; consumers should apply their desired confidence requirements offline.

This feature was implemented by static repository inspection only. It has not been compiled or run under CMSSW in this environment.
