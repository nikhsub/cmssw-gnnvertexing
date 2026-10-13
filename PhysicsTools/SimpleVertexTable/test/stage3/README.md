# Stage 3 HZZc track-only reconstruction feasibility

This is an analysis and dataset-feasibility framework only. It does not train a model, produce graph data for training, cluster jets, or alter the Stage 1/2 producers. It consumes the four existing NanoAOD tables and keeps run/luminosityBlock/event identifiers in all row-level outputs.

## Inputs and verified conventions

The combined CMSSW configuration is the existing [`../quark_reco_track_matching_cfg.py`](../quark_reco_track_matching_cfg.py). It schedules the existing Stage 1 `QuarkFragmentationProducer`, the existing `unpackedTracksAndVertices` collection, and Stage 2 `QuarkRecoTrackMatchProducer`, in that order. Its `VarParsing` interface accepts one or more `inputFiles`, an `outfile`, and `maxEvents`; it writes `QFrag`, `QHadron`, `QGenCharged`, and `QRecoTrack` using `NanoAODOutputModule`. Stage 2 is configured with `maxDeltaR=0.03` and `maxRelPt=0.5`. The NanoAOD `Events` tree supplies `run`, `luminosityBlock`, and `event`.

The analyzer checks every event-local reference and table-column length before use. Its checked producer conventions are:

| Branch | Codes used |
|---|---|
| `QHadron.charmOrigin` | 0 B/not applicable, 1 confirmed prompt D, 2 D from B, 3 unknown D origin |
| `QHadron` / `QGenCharged.associationQuality` | 0 missing, 1 unique, 2 ambiguous, 3 unknown/incomplete (QHadron V2); for QGenCharged, Stage 1's membership/quality values are taken as saved |
| `QGenCharged.membershipType` | 0 unassigned, 1 heavy-hadron decay product, 2 additional ancestry candidate, 3 ambiguous |
| `QRecoTrack.matchStatus` | 0 no candidate, 1 unopposed match, 2 match with competition, 3 lost one-to-one assignment conflict |

A `membershipType=2` QGenCharged row with unique `quarkIdx` is an ancestry candidate, not proven exclusive physical fragmentation membership. Its `fragmentationSourceCertain` value is preserved; Stage 1 currently leaves it 0 for additional candidates. A reconstructed match does not upgrade that truth certainty.

`QRecoTrack` contains `pt`, `eta`, `phi`, `charge`, `isHighPurity`, and match diagnostics. It does **not** contain dxy, dz, their significances, reference-point coordinates, hit counts, or fit-quality values. The default analysis does not require high purity. `--require-high-purity` is an optional robustness pass using the existing `isHighPurity` bit.

## Environment and commands

Build the package in the local CMSSW_15_1_1 source area:

```bash
cd "$CMSSW_BASE/src"
scram b -j 8
```

Run Stage 1 and Stage 2 on the HZZc MiniAODSIM file you supply. Set the shell variable to its real local path; no input path is assumed here:

```bash
export HZZC_MINIAODSIM='/your/actual/local/path/HZZc_MINIAODSIM.root'
cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_matching_cfg.py \
  inputFiles="file:${HZZC_MINIAODSIM}" \
  maxEvents=-1 outfile=HZZc_stage2.root
```

For multiple files, pass a comma-separated list through the same `inputFiles` option:

```bash
export HZZC_MINIAODSIM_1='/your/actual/local/path/HZZc_part1.root'
export HZZC_MINIAODSIM_2='/your/actual/local/path/HZZc_part2.root'
cmsRun PhysicsTools/SimpleVertexTable/test/quark_reco_track_matching_cfg.py \
  inputFiles="file:${HZZC_MINIAODSIM_1},file:${HZZC_MINIAODSIM_2}" \
  maxEvents=-1 outfile=HZZc_stage2.root
```

The config's default `maxEvents` is 100; set `maxEvents=-1` for all events. It uses the existing unpacker and matching thresholds without changing selections. The output NanoAOD file carries event identifiers and all four required flat tables.

Install the standalone Python dependencies in the environment that will read the output:

```bash
python3 -m pip install -r PhysicsTools/SimpleVertexTable/test/stage3/requirements.txt
```

Run synthetic tests without any HZZc input:

```bash
python3 -m unittest discover \
  -s PhysicsTools/SimpleVertexTable/test/stage3 \
  -p 'test_*.py' -v
```

Run the complete analysis on the output file:

```bash
python3 PhysicsTools/SimpleVertexTable/test/stage3/analyze_hzzc.py \
  --input HZZc_stage2.root --output-dir stage3_hzzc --prefix hzzc
```

Analyze a deterministic entry range or multiple outputs (all entries from the first file followed by the next):

```bash
python3 PhysicsTools/SimpleVertexTable/test/stage3/analyze_hzzc.py \
  --input HZZc_stage2_part1.root HZZc_stage2_part2.root \
  --entry-start 0 --max-events 100 --step-size '25 MB' \
  --output-dir stage3_hzzc_first100 --prefix hzzc
```

For the optional high-purity variation, use a separate output directory/prefix and add `--require-high-purity`. The entry range is global across the ordered input list. Processing is chunked by uproot; `--step-size` accepts a byte-size string or an entry count. Outputs are deterministic, uncompressed CSV files.

## Analysis definitions and outputs

Reco acceptance is inclusive: `pt >= 0.4 GeV`, `abs(eta) <= 2.5`. N is the count of distinct accepted QRecoTrack rows that successfully match a QGenCharged row with `membershipType=1` and the target event-local `hadronIdx`. An additional candidate never contributes to N. CSVs include both generator charged decay-descendant counts and accepted matched reco counts. In the row-level `N_category` convenience field, `N3` means exactly three and `N3+` means four or more; `N_ge3_inclusive_category`, the summary, and neighborhood scan use `N3+` to mean three or more, overlapping the exact-N3 population.

For confirmed prompt D hadrons, additional candidates require a valid unique QHadron QFrag association and QGenCharged `membershipType=2`, unique association quality, and the same QFrag row. Per-D rows are retained, including when multiple D hadrons share one QFrag row; the summary separately deduplicates tracks within quark-system counts. Outputs distinguish accepted generator candidates, matched reco tracks, and matched reco tracks passing reco acceptance.

Neighborhoods use a truth D direction for angular diagnostic columns. For neighborhood scans, N1 is centered on its single reconstructed daughter as an explicitly oracle-only diagnostic; N2/N3+ use the vector-summed accepted decay-track momentum direction. These centers are not proposed as inference-time seeds. The scan covers ΔR<0.2, 0.4, 0.6, 0.8, 1.0, and the full event with no angular restriction. The graph size estimates use reconstructed tracks only: all-pairs upper triangle, fixed ΔR radii, and the union of each node's eight nearest angular neighbors. The ΔR<0.4 compatibility diagnostic reports edge ΔR, absolute log pT ratio, and same/opposite-charge counts as reconstructed-only pair features; it does not apply an arbitrary pT-ratio or charge cut. No jet or truth quantity constructs the candidate graph. Approximate memory estimates count 20 raw bytes/node and 12 bytes/undirected edge and omit Python/framework overhead.

Files written under the selected output directory:

| File | Contents |
|---|---|
| `hzzc_truth_hadrons.csv` | One selected B/D object per row, generator and reconstructed decay multiplicities, origin, and QFrag association quality. |
| `hzzc_fragmentation_tracks.csv` | One matched additional candidate per target prompt D, with reco/gen kinematics, ΔR, source-certainty, match status, and residuals; candidate rows remain per-D when systems are shared. |
| `hzzc_feasibility_summary.csv` | N populations, origin/species/kinematic counts, additional-track counts, Stage 2 matching, label confidence, and graph summary statistics. |
| `hzzc_neighborhood_scan.csv` | Aggregated N-category/radius multiplicities, additional-vs-other tracks, pT/charge diagnostics, and D-track neighborhood coverage. |
| `hzzc_training_feasibility.csv` | Separate feasibility counts and limitations for D-pair reconstruction, quark-system reconstruction, and additional-track context. |
| `hzzc_fragmentation_auxiliary.csv` | Per-D observables, raw values, validity masks, N category, additional accepted matched candidate count, and R<0.4 neighborhood counts. Prompt-D target eligibility requires confirmed prompt origin, association quality 1, valid event-local QFrag index with |pdgId|=4, `preHadKinematicsValid=1`, and finite nonzero momentum. |
| `hzzc_fragmentation_auxiliary_summary.csv` | N1, N2, and inclusive N3+ distributions/outlier rates and correlations, including within-D-pT-bin correlations with additional candidate multiplicity. B→D and unknown-origin D are separate diagnostic populations. |
| `hzzc_graph_event_sizes.csv` | Per-event node and edge estimates for all-pairs, angular-radius, and k-nearest-neighbor graphs. |
| `hzzc_report.md` | Concise sample report and limitations. |

The auxiliary study reconstructs 3-momenta from the saved QHadron (`pt,eta,phi`) and pre-hadronization proxy (`preHadPt,preHadEta,preHadPhi`) branches. It computes `zProj=dot(pD,pQ)/|pQ|^2`, `zT=ptD/preHadPt`, `deltaR_DQuark`, and `kTrel=|pD×pQ|/|pQ|`. Values are never clipped. Tukey 1.5-IQR outlier flags are calculated separately by origin category and observable; `zProj` outside [0,1] and `zT>1` also receive explicit physical-boundary flags. These are diagnostics, not universal hard cuts. Prompt-D rows are potential training-only targets; they are never inference inputs. B→D ratios describe D kinematics relative to the saved parent-bottom lineage and must not be called B fragmentation. Unknown-origin rows are diagnostic only. All four targets depend on the unobserved quark reference, so reconstructed-track kinematics can only predict them statistically; angular neighborhoods and candidate multiplicity are measured as possible proxies. The additional-candidate count and candidate separation are truth-assisted evaluation quantities from QGenCharged membership and must not become inference inputs. The R<0.4 local track multiplicity uses the existing N1 oracle or N2+ summed-decay-track center and is therefore a diagnostic, not a proposed inference-time seed.

Pair labels for the D-pair objective require a successful **unopposed** Stage 2 match and known `membershipType=1` ancestry to a selected QHadron; pairs involving unresolved/unmatched tracks are not silently counted as negatives. Quark-system labels have a stricter additional requirement: unique QHadron and QGenCharged quark ancestry, valid QFrag index, and `fragmentationSourceCertain=1`. Unique ancestry candidates with source certainty 0 remain uncertain context. The study does not decide which learning objective is best until the HZZc measurements are available.
