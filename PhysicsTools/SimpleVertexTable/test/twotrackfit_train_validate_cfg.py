"""Dedicated minimal production for prompt-D vertex and momentum studies.

Run with mode=train for compact QTrain tables or mode=validate to also save
the full legacy AssociationGNN validation tables. Both modes retain the same
reconstructed tracks and trackPair rows, including optional two-track fits.
This configuration is for MC MiniAOD with prunedGenParticles and
packedGenParticles available.
"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.register("mode", "train", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Output mode: train or validate")
options.register("outfile", "promptd_train.root", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Output NanoAOD ROOT file")
options.register("maxEvents", 100, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Maximum number of events")
options.register("compressionAlgorithm", "ZSTD", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "ROOT compression algorithm")
options.register("compressionLevel", 5, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "ROOT compression level")
options.parseArguments()

if options.mode not in ("train", "validate"):
    raise ValueError("mode must be 'train' or 'validate'")

compact = options.mode == "train"
process = cms.Process("PROMPTD")
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("106X_upgrade2018_realistic_v16_L1v1")
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

from PhysicsTools.SimpleVertexTable.quark_fragmentation_cff import custom_quark_fragmentation_training
from PhysicsTools.SimpleVertexTable.quark_reco_track_matching_cff import custom_quark_reco_track_matching

# Keep the selected-hadron decay geometry and minimal QFrag p4/origin truth;
# omit pre-hadronization proxy and fragmentation-regression columns.
process = custom_quark_fragmentation_training(
    process,
    saveDecayVertexTruth=True,
    saveFragmentationObservables=False,
    minimalQuarkTruth=True,
    saveDerivedDecayTruth=False,
)
process = custom_quark_reco_track_matching(
    process,
    compact=compact,
    enableTwoTrackVertexFit=True,
    enableVertexFitTiming=True,
    enableTruthMatchedPairFit=True,
    saveFragmentationObservables=False,
)
process.promptDPath = cms.Path(
    process.quark_fragmentation_sequence * process.quark_reco_track_matching_sequence
)

keep = [
    "drop *",
    # Shared reconstructed inputs: the historical 22 node features, all
    # original selected edges and their 14 pair features, with fit columns.
    "keep nanoaodFlatTable_trackTable_*_*",
    "keep nanoaodFlatTable_trackGVFeatures_*_*",
    # Generator matching provenance and decay truth.
    "keep nanoaodFlatTable_quarkFragmentationProducer_QFragTable_*",
    "keep nanoaodFlatTable_quarkFragmentationProducer_QHadronTable_*",
    "keep nanoaodFlatTable_quarkFragmentationProducer_QGenChargedTable_*",
    "keep nanoaodFlatTable_quarkRecoTrackMatchProducer_*_*",
    # Track/hadron/quark relationships and per-hadron decay geometry.
    "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QTrainTrack_*",
    "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QTrainPair_*",
    "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QTrainHadron_*",
    "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QTrainQuark_*",
    "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QTrainLegacyHadron_*",
    # All unique truth-associated same-hadron pairs and the event/PV summary.
    "keep nanoaodFlatTable_truthMatchedPairFitProducer_QTruthMatchedTrack_*",
    "keep nanoaodFlatTable_truthMatchedPairFitProducer_QTruthMatchedHadron_*",
    "keep nanoaodFlatTable_truthMatchedPairFitProducer_QTruthMatchedPairFit_*",
    "keep nanoaodFlatTable_truthMatchedPairFitProducer_QTruthPairFitSummary_*",
]
if not compact:
    # Full legacy validation labels and cross-check tables are useful in
    # validate mode; these are omitted from compact training output.
    keep.extend([
        "keep nanoaodFlatTable_genVertexTable_GVTable_*",
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QGNNTrack_*",
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QGNNPairTruth_*",
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QGNNHadron_*",
    ])

process.out = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string(options.outfile),
    outputCommands=cms.untracked.vstring(*keep),
    compressionAlgorithm=cms.untracked.string(options.compressionAlgorithm),
    compressionLevel=cms.untracked.int32(options.compressionLevel),
)
process.end = cms.EndPath(process.out)
