"""Compact ROOT production for AssociationGNN and QFrag training targets.

The unpacker, historical track/pair feature producer, Stage-2 QFrag matcher,
and target builder are shared with quark_reco_track_matching_cfg.py. This
compact path omits the unused legacy genPartIdx matcher; no reconstructed input
feature or pair is recomputed or filtered.
"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.register("outfile", "qfrag_association_gnn_training.root", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Compact training ROOT file")
options.register("maxEvents", 100, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Maximum events")
options.register("compressionAlgorithm", "ZSTD", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "ROOT lossless compression algorithm")
options.register("compressionLevel", 5, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "ROOT compression level")
options.parseArguments()

process = cms.Process("QFRAGTRAIN")
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

from PhysicsTools.SimpleVertexTable.quark_fragmentation_cff import custom_quark_fragmentation_training
from PhysicsTools.SimpleVertexTable.quark_reco_track_matching_cff import custom_quark_reco_track_matching

process = custom_quark_fragmentation_training(process)
process = custom_quark_reco_track_matching(process, compact=True)
process.qfragTrackPath = cms.Path(
    process.quark_fragmentation_sequence * process.quark_reco_track_matching_sequence
)

process.out = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string(options.outfile),
    outputCommands=cms.untracked.vstring(
        "drop *",
        # Original 22 node features and legacy hadidx/label extension.
        "keep nanoaodFlatTable_trackTable_*_*",
        "keep nanoaodFlatTable_trackGVFeatures_*_*",
        # Original pair features, edge endpoints, and legacy edgeLabel.
        # The corresponding QTrainPair truth rows use this exact row order.
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QTrainTrack_*",
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QTrainPair_*",
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QTrainHadron_*",
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QTrainQuark_*",
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QTrainLegacyHadron_*",
    ),
    compressionAlgorithm=cms.untracked.string(options.compressionAlgorithm),
    compressionLevel=cms.untracked.int32(options.compressionLevel),
)
process.end = cms.EndPath(process.out)
