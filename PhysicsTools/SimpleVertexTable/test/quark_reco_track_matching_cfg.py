import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.register("outfile", "qfrag_association_gnn.root", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Output NanoAOD ROOT file")
options.register("maxEvents", 100, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Maximum events")
options.parseArguments()

process = cms.Process("QFRAGTRACK")
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

from PhysicsTools.SimpleVertexTable.quark_fragmentation_cff import custom_quark_fragmentation_training
from PhysicsTools.SimpleVertexTable.quark_reco_track_matching_cff import custom_quark_reco_track_matching

process = custom_quark_fragmentation_training(process)
process = custom_quark_reco_track_matching(process)
process.qfragTrackPath = cms.Path(
    process.quark_fragmentation_sequence * process.quark_reco_track_matching_sequence
)

process.out = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string(options.outfile),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep nanoaodFlatTable_quarkFragmentationProducer_QFragTable_*",
        "keep nanoaodFlatTable_quarkFragmentationProducer_QHadronTable_*",
        "keep nanoaodFlatTable_quarkFragmentationProducer_QGenChargedTable_*",
        "keep nanoaodFlatTable_quarkRecoTrackMatchProducer_*_*",
        "keep nanoaodFlatTable_trackTable_*_*",
        "keep nanoaodFlatTable_trackGVFeatures_*_*",
        "keep nanoaodFlatTable_genVertexTable_GVTable_*",
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QGNNTrack_*",
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QGNNPairTruth_*",
        "keep nanoaodFlatTable_quarkAssociationGNNTableProducer_QGNNHadron_*",
    ),
)
process.end = cms.EndPath(process.out)
