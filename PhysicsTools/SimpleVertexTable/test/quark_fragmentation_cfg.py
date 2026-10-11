import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.register("outfile", "qfrag_stage1.root", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Output EDM ROOT file")
options.register("maxEvents", 100, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Maximum events")
options.parseArguments()

process = cms.Process("QFRAG")
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

from PhysicsTools.SimpleVertexTable.quark_fragmentation_cff import custom_quark_fragmentation_training
process = custom_quark_fragmentation_training(process)
process.qfragPath = cms.Path(process.quark_fragmentation_sequence)

process.out = cms.OutputModule("NanoAODOutputModule",
    fileName=cms.untracked.string(options.outfile),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep nanoaodFlatTable_quarkFragmentationProducer_QFragTable_*",
        "keep nanoaodFlatTable_quarkFragmentationProducer_QHadronTable_*",
        "keep nanoaodFlatTable_quarkFragmentationProducer_QGenChargedTable_*",
    ),
)
process.end = cms.EndPath(process.out)
