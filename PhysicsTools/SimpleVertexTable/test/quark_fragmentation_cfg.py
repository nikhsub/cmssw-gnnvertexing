import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.register("outfile", "qfrag_stage1.root", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Output EDM ROOT file")
options.register("maxEvents", 100, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Maximum events")
options.register("debugAncestry", False, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Print selected retained generator mother trees")
options.register("debugOnlyAmbiguous", True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Restrict ancestry logs to ambiguous quark associations")
options.register("debugOnlyD", True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Restrict ancestry logs to D hadrons")
options.register("debugV2OldUniqueToV1Unknown", False, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Focus logs on old unique to V1 unknown")
options.register("debugV2LineageReason6", False, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Focus logs on V1 branching-copy reason 6")
options.register("debugV2LineageReason10", False, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Focus logs on V1 unresolved-origin reason 10")
options.register("debugV2InvalidKinematics", False, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Focus logs on V2 invalid momentum references")
options.register("debugMaxEvents", 10, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Maximum processed events eligible for ancestry logs")
options.register("debugMaxHadronsPerEvent", 10, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Maximum logged hadrons per event")
options.register("debugMaxDepth", 30, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Maximum mother-tree depth printed")
options.parseArguments()

process = cms.Process("QFRAG")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.append("QFragAncestry")
process.MessageLogger.cerr.threshold = "INFO"
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

from PhysicsTools.SimpleVertexTable.quark_fragmentation_cff import custom_quark_fragmentation_training
process = custom_quark_fragmentation_training(process)
process.quarkFragmentationProducer.debugAncestry = cms.bool(options.debugAncestry)
process.quarkFragmentationProducer.debugOnlyAmbiguous = cms.bool(options.debugOnlyAmbiguous)
process.quarkFragmentationProducer.debugOnlyD = cms.bool(options.debugOnlyD)
process.quarkFragmentationProducer.debugV2OldUniqueToV1Unknown = cms.bool(options.debugV2OldUniqueToV1Unknown)
process.quarkFragmentationProducer.debugV2LineageReason6 = cms.bool(options.debugV2LineageReason6)
process.quarkFragmentationProducer.debugV2LineageReason10 = cms.bool(options.debugV2LineageReason10)
process.quarkFragmentationProducer.debugV2InvalidKinematics = cms.bool(options.debugV2InvalidKinematics)
process.quarkFragmentationProducer.debugMaxEvents = cms.uint32(options.debugMaxEvents)
process.quarkFragmentationProducer.debugMaxHadronsPerEvent = cms.uint32(options.debugMaxHadronsPerEvent)
process.quarkFragmentationProducer.debugMaxDepth = cms.uint32(options.debugMaxDepth)
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
