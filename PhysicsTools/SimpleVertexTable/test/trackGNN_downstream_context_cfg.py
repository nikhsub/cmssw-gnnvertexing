"""Produce frozen GNN top-K seeds plus jet-independent PF/PV/SV event context."""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.maxEvents = 100
options.register("outfile", "trackGNN_downstream_context.root",
                 VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Output NanoAOD ROOT file")
options.register("isMC", True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Input is simulated MC")
options.register("maxSeeds", 128, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Maximum number of selected physical pairs fitted per event")
options.register("saveEmbeddings", True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Write the 128D GNN track embedding table")
options.register("minChargedPt", 0.1, VarParsing.multiplicity.singleton,
                 VarParsing.varType.float, "Minimum charged PF candidate pT [GeV]")
options.register("minNeutralPt", 0.1, VarParsing.multiplicity.singleton,
                 VarParsing.varType.float, "Minimum neutral PF candidate pT [GeV]")
options.register("maxPFAbsEta", 3.0, VarParsing.multiplicity.singleton,
                 VarParsing.varType.float, "Maximum absolute PF candidate eta")
options.parseArguments()

default_input = (
    "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/"
    "HPlusCharm_3FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_"
    "TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/"
    "106X_upgrade2018_realistic_v16_L1v1-v2/2810000/"
    "22187118-BB4F-074E-B496-C14403E0BAED.root"
)
input_files = options.inputFiles if options.inputFiles else [default_input]

process = cms.Process("GNNCONTEXT")
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if options.isMC:
    process.GlobalTag.globaltag = cms.string("106X_upgrade2018_realistic_v16_L1v1")
else:
    process.GlobalTag.globaltag = cms.string("106X_dataRun2_v35")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(*input_files))

from PhysicsTools.SimpleVertexTable.sv_reco_cff import (
    custom_sv_tracks_training,
    trackPairFeatures,
    trackGNNSeederInference,
    trackGNNSeedProducer,
    trackGNNPFContext,
    trackGNNVertexContext,
    svTable,
)

if options.isMC:
    from PhysicsTools.SimpleVertexTable.mergedGenParticles import customMergedGenParticles
    process.customMergedGenParticles = customMergedGenParticles

# This installs the standard unpacking, filtering, TrackVertexVars, track table,
# and MC matching module definitions. We explicitly schedule them below so seed
# production and context extraction remain independent of generator truth.
process = custom_sv_tracks_training(process, isMC=options.isMC)

# Preserve the exact validated AssociationGNN graph cuts and disable legacy
# all-pair fitting. Only TrackGNNSeedProducer fits its selected top-K pairs.
process.trackPairFeatures = trackPairFeatures.clone(
    maxDR=cms.double(1.5),
    maxDEta=cms.double(9999.0),
    maxDPhi=cms.double(9999.0),
    maxPairEta=cms.double(9999.0),
    maxInvMass=cms.double(5.0),
    maxDca=cms.double(0.1),
    enableTwoTrackVertexFit=cms.bool(False),
    enableVertexFitTiming=cms.bool(False),
)
process.trackGNNSeederInference = trackGNNSeederInference.clone(
    saveEmbeddings=cms.bool(options.saveEmbeddings),
)
process.trackGNNSeedProducer = trackGNNSeedProducer.clone(maxSeeds=cms.uint32(options.maxSeeds))
process.trackGNNPFContext = trackGNNPFContext.clone(
    minChargedPt=cms.double(options.minChargedPt),
    minNeutralPt=cms.double(options.minNeutralPt),
    maxPFAbsEta=cms.double(options.maxPFAbsEta),
)
process.trackGNNVertexContext = trackGNNVertexContext.clone()

# This is the established slimmedSecondaryVertices table with the same 2-sigma
# displacement selection and source-daughter mappings used by the earlier
# validation workflow. CtxSV separately retains the complete source collection.
process.standardSVTable = svTable.clone(
    src=cms.InputTag("unpackedTracksAndVertices", "secondary"),
    pvSrc=cms.InputTag("unpackedTracksAndVertices"),
    filteredTrackIndex=cms.InputTag("filteredTracks", "filteredIndex"),
    sourceSecondaryVertices=cms.InputTag("slimmedSecondaryVertices"),
    unpackedTracks=cms.InputTag("unpackedTracksAndVertices"),
    packedTrackIndex=cms.InputTag("unpackedTracksAndVertices", "packedTrackIndex"),
    lostTrackIndex=cms.InputTag("unpackedTracksAndVertices", "lostTrackIndex"),
    dlenSigMin=cms.double(2.0),
)

if options.isMC:
    from PhysicsTools.SimpleVertexTable.genVertices import custom_GV_producer
    process = custom_GV_producer(process, collection="track")
    from PhysicsTools.SimpleVertexTable.sv_reco_cff import trackPairLabels
    process.trackPairLabels = trackPairLabels.clone()
    process.gvMatchTable.pvSrc = cms.InputTag("unpackedTracksAndVertices")
    process.gvMatchTable.secondaryVertices = cms.InputTag("unpackedTracksAndVertices", "secondary")
    process.gvMatchTable.sourceSecondaryVertices = cms.InputTag("slimmedSecondaryVertices")
    process.gvMatchTable.tracks = cms.InputTag("filteredTracks")
    process.gvMatchTable.dlenSigMin = cms.double(2.0)

modules = (
    process.unpackedTracksAndVertices
    * process.filteredTracks
    * process.trackVertexVars
    * process.trackPairFeatures
    * process.trackGNNSeederInference
    * process.trackGNNSeedProducer
    * process.trackGNNPFContext
    * process.trackGNNVertexContext
    * process.standardSVTable
)
if options.isMC:
    modules = (
        process.customMergedGenParticles
        * modules
        * process.trackGenMatch
        * process.trackTable
        * process.genVertexTable
        * process.trackPairLabels
        * process.gvMatchTable
    )
else:
    modules = modules * process.trackTable
process.contextPath = cms.Path(modules)

output_commands = [
    "drop *",
    "keep nanoaodFlatTable_trackTable_*_*",
    "keep nanoaodFlatTable_trackPairFeatures_*_*",
    "keep nanoaodFlatTable_trackGNNSeederInference_*_*",
    "keep nanoaodFlatTable_trackGNNSeedProducer_*_*",
    "keep nanoaodFlatTable_trackGNNPFContext_*_*",
    "keep nanoaodFlatTable_trackGNNVertexContext_*_*",
    "keep nanoaodFlatTable_standardSVTable_*_*",
]
if options.isMC:
    output_commands += [
        "keep nanoaodFlatTable_genVertexTable_*_*",
        "keep nanoaodFlatTable_trackPairLabels_track_*",
        "keep nanoaodFlatTable_trackPairLabels_trackPair_*",
        "keep nanoaodFlatTable_gvMatchTable_*_*",
    ]

process.nanoOut = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string(options.outfile),
    outputCommands=cms.untracked.vstring(*output_commands),
)
process.end = cms.EndPath(process.nanoOut)
