"""Run frozen AssociationGNN inference and write its full graph outputs.

This config is independent of IVF/modelInference. It reuses the simplelabel
track unpacking/filtering and pair feature producers, retaining the ordered
training graph (including both edge directions) without fits or truth labels.
"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.maxEvents = 100
options.register("outfile", "trackGNN_seeder_validation.root",
                 VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Output NanoAOD ROOT file")
options.register("isMC", True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Input is simulated MC")
options.register("saveEmbeddings", True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Write optional 128D node embeddings")
options.parseArguments()

# A known MiniAODSIM file used by this package. Override with inputFiles=... on
# the command line; data mode uses the same reconstruction path without gen truth.
default_input = (
    "root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/"
    "HPlusCharm_3FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_"
    "TuneCP5_13TeV-amcatnlo-pythia8/MINIAODSIM/"
    "106X_upgrade2018_realistic_v16_L1v1-v2/2810000/"
    "22187118-BB4F-074E-B496-C14403E0BAED.root"
)
input_files = options.inputFiles if options.inputFiles else [default_input]

process = cms.Process("SEEDERVALIDATE")
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
)

if options.isMC:
    from PhysicsTools.SimpleVertexTable.mergedGenParticles import customMergedGenParticles
    process.customMergedGenParticles = customMergedGenParticles

process = custom_sv_tracks_training(process, isMC=options.isMC)

# These are the exact graph cuts used for the frozen model's training input.
# Do not add cuts or collapse reciprocal directed edges.
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

modules = process.sv_track * process.trackPairFeatures * process.trackGNNSeederInference
if options.isMC:
    process.seederPath = cms.Path(process.customMergedGenParticles * modules)
else:
    process.seederPath = cms.Path(modules)

process.nanoOut = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string(options.outfile),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep nanoaodFlatTable_trackTable_*_*",
        "keep nanoaodFlatTable_trackPairFeatures_*_*",
        "keep nanoaodFlatTable_trackGNNSeederInference_*_*",
    ),
)
process.end = cms.EndPath(process.nanoOut)
