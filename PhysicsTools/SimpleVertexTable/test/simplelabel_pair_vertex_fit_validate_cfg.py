"""Validation NanoAOD with simplelabel and standard MiniAOD IVF SVs.

This is the validation counterpart to simplelabel_pair_vertex_fit_cfg.py. It
keeps the simplelabel track/pair truth and pair-fit outputs, and adds SV
properties, constituent tracks, and GV matching for the standard
slimmedSecondaryVertices collection. It does not rerun the score-based IVF.
"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.maxEvents = 100
options.register(
    "outfile",
    "simplelabel_pair_vertex_fit_validation.root",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Validation NanoAOD ROOT file",
)
options.parseArguments()

default_input = (
    "root://cms-xrd-global.cern.ch//store/mc/"
    "RunIISummer20UL18MiniAODv2/"
    "HPlusCharm_3FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_"
    "TuneCP5_13TeV-amcatnlo-pythia8/"
    "MINIAODSIM/"
    "106X_upgrade2018_realistic_v16_L1v1-v2/"
    "2810000/"
    "22187118-BB4F-074E-B496-C14403E0BAED.root"
)
input_files = options.inputFiles if options.inputFiles else [default_input]

process = cms.Process("SIMPLELABELVALIDATE")
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("106X_upgrade2018_realistic_v16_L1v1")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*input_files),
)

# The legacy track matcher and GV producer both consume this merged gen view.
from PhysicsTools.SimpleVertexTable.mergedGenParticles import customMergedGenParticles
process.customMergedGenParticles = customMergedGenParticles

# Use the same simplelabel track/table sequence as the training config. The
# unpacker also converts slimmedSecondaryVertices into a reco::VertexCollection
# whose daughters reference the unpacked track collection.
from PhysicsTools.SimpleVertexTable.sv_reco_cff import (
    custom_sv_tracks_training,
    trackPairFeatures,
    trackPairLabels,
    svTable,
)

process = custom_sv_tracks_training(process, isMC=True)

# Simplelabel generator-vertex truth, with the same module labels as training.
from PhysicsTools.SimpleVertexTable.genVertices import custom_GV_producer
process = custom_GV_producer(process, collection="track")

process.trackPairFeatures = trackPairFeatures.clone(
    enableTwoTrackVertexFit=cms.bool(True),
    enableVertexFitTiming=cms.bool(False),
)
process.trackPairLabels = trackPairLabels.clone()

# Save the standard slimmed IVF SVs and their unpacked track constituents.
# Match the 2-sigma SV selection in GVMatchProducer so its extension rows have
# the same order and row count as SVTableProducer's SV table.
process.standardSVTable = svTable.clone(
    src=cms.InputTag("unpackedTracksAndVertices", "secondary"),
    pvSrc=cms.InputTag("unpackedTracksAndVertices"),
    filteredTrackIndex=cms.InputTag("filteredTracks", "filteredIndex"),
    dlenSigMin=cms.double(2.0),
)
process.gvMatchTable.pvSrc = cms.InputTag("unpackedTracksAndVertices")
process.gvMatchTable.secondaryVertices = cms.InputTag(
    "unpackedTracksAndVertices", "secondary"
)
process.gvMatchTable.tracks = cms.InputTag("filteredTracks")
process.gvMatchTable.dlenSigMin = cms.double(2.0)

process.validationPath = cms.Path(
    process.customMergedGenParticles
    * process.sv_track
    * process.genVertexTable
    * process.trackPairFeatures
    * process.trackPairLabels
    * process.standardSVTable
    * process.gvMatchTable
)

process.nanoOut = cms.OutputModule(
    "NanoAODOutputModule",
    fileName=cms.untracked.string(options.outfile),
    outputCommands=cms.untracked.vstring(
        "drop *",
        # Generator hadrons and charged daughters used by simplelabel truth.
        "keep nanoaodFlatTable_genVertexTable_*_*",
        # Reconstructed track features and simplelabel truth labels.
        "keep nanoaodFlatTable_trackTable_*_*",
        "keep nanoaodFlatTable_trackPairLabels_track_*",
        # Selected ordered pairs, fit observables, and same-order edge labels.
        "keep nanoaodFlatTable_trackPairFeatures_*_*",
        "keep nanoaodFlatTable_trackPairLabels_trackPair_*",
        # Standard slimmed IVF SV properties, PV properties, and constituents.
        "keep nanoaodFlatTable_standardSVTable_*_*",
        # GV-SV matches, matched SV constituents, and daughter-to-track matches.
        "keep nanoaodFlatTable_gvMatchTable_*_*",
    ),
)
process.end = cms.EndPath(process.nanoOut)
