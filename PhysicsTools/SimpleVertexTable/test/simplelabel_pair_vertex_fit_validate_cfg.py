"""Validation NanoAOD with simplelabel and standard MiniAOD IVF SVs.

This is the validation counterpart to simplelabel_pair_vertex_fit_cfg.py. It
keeps the simplelabel track/pair truth and top-ranked seed-fit outputs, and adds SV
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
options.register(
    "saveEmbeddings",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Write the optional 128-component seeder embedding table",
)
options.register(
    "maxSeeds",
    128,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Maximum number of top-ranked physical pairs to fit per event",
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
    trackGNNSeederInference,
    trackGNNSeedProducer,
    svTable,
)

process = custom_sv_tracks_training(process, isMC=True)

# Simplelabel generator-vertex truth, with the same module labels as training.
from PhysicsTools.SimpleVertexTable.genVertices import custom_GV_producer
process = custom_GV_producer(process, collection="track")

process.trackPairFeatures = trackPairFeatures.clone(
    enableTwoTrackVertexFit=cms.bool(False),
    enableVertexFitTiming=cms.bool(False),
)
# Inference consumes the same directed trackPair table, including its feature
# columns and pair-row order. The seed producer performs fits after ranking.
process.trackGNNSeederInference = trackGNNSeederInference.clone(
    saveEmbeddings=cms.bool(options.saveEmbeddings),
)
process.trackGNNSeedProducer = trackGNNSeedProducer.clone(
    maxSeeds=cms.uint32(options.maxSeeds),
)
process.trackPairLabels = trackPairLabels.clone()

# Save the standard slimmed IVF SVs and their unpacked track constituents.
# Match the 2-sigma SV selection in GVMatchProducer so its extension rows have
# the same order and row count as SVTableProducer's SV table.
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
process.gvMatchTable.pvSrc = cms.InputTag("unpackedTracksAndVertices")
process.gvMatchTable.secondaryVertices = cms.InputTag(
    "unpackedTracksAndVertices", "secondary"
)
process.gvMatchTable.sourceSecondaryVertices = cms.InputTag("slimmedSecondaryVertices")
process.gvMatchTable.tracks = cms.InputTag("filteredTracks")
process.gvMatchTable.dlenSigMin = cms.double(2.0)

process.validationPath = cms.Path(
    process.customMergedGenParticles
    * process.unpackedTracksAndVertices
    * process.filteredTracks
    * process.trackVertexVars
    * process.trackPairFeatures
    * process.trackGNNSeederInference
    * process.trackGNNSeedProducer
    # MC truth matching and table outputs follow the truth-independent ranking/fits.
    * process.trackGenMatch
    * process.trackTable
    * process.genVertexTable
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
        # Original directed pairs and GNN predictions, without all-pair fits.
        "keep nanoaodFlatTable_trackPairFeatures_*_*",
        # GNN node outputs, optional embeddings, and the same-row-order
        # trackPair score extension (merged with the original trackPair table).
        "keep nanoaodFlatTable_trackGNNSeederInference_*_*",
        "keep nanoaodFlatTable_trackPairLabels_trackPair_*",
        # Top-K unique pairs and fits performed only for selected seeds.
        "keep nanoaodFlatTable_trackGNNSeedProducer_*_*",
        # Standard slimmed IVF SV properties, PV properties, and constituents.
        "keep nanoaodFlatTable_standardSVTable_*_*",
        # GV-SV matches, matched SV constituents, and daughter-to-track matches.
        "keep nanoaodFlatTable_gvMatchTable_*_*",
    ),
)
process.end = cms.EndPath(process.nanoOut)
