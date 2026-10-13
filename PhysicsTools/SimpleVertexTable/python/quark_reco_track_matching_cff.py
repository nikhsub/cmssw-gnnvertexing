import FWCore.ParameterSet.Config as cms
from PhysicsTools.SimpleVertexTable.sv_reco_cff import (
    trackGenMatch as baseTrackGenMatch,
    trackTable as baseTrackTable,
    trackVertexVars as baseTrackVertexVars,
    custom_train_features,
)
from PhysicsTools.SimpleVertexTable.mergedGenParticles import mergedGenParticles as mergedGenParticlesConfig
from PhysicsTools.SimpleVertexTable.genVertices import genVertexTable as genVertexTableConfig


# Keep this configuration identical to unpackedTracksAndVertices in
# sv_reco_cff.py. That producer is the graph-vertex analysis input collection:
# unpacked charged PF candidates plus lost tracks with the package's 0.4 GeV
# unpacker threshold. It produces the reco::TrackCollection consumed below.
unpackedTracksAndVertices = cms.EDProducer(
    "PATTrackAndVertexUnpacker",
    slimmedVertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    slimmedSecondaryVertices=cms.InputTag("slimmedSecondaryVertices"),
    additionalTracks=cms.InputTag("lostTracks"),
    packedCandidates=cms.InputTag("packedPFCandidates"),
    minTrackPt=cms.double(0.4),
)

quarkRecoTrackMatchProducer = cms.EDProducer(
    "QuarkRecoTrackMatchProducer",
    tracks=cms.InputTag("unpackedTracksAndVertices"),
    genChargedTable=cms.InputTag("quarkFragmentationProducer", "QGenChargedTable"),
    maxDeltaR=cms.double(0.03),
    maxRelPt=cms.double(0.5),
)

# Reuse the original AssociationGNN input producers and feature definitions.
qfragTrackVertexVars = baseTrackVertexVars.clone(
    tracks=cms.InputTag("unpackedTracksAndVertices"),
    primaryVertices=cms.InputTag("unpackedTracksAndVertices"),
)
qfragTrackGenMatch = baseTrackGenMatch.clone(
    tracks=cms.InputTag("unpackedTracksAndVertices"),
    genParticles=cms.InputTag("mergedGenParticles"),
)
qfragTrackTable = baseTrackTable.clone(src=cms.InputTag("unpackedTracksAndVertices"))
qfragMergedGenParticles = mergedGenParticlesConfig.clone()
qfragGenVertexTable = genVertexTableConfig.clone(genParticles=cms.InputTag("mergedGenParticles"))

quarkAssociationGNNTableProducer = cms.EDProducer(
    "QuarkAssociationGNNTableProducer",
    recoTrackTable=cms.InputTag("quarkRecoTrackMatchProducer"),
    genChargedTable=cms.InputTag("quarkFragmentationProducer", "QGenChargedTable"),
    hadronTable=cms.InputTag("quarkFragmentationProducer", "QHadronTable"),
    pairTable=cms.InputTag("trackGVFeatures", "trackPair"),
)


def custom_quark_reco_track_matching(process):
    """Schedule legacy AssociationGNN features and Stage-2 QFrag associations.

    The caller schedules Stage 1 first. Historical track/pair calculations and
    labels come directly from TrackGVFeaturesProducer; QFrag truth is appended
    by QuarkAssociationGNNTableProducer without changing Stage 2 assignments.
    """
    process.unpackedTracksAndVertices = unpackedTracksAndVertices.clone()
    process.mergedGenParticles = qfragMergedGenParticles.clone()
    process.trackVertexVars = qfragTrackVertexVars.clone()
    process.trackGenMatch = qfragTrackGenMatch.clone()
    process.trackTable = qfragTrackTable.clone()
    process.genVertexTable = qfragGenVertexTable.clone()
    process = custom_train_features(process, gvProducerLabel="genVertexTable")
    # TrackGVFeaturesProducer declares this legacy token but does not read it;
    # bind it to an available reco::VertexCollection in this training-only path.
    process.trackGVFeatures.secondaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
    process.quarkRecoTrackMatchProducer = quarkRecoTrackMatchProducer.clone()
    process.quarkAssociationGNNTableProducer = quarkAssociationGNNTableProducer.clone()
    process.quark_reco_track_matching_sequence = cms.Sequence(
        process.mergedGenParticles * process.unpackedTracksAndVertices * process.trackVertexVars *
        process.trackGenMatch * process.trackTable * process.genVertexTable * process.trackGVFeatures *
        process.quarkRecoTrackMatchProducer * process.quarkAssociationGNNTableProducer
    )
    return process
