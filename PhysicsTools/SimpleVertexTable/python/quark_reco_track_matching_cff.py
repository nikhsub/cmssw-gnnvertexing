import FWCore.ParameterSet.Config as cms


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


def custom_quark_reco_track_matching(process):
    """Schedule the existing unpacked tracks and Stage 2 table producer.

    The caller must schedule the Stage 1 QuarkFragmentationProducer before
    this sequence so that its QGenChargedTable product is available.
    """
    process.unpackedTracksAndVertices = unpackedTracksAndVertices.clone()
    process.quarkRecoTrackMatchProducer = quarkRecoTrackMatchProducer.clone()
    process.quark_reco_track_matching_sequence = cms.Sequence(
        process.unpackedTracksAndVertices * process.quarkRecoTrackMatchProducer
    )
    return process
