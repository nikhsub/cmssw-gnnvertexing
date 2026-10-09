import FWCore.ParameterSet.Config as cms

quarkFragmentationProducer = cms.EDProducer(
    "QuarkFragmentationProducer",
    prunedGenParticles=cms.InputTag("prunedGenParticles"),
    packedGenParticles=cms.InputTag("packedGenParticles"),
    debugAncestry=cms.bool(False),
    debugOnlyAmbiguous=cms.bool(True),
    debugOnlyD=cms.bool(True),
    debugMaxEvents=cms.uint32(10),
    debugMaxHadronsPerEvent=cms.uint32(10),
    debugMaxDepth=cms.uint32(30),
)

def custom_quark_fragmentation_training(process):
    """Add the generator-only Stage-1 tables without GV or track producers."""
    process.quarkFragmentationProducer = quarkFragmentationProducer.clone()
    process.quark_fragmentation_sequence = cms.Sequence(process.quarkFragmentationProducer)
    return process
