import FWCore.ParameterSet.Config as cms

quarkFragmentationProducer = cms.EDProducer(
    "QuarkFragmentationProducer",
    prunedGenParticles=cms.InputTag("prunedGenParticles"),
    packedGenParticles=cms.InputTag("packedGenParticles"),
)

def custom_quark_fragmentation_training(process):
    """Add the generator-only Stage-1 tables without GV or track producers."""
    process.quarkFragmentationProducer = quarkFragmentationProducer.clone()
    process.quark_fragmentation_sequence = cms.Sequence(process.quarkFragmentationProducer)
    return process
