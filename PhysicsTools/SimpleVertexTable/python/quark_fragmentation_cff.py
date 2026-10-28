import FWCore.ParameterSet.Config as cms

quarkFragmentationProducer = cms.EDProducer(
    "QuarkFragmentationProducer",
    prunedGenParticles=cms.InputTag("prunedGenParticles"),
    packedGenParticles=cms.InputTag("packedGenParticles"),
)

def custom_quark_fragmentation_training(process, saveDecayVertexTruth=False,
                                       saveFragmentationObservables=True,
                                       minimalQuarkTruth=False,
                                       saveDerivedDecayTruth=True):
    """Add the generator-only Stage-1 tables without GV or track producers."""
    process.quarkFragmentationProducer = quarkFragmentationProducer.clone(
        saveDecayVertexTruth=cms.untracked.bool(saveDecayVertexTruth),
        saveFragmentationObservables=cms.untracked.bool(saveFragmentationObservables),
        minimalQuarkTruth=cms.untracked.bool(minimalQuarkTruth),
        saveDerivedDecayTruth=cms.untracked.bool(saveDerivedDecayTruth),
    )
    process.quark_fragmentation_sequence = cms.Sequence(process.quarkFragmentationProducer)
    return process
