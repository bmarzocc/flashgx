import FWCore.ParameterSet.Config as cms

flashgxtruepileup = cms.EDAnalyzer("FlashGXTruePileup",

    pileupInfo                        = cms.InputTag("slimmedAddPileupInfo")

)
