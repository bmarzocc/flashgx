
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("MicroAODAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff") # gives deprecated message in 80X but still runs
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'94X_mc2017_realistic_v14','')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 100 )
                                                                       
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),                       
    fileNames = cms.untracked.vstring(
    LISTOFFILES)
)

process.load('NAMECfi')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('OUTPUTNAME.root')
)

import HLTrigger.HLTfilters.hltHighLevel_cfi 
process.HLTfilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    HLTPaths = ["HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_v*",
                "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_v*",
                "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_v*",
                #"HLT_DiJet110_35_Mjj650_PFMET110_v*",
                #"HLT_DiJet110_35_Mjj650_PFMET120_v*",
                #"HLT_DiJet110_35_Mjj650_PFMET130_v*"
               ],
    #throw = True,
    andOr = cms.bool(True)
)

process.p = cms.Path(
    process.HLTfilter*
    process.flashggUnpackedJets*
    process.flashgxanalysis
)
