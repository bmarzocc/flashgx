
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
    '/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd//W1JetsToLNu_LHEWpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180606_172207/0000/myMicroAODOutputFile_220.root')
)

process.load('flashgx.MicroAODAnalysis.FlashGXAnalysis_cfi')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('FlashGXAnalysis_VBFH.root')
)

process.p = cms.Path(
    process.flashggUnpackedJets*
    process.flashgxanalysis
)
