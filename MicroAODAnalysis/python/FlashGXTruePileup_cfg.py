
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
    '/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_1.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_10.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_100.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_11.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_12.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_13.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_14.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_15.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_16.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_17.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_18.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_19.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_2.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_20.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_21.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_22.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_23.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_24.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_25.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_26.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_27.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_28.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_29.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_3.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_30.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_31.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_32.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_33.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_34.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_35.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_36.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_37.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_38.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_39.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_4.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_40.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_41.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_42.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_43.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_44.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_45.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_46.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_47.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_48.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_49.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_5.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_50.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_51.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_52.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_53.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_54.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_55.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_56.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_57.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_58.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_59.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_6.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_60.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_61.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_62.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_63.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_64.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_65.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_66.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_67.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_68.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_69.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_7.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_70.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_71.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_72.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_73.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_74.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_75.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_76.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_77.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_78.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_79.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_8.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_80.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_81.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_82.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_83.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_84.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_85.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_86.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_87.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_88.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_89.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_9.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_90.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_91.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_92.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_93.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_94.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_95.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_96.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_97.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_98.root','/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD/RunIIFall17_DarkPhoton_v2-prod-uAOD-300-65-g97bd5bfd-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/180615_173327/0000/myMicroAODOutputFile_99.root')
)

process.load('flashgx.MicroAODAnalysis.FlashGXTruePileup_cfi')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('FlashGXTruePileup.root')
)

process.p = cms.Path(
    process.flashgxtruepileup
)