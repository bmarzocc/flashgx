import FWCore.ParameterSet.Config as cms
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections

flashggUnpackedJets = cms.EDProducer("FlashggVectorVectorJetUnpacker",
                                     JetsTag = cms.InputTag("flashggFinalJets"),
                                     NCollections = cms.uint32(maxJetCollections)
)

UnpackedJetCollectionVInputTag = cms.VInputTag()
for i in range(0,maxJetCollections):
   UnpackedJetCollectionVInputTag.append(cms.InputTag('flashggUnpackedJets',str(i)))


flashgxanalysis = cms.EDAnalyzer("FlashGXAnalysis",

    triggerResults                    = cms.InputTag("TriggerResults","","HLT"),
    #rhoCollection                     = cms.InputTag("fixedGridRhoAll"),
    rhoCollection                     = cms.InputTag("fixedGridRhoFastjetAll"),
    genParticleCollection             = cms.InputTag("flashggPrunedGenParticles"),
    vertexCollection                  = cms.InputTag("offlineSlimmedPrimaryVertices"),
    photonCollection                  = cms.InputTag("flashggRandomizedPhotons"),
    muonCollection                    = cms.InputTag("flashggSelectedMuons"),
    electronCollection                = cms.InputTag("flashggSelectedElectrons"),
    metCollection                     = cms.InputTag("flashggMets"),
    jetCollections                    = cms.InputTag("flashggFinalJets"),
    DiPhotonTag                       = cms.InputTag("flashggDiPhotons"),
    
    #Higgs_pt reweighting
    #higgsPtReweighting                = cms.string("data/GluGluHToGX_pt_reweighting.root"),
    higgsPtReweighting                = cms.string("data/VBFHToGX_pt_reweighting.root"),

    #triggerPaths                      = cms.vstring("NONE"), #No trigger studies
    #triggerPaths                      = cms.vstring("ALL"), #Make trigger studies
    triggerPaths                      = cms.vstring("HLT_Photon60_R9Id90_CaloIdL_IsoL",
                                                    "HLT_Photon20_HoverELoose",
                                                    "HLT_Photon30_HoverELoose",
                                                    "HLT_Photon40_HoverELoose",
                                                    "HLT_Photon50_HoverELoose",
                                                    "HLT_Photon60_HoverELoose",
                                                    "HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50",
                                                    "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3",
                                                    "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3"
                                                    "HLT_DiJet110_35_Mjj650_PFMET110",
                                                    "HLT_DiJet110_35_Mjj650_PFMET120",
                                                    "HLT_DiJet110_35_Mjj650_PFMET130",
                                                    "HLT_TripleJet110_35_35_Mjj650_PFMET110",
                                                    "HLT_TripleJet110_35_35_Mjj650_PFMET120",
                                                    "HLT_TripleJet110_35_35_Mjj650_PFMET130"
                                                   ), 

    photonPtThres                     = cms.double(20.),
    photonR9Thres                     = cms.double(0.8),

    applyCuBasedPhotonID              = cms.bool(False),
    applyMVAPhotonID                  = cms.bool(True),
    applyHggPhotonID                  = cms.bool(False),

    #80X MVAphotonID 
    #phoMvaThres                       = cms.vdouble(0.2,0.2), #tight: EB,EE
    #phoMvaThres                       = cms.vdouble(0.68,0.60), #loose: EB,EE
    #passElectronVetoThres             = cms.vdouble(0.5,0.5), #EB,EE
    
    #92X MVAphotonID 
    phoMvaThres                       = cms.vdouble(0.27,0.14), #tight: EB,EE
    #phoMvaThres                       = cms.vdouble(0.67,0.54), #loose: EB,EE
    passElectronVetoThres             = cms.vdouble(0.5,0.5), #EB,EE

    #isoEffArea 80X
    #isoCorrEtaBinsCh                  = cms.vdouble(0.0000,1.0000,1.4790,2.0000,2.2000,2.3000,2.4000,999.),  
    #isoCorrAreaValCh                  = cms.vdouble(0.0360,0.0377,0.0306,0.0283,0.0254,0.0217,0.0167),
    #isoCorrEtaBinsNh                  = cms.vdouble(0.0000,1.0000,1.4790,2.0000,2.2000,2.3000,2.4000,999.),  
    #isoCorrAreaValNh                  = cms.vdouble(0.0597,0.0807,0.0629,0.0197,0.0184,0.0284,0.0591),
    #isoCorrEtaBinsPh                  = cms.vdouble(0.0000,1.0000,1.4790,2.0000,2.2000,2.3000,2.4000,999.),  
    #isoCorrAreaValPh                  = cms.vdouble(0.1210,0.1107,0.0699,0.1056,0.1457,0.1719,0.1998),
     
    #Tight 80X CutBased photonID 
    #hadronicOverEmThres               = cms.vdouble(0.02690,0.02130), #EB,EE
    #sigmaIetaIetaThres                = cms.vdouble(0.00994,0.03000), #EB,EE
    #egChargedHadronIsoThres           = cms.vdouble(0.20200,0.03400), #EB,EE
    #egNeutralHadronIsoEBThres         = cms.vdouble(0.26400,0.01480,0.000017), #p0,p1,p2
    #egNeutralHadronIsoEEThres         = cms.vdouble(0.58600,0.01630,0.000014), #p0,p1,p2
    #egPhotonIsoEBThres                = cms.vdouble(2.36200,0.00400), #p0,p1
    #egPhotonIsoEEThres                = cms.vdouble(2.61700,0.00340), #p0,p1

    #Medium 80X CutBased photonID 
    #hadronicOverEmThres               = cms.vdouble(0.03960,0.02190), #EB,EE
    #sigmaIetaIetaThres                = cms.vdouble(0.01022,0.03001), #EB,EE
    #egChargedHadronIsoThres           = cms.vdouble(0.44100,0.44200), #EB,EE
    #egNeutralHadronIsoEBThres         = cms.vdouble(2.72500,0.01480,0.000017), #p0,p1,p2
    #egNeutralHadronIsoEEThres         = cms.vdouble(1.71500,0.01630,0.000014), #p0,p1,p2
    #egPhotonIsoEBThres                = cms.vdouble(2.57100,0.00470), #p0,p1
    #egPhotonIsoEEThres                = cms.vdouble(3.86300,0.00340), #p0,p1

    #Loose 80X CutBased photonID 
    #hadronicOverEmThres               = cms.vdouble(0.05970,0.04810), #EB,EE
    #sigmaIetaIetaThres                = cms.vdouble(0.01031,0.03013), #EB,EE
    #egChargedHadronIsoThres           = cms.vdouble(1.29500,1.01100), #EB,EE
    #egNeutralHadronIsoEBThres         = cms.vdouble(10.9100,0.01480,0.000017), #p0,p1,p2
    #egNeutralHadronIsoEEThres         = cms.vdouble(5.93100,0.01630,0.000014), #p0,p1,p2
    #egPhotonIsoEBThres                = cms.vdouble(3.63000,0.00470), #p0,p1
    #egPhotonIsoEEThres                = cms.vdouble(6.64100,0.00340), #p0,p1

    #isoEffArea 92X
    isoCorrEtaBinsCh                  = cms.vdouble(0.0000,1.0000,1.4790,2.0000,2.2000,2.3000,2.4000,999.),  
    isoCorrAreaValCh                  = cms.vdouble(0.0385,0.0468,0.0435,0.0378,0.0338,0.0314,0.0269),
    isoCorrEtaBinsNh                  = cms.vdouble(0.0000,1.0000,1.4790,2.0000,2.2000,2.3000,2.4000,999.),  
    isoCorrAreaValNh                  = cms.vdouble(0.0636,0.1103,0.0759,0.0236,0.0151,0.00007,0.0132),
    isoCorrEtaBinsPh                  = cms.vdouble(0.0000,1.0000,1.4790,2.0000,2.2000,2.3000,2.4000,999.),  
    isoCorrAreaValPh                  = cms.vdouble(0.1240,0.1093,0.0631,0.0779,0.0999,0.1155,0.1373),

    #Tight 92X CutBased photonID 
    #hadronicOverEmThres               = cms.vdouble(0.0200,0.0250), #EB,EE
    #sigmaIetaIetaThres                = cms.vdouble(0.0103,0.0271), #EB,EE
    #egChargedHadronIsoThres           = cms.vdouble(1.1580,0.5750), #EB,EE
    #egNeutralHadronIsoEBThres         = cms.vdouble(1.2670,0.0126,0.000026), #p0,p1,p2
    #egNeutralHadronIsoEEThres         = cms.vdouble(8.9160,0.0119,0.000025), #p0,p1,p2
    #egPhotonIsoEBThres                = cms.vdouble(2.0650,0.0035), #p0,p1
    #egPhotonIsoEEThres                = cms.vdouble(3.2720,0.0040), #p0,p1

    #Medium 92X CutBased photonID 
    hadronicOverEmThres               = cms.vdouble(0.0350,0.0270), #EB,EE
    sigmaIetaIetaThres                = cms.vdouble(0.0103,0.0271), #EB,EE
    egChargedHadronIsoThres           = cms.vdouble(1.4160,1.0120), #EB,EE
    egNeutralHadronIsoEBThres         = cms.vdouble(2.4910,0.0126,0.000026), #p0,p1,p2
    egNeutralHadronIsoEEThres         = cms.vdouble(9.1310,0.0119,0.000025), #p0,p1,p2
    egPhotonIsoEBThres                = cms.vdouble(2.9520,0.0040), #p0,p1
    egPhotonIsoEEThres                = cms.vdouble(4.0950,0.0040), #p0,p1

    #Loose 92X CutBased photonID 
    #hadronicOverEmThres               = cms.vdouble(0.1050,0.0290), #EB,EE
    #sigmaIetaIetaThres                = cms.vdouble(0.0103,0.0276), #EB,EE
    #egChargedHadronIsoThres           = cms.vdouble(2.8390,2.1500), #EB,EE
    #egNeutralHadronIsoEBThres         = cms.vdouble(9.1880,0.0126,0.000026), #p0,p1,p2
    #egNeutralHadronIsoEEThres         = cms.vdouble(10.471,0.0119,0.000025), #p0,p1,p2
    #egPhotonIsoEBThres                = cms.vdouble(2.9560,0.0035), #p0,p1
    #egPhotonIsoEEThres                = cms.vdouble(4.8950,0.0040), #p0,p1


    #photonHggID
    egChargedHadronIsoHggThres        = cms.double(20.00),
    egChargedHadronIsoOverPtHggThres  = cms.double(0.30),
    hadronicOverEmHggThres            = cms.double(0.08),
    photonIdMvaDHggThres              = cms.double(-0.90),
    photonPtHgg                       = cms.double(20.),
    passElectronVetoHggThres          = cms.double(0.50),
    r9ThresMinHgg                     = cms.vdouble(0.85,0.80,0.90,0.80), #highEB,lowEB,highEE,lowEE   
    r9ThresMaxHgg                     = cms.vdouble(1.00,0.85,1.00,0.90), #highEB,lowEB,highEE,lowEE
    pfPhoIso03HggThres                = cms.vdouble(100000.00,4.00,100000.00,4.00), #highEB,lowEB,highEE,lowEE    
    trkSumPtHollowConeDR03HggThres    = cms.vdouble(100000.00,6.00,100000.00,6.00), #highEB,lowEB,highEE,lowEE   
    sigmaIetaIetaHggThres             = cms.vdouble(100000.00,0.015,100000.00,0.035), #highEB,lowEB,highEE,lowEE  
    isoCorrEtaBinsHgg                 = cms.vdouble(0.,0.9,1.5,2,2.2,3),  
    isoCorrAreaValHgg                 = cms.vdouble(0.16544,0.16544,0.13212,0.13212,0.13212),

    #data metPhi corrections and selections
    #metPhiCorr0                   = cms.vdouble(2.74771,1.46669), #x,y
    #metPhiCorr1                   = cms.vdouble(0.00842948,0.0295922), #x,y

    #MC metPhi corrections and selections
    #metPhiCorr0                    = cms.vdouble(-3.17373,0.480606), #x,y
    #metPhiCorr1                    = cms.vdouble(-0.024752,0.0297115), #x,y 
 
    #No MET-phi corrections   
    metPhiCorr0                    = cms.vdouble(0.,0.), #x,y
    metPhiCorr1                    = cms.vdouble(0.,0.), #x,y

    metThres                       = cms.double(40.00),
    dPhiPhotonMetThres             = cms.double(0.),

    #electrons and muons selection
    useStdMuonID                   = cms.bool(False),
    useStdElectronID               = cms.bool(True),
    leptonPtThres                  = cms.double(20),
    muonEtaThres                   = cms.double(2.4), 
    muPFIsoSumRelThres             = cms.double(0.25),
    muMiniIsoSumRelThres           = cms.double(0.06),	 
    electronEtaThres               = cms.vdouble(1.4442,1.566,2.5),
    electronIsoThres               = cms.double(0.15),
    elMiniIsoThres                 = cms.vdouble(0.045,0.08), #EB, EE
    impactParam                    = cms.vdouble(0.0261,0.41,0.118,0.822), #transEB, longEB, transEE, longEE
    useElectronMVARecipe           = cms.bool(True),
    useElectronLooseID             = cms.bool(True),

    #vbfjets selections
    inputTagJets                   = UnpackedJetCollectionVInputTag,
    jetPtThres                     = cms.double(20.),
    jetEtaThres                    = cms.double(4.7),
    drJetPhoCut                    = cms.double(0.4),
    usePuJetID                     = cms.bool(True),
    useJetID                       = cms.bool(True),
    jetIDLevel                     = cms.string("Loose"),
    merge3rdJet                    = cms.bool(True),
    #merge3rdJet                    = cms.bool(False),
    thirdJetDRCut                  = cms.double(0.8),
    pujidEtaBins                   = cms.vdouble(0.00,2.50,2.75,3.00,5.00),  
    #pujidWpPtBin1                  = cms.vdouble(0.69, -0.35, -0.26, -0.21), # 94X tight
    #pujidWpPtBin2                  = cms.vdouble(0.69, -0.35, -0.26, -0.21), # 94X tight
    #pujidWpPtBin3                  = cms.vdouble(0.69, -0.35, -0.26, -0.21), # 94X tight
    #pujidWpPtBin4                  = cms.vdouble(0.86, -0.10, -0.05, -0.01), # 94X tight
    pujidWpPtBin1                  = cms.vdouble(0.18, -0.55, -0.42, -0.36), # 94X medium
    pujidWpPtBin2                  = cms.vdouble(0.18, -0.55, -0.42, -0.36), # 94X medium
    pujidWpPtBin3                  = cms.vdouble(0.18, -0.55, -0.42, -0.36), # 94X medium
    pujidWpPtBin4                  = cms.vdouble(0.61, -0.35, -0.23, -0.17), # 94X medium
    #pujidWpPtBin1                  = cms.vdouble(-0.97, -0.68, -0.53, -0.47), # 94X loose
    #pujidWpPtBin2                  = cms.vdouble(-0.97, -0.68, -0.53, -0.47), # 94X loose
    #pujidWpPtBin3                  = cms.vdouble(-0.97, -0.68, -0.53, -0.47), # 94X loose
    #pujidWpPtBin4                  = cms.vdouble(-0.89, -0.52, -0.38, -0.30), # 94X loose
    rmsforwardCut                  = cms.double(3.0), # default was 0.03 , running on loose pujid
    
    deltaEtaJets                   = cms.double(1.),
    invMassJets                    = cms.double(300.),
)
