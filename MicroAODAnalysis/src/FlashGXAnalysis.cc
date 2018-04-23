// system include files
#include <memory>
#include <iostream>
#include <fstream>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/LooperFactory.h"
#include "FWCore/Framework/interface/ESProducerLooper.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
//#include "PhysicsTools/Utilities/macros/setTDRStyle.C"

// flashgg includes
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "flashgg/DataFormats/interface/SingleVertexView.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

#include "flashgx/MicroAODAnalysis/interface/FlashGXAnalysis.h"

#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TRegexp.h"

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <set>
#include <assert.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include <TMath.h>
#include <Math/VectorUtil.h>
//#include <boost/tokenizer.hpp>

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;
using namespace flashgg;

//
// constructors and destructor
//
FlashGXAnalysis::FlashGXAnalysis(const edm::ParameterSet& iConfig)
{

  triggerToken_          = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"));
  rhoToken_              = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"));
  genParticleToken_      = consumes<View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
  vertexToken_           = consumes<View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  photonToken_           = consumes<View<flashgg::Photon> >(iConfig.getParameter<edm::InputTag>("photonCollection"));
  muonToken_             = consumes<View<flashgg::Muon> >(iConfig.getParameter<edm::InputTag>("muonCollection"));
  electronToken_         = consumes<View<flashgg::Electron> >(iConfig.getParameter<edm::InputTag>("electronCollection"));
  metToken_              = consumes<View<flashgg::Met> >(iConfig.getParameter<edm::InputTag>("metCollection"));

  inputTagJets_  = iConfig.getParameter<std::vector<edm::InputTag> >("inputTagJets");
  for (unsigned i=0;i<inputTagJets_.size();i++) 
  {
       auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
       jetsToken_.push_back(token);
  }

  higgsPtReweighting_               = iConfig.getParameter<string>("higgsPtReweighting");

  triggerPaths_                     = iConfig.getParameter<vector<string> >("triggerPaths");

  photonPtThres_                    = iConfig.getParameter<double>("photonPtThres");
  photonR9Thres_                    = iConfig.getParameter<double>("photonR9Thres");
  deltaEtaJets_                     = iConfig.getParameter<double>("deltaEtaJets");
  invMassJets_                      = iConfig.getParameter<double>("invMassJets");
  applyCuBasedPhotonID_             = iConfig.getParameter<bool>("applyCuBasedPhotonID");
  applyMVAPhotonID_                 = iConfig.getParameter<bool>("applyMVAPhotonID");
  applyHggPhotonID_                 = iConfig.getParameter<bool>("applyHggPhotonID");

  phoMvaThres_                      = iConfig.getParameter<vector<double> >("phoMvaThres");
  passElectronVetoThres_            = iConfig.getParameter<vector<double> >("passElectronVetoThres");

  isoCorrEtaBinsCh_                 = iConfig.getParameter<vector<double> >("isoCorrEtaBinsCh");
  isoCorrAreaValCh_                 = iConfig.getParameter<vector<double> >("isoCorrAreaValCh");
  isoCorrEtaBinsNh_                 = iConfig.getParameter<vector<double> >("isoCorrEtaBinsNh");
  isoCorrAreaValNh_                 = iConfig.getParameter<vector<double> >("isoCorrAreaValNh");
  isoCorrEtaBinsPh_                 = iConfig.getParameter<vector<double> >("isoCorrEtaBinsPh");
  isoCorrAreaValPh_                 = iConfig.getParameter<vector<double> >("isoCorrAreaValPh");
  hadronicOverEmThres_              = iConfig.getParameter<vector<double> >("hadronicOverEmThres");
  sigmaIetaIetaThres_               = iConfig.getParameter<vector<double> >("sigmaIetaIetaThres");
  egChargedHadronIsoThres_          = iConfig.getParameter<vector<double> >("egChargedHadronIsoThres");
  egNeutralHadronIsoEBThres_        = iConfig.getParameter<vector<double> >("egNeutralHadronIsoEBThres");
  egNeutralHadronIsoEEThres_        = iConfig.getParameter<vector<double> >("egNeutralHadronIsoEEThres");
  egPhotonIsoEBThres_               = iConfig.getParameter<vector<double> >("egPhotonIsoEBThres");
  egPhotonIsoEEThres_               = iConfig.getParameter<vector<double> >("egPhotonIsoEEThres");
  
  isoCorrEtaBinsHgg_                = iConfig.getParameter<vector<double> >("isoCorrEtaBinsHgg");
  isoCorrAreaValHgg_                = iConfig.getParameter<vector<double> >("isoCorrAreaValHgg");
  egChargedHadronIsoHggThres_       = iConfig.getParameter<double>("egChargedHadronIsoHggThres"); 
  egChargedHadronIsoOverPtHggThres_ = iConfig.getParameter<double>("egChargedHadronIsoOverPtHggThres"); 
  hadronicOverEmHggThres_           = iConfig.getParameter<double>("hadronicOverEmHggThres"); 
  photonPtHggThres_                 = iConfig.getParameter<double>("photonPtHgg"); 
  photonIdMvaDHggThres_             = iConfig.getParameter<double>("photonIdMvaDHggThres"); 
  r9ThresMinHgg_                    = iConfig.getParameter<vector<double> >("r9ThresMinHgg");
  r9ThresMaxHgg_                    = iConfig.getParameter<vector<double> >("r9ThresMaxHgg");
  passElectronVetoHggThres_         = iConfig.getParameter<double>("passElectronVetoHggThres"); 
  pfPhoIso03HggThres_               = iConfig.getParameter<vector<double> >("pfPhoIso03HggThres");
  trkSumPtHollowConeDR03HggThres_   = iConfig.getParameter<vector<double> >("trkSumPtHollowConeDR03HggThres"); 
  sigmaIetaIetaHggThres_            = iConfig.getParameter<vector<double> >("sigmaIetaIetaHggThres"); 

  metPhiCorr0_                      = iConfig.getParameter<vector<double> >("metPhiCorr0");
  metPhiCorr1_                      = iConfig.getParameter<vector<double> >("metPhiCorr1");
  metThres_                         = iConfig.getParameter<double>("metThres");
  dPhiPhotonMetThres_               = iConfig.getParameter<double>("dPhiPhotonMetThres");

  useStdMuonID_                     = iConfig.getParameter<bool>("useStdMuonID");
  muonEtaThres_                     = iConfig.getParameter<double>("muonEtaThres");
  leptonPtThres_                    = iConfig.getParameter<double>("leptonPtThres");
  muMiniIsoSumRelThres_             = iConfig.getParameter<double>("muMiniIsoSumRelThres");
  muPFIsoSumRelThres_               = iConfig.getParameter<double>("muPFIsoSumRelThres");

  useStdElectronID_                 = iConfig.getParameter<bool>("useStdElectronID");
  electronEtaThres_                 = iConfig.getParameter<vector<double> >("electronEtaThres");    
  elMiniIsoThres_                   = iConfig.getParameter<vector<double> >("elMiniIsoThres");
  impactParam_                      = iConfig.getParameter<vector<double> >("impactParam");
  useElectronMVARecipe_             = iConfig.getParameter<bool>("useElectronMVARecipe");
  useElectronLooseID_               = iConfig.getParameter<bool>("useElectronLooseID");

  jetPtThres_                       = iConfig.getParameter<double>("jetPtThres");
  jetEtaThres_                      = iConfig.getParameter<double>("jetEtaThres");             
  drJetPhoCut_                      = iConfig.getParameter<double>("drJetPhoCut");

  usePuJetID_                       = iConfig.getParameter<bool>("usePuJetID");
  useJetID_                         = iConfig.getParameter<bool>("useJetID");
  merge3rdJet_                      = iConfig.getParameter<bool>("merge3rdJet");
  thirdJetDRCut_                    = iConfig.getParameter<double>("thirdJetDRCut");
  rmsforwardCut_                    = iConfig.getParameter<double>("rmsforwardCut");
  jetIDLevel_                       = iConfig.getParameter<string>("jetIDLevel");
  drJetPhoVBFCut_                   = iConfig.getParameter<double>("drJetPhoVBFCut");
  pujid_wp_pt_bin_1_                = iConfig.getParameter<std::vector<double> >("pujidWpPtBin1");
  pujid_wp_pt_bin_2_                = iConfig.getParameter<std::vector<double> >("pujidWpPtBin2");
  pujid_wp_pt_bin_3_                = iConfig.getParameter<std::vector<double> >("pujidWpPtBin3");

  //output file and historgrams
  h_numberOfEvents        = iFile->make<TH1D>("h_numberOfEvents","h_numberOfEvents",10,-0.5,9.5);
  h_Efficiency            = iFile->make<TH1D>("h_Efficiency","h_Efficiency",10,-0.5,9.5);

  h_Pho_Energy_noCut      = iFile->make<TH1D>("h_Pho_Energy_noCut","h_Pho_Energy_noCut",100,0.,600.);
  h_Pho_Pt_noCut          = iFile->make<TH1D>("h_Pho_Pt_noCut","h_Pho_Pt_noCut",100,0.,300.);
  h_Pho_Eta_noCut         = iFile->make<TH1D>("h_Pho_Eta_noCut","h_Pho_Eta_noCut",100,-3.,3.);
  h_Pho_Phi_noCut         = iFile->make<TH1D>("h_Pho_Phi_noCut","h_Pho_Phi_noCut",100,-3.5,3.5);
  h_Pho_Energy            = iFile->make<TH1D>("h_Pho_Energy","h_Pho_Energy",100,0.,600.);
  h_Pho_Pt                = iFile->make<TH1D>("h_Pho_Pt","h_Pho_Pt",100,0.,300.);
  h_Pho_Eta               = iFile->make<TH1D>("h_Pho_Eta","h_Pho_Eta",100,-3.,3.);
  h_Pho_Phi               = iFile->make<TH1D>("h_Pho_Phi","h_Pho_Phi",100,-3.5,3.5);

  h_MET_Pt_noCut          = iFile->make<TH1D>("h_MET_Pt_noCut","h_MET_Pt_noCut",100,0.,300.);
  h_MET_Phi_noCut         = iFile->make<TH1D>("h_MET_Phi_noCut","h_MET_Phi_noCut",100,-3.5,3.5);
  h_MET_Pt                = iFile->make<TH1D>("h_MET_Pt","h_MET_Pt",100,0.,300.);
  h_MET_Phi               = iFile->make<TH1D>("h_MET_Phi","h_MET_Phi",100,-3.5,3.5);
 
  h_dPhiMETPho_noCut      = iFile->make<TH1D>("h_dPhiMETPho_noCut","h_dPhiMETPho_noCut",100,-3.5,3.5);
  h_MtMETPho_noCut        = iFile->make<TH1D>("h_MtMETPho_noCut","h_MtMETPho_noCut",100,0.,200.);
  h_dPhiMETPho            = iFile->make<TH1D>("h_dPhiMETPho","h_dPhiMETPho",100,-3.5,3.5);
  h_MtMETPho              = iFile->make<TH1D>("h_MtMETPho","h_MtMETPho",100,0.,200.);

  h_Jet_Energy            = iFile->make<TH1D>("h_Jet_Energy","h_Jet_Energy",100,0.,600.);
  h_Jet_Pt                = iFile->make<TH1D>("h_Jet_Pt","h_Jet_Pt",100,0.,300.);
  h_Jet_Eta               = iFile->make<TH1D>("h_Jet_Eta","h_Jet_Eta",100,-3.,3.);
  h_Jet_Phi               = iFile->make<TH1D>("h_Jet_Phi","h_Jet_Phi",100,-3.5,3.5);

  h_dPhiJetPho            = iFile->make<TH1D>("h_dPhiJetPho","h_dPhiJetPho",100,-3.5,3.5);
  h_dEtaJetPho            = iFile->make<TH1D>("h_dEtaJetPho","h_dEtaJetPho",100,-3.5,3.5);
  h_dRJetPho              = iFile->make<TH1D>("h_dRJetPho","h_dRJetPho",100,0.,6.);

  h_nJets_noCut           = iFile->make<TH1D>("h_nJets_noCut","h_nJets_noCut",8,-0.5,7.5);
  h_nJets                 = iFile->make<TH1D>("h_nJets","h_nJets",8,-0.5,7.5);
  h_nVBFJets              = iFile->make<TH1D>("h_nVBFJets","h_nVBFJets",8,-0.5,7.5);

  if(triggerPaths_.at(0)=="ALL")
  {
     unsigned int maxSize = 1000;
     HLT_Names.resize(maxSize);
     HLT_PhotonPt_efficiency_tmp.resize(maxSize);
     HLT_MET_efficiency_tmp.resize(maxSize);
     HLT_PhotonPt_efficiency.resize(maxSize);
     HLT_MET_efficiency.resize(maxSize);
     for(unsigned int iPath=0; iPath<maxSize; iPath++) 
     {
         char  photonName[500];
         sprintf (photonName,"PhotonPt_efficiency_%d",iPath);
         HLT_PhotonPt_efficiency_tmp[iPath] = new TH1D(photonName,photonName,20,0.,100.);
         HLT_PhotonPt_efficiency_tmp[iPath]->Sumw2();
         char  metName[500];
         sprintf (metName,"MET_efficiency_%d",iPath);
         HLT_MET_efficiency_tmp[iPath] = new TH1D(metName,metName,10,0.,200.);
         HLT_MET_efficiency_tmp[iPath]->Sumw2();
     }
  }else if(triggerPaths_.at(0)!="ALL" && triggerPaths_.at(0)!="NONE")
  {
     unsigned int maxSize = triggerPaths_.size();
     HLT_Names.resize(maxSize);
     HLT_PhotonPt_efficiency_tmp.resize(maxSize);
     HLT_MET_efficiency_tmp.resize(maxSize);
     HLT_PhotonPt_efficiency.resize(maxSize);
     HLT_MET_efficiency.resize(maxSize);
     for(unsigned int iPath=0; iPath<maxSize; iPath++) 
     {
         HLT_PhotonPt_efficiency_tmp[iPath] = new TH1D(std::string("PhotonPt_efficiency_"+triggerPaths_.at(iPath)).c_str(),std::string("PhotonPt_efficiency_"+triggerPaths_.at(iPath)).c_str(),20,0.,100.);
         HLT_PhotonPt_efficiency_tmp[iPath]->Sumw2();
         HLT_MET_efficiency_tmp[iPath] = new TH1D(std::string("MET_efficiency_"+triggerPaths_.at(iPath)).c_str(),std::string("MET_efficiency_"+triggerPaths_.at(iPath)).c_str(),10,0.,200.);
         HLT_MET_efficiency_tmp[iPath]->Sumw2();
     }
  }
  HLT_PhotonPt_denum = new TH1D("HLT_PhotonPt_denum","",20,0.,100.);
  HLT_PhotonPt_denum->Sumw2();
  HLT_MET_denum = new TH1D("HLT_MET_denum","",10,0.,200.);
  HLT_MET_denum->Sumw2();
  
  nTot=0.;
  nTot_step0=0.;
  nTotSelected_step0=0.;
  nTotSelected_step1=0.;
  nTotSelected_step2=0.;
  nTotSelected_step3=0.;
  nTotSelected_step4=0.; 
  nTotSelected_step5=0.;
  nTotSelected_step6=0.;
  nTotSelected_step7=0.;
  nTotSelected_step8=0.;

  TFile* higgsPtFile = TFile::Open(higgsPtReweighting_.c_str()); 
  h_higgsPtReweighting = (TH1D*)higgsPtFile->Get("h_Higgs_pt_reweighting");
}

FlashGXAnalysis::~FlashGXAnalysis()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void FlashGXAnalysis::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{

  edm::Handle<edm::TriggerResults> trigResults;
  ev.getByToken(triggerToken_,trigResults);

  edm::Handle<double>  rho_;
  ev.getByToken(rhoToken_,rho_);
  double rho = *rho_;

  Handle<View<reco::GenParticle> > genParticles;
  ev.getByToken(genParticleToken_, genParticles );

  Handle<View<reco::Vertex> > vertices;
  ev.getByToken(vertexToken_, vertices); 
  edm::Ptr<reco::Vertex> vtx = vertices->ptrAt(0); 

  Handle<View<flashgg::Photon> > photons;
  ev.getByToken(photonToken_, photons);  
  
  Handle<View<flashgg::Muon> > muons;
  ev.getByToken(muonToken_, muons);

  Handle<View<flashgg::Electron> > electrons;
  ev.getByToken(electronToken_, electrons);

  Handle<View<flashgg::Met> > mets;
  ev.getByToken(metToken_, mets);  
  if(mets->size() != 1) std::cout << "WARNING number of MET is not equal to 1" << std::endl; 
  edm::Ptr<flashgg::Met> theMet = mets->ptrAt(0);

  JetCollectionVector Jets(inputTagJets_.size());
  for(size_t j=0;j<inputTagJets_.size();++j)
      ev.getByToken(jetsToken_[j],Jets[j]);

  //compute the higgsPt weight
  double wPt = ptWeight(genParticles);
  //std::cout << "wPt = " << wPt << std::endl;

  nTot_step0+=wPt;
  if(photons->size() == 0) nTotSelected_step0+=wPt;
  if(photons->size() == 0) return;

  //apply global event selection
  nTot+=wPt;

  //compute trigger efficiencies
  int selectedPhoton=0;
  double maxPt=0.;
  for(unsigned int phoIndex=0; phoIndex<photons->size(); phoIndex++) 
  {
     edm::Ptr<flashgg::Photon> ipho = photons->ptrAt(phoIndex); 
     if(ipho->pt()>maxPt)
     {
        maxPt = ipho->pt();
        selectedPhoton = phoIndex;
     }        
  }  

  edm::Ptr<flashgg::Photon> thePhoton = photons->ptrAt(selectedPhoton);
  int n_triggers=0;
  const edm::TriggerNames& trigNames = ev.triggerNames(*trigResults);   
  if(triggerPaths_.at(0)=="NONE") n_triggers = 1;
  else if(triggerPaths_.at(0)!="NONE" && triggerPaths_.at(0)!="ALL")
  {
     for(unsigned iPath=0; iPath<triggerPaths_.size(); iPath++)
     {
         bool passTrig=getHLTResults(trigResults,trigNames,triggerPaths_.at(iPath));
         if(passTrig==true)
         {
            HLT_MET_efficiency_tmp.at(iPath)->Fill(theMet->getCorPt(),wPt); 
            HLT_PhotonPt_efficiency_tmp.at(iPath)->Fill(thePhoton->pt(),wPt); 
            n_triggers++;
         }
     }
  }
  else if(triggerPaths_.at(0)=="ALL")
  {
     if(HLT_pathsSize==0) HLT_pathsSize = trigResults->size();
     n_triggers=(int)trigResults->size();
     if(trigResults->size()>1000) HLT_Names.resize(trigResults->size());
     for(unsigned iPath=0; iPath<trigResults->size(); iPath++)
     {
         if(HLT_Names.at(iPath)=="") HLT_Names.at(iPath) = string(trigNames.triggerName(iPath));
         bool passTrig=getHLTResults(trigResults,trigNames,trigNames.triggerName(iPath));
         if(passTrig==true)
         {
            HLT_MET_efficiency_tmp.at(iPath)->Fill(theMet->getCorPt(),wPt); 
            HLT_PhotonPt_efficiency_tmp.at(iPath)->Fill(thePhoton->pt(),wPt); 
         }
     }
  }
  //apply trigger selections
  if(n_triggers<1) return;
  nTotSelected_step1+=wPt;

  //apply photonID and select the best photon
  selectedPhoton=0;
  int nSelectedPhotons=0;
  maxPt=0.;
  for(unsigned int phoIndex=0; phoIndex<photons->size(); phoIndex++) 
  {
     edm::Ptr<flashgg::Photon> ipho = photons->ptrAt(phoIndex); 
     
     if(applyCuBasedPhotonID_)
     {
        if(isGoodPhotonCutBased(ipho,rho) && ipho->pt()>photonPtThres_ && ipho->full5x5_r9()>photonR9Thres_ && isInEcal(ipho)) nSelectedPhotons++;  
        else continue;
     }
     else if(applyMVAPhotonID_)
     {
        if(isGoodPhotonMVA(ipho) && ipho->pt()>photonPtThres_ && ipho->full5x5_r9()>photonR9Thres_ && isInEcal(ipho)) nSelectedPhotons++;  
        else continue;
     }
     else if(applyHggPhotonID_)
     {
        if(isGoodPhotonHgg(ipho,vtx,rho) && ipho->pt()>photonPtThres_ && ipho->full5x5_r9()>photonR9Thres_ && isInEcal(ipho)) nSelectedPhotons++;  
        else continue;
     }
   
     if(ipho->pt()>maxPt)
     {
        maxPt = ipho->pt();
        selectedPhoton = phoIndex;
     }        
  }   

  //apply preselection: photonID, photon pt and MET 
  if(nSelectedPhotons!=1) return;
  nTotSelected_step2+=wPt;
  thePhoton = photons->ptrAt(selectedPhoton);

  if(theMet->getCorPt()<metThres_) return; 
  //deltaPhi(MET-Photon) 
  //double deltaPhiPhotonMet = metPhiCorr(theMet)-thePhoton->phi();
  double deltaPhiPhotonMet = deltaPhi(theMet->corPhi(),thePhoton->phi());
  //if(fabs(deltaPhiPhotonMet)<dPhiPhotonMetThres_)  //skip if close
  //   if(fabs(deltaPhiPhotonMet)-TMath::Pi()<dPhiPhotonMetThres_) //skip if close but on other side of phi=0
  if(fabs(deltaPhiPhotonMet)<dPhiPhotonMetThres_) return;
  nTotSelected_step3+=wPt;

  vector<edm::Ptr<flashgg::Jet> > gJets = goodJets(Jets,thePhoton,0);
 
  if(gJets.size()<2) return;
  nTotSelected_step4+=wPt;

  unsigned selectedJet1=-1;
  unsigned selectedJet2=-1;

  maxPt=0.;
  int selectedJet=-1;
  for(unsigned int jetIndex=0;jetIndex<gJets.size();jetIndex++) 
  {
      edm::Ptr<flashgg::Jet> theJet = gJets[jetIndex];
      
      if(theJet->pt()>maxPt)
      {
         maxPt = theJet->pt();
         selectedJet = jetIndex;  
      }   
  }
  selectedJet1=selectedJet;

  maxPt=0.;
  selectedJet=-1;
  for(unsigned int jetIndex=0;jetIndex<gJets.size();jetIndex++) 
  {
      if(jetIndex==selectedJet1) continue;

      edm::Ptr<flashgg::Jet> theJet = gJets[jetIndex];
      
      if(theJet->pt()>maxPt)
      {
         maxPt = theJet->pt();
         selectedJet = jetIndex;  
      }   
  }
  selectedJet2=selectedJet;

  if((int)selectedJet1==-1 || (int)selectedJet2==-1) return;

  edm::Ptr<flashgg::Jet> theJet1 = gJets[selectedJet1];
  edm::Ptr<flashgg::Jet> theJet2 = gJets[selectedJet2];
  float invMassJets = sqrt((theJet1->energy()+theJet2->energy())*(theJet1->energy()+theJet2->energy()) -
                           (theJet1->px()+theJet2->px())*(theJet1->px()+theJet2->px()) -
                           (theJet1->py()+theJet2->py())*(theJet1->py()+theJet2->py()) -
                           (theJet1->pz()+theJet2->pz())*(theJet1->pz()+theJet2->pz()));

  if(invMassJets<invMassJets_) return;
  nTotSelected_step5+=wPt;
  if(fabs(theJet1->eta()-theJet2->eta())<deltaEtaJets_) return;
  nTotSelected_step6+=wPt;

  HLT_PhotonPt_denum->Fill(thePhoton->pt(),wPt);
  HLT_MET_denum->Fill(theMet->getCorPt(),wPt);

  h_Pho_Energy_noCut->Fill(thePhoton->energy(),wPt);
  h_Pho_Pt_noCut->Fill(thePhoton->pt(),wPt);
  h_Pho_Eta_noCut->Fill(thePhoton->eta(),wPt);
  h_Pho_Phi_noCut->Fill(thePhoton->phi(),wPt);

  h_MET_Pt_noCut->Fill(theMet->getCorPt(),wPt);
  h_MET_Phi_noCut->Fill(metPhiCorr(theMet),wPt);

  h_dPhiMETPho_noCut->Fill(deltaPhiPhotonMet,wPt);
  h_MtMETPho_noCut->Fill(computeMtMETPhoton(thePhoton,theMet),wPt);

  //muon selections
  std::vector<edm::Ptr<flashgg::Muon> > goodMuons;
  if(!useStdMuonID_)
  {
     goodMuons = selectAllMuonsSum16(muons->ptrs(),vertices->ptrs(),muonEtaThres_,leptonPtThres_,muMiniIsoSumRelThres_);
  }else{
     goodMuons = selectAllMuons(muons->ptrs(),vertices->ptrs(),muonEtaThres_,leptonPtThres_,muPFIsoSumRelThres_);
  }

  //electron selections
  std::vector<edm::Ptr<Electron> > goodElectrons;
  if(!useStdElectronID_) 
  {
     goodElectrons = selectAllElectronsSum16(electrons->ptrs(),vertices->ptrs(),leptonPtThres_,electronEtaThres_,true,true,elMiniIsoThres_[0],
                                             elMiniIsoThres_[1],impactParam_[0],impactParam_[1],impactParam_[2],impactParam_[3],
                                             rho,ev.isRealData());
  }else{
     goodElectrons = selectStdAllElectrons(electrons->ptrs(),vertices->ptrs(),leptonPtThres_,electronEtaThres_,useElectronMVARecipe_,useElectronLooseID_,
                                           rho,ev.isRealData());
  }

  //no good muons and electrons
  if(goodMuons.size() != 0 || goodElectrons.size() != 0) return;
  nTotSelected_step7+=wPt;
 
  h_Pho_Energy->Fill(thePhoton->energy(),wPt);
  h_Pho_Pt->Fill(thePhoton->pt(),wPt);
  h_Pho_Eta->Fill(thePhoton->eta(),wPt);
  h_Pho_Phi->Fill(thePhoton->phi(),wPt);

  h_MET_Pt->Fill(theMet->getCorPt(),wPt);
  h_MET_Phi->Fill(metPhiCorr(theMet),wPt);

  h_dPhiMETPho->Fill(deltaPhiPhotonMet,wPt);
  h_MtMETPho->Fill(computeMtMETPhoton(thePhoton,theMet),wPt);

  //jet selections
  h_nJets_noCut->Fill(Jets[0]->size(),wPt);
  
  for(unsigned int jetIndex = 0;jetIndex<gJets.size();jetIndex++) 
  {
      edm::Ptr<flashgg::Jet> theJet = gJets[jetIndex];

      h_Jet_Energy->Fill(theJet->energy(),wPt);
      h_Jet_Pt->Fill(theJet->pt(),wPt);
      h_Jet_Eta->Fill(theJet->eta(),wPt);
      h_Jet_Phi->Fill(theJet->phi(),wPt);

      h_dPhiJetPho->Fill(deltaPhi(theJet->phi(),thePhoton->phi()),wPt);
      h_dEtaJetPho->Fill(theJet->eta()-thePhoton->eta(),wPt);
      h_dRJetPho->Fill(deltaR(theJet->eta(),theJet->phi(),thePhoton->eta(),thePhoton->phi()),wPt);
     
  }
  h_nJets->Fill((int)gJets.size(),wPt);

  vector<edm::Ptr<flashgg::Jet> > vbfjets = vbfJets(Jets,thePhoton,vtx,0);
  //std::cout << "nVbfjets = " << vbfjets.size() << std::endl;
  h_nVBFJets->Fill((int)vbfjets.size(),wPt);

  if(vbfjets.size() != 0) nTotSelected_step8+=wPt; 
 
}

        
void FlashGXAnalysis::beginJob()
{

}

void FlashGXAnalysis::endJob() 
{

   if(triggerPaths_.at(0)!="NONE" && triggerPaths_.at(0)!="ALL")
   {
      for(unsigned int iPath=0; iPath<triggerPaths_.size(); iPath++) 
      {
          HLT_PhotonPt_efficiency_tmp[iPath]->Divide(HLT_PhotonPt_denum);
          HLT_MET_efficiency_tmp[iPath]->Divide(HLT_MET_denum);
          if(isGoodTrigger(HLT_PhotonPt_efficiency_tmp[iPath],0.6) && isGoodTrigger(HLT_MET_efficiency_tmp[iPath],0.6))
          {
             HLT_PhotonPt_efficiency[iPath] = iFile->make<TH1D>(std::string("PhotonPt_efficiency_"+triggerPaths_.at(iPath)).c_str(),std::string("PhotonPt_efficiency_"+triggerPaths_.at(iPath)).c_str(),20,0.,100.);
             saveTrigger(HLT_PhotonPt_efficiency_tmp[iPath],HLT_PhotonPt_efficiency[iPath]);
             HLT_MET_efficiency[iPath] = iFile->make<TH1D>(std::string("MET_efficiency_"+triggerPaths_.at(iPath)).c_str(),std::string("MET_efficiency_"+triggerPaths_.at(iPath)).c_str(),10,0.,200.);
             saveTrigger(HLT_MET_efficiency_tmp[iPath],HLT_MET_efficiency[iPath]);
          }
      }
   }
   else if(triggerPaths_.at(0)=="ALL")
   {
      for(unsigned int iPath=0; iPath<HLT_pathsSize; iPath++) 
      {
          HLT_PhotonPt_efficiency_tmp[iPath]->Divide(HLT_PhotonPt_denum);
          HLT_MET_efficiency_tmp[iPath]->Divide(HLT_MET_denum);
          if(isGoodTrigger(HLT_PhotonPt_efficiency_tmp[iPath],0.6) && isGoodTrigger(HLT_MET_efficiency_tmp[iPath],0.6))
          {
             if (HLT_Names.at(iPath).find("HLT")==std::string::npos) continue;

             HLT_PhotonPt_efficiency[iPath] = iFile->make<TH1D>(std::string("PhotonPt_efficiency_"+HLT_Names.at(iPath)).c_str(),std::string("PhotonPt_efficiency_"+HLT_Names.at(iPath)).c_str(),20,0.,100.);
             saveTrigger(HLT_PhotonPt_efficiency_tmp[iPath],HLT_PhotonPt_efficiency[iPath]);
             HLT_MET_efficiency[iPath] = iFile->make<TH1D>(std::string("MET_efficiency_"+HLT_Names.at(iPath)).c_str(),std::string("MET_efficiency_"+HLT_Names.at(iPath)).c_str(),10,0.,200.);
             saveTrigger(HLT_MET_efficiency_tmp[iPath],HLT_MET_efficiency[iPath]);
          }
      }
   }

   double eff=-1.;
   double eff_error=0.;
   h_numberOfEvents->SetBinContent(1,nTot); 
   
   eff=(double)nTotSelected_step1/(double)nTot;
   eff_error=sqrt(eff*(1-eff)/(double)nTot);
   h_numberOfEvents->SetBinContent(2,nTotSelected_step1); 
   h_Efficiency->SetBinContent(2,eff); 
   h_Efficiency->SetBinError(2,eff_error); 

   eff=(double)nTotSelected_step2/(double)nTot;
   eff_error=sqrt(eff*(1-eff)/(double)nTot);
   h_numberOfEvents->SetBinContent(3,nTotSelected_step2); 
   h_Efficiency->SetBinContent(3,eff); 
   h_Efficiency->SetBinError(3,eff_error); 
 
   eff=(double)nTotSelected_step3/(double)nTot;
   eff_error=sqrt(eff*(1-eff)/(double)nTot);
   h_numberOfEvents->SetBinContent(4,nTotSelected_step3); 
   h_Efficiency->SetBinContent(4,eff); 
   h_Efficiency->SetBinError(4,eff_error);

   eff=(double)nTotSelected_step4/(double)nTot;
   eff_error=sqrt(eff*(1-eff)/(double)nTot);
   h_numberOfEvents->SetBinContent(5,nTotSelected_step4); 
   h_Efficiency->SetBinContent(5,eff); 
   h_Efficiency->SetBinError(5,eff_error);

   eff=(double)nTotSelected_step5/(double)nTot;
   eff_error=sqrt(eff*(1-eff)/(double)nTot);
   h_numberOfEvents->SetBinContent(6,nTotSelected_step5); 
   h_Efficiency->SetBinContent(6,eff); 
   h_Efficiency->SetBinError(6,eff_error);
   
   eff=(double)nTotSelected_step6/(double)nTot;
   eff_error=sqrt(eff*(1-eff)/(double)nTot);
   h_numberOfEvents->SetBinContent(7,nTotSelected_step6); 
   h_Efficiency->SetBinContent(7,eff); 
   h_Efficiency->SetBinError(7,eff_error);

   eff=(double)nTotSelected_step7/(double)nTot;
   eff_error=sqrt(eff*(1-eff)/(double)nTot);
   h_numberOfEvents->SetBinContent(8,nTotSelected_step7); 
   h_Efficiency->SetBinContent(8,eff); 
   h_Efficiency->SetBinError(8,eff_error);
   
   eff=(double)nTotSelected_step8/(double)nTot;
   eff_error=sqrt(eff*(1-eff)/(double)nTot);
   h_numberOfEvents->SetBinContent(9,nTotSelected_step8); 
   h_Efficiency->SetBinContent(9,eff); 
   h_Efficiency->SetBinError(9,eff_error);
   

   cout << "------------- Preselections & HLT paths -------------" << endl; 
   cout << "Events with 0 photons  : " << (double)nTotSelected_step0/(double)nTot_step0 << endl; 
   cout << "HLT selections         : " << (double)nTotSelected_step1/(double)nTot << ", wrt HLT: " << (double)nTotSelected_step1/(double)nTotSelected_step1 << endl; 
   cout << "Photon selections      : " << (double)nTotSelected_step2/(double)nTot << ", wrt HLT: " << (double)nTotSelected_step2/(double)nTotSelected_step1 << endl;    
   cout << "MET selections         : " << (double)nTotSelected_step3/(double)nTot << ", wrt HLT: " << (double)nTotSelected_step3/(double)nTotSelected_step1 << endl;    
   cout << "Good jets selections   : " << (double)nTotSelected_step4/(double)nTot << ", wrt HLT: " << (double)nTotSelected_step4/(double)nTotSelected_step1 << endl;
   cout << "VBF invMass selection  : " << (double)nTotSelected_step5/(double)nTot << ", wrt HLT: " << (double)nTotSelected_step5/(double)nTotSelected_step1 << endl;
   cout << "VBF dEta selection     : " << (double)nTotSelected_step6/(double)nTot << ", wrt HLT: " << (double)nTotSelected_step6/(double)nTotSelected_step1 << endl;
   cout << "No good electrons/muons: " << (double)nTotSelected_step7/(double)nTot << ", wrt HLT: " << (double)nTotSelected_step7/(double)nTotSelected_step1 << endl; 
   cout << "VBF final selections   : " << (double)nTotSelected_step8/(double)nTot << ", wrt HLT: " << (double)nTotSelected_step8/(double)nTotSelected_step1 << endl;  
   
}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double FlashGXAnalysis::ptWeight(Handle<View<reco::GenParticle> > genParticles)
{
   double ptWeight=1.;
   int higgsIndex=-1;
   
   for(unsigned int genIndex = 0; genIndex<genParticles->size(); genIndex++) 
   {
     edm::Ptr<reco::GenParticle> genPart = genParticles->ptrAt(genIndex); 
     if(genPart->pdgId()==25 && genPart->status()==62 && genPart->mother()->pdgId()==25) higgsIndex=genIndex;
   }  

   
   if(higgsIndex!=-1)
   {
      edm::Ptr<reco::GenParticle> genPart = genParticles->ptrAt(higgsIndex);  
      if(h_higgsPtReweighting->FindBin(genPart->pt())<=h_higgsPtReweighting->GetNbinsX())
      {
         ptWeight=h_higgsPtReweighting->GetBinContent(h_higgsPtReweighting->FindBin(genPart->pt()));
      }else{
       ptWeight=h_higgsPtReweighting->GetBinContent(h_higgsPtReweighting->GetNbinsX()); 
      }
      //std::cout << "Higgs pt:" << genPart->pt() << " " << ptWeight << std::endl; 
   }  

   return ptWeight;
}
///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool FlashGXAnalysis::getHLTResults(edm::Handle<edm::TriggerResults> trigResults,const edm::TriggerNames& trigNames,std::string path)
{

  int hltCount = trigResults->size();
  TRegexp reg(TString( path.c_str()) );
  for (int ii = 0 ; ii != hltCount; ++ii) {
    TString hltName_tstr(trigNames.triggerName(ii));
    if ( hltName_tstr.Contains(reg) )
      return trigResults->accept(ii); // False or True depending if it fired.
  }
  return false;
}

bool FlashGXAnalysis::isGoodTrigger(TH1D* efficiency, double threshold)
{
  bool isGood = false;
  for(int iBin=1; iBin<=efficiency->GetNbinsX(); iBin++)
      if(efficiency->GetBinContent(iBin)>=threshold) isGood = true;

  return isGood; 
}

void FlashGXAnalysis::saveTrigger(TH1D* efficiency_tmp, TH1D* efficiency)
{
  for(int iBin=1; iBin<=efficiency_tmp->GetNbinsX(); iBin++)
      efficiency->SetBinContent(iBin,efficiency_tmp->GetBinContent(iBin));
}
///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool FlashGXAnalysis::isGoodPhotonHgg(edm::Ptr<flashgg::Photon> ipho, edm::Ptr<reco::Vertex> vtx, double rho)
{

  bool isGood = false;

  double effArea = 0.;
  for(unsigned int iBin=0; iBin<isoCorrEtaBinsHgg_.size()-1; iBin++)
      if(fabs(ipho->superCluster()->eta())>isoCorrEtaBinsHgg_.at(iBin) && fabs(ipho->superCluster()->eta())<=isoCorrEtaBinsHgg_.at(iBin+1)) effArea = isoCorrAreaValHgg_.at(iBin);
  double pfPhoIso03 = std::max(ipho->pfPhoIso03()-rho*effArea,0.);

  if(ipho->egChargedHadronIso()<egChargedHadronIsoHggThres_ || ipho->egChargedHadronIso()/ipho->pt()<egChargedHadronIsoOverPtHggThres_)
  { 
     if(ipho->hadronicOverEm()<hadronicOverEmHggThres_ && ipho->pt()>photonPtHggThres_ && ipho->phoIdMvaDWrtVtx(vtx)>photonIdMvaDHggThres_ && ipho->passElectronVeto()>passElectronVetoHggThres_ && isInEcal(ipho))   
     {
        if(ipho->isEB() && ipho->full5x5_r9()>r9ThresMinHgg_[0] && ipho->full5x5_r9()<=r9ThresMaxHgg_[0])
        {
           if(pfPhoIso03<pfPhoIso03HggThres_[0] && ipho->trkSumPtHollowConeDR03()<trkSumPtHollowConeDR03HggThres_[0] && ipho->full5x5_sigmaIetaIeta()<sigmaIetaIetaHggThres_[0]) 
              isGood = true;
        }    
        else if(ipho->isEB() && ipho->full5x5_r9()>r9ThresMinHgg_[1] && ipho->full5x5_r9()<=r9ThresMaxHgg_[1])
        {
           if(pfPhoIso03<pfPhoIso03HggThres_[1] && ipho->trkSumPtHollowConeDR03()<trkSumPtHollowConeDR03HggThres_[1] && ipho->full5x5_sigmaIetaIeta()<sigmaIetaIetaHggThres_[1]) 
              isGood = true;
        }  
        else if(ipho->isEB() && ipho->full5x5_r9()>r9ThresMinHgg_[2] && ipho->full5x5_r9()<=r9ThresMaxHgg_[2])
        {
           if(pfPhoIso03<pfPhoIso03HggThres_[2] && ipho->trkSumPtHollowConeDR03()<trkSumPtHollowConeDR03HggThres_[2] && ipho->full5x5_sigmaIetaIeta()<sigmaIetaIetaHggThres_[2]) 
              isGood = true; 
        }    
        else if(ipho->isEB() && ipho->full5x5_r9()>r9ThresMinHgg_[3] && ipho->full5x5_r9()<=r9ThresMaxHgg_[3])
        {
           if(pfPhoIso03<pfPhoIso03HggThres_[3] && ipho->trkSumPtHollowConeDR03()<trkSumPtHollowConeDR03HggThres_[3] && ipho->full5x5_sigmaIetaIeta()<sigmaIetaIetaHggThres_[3]) 
              isGood = true;
        }   
     }
  }      
  return isGood;

}

bool FlashGXAnalysis::isGoodPhotonCutBased(edm::Ptr<flashgg::Photon> ipho, double rho)
{

  bool isGood = false;

  double effAreaCh = 0.;
  double effAreaNh = 0.;
  double effAreaPh = 0.;

  for(unsigned int iBin=0; iBin<isoCorrEtaBinsCh_.size()-1; iBin++)
       if(fabs(ipho->superCluster()->eta())>isoCorrEtaBinsCh_.at(iBin) && fabs(ipho->superCluster()->eta())<=isoCorrEtaBinsCh_.at(iBin+1)) effAreaCh=isoCorrAreaValCh_.at(iBin);       
  for(unsigned int iBin=0; iBin<isoCorrEtaBinsNh_.size()-1; iBin++)
       if(fabs(ipho->superCluster()->eta())>isoCorrEtaBinsNh_.at(iBin) && fabs(ipho->superCluster()->eta())<=isoCorrEtaBinsNh_.at(iBin+1)) effAreaNh=isoCorrAreaValNh_.at(iBin);
  for(unsigned int iBin=0; iBin<isoCorrEtaBinsPh_.size()-1; iBin++)
       if(fabs(ipho->superCluster()->eta())>isoCorrEtaBinsPh_.at(iBin) && fabs(ipho->superCluster()->eta())<=isoCorrEtaBinsPh_.at(iBin+1)) effAreaPh=isoCorrAreaValPh_.at(iBin);

  double egChargedHadronIso = std::max(ipho->egChargedHadronIso()-rho*effAreaCh,0.);
  double egNeutralHadronIso = std::max(ipho->egNeutralHadronIso()-rho*effAreaNh,0.);
  double egPhotonIso = std::max(ipho->egPhotonIso()-rho*effAreaPh,0.);
  //double egPhotonIso = std::max(ipho->pfPhoIso03()-rho*effAreaPh,0.);
 
  double egNeutralHadronIsoEBThres = egNeutralHadronIsoEBThres_[0] + egNeutralHadronIsoEBThres_[1]*ipho->pt() + egNeutralHadronIsoEBThres_[2]*ipho->pt()*ipho->pt();
  double egNeutralHadronIsoEEThres = egNeutralHadronIsoEEThres_[0] + egNeutralHadronIsoEEThres_[1]*ipho->pt() + egNeutralHadronIsoEEThres_[2]*ipho->pt()*ipho->pt();
  double egPhotonIsoEBThres = egPhotonIsoEBThres_[0] + egPhotonIsoEBThres_[1]*ipho->pt();
  double egPhotonIsoEEThres = egPhotonIsoEEThres_[0] + egPhotonIsoEEThres_[1]*ipho->pt();
  
  if(ipho->isEB())
  {
     if(ipho->hadronicOverEm() < hadronicOverEmThres_[0] && ipho->sigmaIetaIeta() < sigmaIetaIetaThres_[0] && egChargedHadronIso<egChargedHadronIsoThres_[0] && egNeutralHadronIso<egNeutralHadronIsoEBThres && egPhotonIso<egPhotonIsoEBThres && isInEcal(ipho))    
        isGood = true;
  }
  else if(ipho->isEE())
  {
     if(ipho->hadronicOverEm() < hadronicOverEmThres_[1] && ipho->sigmaIetaIeta() < sigmaIetaIetaThres_[1] && egChargedHadronIso<egChargedHadronIsoThres_[1] && egNeutralHadronIso<egNeutralHadronIsoEEThres && egPhotonIso<egPhotonIsoEEThres && isInEcal(ipho))    
        isGood = true;
  }
  return isGood;

}

bool FlashGXAnalysis::isGoodPhotonMVA(edm::Ptr<flashgg::Photon> ipho)
{

  bool isGood = false;
  if(ipho->isEB())
  {
     if(ipho->egmMVA()>phoMvaThres_[0] && ipho->passElectronVeto()>passElectronVetoThres_[0]) isGood = true;
  }
  else if(ipho->isEE())
  {
     if(ipho->egmMVA()>phoMvaThres_[1] && ipho->passElectronVeto()>passElectronVetoThres_[1]) isGood = true;
  }
  return isGood; 

}
///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool FlashGXAnalysis::isInEcal(edm::Ptr<flashgg::Photon> ipho)
{
   bool isGood = false;
   if(fabs(ipho->superCluster()->eta())<2.5 && (fabs(ipho->superCluster()->eta())<1.4442 || fabs(ipho->superCluster()->eta())>1.566)) 
      isGood = true;
   
   return isGood;   
}
///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double FlashGXAnalysis::metPhiCorr(edm::Ptr<flashgg::Met> met)
{
   double newPhi=met->getCorPhi();
   
   double oldPx=met->getCorPx();
   double oldPy=met->getCorPy();
   double newPx=oldPx - (metPhiCorr0_[0]+metPhiCorr1_[0]*met->corSumEt());
   double newPy=oldPy - (metPhiCorr0_[1]+metPhiCorr1_[1]*met->corSumEt());

   newPhi = atan(newPy/newPx);
   if(newPx<0. && newPy<0.)
      newPhi = -TMath::Pi() + newPhi;
   else if(newPx<0. && newPy>=0.)
      newPhi = TMath::Pi() + newPhi;

   return newPhi;
}

double FlashGXAnalysis::computeMtMETPhoton(edm::Ptr<flashgg::Photon> ipho, edm::Ptr<flashgg::Met> met)
{

     double Mt = 0.;

     double oldPx=met->getCorPx();
     double oldPy=met->getCorPy();
     double newPx=oldPx - (metPhiCorr0_[0]+metPhiCorr1_[0]*met->corSumEt());
     double newPy=oldPy - (metPhiCorr0_[1]+metPhiCorr1_[1]*met->corSumEt());

     TLorentzVector* pho = new TLorentzVector();
     pho->SetPtEtaPhiE(ipho->pt(),ipho->eta(),ipho->phi(),ipho->energy());

     double m1 = 0.;
     double m2 = 0.;
     double Et1 = sqrt(m1*m1 + pho->Px()*pho->Px() + pho->Py()*pho->Py());
     double Et2 = sqrt(m2*m2 + newPx*newPx + newPy*newPy);
     double pt2 = pho->Px()*newPx + pho->Py()*newPy;
     Mt = sqrt(m1*m1 + m2*m2 + 2*(Et1*Et2 - pt2));
     return Mt;
}
///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
vector<edm::Ptr<flashgg::Jet> > FlashGXAnalysis::goodJets(JetCollectionVector Jets, edm::Ptr<flashgg::Photon> pho, int iColl)
{
   vector<edm::Ptr<flashgg::Jet> > goodJets;

   for(unsigned int jetIndex = 0;jetIndex<Jets[iColl]->size();jetIndex++) 
   {
      edm::Ptr<flashgg::Jet> jet = Jets[iColl]->ptrAt(jetIndex);
      float dRPhoJet = deltaR(jet->eta(),jet->phi(),pho->superCluster()->eta(),pho->superCluster()->phi());

      if(fabs(jet->eta())>jetEtaThres_) { continue; }
      if(!jet->passesJetID(flashgg::Loose)) { continue; }
      if(dRPhoJet<drJetPhoCut_) { continue; }
      if(jet->pt()<jetPtThres_) { continue; }

      goodJets.push_back(jet);
   }
   return goodJets;
}

vector<edm::Ptr<flashgg::Jet> > FlashGXAnalysis::vbfJets(JetCollectionVector Jets, edm::Ptr<flashgg::Photon> pho, edm::Ptr<reco::Vertex> vtx, int iColl)
{
   std::vector<edm::Ptr<flashgg::Jet> > vbfJets;

   std::vector<std::pair<double,double> > eta_cuts_(4);
   eta_cuts_[0] = std::make_pair(0,   2.50);
   eta_cuts_[1] = std::make_pair(2.50,2.75);
   eta_cuts_[2] = std::make_pair(2.75,3.00);
   eta_cuts_[3] = std::make_pair(3.00,10.);
   
   std::pair <int, int>  dijet_indices(-1,-1);
   std::pair <float, float> dijet_pts(-1.,-1.);
   int jet_3_index = -1;
   float jet_3_pt = -1.;

   bool hasValidVBFDiJet  = 0;
   bool hasValidVBFTriJet = 0;
            
   int n_jets_count = 0;

   for(unsigned int jetIndex = 0;jetIndex<Jets[iColl]->size();jetIndex++) 
   {
      edm::Ptr<flashgg::Jet> jet = Jets[iColl]->ptrAt(jetIndex);

      if(usePuJetID_ && !jet->passesPuJetId(vtx)){ continue;}
      if(useJetID_)
      {
          if(jetIDLevel_ == "Loose" && !jet->passesJetID(flashgg::Loose )) continue;
          if(jetIDLevel_ == "Tight" && !jet->passesJetID(flashgg::Tight )) continue;
      }
  
      if(fabs(jet->eta())>2.5 && jet->rms()>rmsforwardCut_){ continue; } 

                
      if((!pujid_wp_pt_bin_1_.empty()) && (!pujid_wp_pt_bin_2_.empty()) && (!pujid_wp_pt_bin_3_.empty()))
      {
          bool pass=false;
          for(UInt_t eta_bin=0; eta_bin < pujid_wp_pt_bin_1_.size(); eta_bin++)
          {
              if(fabs(jet->eta())> eta_cuts_[eta_bin].first && fabs(jet->eta())<=eta_cuts_[eta_bin].second)
              {
                 if(jet->pt()>20. && jet->pt()<=30.  && jet->puJetIdMVA()>pujid_wp_pt_bin_1_[eta_bin]) pass=true;                                
                 if(jet->pt()>30. && jet->pt()<=50.  && jet->puJetIdMVA()>pujid_wp_pt_bin_2_[eta_bin]) pass=true;
                 if(jet->pt()>50. && jet->pt()<=100. && jet->puJetIdMVA()>pujid_wp_pt_bin_3_[eta_bin]) pass=true;
                 if(jet->pt()>100.) pass = true;          
              }
          }                   
          if (!pass) continue;
      }

      if(fabs(jet->eta())> 4.7) { continue; }
      if(deltaR(jet->eta(),pho->eta(),jet->phi(),pho->phi())<drJetPhoVBFCut_) { continue; }

      if(jet->pt()>dijet_pts.first)
      {
         dijet_indices.second = dijet_indices.first;
         dijet_pts.second = dijet_pts.first;
 
         dijet_indices.first = jetIndex;
         dijet_pts.first = jet->pt();
      }
      else if(jet->pt()>dijet_pts.second)
      { 
         // for the 3rd jets
         jet_3_index = dijet_indices.second;
         jet_3_pt = dijet_pts.second;
                    
         dijet_indices.second = jetIndex;
         dijet_pts.second = jet->pt();         
      }
      else if(jet->pt()>jet_3_pt)
      {
         jet_3_index = jetIndex;
         jet_3_pt = jet->pt();
      }

      n_jets_count++;
      // if the jet's pt is neither higher than the lead jet or sublead jet, then forget it!
      if(dijet_indices.first != -1 && dijet_indices.second != -1) {hasValidVBFDiJet  = 1;}
      if(hasValidVBFDiJet && jet_3_index != -1) {hasValidVBFTriJet = 1;}
   }

   if(n_jets_count == 0) vbfJets.resize(0);
   else if(n_jets_count == 1)
   {
      vbfJets.resize(1);
      if(dijet_indices.first != -1 && dijet_indices.second == -1) vbfJets.push_back(Jets[iColl]->ptrAt(dijet_indices.first));
      if(dijet_indices.first == -1 && dijet_indices.second == -1) vbfJets.push_back(Jets[iColl]->ptrAt(dijet_indices.second));
   }  
   else if(n_jets_count == 2 && hasValidVBFDiJet)
   {
      vbfJets.resize(2);
      vbfJets.push_back(Jets[iColl]->ptrAt(dijet_indices.first));
      vbfJets.push_back(Jets[iColl]->ptrAt(dijet_indices.second));
   }  
   else if(n_jets_count == 3 && hasValidVBFTriJet)
   {
      vbfJets.resize(3);
      vbfJets.push_back(Jets[iColl]->ptrAt(dijet_indices.first));
      vbfJets.push_back(Jets[iColl]->ptrAt(jet_3_index));
      vbfJets.push_back(Jets[iColl]->ptrAt(dijet_indices.second));
   }  

   return vbfJets;
}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FlashGXAnalysis);











