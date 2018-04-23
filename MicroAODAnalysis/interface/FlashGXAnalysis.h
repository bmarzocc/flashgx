#ifndef FlashGXAnalysis_H
#define FlashGXAnalysis_H

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
#include "PhysicsTools/Utilities/macros/setTDRStyle.C"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

// flashgg includes
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "flashgg/DataFormats/interface/SingleVertexView.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

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

#include <TMath.h>
#include <Math/VectorUtil.h>
//#include <boost/tokenizer.hpp>

typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

class FlashGXAnalysis : public edm::EDAnalyzer
{
      public:
         explicit FlashGXAnalysis(const edm::ParameterSet&);
	 ~FlashGXAnalysis();
  
  
      private:
	 virtual void beginJob() ;
	 virtual void analyze(const edm::Event&, const edm::EventSetup&);
         virtual void endJob() ;

      // ----------additional functions-------------------
      double ptWeight(Handle<View<reco::GenParticle> > genParticles);
      bool getHLTResults(edm::Handle<edm::TriggerResults> trigResults,const edm::TriggerNames& trigNames,std::string path);
      bool isGoodPhotonHgg(edm::Ptr<flashgg::Photon> pho, edm::Ptr<reco::Vertex> vtx, double rho); 
      bool isGoodPhotonCutBased(edm::Ptr<flashgg::Photon> ipho, double rho);
      bool isGoodPhotonMVA(edm::Ptr<flashgg::Photon> ipho);
      bool isInEcal(edm::Ptr<flashgg::Photon> ipho);
      bool isGoodTrigger(TH1D* efficiency, double threshold);
      void saveTrigger(TH1D* efficiency_tmp, TH1D* efficiency);
      double phoIsoCorr(edm::Ptr<flashgg::Photon> ipho, double rho);
      double metPhiCorr(edm::Ptr<flashgg::Met> met);
      double computeMtMETPhoton(edm::Ptr<flashgg::Photon> ipho, edm::Ptr<flashgg::Met> met);
      vector<edm::Ptr<flashgg::Jet> > goodJets(JetCollectionVector Jets, edm::Ptr<flashgg::Photon> pho, int iColl);
      vector<edm::Ptr<flashgg::Jet> > vbfJets(JetCollectionVector Jets, edm::Ptr<flashgg::Photon> pho, edm::Ptr<reco::Vertex> vtx, int iColl);
      
      edm::EDGetTokenT<edm::TriggerResults> triggerToken_;  
      edm::EDGetTokenT<double> rhoToken_;  
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticleToken_;  
      edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;  
      edm::EDGetTokenT<edm::View<flashgg::Photon> > photonToken_; 
      edm::EDGetTokenT<edm::View<flashgg::Muon> > muonToken_; 
      edm::EDGetTokenT<edm::View<flashgg::Electron> > electronToken_; 
      edm::EDGetTokenT<edm::View<flashgg::Met> > metToken_;
      std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > jetsToken_;

      std::string higgsPtReweighting_;
      TH1D* h_higgsPtReweighting;

      std::vector<std::string> triggerPaths_;

      double photonPtThres_;
      double photonR9Thres_;
      bool applyCuBasedPhotonID_;   
      bool applyMVAPhotonID_;
      bool applyHggPhotonID_;
      double deltaEtaJets_;
      double invMassJets_;

      std::vector<double> phoMvaThres_;
      std::vector<double> passElectronVetoThres_;
      std::vector<double> isoCorrEtaBinsCh_;
      std::vector<double> isoCorrAreaValCh_;
      std::vector<double> isoCorrEtaBinsNh_;
      std::vector<double> isoCorrAreaValNh_;
      std::vector<double> isoCorrEtaBinsPh_;
      std::vector<double> isoCorrAreaValPh_;
      std::vector<double> hadronicOverEmThres_;
      std::vector<double> sigmaIetaIetaThres_;
      std::vector<double> egChargedHadronIsoThres_;
      std::vector<double> egNeutralHadronIsoEBThres_;
      std::vector<double> egNeutralHadronIsoEEThres_;
      std::vector<double> egPhotonIsoEBThres_;
      std::vector<double> egPhotonIsoEEThres_;

      std::vector<double> isoCorrEtaBinsHgg_;
      std::vector<double> isoCorrAreaValHgg_;
      double egChargedHadronIsoHggThres_;
      double egChargedHadronIsoOverPtHggThres_;
      double hadronicOverEmHggThres_;
      double photonPtHggThres_;
      double photonIdMvaDHggThres_;
      std::vector<double> r9ThresMinHgg_;
      std::vector<double> r9ThresMaxHgg_;
      double passElectronVetoHggThres_;
      std::vector<double> pfPhoIso03HggThres_;
      std::vector<double> trkSumPtHollowConeDR03HggThres_;
      std::vector<double> sigmaIetaIetaHggThres_;

      double metThres_;
      std::vector<double> metPhiCorr0_;
      std::vector<double> metPhiCorr1_;
      double dPhiPhotonMetThres_;

      bool useStdMuonID_;
      double muonEtaThres_;
      double leptonPtThres_;
      double muMiniIsoSumRelThres_;
      double muPFIsoSumRelThres_;
      
      bool useStdElectronID_;
      std::vector<double> electronEtaThres_;
      std::vector<double> elMiniIsoThres_;
      std::vector<double> impactParam_;
      bool useElectronMVARecipe_;
      bool useElectronLooseID_;

      std::vector<edm::InputTag> inputTagJets_;

      double jetPtThres_;
      double jetEtaThres_;
      double drJetPhoCut_;

      bool usePuJetID_;
      bool useJetID_;
      bool merge3rdJet_;
      double thirdJetDRCut_;
      double rmsforwardCut_;
      std::string jetIDLevel_;
      double drJetPhoVBFCut_;
      std::vector<double> pujid_wp_pt_bin_1_;
      std::vector<double> pujid_wp_pt_bin_2_;
      std::vector<double> pujid_wp_pt_bin_3_;
      
      edm::Service<TFileService> iFile;
      
      TH1D *h_numberOfEvents;
      TH1D *h_Efficiency;

      TH1D *h_Pho_Energy_noCut;
      TH1D *h_Pho_Pt_noCut;
      TH1D *h_Pho_Eta_noCut;        
      TH1D *h_Pho_Phi_noCut;            
      TH1D *h_Pho_Energy;      
      TH1D *h_Pho_Pt;          
      TH1D *h_Pho_Eta;         
      TH1D *h_Pho_Phi;         

      TH1D *h_MET_Pt_noCut;       
      TH1D *h_MET_Phi_noCut;   
      TH1D *h_MET_Pt;         
      TH1D *h_MET_Phi;      

      TH1D *h_dPhiMETPho_noCut;            
      TH1D *h_MtMETPho_noCut;   
      TH1D *h_dPhiMETPho;            
      TH1D *h_MtMETPho;             

      TH1D *h_Jet_Energy;            
      TH1D *h_Jet_Pt;                
      TH1D *h_Jet_Eta;               
      TH1D *h_Jet_Phi;               

      TH1D *h_dPhiJetPho;            
      TH1D *h_dEtaJetPho;            
      TH1D *h_dRJetPho;    

      TH1D *h_nJets_noCut; 
      TH1D *h_nJets;  
      TH1D *h_nVBFJets;       

      unsigned int HLT_pathsSize = 0;
      std::vector<std::string> HLT_Names; 
      std::vector<TH1D*> HLT_PhotonPt_efficiency_tmp; 
      std::vector<TH1D*> HLT_MET_efficiency_tmp; 
      std::vector<TH1D*> HLT_PhotonPt_efficiency; 
      std::vector<TH1D*> HLT_MET_efficiency; 
      TH1D* HLT_PhotonPt_denum;
      TH1D* HLT_MET_denum;

      double nTot; 
      double nTot_step0;
      double nTotSelected_step0;
      double nTotSelected_step1;
      double nTotSelected_step2;
      double nTotSelected_step3;
      double nTotSelected_step4;
      double nTotSelected_step5;
      double nTotSelected_step6;
      double nTotSelected_step7;
      double nTotSelected_step8;
};

#endif


