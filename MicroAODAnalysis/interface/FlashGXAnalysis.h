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
//#include "PhysicsTools/Utilities/macros/setTDRStyle.C"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// flashgg includes
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Photon.h"
#include "flashgg/DataFormats/interface/SingleVertexView.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"

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

#define MAXVEC 20

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
      std::map<std::string,std::vector<double> > mcInfos(string mcInfos_);
      double ptWeight(Handle<View<reco::GenParticle> > genParticles);
      double puWeight(edm::Handle<std::vector<PileupSummaryInfo> > puInfo, TH1D* mcPuHist, TH1D* dataPuHist);
      bool getHLTResults(edm::Handle<edm::TriggerResults> trigResults,const edm::TriggerNames& trigNames,std::string path);
      bool isGoodPhotonHgg(edm::Ptr<flashgg::Photon> pho, edm::Ptr<reco::Vertex> vtx, double rho); 
      bool isGoodPhotonCutBased(edm::Ptr<flashgg::Photon> ipho, double rho);
      bool isGoodPhotonMVA(edm::Ptr<flashgg::Photon> ipho);
      std::pair<int,int> selectedPhotonIndex(Handle<View<flashgg::Photon> > photons,edm::Ptr<reco::Vertex> vtx,double rho);
      bool isInEcal(edm::Ptr<flashgg::Photon> ipho);
      bool isGoodTrigger(TGraphAsymmErrors* efficiency, double threshold);
      TH1D* fixBins(TH1D* h);
      void saveTrigger(TGraphAsymmErrors* efficiency_tmp, TGraphAsymmErrors* efficiency);
      double phoIsoCorr(edm::Ptr<flashgg::Photon> ipho, double rho);
      double metPhiCorr(edm::Ptr<flashgg::Met> met);
      double computeMtMETPhoton(edm::Ptr<flashgg::Photon> ipho, edm::Ptr<flashgg::Met> met);
      std::vector<reco::Candidate::LorentzVector> goodJets(JetCollectionVector Jets, edm::Ptr<flashgg::Photon> pho, int iColl);
      std::vector<reco::Candidate::LorentzVector> vbfJets(JetCollectionVector Jets, edm::Ptr<flashgg::Photon> pho, edm::Ptr<reco::Vertex> vtx, int iColl);
      void setTree(TTree *outTree);
      void fillTree(TTree *outTree,int irun,int ilumi,long int ievent,long int ibx,double wTot,double rho,Handle<View<reco::Vertex> > vertices,Handle<View<flashgg::Met> > mets,Handle<View<flashgg::Photon> > photons, int selPhoton_pos, JetCollectionVector Jets, std::vector<reco::Candidate::LorentzVector>* vbfJets, int iColl);

      edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;  
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerToken_;  
      edm::EDGetTokenT<double> rhoToken_;  
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticleToken_;  
      edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;  
      edm::EDGetTokenT<edm::View<flashgg::Photon> > photonToken_; 
      edm::EDGetTokenT<edm::View<flashgg::Muon> > muonToken_; 
      edm::EDGetTokenT<edm::View<flashgg::Electron> > electronToken_; 
      edm::EDGetTokenT<edm::View<flashgg::Met> > metToken_;
      std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > jetsToken_;

      std::string sampleType_;
      std::string sampleName_;
      std::string mcInfos_;
      std::map<std::string,std::vector<double> > infos_;

      TFile* mcPuFile_;
      std::string mcPuPath_;
      TH1D* mcPuHist_;

      TFile* dataPuFile_;
      std::string dataPuPath_;
      TH1D* dataPuHist_;

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
      vector<reco::Candidate::LorentzVector> vbfjets;

      double jetPtThres_;
      double jetEtaThres_;
      double drJetPhoCut_;
      bool usePuJetID_;
      bool useJetID_;
      bool merge3rdJet_;
      double thirdJetDRCut_;
      double rmsforwardCut_;
      std::string jetIDLevel_;
      std::vector<double> pujidEtaBins_;
      std::vector<double> pujid_wp_pt_bin_1_;
      std::vector<double> pujid_wp_pt_bin_2_;
      std::vector<double> pujid_wp_pt_bin_3_;
      std::vector<double> pujid_wp_pt_bin_4_;
      
      edm::Service<TFileService> iFile;
      
      TH1D *h_numberOfEvents;
      TH1D *h_Efficiency;
      TH1D *h_Efficiency_num;
      TH1D *h_Efficiency_denum;

      TH1D *h_nVtx;
      TH1D *h_dataPileup;
      TH1D *h_mcPileup;

      TH1D *h_Pho_Energy_noCut;
      TH1D *h_Pho_Pt_noCut;
      TH1D *h_Pho_Eta_noCut;        
      TH1D *h_Pho_Phi_noCut;            
      TH1D *h_Pho_Energy;      
      TH1D *h_Pho_Pt;          
      TH1D *h_Pho_Eta;         
      TH1D *h_Pho_Phi; 
      TH1D *h_Pho_Energy_Final;
      TH1D *h_Pho_Pt_Final;
      TH1D *h_Pho_Eta_Final;
      TH1D *h_Pho_Phi_Final;    

      TH1D *h_MET_Pt_noCut;       
      TH1D *h_MET_Phi_noCut;   
      TH1D *h_MET_Pt;         
      TH1D *h_MET_Phi;      
      TH1D *h_MET_Pt_Final;         
      TH1D *h_MET_Phi_Final;   
  
      TH1D *h_dPhiMETPho_noCut;            
      TH1D *h_MtMETPho_noCut;   
      TH1D *h_dPhiMETPho;            
      TH1D *h_MtMETPho;    
      TH1D *h_dPhiMETPho_Final;            
      TH1D *h_MtMETPho_Final;             

      TH1D *h_nJets_noCut;
      TH1D *h_Jet_Energy_noCut;            
      TH1D *h_Jet_Pt_noCut;                
      TH1D *h_Jet_Eta_noCut;               
      TH1D *h_Jet_Phi_noCut;               
      TH1D *h_dPhiJetPho_noCut;            
      TH1D *h_dEtaJetPho_noCut;            
      TH1D *h_dRJetPho_noCut;    

      TH1D *h_nJets; 
      TH1D *h_Jet_Energy;            
      TH1D *h_Jet_Pt;                
      TH1D *h_Jet_Eta;               
      TH1D *h_Jet_Phi;               
      TH1D *h_dPhiJetPho;            
      TH1D *h_dEtaJetPho;            
      TH1D *h_dRJetPho;      

      TH1D *h_VBFjet1_Energy;  
      TH1D *h_VBFjet1_Pt;
      TH1D *h_VBFjet1_Eta;
      TH1D *h_VBFjet1_Phi;
      TH1D *h_VBFjet1_Energy_Final;  
      TH1D *h_VBFjet1_Pt_Final;
      TH1D *h_VBFjet1_Eta_Final;
      TH1D *h_VBFjet1_Phi_Final;
     
      
      TH1D *h_VBFjet2_Energy;  
      TH1D *h_VBFjet2_Pt;
      TH1D *h_VBFjet2_Eta;
      TH1D *h_VBFjet2_Phi;
      TH1D *h_VBFjet2_Energy_Final;  
      TH1D *h_VBFjet2_Pt_Final;
      TH1D *h_VBFjet2_Eta_Final;
      TH1D *h_VBFjet2_Phi_Final;

      TH1D *h_VBFjets_invMass;
      TH1D *h_VBFjets_dEta;
      TH1D *h_VBFjets_invMass_Final;
      TH1D *h_VBFjets_dEta_Final;
      
      unsigned int HLT_pathsSize = 0;
      std::vector<std::string> HLT_Names; 
      std::vector<TH1D*> HLT_PhotonPt_efficiency_tmp;
      TH1D* HLT_PhotonPt_efficiency_AllOR_tmp; 
      std::vector<TH1D*> HLT_MET_efficiency_tmp; 
      TH1D* HLT_MET_efficiency_AllOR_tmp; 
      std::vector<TGraphAsymmErrors*> HLT_PhotonPt_efficiency; 
      std::vector<TGraphAsymmErrors*> HLT_MET_efficiency; 
      TGraphAsymmErrors* HLT_PhotonPt_efficiency_AllOR; 
      TGraphAsymmErrors* HLT_MET_efficiency_AllOR; 
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

      TTree *outTree; 
      int run;
      int lumi;
      long int event; 
      long int bx;
      double weight;
      int nVtx;
      double MET_pt; 
      double MET_phi; 
      double deltaPhi_MET_photon[MAXVEC]; 
      double Mt_MET_photon[MAXVEC]; 
      double photon_pt[MAXVEC]; 
      double photon_energy[MAXVEC]; 
      double photon_phi[MAXVEC]; 
      double photon_r9[MAXVEC]; 
      double photon_eta[MAXVEC]; 
      double photon_egMVA[MAXVEC]; 
      bool photon_passElectronVeto[MAXVEC]; 
      bool photon_passCutID[MAXVEC];
      bool photon_passHggID[MAXVEC]; 
      bool photon_passMVAID[MAXVEC]; 
      int selPhoton_index;
      double jet_pt[MAXVEC]; 
      double jet_energy[MAXVEC]; 
      double jet_phi[MAXVEC]; 
      double jet_eta[MAXVEC];
      double deltaPhi_jet_photon[MAXVEC]; 
      double deltaEta_jet_photon[MAXVEC]; 
      double deltaR_jet_photon[MAXVEC];
      int n_jets; 
      double vbfjet_pt[MAXVEC]; 
      double vbfjet_energy[MAXVEC]; 
      double vbfjet_phi[MAXVEC]; 
      double vbfjet_eta[MAXVEC];
      double deltaPhi_vbfjet_photon[MAXVEC]; 
      double deltaEta_vbfjet_photon[MAXVEC]; 
      double deltaR_vbfjet_photon[MAXVEC];
      int n_vbfjets; 
      double deltaEta_selVbfJets;  
      double invMass_selVbfJets;  
};

#endif


