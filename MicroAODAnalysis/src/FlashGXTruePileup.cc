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

#include "flashgx/MicroAODAnalysis/interface/FlashGXTruePileup.h"

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
FlashGXTruePileup::FlashGXTruePileup(const edm::ParameterSet& iConfig)
{

  puInfoToken_           = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupInfo"));
  h_TruePileup           = iFile->make<TH1D>("h_TruePileup","h_TruePileup",100,-0.5,99.5);

}

FlashGXTruePileup::~FlashGXTruePileup()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void FlashGXTruePileup::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{

  edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
  ev.getByToken(puInfoToken_,puInfo);
      
  for(std::vector<PileupSummaryInfo>::const_iterator pu = puInfo->begin(); pu != puInfo->end() ; ++pu) 
  {     
      if(pu->getBunchCrossing()==0)
      {
          h_TruePileup->Fill(pu->getTrueNumInteractions());
          break;
      }
       
  }

}

        
void FlashGXTruePileup::beginJob()
{

}

void FlashGXTruePileup::endJob() 
{

}

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FlashGXTruePileup);

