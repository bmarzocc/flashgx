#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TKey.h"

#include "FPCanvasStyle.C"
#include "setStyle.C"

#include<iostream>
#include<string>
#include<fstream>

using namespace std;

void check_events() {
  
  // inputs
  TFile* inFile = TFile::Open("/eos/cms//store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v4/prod-uAOD-300-66-g644fd526/SinglePhoton_ntuples.root");
  TTree* tree = (TTree*)inFile->Get("flashgxanalysis/FlashGXTree");
  
  Int_t           run;
  Int_t           lumi;
  Long64_t        event;
  Long64_t        bx;

  TBranch        *b_run;   //!
  TBranch        *b_lumi;   //!
  TBranch        *b_event;   //!
  TBranch        *b_bx;   //!

  tree->SetBranchAddress("run", &run, &b_run);
  tree->SetBranchAddress("lumi", &lumi, &b_lumi);
  tree->SetBranchAddress("event", &event, &b_event);
  tree->SetBranchAddress("bx", &bx, &b_bx);

  Long64_t nentries = tree->GetEntries();
  for(Long64_t i=0;i<nentries;i++) {
     tree->GetEntry(i);
     std::cout << run << " " << lumi << " " << event << std::endl; 
  }
  
}
