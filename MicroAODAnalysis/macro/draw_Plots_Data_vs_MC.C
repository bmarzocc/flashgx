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

//#include "FPCanvasStyle.C"
//#include "setStyle.C"
#include "centerOfMassEnergy.h"
#include "CMS_lumi.h"

#include<iostream>
#include<string>
#include<fstream>

using namespace std;

void drawTH1dataMCstack(TH1D* h1,TH1D* h2, vector<TH1D*> vecMC, 
			const string& xAxisNameTmp, const string& yAxisName, const string& canvasName, 
			const string& outputDIR, 
			const string& legEntry1, const vector<string>& vecLegEntryMC, 
			const string& ratioPadYaxisName, const Double_t lumi, const Int_t rebinFactor, 
			const Bool_t normalizeMCToData,
			const Int_t draw_both0_noLog1_onlyLog2,
			const Double_t minFractionToBeInLegend,
			const Int_t fillStyle
			);

Bool_t getAxisRangeFromUser(string& axisName, Double_t& min, Double_t& max, 
			    const string& axisNameTmp, 
			    const string& separator, 
			    const string& rangeSeparator
			    ); 

void myRebinHisto(TH1 *h, const Int_t rebinFactor);

void draw_Plots_Data_vs_MC() {
  
  gStyle->SetOptStat(0);

  TH1D* h_tmp;
  TKey* key ;
  TObject* obj ;
  TList* list;

  // inputs
  TFile* inFile_data = TFile::Open("/eos/cms/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v4/prod-uAOD-300-66-g644fd526/SinglePhoton_ntuples.root");
  TDirectory* dir_data = inFile_data->GetDirectory("flashgxanalysis");
  vector<TH1D*> histos_data;
 
  list = dir_data->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next_data(list) ;
  
  while ( (key = (TKey*)next_data()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1D*)inFile_data->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_data.push_back(h_tmp);
    }
  }

  TFile* inFile_sig = TFile::Open("/eos/cms/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VBFHToGX_M125_13TeV_amcatnlo_pythia8_RunIIFall17_PreMixing_MicroAOD_ntuples.root");
  TDirectory* dir_sig = inFile_sig->GetDirectory("flashgxanalysis");
  vector<TH1D*> histos_sig;
 
  list = dir_sig->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next_sig(list) ;
  
  while ( (key = (TKey*)next_sig()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1D*)inFile_sig->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_sig.push_back(h_tmp);
    }
  }

  TFile* inFile_WJets = TFile::Open("/eos/cms/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/WJetsToLNu_ntuples.root");
  TDirectory* dir_WJets = inFile_WJets->GetDirectory("flashgxanalysis");
  vector<TH1D*> histos_WJets;
 
  list = dir_WJets->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next_WJets(list) ;
  
  while ( (key = (TKey*)next_WJets()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1D*)inFile_WJets->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_WJets.push_back(h_tmp);
    }
  }

  TFile* inFile_WToENu = TFile::Open("/eos/cms/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/WToENu_ntuples.root");
  TDirectory* dir_WToENu = inFile_WToENu->GetDirectory("flashgxanalysis");
  vector<TH1D*> histos_WToENu;
 
  list = dir_WToENu->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next_WToENu(list) ;
  
  while ( (key = (TKey*)next_WToENu()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1D*)inFile_WToENu->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_WToENu.push_back(h_tmp);
    }
  }

  TFile* inFile_ZJets = TFile::Open("/eos/cms/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/ZJetsToNuNu_ntuples.root");
  TDirectory* dir_ZJets = inFile_ZJets->GetDirectory("flashgxanalysis");
  vector<TH1D*> histos_ZJets;

  list = dir_ZJets->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next_ZJets(list) ;
  
  while ( (key = (TKey*)next_ZJets()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1D*)inFile_ZJets->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_ZJets.push_back(h_tmp);
    }
  }

  TFile* inFile_VGGJets = TFile::Open("/eos/cms/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v2/prod-uAOD-300-65-g97bd5bfd/VGGJets_ntuples.root");
  TDirectory* dir_VGGJets = inFile_VGGJets->GetDirectory("flashgxanalysis");
  vector<TH1D*> histos_VGGJets;

  list = dir_VGGJets->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next_VGGJets(list) ;
  
  while ( (key = (TKey*)next_VGGJets()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1D*)inFile_VGGJets->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_VGGJets.push_back(h_tmp);
    }
  }

  TFile* inFile_DYJets = TFile::Open("/eos/cms/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v3/prod-uAOD-300-65-g97bd5bfd/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_ntuples.root");
  TDirectory* dir_DYJets = inFile_DYJets->GetDirectory("flashgxanalysis");
  vector<TH1D*> histos_DYJets;

  list = dir_DYJets->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next_DYJets(list) ;
  
  while ( (key = (TKey*)next_DYJets()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1D*)inFile_DYJets->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_DYJets.push_back(h_tmp);
    }
  }
  
  TFile* inFile_DiPhotonJetsBox = TFile::Open("/eos/cms/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v3/prod-uAOD-300-65-g97bd5bfd/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_ntuples.root");
  TDirectory* dir_DiPhotonJetsBox = inFile_DiPhotonJetsBox->GetDirectory("flashgxanalysis");
  vector<TH1D*> histos_DiPhotonJetsBox;

  list = dir_DiPhotonJetsBox->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next_DiPhotonJetsBox(list) ;
  
  while ( (key = (TKey*)next_DiPhotonJetsBox()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1D*)inFile_DiPhotonJetsBox->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_DiPhotonJetsBox.push_back(h_tmp);
    }
  }
  
  TFile* inFile_QCD = TFile::Open("/eos/cms/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v3/prod-uAOD-300-65-g97bd5bfd/QCD_ntuples.root");
  TDirectory* dir_QCD = inFile_QCD->GetDirectory("flashgxanalysis");
  vector<TH1D*> histos_QCD;

  list = dir_QCD->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next_QCD(list) ;
  
  while ( (key = (TKey*)next_QCD()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1D*)inFile_QCD->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_QCD.push_back(h_tmp);
    }
  }

  TFile* inFile_GJet = TFile::Open("/eos/cms/store/group/phys_higgs/HiggsExo/HToGX/bmarzocc/flashgg/RunIIFall17_DarkPhoton_v3/prod-uAOD-300-65-g97bd5bfd/GJet_ntuples.root");
  TDirectory* dir_GJet = inFile_GJet->GetDirectory("flashgxanalysis");
  vector<TH1D*> histos_GJet;

  list = dir_GJet->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next_GJet(list) ;
  
  while ( (key = (TKey*)next_GJet()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1D*)inFile_GJet->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_GJet.push_back(h_tmp);
    }
  }
  
  gROOT->SetBatch(kTRUE);
  
  vector<string> legendEntriesMC;
  legendEntriesMC.resize(7);
  legendEntriesMC[0] = "WJets"; 
  legendEntriesMC[1] = "W#rightarrow e#nu";
  legendEntriesMC[2] = "ZJets";
  legendEntriesMC[3] = "W/Z + GGJets";
  legendEntriesMC[4] = "DYJets";
  legendEntriesMC[5] = "DiPhotonJets";
  legendEntriesMC[6] = "GJets";
  //legendEntriesMC[7] = "QCD";

  vector<TH1D*> vecMC;
 
  for(unsigned int ii=0; ii<histos_sig.size(); ii++)
  {

      if(string(histos_WJets.at(ii)->GetName()).find("h_Efficiency")!=std::string::npos) continue;

      string xTitle = histos_WJets.at(ii)->GetName();
      xTitle.erase(0,2);
      if(xTitle.find("_Final")!=std::string::npos) xTitle.erase(xTitle.end()-6,xTitle.end());
      if(xTitle.find("_noCut")!=std::string::npos) xTitle.erase(xTitle.end()-6,xTitle.end());
      if(xTitle.find("Pt")!=std::string::npos) xTitle = xTitle + string(" (GeV)");
      if(xTitle.find("MtMETPho")!=std::string::npos) xTitle = xTitle + string(" (GeV)");
      if(xTitle.find("invMass")!=std::string::npos) xTitle = xTitle + string(" (GeV)");

      vecMC.resize(7);

      vecMC[0] = histos_WJets.at(ii);
      vecMC[0]->SetFillColor(kRed);
      vecMC[1] = histos_WToENu.at(ii);
      vecMC[1]->SetFillColor(kBlue);
      vecMC[2] = histos_ZJets.at(ii);
      vecMC[2]->SetFillColor(kCyan+1);
      vecMC[3] = histos_VGGJets.at(ii);
      vecMC[3]->SetFillColor(kYellow);
      vecMC[4] = histos_DYJets.at(ii);
      vecMC[4]->SetFillColor(kGreen);
      vecMC[5] = histos_DiPhotonJetsBox.at(ii);
      vecMC[5]->SetFillColor(kOrange+1);
      vecMC[6] = histos_GJet.at(ii);
      vecMC[6]->SetFillColor(kViolet);
      //vecMC[7] = histos_QCD.at(ii);
      //vecMC[7]->SetFillColor(kGreen+2);
 
      TH1D* h1 = (TH1D*)histos_data.at(ii)->Clone();   
      TH1D* h2 = (TH1D*)histos_sig.at(ii)->Clone();
      
 
      drawTH1dataMCstack(h1,h2,vecMC,
			 xTitle.c_str(),vecMC.at(0)->GetYaxis()->GetTitle(),histos_WJets.at(ii)->GetName(),
                         "","Data",legendEntriesMC,"Data/MC",7.7,2,true,1,0.0,3001);

      drawTH1dataMCstack(h1,h2,vecMC,
			 xTitle.c_str(),vecMC.at(0)->GetYaxis()->GetTitle(),histos_WJets.at(ii)->GetName(),
                         "","Data",legendEntriesMC,"Data/MC",7.7,2,true,0,0.0,3001);

      vecMC.clear();      
  }
  
}

void drawTH1dataMCstack(TH1D* h1 = NULL, TH1D* h2 = NULL, vector<TH1D*> vecMC = {}, 
			const string& xAxisNameTmp = "", const string& yAxisName = "Events", const string& canvasName = "default", 
			const string& outputDIR = "./", 
			const string& legEntry1 = "data", const vector<string>& vecLegEntryMC = {""}, 
			const string& ratioPadYaxisName = "data/MC", const Double_t lumi = -1.0, const Int_t rebinFactor = 1, 
			const Bool_t normalizeMCToData = false,
			const Int_t draw_both0_noLog1_onlyLog2 = 0,
			const Double_t minFractionToBeInLegend = 0.001,
			const Int_t fillStyle = 3001
			)
{

  string xAxisName = "";
  Double_t xmin = 0;
  Double_t xmax = 0;
  Bool_t setXAxisRangeFromUser = getAxisRangeFromUser(xAxisName, xmin, xmax, xAxisNameTmp,"::",",");

  string yAxisNameRatio = "";
  Double_t yminRatio = 0.5;
  Double_t ymaxRatio = 1.5;
  Bool_t setYAxisRatioRangeFromUser = getAxisRangeFromUser(yAxisNameRatio, yminRatio, ymaxRatio, ratioPadYaxisName,"::",",");

  // cout << "xAxisName = " << xAxisName << "   xmin = " << xmin << "  xmax = " << xmax << endl;

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 

  if (vecMC.size() != vecLegEntryMC.size()) {
    cout << "Error: legend has different number of entries than stack elements, please check. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  for (UInt_t i = 0; i < vecMC.size(); i++) {
    vecMC[i]->SetStats(0);
  }
  
  //Int_t colorList[] = {kCyan, kViolet, kBlue, kRed, kYellow, kGreen, kOrange+1, kCyan+2, kGreen+2, kGray}; 
  // the first color is for the main object. This array may contain more values than vecMC.size()
  
  Double_t dataNorm = h1->Integral();
  Double_t stackNorm = 0.0; 

  if (normalizeMCToData) {
    for (UInt_t ibin = 0; ibin < vecMC.size(); ibin++) {
      stackNorm += vecMC[ibin]->Integral();
    }
  }

  THStack* hMCstack = new THStack("hMCstack","");
  for (UInt_t j = 0; j < vecMC.size(); j++) {
    float norm = 7742.163000/40000.;
    if(!normalizeMCToData) vecMC[j]->Scale(norm);
    vecMC[j]->SetFillStyle(fillStyle);
    if(string(vecMC[j]->GetName()).find("nJets")==std::string::npos && string(vecMC[j]->GetName()).find("nVtx")==std::string::npos && string(vecMC[j]->GetName()).find("numberOfEvents")==std::string::npos && string(vecMC[j]->GetName()).find("h_Efficiency_num")==std::string::npos && string(vecMC[j]->GetName()).find("h_Efficiency_denum")==std::string::npos) myRebinHisto(vecMC[j],rebinFactor);
    if(normalizeMCToData) vecMC[j]->Scale(dataNorm/stackNorm);
    hMCstack->Add(vecMC[int(vecMC.size())-j-1]);  // add last element as the first one (last element added in stack goes on top)
  }

  cout << "N-Histos: " << hMCstack->GetNhists() << endl;

  THStack* hMCstack2 = (THStack*)hMCstack->Clone();
  TH1D* stackCopy = (TH1D*) hMCstack2->GetStack()->Last(); // used to make ratioplot without affecting the plot and setting maximum

  if (h1 == NULL) h1 = (TH1D*) stackCopy->Clone("pseudoData");
  if(string(h1->GetName()).find("nJets")==std::string::npos && string(h1->GetName()).find("nVtx")==std::string::npos && string(h1->GetName()).find("numberOfEvents")==std::string::npos && string(h1->GetName()).find("h_Efficiency_num")==std::string::npos && string(h1->GetName()).find("h_Efficiency_denum")==std::string::npos) myRebinHisto(h1,rebinFactor);

  if (h2 == NULL) h2 = (TH1D*) stackCopy->Clone("pseudoData");
  if(string(h2->GetName()).find("nJets")==std::string::npos && string(h2->GetName()).find("nVtx")==std::string::npos && string(h2->GetName()).find("numberOfEvents")==std::string::npos && string(h2->GetName()).find("h_Efficiency_num")==std::string::npos && string(h2->GetName()).find("h_Efficiency_denum")==std::string::npos) myRebinHisto(h2,rebinFactor);

  if (normalizeMCToData) h2->Scale(dataNorm/h2->Integral());

  h1->SetStats(0);

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);
  canvas->SetLeftMargin(0.16);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetLeftMargin(0.16);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);


  TH1* frame =  (TH1*) h1->Clone("frame");
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->SetStats(0);

  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);

  h2->SetLineColor(kRed);
  h2->SetMarkerColor(kRed);
  h2->SetMarkerStyle(20);
  //h2->SetMarkerSize(2);
  h2->SetLineWidth(2);

  stackCopy->SetFillColor(kGray+1);
  stackCopy->SetMarkerStyle(0);

  h1->GetXaxis()->SetLabelSize(0);
  h1->GetXaxis()->SetTitle(0);
  h1->GetYaxis()->SetTitle(yAxisName.c_str());
  h1->GetYaxis()->SetTitleOffset(1.2);
  // h1->GetYaxis()->SetTitleOffset(0.8);  // was 1.03 without setting also the size
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetLabelSize(0.04);
  h1->GetYaxis()->SetRangeUser(0.0, max(h1->GetMaximum(),stackCopy->GetMaximum()) * 1.2);
  if (setXAxisRangeFromUser) h1->GetXaxis()->SetRangeUser(xmin,xmax);
  hMCstack->Draw("HIST");
  if(string(h1->GetName()).find("Final")==std::string::npos) h1->Draw("EP");
  hMCstack->Draw("HIST SAME");
  stackCopy->Draw("e2same");
  //h2->Draw("HIST SAME");
  if(string(h1->GetName()).find("Final")==std::string::npos) h1->Draw("EP SAME");

  TLegend leg (0.65,0.7,0.95,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(h1,legEntry1.c_str(),"PLE");
  leg.AddEntry(h2,"Signal","L");
  for (UInt_t i = 0; i < vecMC.size(); i++) {
    if (vecMC[i]->Integral()/stackCopy->Integral() > minFractionToBeInLegend)
      leg.AddEntry(vecMC[i],vecLegEntryMC[i].c_str(),"F");
  }
  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  //  CMS_lumi(canvas,Form("%.1f",lumi));
  if (lumi < 0) CMS_lumi(canvas,"",false,false);
  else CMS_lumi(canvas,Form("%.1f",lumi),false,false);
  setTDRStyle();

  pad2->Draw();
  pad2->cd();

  frame->Reset("ICES");
  frame->GetYaxis()->SetRangeUser(yminRatio,ymaxRatio);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle(yAxisNameRatio.c_str());
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->CenterTitle();
  if (setXAxisRangeFromUser) frame->GetXaxis()->SetRangeUser(xmin,xmax);
  frame->GetXaxis()->SetTitle(xAxisName.c_str());
  // frame->GetXaxis()->SetTitleOffset(0.8);
  frame->GetXaxis()->SetTitleSize(0.05);

  TH1D* ratio = (TH1D*) h1->Clone("ratio");
  TH1D* den_noerr = (TH1D*) stackCopy->Clone("den_noerr");
  TH1D* den = (TH1D*) stackCopy->Clone("den");
  for(int iBin = 1; iBin < den->GetNbinsX()+1; iBin++)
    den_noerr->SetBinError(iBin,0.);

  ratio->Divide(den_noerr);
  den->Divide(den_noerr);
  den->SetFillColor(kGray+1);
  den->SetMarkerStyle(0);
  frame->Draw();
  ratio->SetMarkerSize(0.85);
  if(string(h1->GetName()).find("Final")==std::string::npos) ratio->Draw("EPsame");
  den->Draw("E2same");

  TF1* line = new TF1("horiz_line","1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  if(string(h1->GetName()).find("Final")==std::string::npos) ratio->Draw("EPsame");
  pad2->RedrawAxis("sameaxis");

  // Calculate chi2                                                                                        
  /* double chi2 = h1->Chi2Test(stackCopy,"CHI2/NDF WW"); */
  /* TLegend leg2 (0.14,0.25,0.32,0.28,NULL,"brNDC"); */
  /* leg2.SetFillColor(0); */
  /* leg2.SetFillStyle(1); */
  /* leg2.SetBorderSize(0); */
  /* leg2.SetLineColor(0); */
  /* leg2.AddEntry((TObject*)0,Form("#chi^{2}/ndf = %.2f",chi2),""); */
  /* leg2.Draw("same"); */

  if (draw_both0_noLog1_onlyLog2 != 2) {
    canvas->SaveAs((outputDIR + canvasName + ".png").c_str());
    canvas->SaveAs((outputDIR + canvasName + ".pdf").c_str());
  }

  /* Double_t minY_log = max(0.001, min(h1->GetBinContent(h1->GetMinimumBin()),stackCopy->GetBinContent(stackCopy->GetMinimumBin())) ) * 0.01; */
  /* if ( fabs(minY_log) < 0.00001) minY_log = 0.001;  // avoid zero */
  /* h1->GetYaxis()->SetRangeUser(minY_log * 0.8, max(h1->GetMaximum(),stackCopy->GetMaximum())*100); */
  //Double_t minY_log = max(0.001,min(h1->GetMinimum(),stackCopy->GetMinimum())*0.8);
  //if (min(h1->GetBinContent(h1->GetMinimumBin()),stackCopy->GetBinContent(stackCopy->GetMinimumBin()) > 0) minY_log = 0.5;

  /* Double_t minY_log = -1.0; */
  /* TH1D* lowerStackObject = (TH1D*) hMCstack->GetStack()->First(); // */

  if (draw_both0_noLog1_onlyLog2 != 1) {

    Double_t minY_log = 1e20;
    TH1D* lowerStackObject = (TH1D*) hMCstack->GetStack()->First();
    /* if (lowerStackObject->GetBinContent(stackCopy->GetMinimumBin()) > 0.000001 ) { */
    /*   // if no empty bin in stack */
    /*   minY_log = 0.05 * lowerStackObject->GetBinContent(lowerStackObject->GetMinimumBin()); */
    /* } else { */
    /*   minY_log = 0.01 * stackCopy->GetBinContent(stackCopy->GetMinimumBin()); */
    /* } */

    for (Int_t ibin = 0; ibin <= lowerStackObject->GetNbinsX(); ibin++ ) {
      if (lowerStackObject->GetBinContent(ibin) > 0.0000001 && minY_log > lowerStackObject->GetBinContent(ibin)) minY_log = lowerStackObject->GetBinContent(ibin);
    }

    if (minY_log < 0.000001) minY_log = 0.1;
    minY_log = 0.05 * minY_log;

    h1->GetYaxis()->SetRangeUser(minY_log,max(h1->GetMaximum(),stackCopy->GetMaximum())*100);
    canvas->SetLogy();
    pad2->RedrawAxis("sameaxis");
    /* if (lumi < 0) CMS_lumi(canvas,"",true,false); */
    /* else CMS_lumi(canvas,Form("%.1f",lumi),true,false); */
    canvas->SaveAs((outputDIR + canvasName + "_logY.png").c_str());
    canvas->SaveAs((outputDIR + canvasName + "_logY.pdf").c_str());
    canvas->SetLogy(0);
  }

  delete hMCstack;
  delete line;
  delete pad2;
  delete canvas;

}

Bool_t getAxisRangeFromUser(string& axisName, Double_t& min, Double_t& max, 
			    const string& axisNameTmp = "", 
			    const string& separator = "::", 
			    const string& rangeSeparator = ","
			    ) {
  
  Bool_t setXAxisRangeFromUser = false;
  size_t pos = axisNameTmp.find(separator);
    
  if (pos != string::npos) {
    string xrange = "";
    setXAxisRangeFromUser = true;
    axisName.assign(axisNameTmp, 0, pos);
    xrange.assign(axisNameTmp, pos + separator.size(), string::npos);
    pos = xrange.find(rangeSeparator);
    string numString = "";
    numString.assign(xrange,0,pos);
    min = std::stod(numString);
    numString.assign(xrange,pos + rangeSeparator.size(), string::npos);
    max = std::stod(numString);
  } else {
    axisName = axisNameTmp;
  }

  return setXAxisRangeFromUser;

}

void myRebinHisto(TH1 *h, const Int_t rebinFactor = 1) {

  if (rebinFactor != 1) {
    h->Rebin(rebinFactor);
    //if ( (h->GetNbinsX() % rebinFactor) != 0) myAddOverflowInLastBin(h);  // better not to add overflow to last bin inside this function
  }

}

