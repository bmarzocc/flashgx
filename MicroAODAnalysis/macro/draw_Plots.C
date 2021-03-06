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

void drawHisto(TH1F* h, string outputDir, string unit, string label1);
void drawEfficiency(TGraphAsymmErrors* eff, string outputDir, string unit, string label1);

void draw_Plots() {
  
  gStyle->SetOptStat(0);

  // inputs
  TFile* inFile1 = TFile::Open("../FlashGXAnalysis.root");
  
  TDirectory* dir1 = inFile1->GetDirectory("flashgxanalysis");
  
  vector<TH1F*> histos_f1;
  vector<TGraphAsymmErrors*> graph_f1;
 
  TH1F* h_tmp;
  TGraphAsymmErrors* g_tmp;

  TList* list = dir1->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list) ;
  TKey* key ;
  TObject* obj ;
      
  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")==0)
         || (obj->InheritsFrom("TH2"))
	 || (obj->InheritsFrom("TH1")) 
       ) {
      h_tmp = (TH1F*)inFile1->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      histos_f1.push_back(h_tmp);
    }
    if ((obj->InheritsFrom("TGraphAsymmErrors"))
       ) {
      g_tmp = (TGraphAsymmErrors*)inFile1->Get((string("flashgxanalysis/")+string(obj->GetName())).c_str());
      graph_f1.push_back(g_tmp);
    }
    
  }

  gROOT->SetBatch(kTRUE);
  for(unsigned int ii=0; ii<histos_f1.size(); ii++)
  {
      if(std::string(histos_f1.at(ii)->GetName()).find("h_n")==std::string::npos && std::string(histos_f1.at(ii)->GetName()).find("efficiency")==std::string::npos) histos_f1.at(ii)->Rebin(2);
      if(std::string(histos_f1.at(ii)->GetName()).find("efficiency")==std::string::npos) drawHisto(histos_f1.at(ii), string("../output"), string("GeV"), string("VBFHToGX"));
  }
  for(unsigned int ii=0; ii<graph_f1.size(); ii++)
  {
      drawEfficiency(graph_f1.at(ii), string("../output"), string("GeV"), string("VBFHToGX"));
  }
  
}

void drawHisto(TH1F* h, string outputDir = "../output", string unit = "GeV", string label1 = "file1")
{
     h->SetLineColor(kBlack);
     h->SetLineWidth(2);
     
     // plotting
     setStyle();

     int nbins = h->GetNbinsX();
     float xmin = h->GetBinCenter(1)-h->GetBinWidth(1)/2;
     float xmax = h->GetBinCenter(nbins)+h->GetBinWidth(nbins)/2;
     float ymax = 1.1*h->GetBinContent(h->GetMaximumBin());
     TH2F* H2 = new TH2F("H2","",nbins,xmin,xmax,1000,0.,ymax);
     if(string(h->GetName()).find("eta")!=std::string::npos || string(h->GetName()).find("phi")!=std::string::npos || string(h->GetName()).find("delta")!=std::string::npos || string(h->GetName()).find("Fraction")!=std::string::npos) H2->GetXaxis()->SetTitle((string(h->GetName())).c_str());
     else H2->GetXaxis()->SetTitle((string(h->GetName())+" ("+unit+")").c_str());
     H2->GetYaxis()->SetTitle("events");

     float maximum = h->GetMaximum();
     if(maximum < h->GetMaximum()) maximum = h->GetMaximum();
     H2->GetYaxis()->SetRangeUser(0.9,1.5*maximum);

     TLegend *leg;
     leg = new TLegend(0.68,0.70,0.88,0.90);
     leg->SetFillStyle(0);
     leg->SetBorderSize(0);
     leg->SetTextSize(0.03);
     leg->SetFillColor(0);
     leg->AddEntry(h, label1.c_str(), "l");
     
     TCanvas* c1 = new TCanvas("c1","c1",1);
     FPCanvasStyle(c1);
     H2->Draw();
     if(string(h->GetName()).find("HLT")!=std::string::npos) h->Draw("P,same");
     else h->Draw("H,same");
     //leg->Draw("same");
     c1->SaveAs((outputDir+"/"+string(h->GetName())+".png").c_str());
     c1->SaveAs((outputDir+"/"+string(h->GetName())+".pdf").c_str());

     TCanvas* c2 = new TCanvas("c2","c2",1);
     FPCanvasStyle(c2);
     c2->SetLogy();
     H2->Draw();
     if(string(h->GetName()).find("HLT")!=std::string::npos) h->Draw("P,same");
     else h->Draw("H,same");
     //leg->Draw("same");
     c2->SaveAs((outputDir+"/"+string(h->GetName())+"_log.png").c_str());
     c2->SaveAs((outputDir+"/"+string(h->GetName())+"_log.pdf").c_str());

     delete c1;
     delete c2;
     delete H2;
     delete leg;
}

void drawEfficiency(TGraphAsymmErrors* eff, string outputDir = "../output", string unit = "GeV", string label1 = "file1")
{
     eff->SetLineColor(kBlack);
     eff->SetLineWidth(2);
     eff->GetYaxis()->SetRangeUser(0.,1.01);

     // plotting
     setStyle();

     if(string(eff->GetName()).find("HLT")!=std::string::npos && string(eff->GetName()).find("PhotonPt_efficiency")!=std::string::npos)
     {
        eff->GetXaxis()->SetTitle("photon_pt (GeV)");
        eff->GetYaxis()->SetTitle("efficiency");
        eff->SetMarkerSize(1.);
        eff->SetMarkerStyle(20);
     }
     if(string(eff->GetName()).find("HLT")!=std::string::npos && string(eff->GetName()).find("MET_efficiency")!=std::string::npos)
     {
        eff->GetXaxis()->SetTitle("MET (GeV)");
        eff->GetYaxis()->SetTitle("efficiency");
        eff->SetMarkerSize(1.);
        eff->SetMarkerStyle(20);
     }

     TLegend *leg;
     leg = new TLegend(0.68,0.70,0.88,0.90);
     leg->SetFillStyle(0);
     leg->SetBorderSize(0);
     leg->SetTextSize(0.03);
     leg->SetFillColor(0);
     leg->AddEntry(eff, label1.c_str(), "l");
     
     TCanvas* c1 = new TCanvas("c1","c1",1);
     FPCanvasStyle(c1);
     eff->Draw("AP");
     //leg->Draw("same");
     c1->SaveAs((outputDir+"/"+string(eff->GetName())+".png").c_str());
     c1->SaveAs((outputDir+"/"+string(eff->GetName())+".pdf").c_str());

     /*TCanvas* c2 = new TCanvas("c2","c2",1);
     FPCanvasStyle(c2);
     c2->SetLogy();
     eff->Draw("AP");
     //leg->Draw("same");
     c2->SaveAs((outputDir+"/"+string(eff->GetName())+"_log.png").c_str());
     c2->SaveAs((outputDir+"/"+string(eff->GetName())+"_log.pdf").c_str());*/

     delete c1;
     //delete c2;
     delete leg;
}


