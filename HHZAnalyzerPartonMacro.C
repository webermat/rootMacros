#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"
#include "THStack.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGaxis.h"
#include "TLorentzVector.h"
#include <TVirtualFitter.h>
#include "TProfile.h"
#include "TColor.h"
#include <vector>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
//#include "RooNovosibirsk.h"

using namespace std;



float DeltaPhi(float Phi1,float Phi2){
  float deltaphi=fabs(Phi1-Phi2);
  if(deltaphi>M_PI){
    deltaphi=2*M_PI-deltaphi;
  }
  return deltaphi;
}

float DeltaPhiDir(float Phi1,float Phi2){
  float deltaphi=Phi1-Phi2;
  if(deltaphi>M_PI){
    deltaphi=deltaphi-2*M_PI;
  }
  if(deltaphi<(-M_PI)){
    deltaphi=2*M_PI+deltaphi;
  }
  return deltaphi;
}


TCanvas* setUpperCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,10,50,600,500);
    c1->cd();
    gPad->SetTopMargin(0.06);
    gPad->SetRightMargin(0.15);
    return c1;
}


void CLICdpStyle()
{
  gROOT->SetStyle("Plain"); /*Default white background for all plots*/
  /* set bkg color of all to kWhite: white, but not 0*/
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetStatColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFillColor(10);
  gStyle->SetTitleFillColor(kWhite);
  
  
   /* SetPaperSize wants width & height in cm: A4 is 20,26 & US is 20,24*/
   gStyle->SetPaperSize(20, 26); 
   /* No yellow border around histogram*/
   gStyle->SetDrawBorder(0);
   /* remove border of canvas*/
   gStyle->SetCanvasBorderMode(0);
   /* remove border of pads*/
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetLegendBorderSize(0);
  
   /* default text size*/
   gStyle->SetTextSize(0.05);
   gStyle->SetTitleSize(0.06,"xyz");
   gStyle->SetLabelSize(0.06,"xyz");
   /* title offset: distance between given text and axis, here x,y,z*/
   gStyle->SetLabelOffset(0.015,"xyz");
   gStyle->SetTitleOffset(1.2,"yz"); //equivalent to: gStyle->SetTitleYOffset(1.2);
   gStyle->SetTitleOffset(1.17,"x");



   /* Use visible font for all text*/
   int font = 42; 
   gStyle->SetTitleFont(font);
   gStyle->SetTitleFontSize(0.05);
   gStyle->SetStatFont(font);
   gStyle->SetStatFontSize(0.07);
   gStyle->SetTextFont(font);
   gStyle->SetLabelFont(font,"xyz");
   gStyle->SetTitleFont(font,"xyz");
   gStyle->SetTitleBorderSize(0);
   gStyle->SetStatBorderSize(1);
   //ROSA
   //gStyle->SetLegendFont(font);

   /* big marker points*/
   gStyle->SetMarkerStyle(1);
   gStyle->SetLineWidth(2);  
   gStyle->SetMarkerSize(1.2);
   /*set palette in 2d histogram to nice and colorful one*/
   gStyle->SetPalette(1,0); 

   /*No title for histograms*/
   gStyle->SetOptTitle(0);
   /* show the errors on the stat box */
   gStyle->SetOptStat(0); 
   /* show errors on fitted parameters*/
   gStyle->SetOptFit(0); 
   /* number of decimals used for errors*/
   gStyle->SetEndErrorSize(5);   

   /* set line width to 2 by default so that histograms are visible when printed small
      idea: emphasize the data, not the frame around*/
   gStyle->SetHistLineWidth(2);
   gStyle->SetFrameLineWidth(2);
   gStyle->SetFuncWidth(2);
   gStyle->SetHistLineColor(kBlack);
   gStyle->SetFuncColor(kBlack);
   gStyle->SetLabelColor(kBlack,"xyz");

   //set the margins
   gStyle->SetPadBottomMargin(0.18);
   gStyle->SetPadTopMargin(0.11);
   gStyle->SetPadRightMargin(0.08);
   gStyle->SetPadLeftMargin(0.17);
   
   //set the number of divisions to show
   gStyle->SetNdivisions(506, "xy");
   
   //turn off xy grids
   gStyle->SetPadGridX(0);
   gStyle->SetPadGridY(0);
   
   //set the tick mark style
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);

   gStyle->SetCanvasDefW(800);
   gStyle->SetCanvasDefH(700);

   gROOT->ForceStyle();
}

void plot_trips(){

  CLICdpStyle();

  int nbins=12;
  float n_limit_low=0.5;
  float n_limit_high =12.5; 

  TH1F* h_multi_trips_per_months=new TH1F("h_multi_trips_per_months","",nbins,n_limit_low,n_limit_high);
  h_multi_trips_per_months->Sumw2();
  h_multi_trips_per_months->SetFillColor(kBlue);
  h_multi_trips_per_months->SetLineColor(kBlue);
  h_multi_trips_per_months->GetYaxis()->SetTitle("#trips");
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(1,"January");
  h_multi_trips_per_months->SetBinContent(1,7);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(2,"February");
  h_multi_trips_per_months->SetBinContent(2,6);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(3,"March");
  h_multi_trips_per_months->SetBinContent(3,10);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(4,"April");
  h_multi_trips_per_months->SetBinContent(4,11);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(5,"May");
  h_multi_trips_per_months->SetBinContent(5,12);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(6,"June");
  h_multi_trips_per_months->SetBinContent(6,10);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(7,"July");
  h_multi_trips_per_months->SetBinContent(7,14);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(8,"August");
  h_multi_trips_per_months->SetBinContent(8,8);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(9,"September");
  h_multi_trips_per_months->SetBinContent(9,17);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(10,"October");
  h_multi_trips_per_months->SetBinContent(10,7);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(11,"November");
  h_multi_trips_per_months->SetBinContent(11,6);
  h_multi_trips_per_months->GetXaxis()->SetBinLabel(12,"December");
  h_multi_trips_per_months->SetBinContent(12,7);

  
  TH1F* h_single_trips_per_months=new TH1F("h_single_trips_per_months","",nbins,n_limit_low,n_limit_high);
  h_single_trips_per_months->Sumw2();
  h_single_trips_per_months->SetLineColor(kRed);
  h_single_trips_per_months->SetFillColor(kRed);
  h_single_trips_per_months->GetYaxis()->SetTitle("#trips");
  h_single_trips_per_months->GetXaxis()->SetBinLabel(1,"January");
  h_single_trips_per_months->SetBinContent(1,9);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(2,"February");
  h_single_trips_per_months->SetBinContent(2,4);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(3,"March");
  h_single_trips_per_months->SetBinContent(3,3);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(4,"April");
  h_single_trips_per_months->SetBinContent(4,7);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(5,"May");
  h_single_trips_per_months->SetBinContent(5,13);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(6,"June");
  h_single_trips_per_months->SetBinContent(6,14);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(7,"July");
  h_single_trips_per_months->SetBinContent(7,13);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(8,"August");
  h_single_trips_per_months->SetBinContent(8,18);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(9,"September");
  h_single_trips_per_months->SetBinContent(9,9);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(10,"October");
  h_single_trips_per_months->SetBinContent(10,9);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(11,"November");
  h_single_trips_per_months->SetBinContent(11,10);
  h_single_trips_per_months->GetXaxis()->SetBinLabel(12,"December");
  h_single_trips_per_months->SetBinContent(12,8);


  TCanvas* canvas_trip_per_month =setUpperCanvas("canvas_trip_per_month");
  THStack *hs_trips_per_months = new THStack("hs_trips_per_months","");
  //hs_trips_per_months->GetYaxis()->SetTitle("#trips");
  hs_trips_per_months->Add(h_multi_trips_per_months);
  hs_trips_per_months->Add(h_single_trips_per_months);
  hs_trips_per_months->Draw();
  hs_trips_per_months->GetYaxis()->SetTitle("# trips");
  canvas_trip_per_month->Modified();

  TLegend *leg_trip_per_month= new TLegend(0.20,0.735,0.37,0.90);
  leg_trip_per_month->SetBorderSize(0);
  leg_trip_per_month->SetTextAlign(12);
  leg_trip_per_month->SetTextSize(0.050);
  leg_trip_per_month->SetTextFont(42);
  leg_trip_per_month->SetMargin(0.15);
  leg_trip_per_month->SetLineColor(1);
  leg_trip_per_month->SetLineStyle(1);
  leg_trip_per_month->SetLineWidth(1);
  leg_trip_per_month->SetFillColor(0);
  leg_trip_per_month->SetFillStyle(1001);
  leg_trip_per_month->AddEntry(h_multi_trips_per_months,"multi day");
  leg_trip_per_month->AddEntry(h_single_trips_per_months,"single day");
  leg_trip_per_month->Draw();


  TH1F* h_multi_trips_duration_per_months=new TH1F("h_multi_trips_duration_per_months","",nbins,n_limit_low,n_limit_high);
  h_multi_trips_duration_per_months->Sumw2();
  h_multi_trips_duration_per_months->SetFillColor(kBlue);
  h_multi_trips_duration_per_months->SetLineColor(kBlue);
  h_multi_trips_duration_per_months->GetYaxis()->SetTitle("#days on trips");
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(1,"January");
  h_multi_trips_duration_per_months->SetBinContent(1,31);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(2,"February");
  h_multi_trips_duration_per_months->SetBinContent(2,13);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(3,"March");
  h_multi_trips_duration_per_months->SetBinContent(3,46);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(4,"April");
  h_multi_trips_duration_per_months->SetBinContent(4,42);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(5,"May");
  h_multi_trips_duration_per_months->SetBinContent(5,31);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(6,"June");
  h_multi_trips_duration_per_months->SetBinContent(6,70);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(7,"July");
  h_multi_trips_duration_per_months->SetBinContent(7,68);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(8,"August");
  h_multi_trips_duration_per_months->SetBinContent(8,40);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(9,"September");
  h_multi_trips_duration_per_months->SetBinContent(9,60);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(10,"October");
  h_multi_trips_duration_per_months->SetBinContent(10,21);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(11,"November");
  h_multi_trips_duration_per_months->SetBinContent(11,22);
  h_multi_trips_duration_per_months->GetXaxis()->SetBinLabel(12,"December");
  h_multi_trips_duration_per_months->SetBinContent(12,28);
  
  TH1F* h_single_trips_duration_per_months=new TH1F("h_single_trips_duration_per_months","",nbins,n_limit_low,n_limit_high);
  h_single_trips_duration_per_months->Sumw2();
  h_single_trips_duration_per_months->SetFillColor(kRed);
  h_single_trips_duration_per_months->SetLineColor(kRed);
  h_single_trips_duration_per_months->GetYaxis()->SetTitle("#days on trips");
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(1,"January");
  h_single_trips_duration_per_months->SetBinContent(1,9);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(2,"February");
  h_single_trips_duration_per_months->SetBinContent(2,4);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(3,"March");
  h_single_trips_duration_per_months->SetBinContent(3,3);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(4,"April");
  h_single_trips_duration_per_months->SetBinContent(4,7);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(5,"May");
  h_single_trips_duration_per_months->SetBinContent(5,13);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(6,"June");
  h_single_trips_duration_per_months->SetBinContent(6,14);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(7,"July");
  h_single_trips_duration_per_months->SetBinContent(7,13);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(8,"August");
  h_single_trips_duration_per_months->SetBinContent(8,18);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(9,"September");
  h_single_trips_duration_per_months->SetBinContent(9,9);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(10,"October");
  h_single_trips_duration_per_months->SetBinContent(10,9);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(11,"November");
  h_single_trips_duration_per_months->SetBinContent(11,10);
  h_single_trips_duration_per_months->GetXaxis()->SetBinLabel(12,"December");
  h_single_trips_duration_per_months->SetBinContent(12,8);

  TCanvas* canvas_trips_duration_per_month =setUpperCanvas("canvas_trips_duration_per_month");
  THStack *hs_trips_duration_per_months = new THStack("hs_trips_duration_per_months","");
  //hs_trips_duration_per_months->GetYaxis()->SetTitle("#trips_duration");
  hs_trips_duration_per_months->Add(h_multi_trips_duration_per_months);
  hs_trips_duration_per_months->Add(h_single_trips_duration_per_months);
  hs_trips_duration_per_months->Draw();
  hs_trips_duration_per_months->GetYaxis()->SetTitle("# total days spend on trips");
  canvas_trips_duration_per_month->Modified();

  TLegend *leg_trips_duration_per_month= new TLegend(0.20,0.735,0.37,0.90);
  leg_trips_duration_per_month->SetBorderSize(0);
  leg_trips_duration_per_month->SetTextAlign(12);
  leg_trips_duration_per_month->SetTextSize(0.050);
  leg_trips_duration_per_month->SetTextFont(42);
  leg_trips_duration_per_month->SetMargin(0.15);
  leg_trips_duration_per_month->SetLineColor(1);
  leg_trips_duration_per_month->SetLineStyle(1);
  leg_trips_duration_per_month->SetLineWidth(1);
  leg_trips_duration_per_month->SetFillColor(0);
  leg_trips_duration_per_month->SetFillStyle(1001);
  leg_trips_duration_per_month->AddEntry(h_multi_trips_duration_per_months,"multi day");
  leg_trips_duration_per_month->AddEntry(h_single_trips_duration_per_months,"single day");
  leg_trips_duration_per_month->Draw();


}


TCanvas* setRatioCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,0,50,800,300);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}


void fill_HZ_histograms(TFile* file, std::vector<TH1F*> h_hist_vec, std::vector<TH2F*> h_hist_vec_2D){

  TTree* tree = (TTree*)file->Get("showerData");

  //totPFO includes ALL PFO's -->also the isolated particles
  float totPFO_Px=0;
  float totPFO_Py=0;
  float totPFO_Pz=0;
  float totPFO_E=0;
  //includes ALL stable visible Gen particles --> also isolated particles
  float true_Px=0;
  float true_Py=0;
  float true_Pz=0;
  float true_E=0;

 //includes ALL stable visible Gen particles
  float true_inv_Px=0;
  float true_inv_Py=0;
  float true_inv_Pz=0;
  float true_inv_E=0;


  vector<float> *genTrueLepPh_Px=0;
  vector<float> *genTrueLepPh_Py=0;
  vector<float> *genTrueLepPh_Pz=0;
  vector<float> *genTrueLepPh_E=0;
  vector<float> *genTrueLepPh_relIso=0;
  vector<int> *genTrueLepPh_PDGID=0;
  vector<int> *genTrueLepPh_ANC_PDGID=0;

  vector<float> *isoPartGenDR10_Px=0;
  vector<float> *isoPartGenDR10_Py=0;
  vector<float> *isoPartGenDR10_Pz=0;
  vector<float> *isoPartGenDR10_E=0;
  vector<float> *isoPartGenDR10_relIso=0;
  vector<int> *isoPartGenDR10_PDGID=0;

  vector<float> *isoPartRecoDR10_Px=0;
  vector<float> *isoPartRecoDR10_Py=0;
  vector<float> *isoPartRecoDR10_Pz=0;
  vector<float> *isoPartRecoDR10_E=0;
  vector<float> *isoPartRecoDR10_relIso=0;
  vector<int> *isoPartRecoDR10_PDGID=0;
 
  vector<float> *recojet_Px=0;
  vector<float> *recojet_Py=0;
  vector<float> *recojet_Pz=0;
  vector<float> *recojet_E=0;

  vector<float> *genjet_Px=0;
  vector<float> *genjet_Py=0;
  vector<float> *genjet_Pz=0;
  vector<float> *genjet_E=0;

  vector<float> *recojet_subjet_Px=0;
  vector<float> *recojet_subjet_Py=0;
  vector<float> *recojet_subjet_Pz=0;
  vector<float> *recojet_subjet_E=0;
  vector<int> *recojet_subjet_jetindex=0;


  vector<float> *genjet_subjet_Px=0;
  vector<float> *genjet_subjet_Py=0;
  vector<float> *genjet_subjet_Pz=0;
  vector<float> *genjet_subjet_E=0;
  vector<int> *genjet_subjet_jetindex=0;



  vector<float> *trueME_Px=0;
  vector<float> *trueME_Py=0;
  vector<float> *trueME_Pz=0;
  vector<float> *trueME_E=0;
  vector<int> *trueME_PDGID=0;

  tree->SetBranchAddress("trueME_E", &trueME_E);
  tree->SetBranchAddress("trueME_Px", &trueME_Px);
  tree->SetBranchAddress("trueME_Py", &trueME_Py);
  tree->SetBranchAddress("trueME_Pz", &trueME_Pz);
  tree->SetBranchAddress("trueME_PDGID", &trueME_PDGID);

  tree->SetBranchAddress("totPFO_E", &totPFO_E);
  tree->SetBranchAddress("totPFO_Px", &totPFO_Px);
  tree->SetBranchAddress("totPFO_Py", &totPFO_Py);
  tree->SetBranchAddress("totPFO_Pz", &totPFO_Pz);

  tree->SetBranchAddress("true_E", &true_E);
  tree->SetBranchAddress("true_Px", &true_Px);
  tree->SetBranchAddress("true_Py", &true_Py);
  tree->SetBranchAddress("true_Pz", &true_Pz);

  tree->SetBranchAddress("true_inv_E", &true_inv_E);
  tree->SetBranchAddress("true_inv_Px", &true_inv_Px);
  tree->SetBranchAddress("true_inv_Py", &true_inv_Py);
  tree->SetBranchAddress("true_inv_Pz", &true_inv_Pz);

  tree->SetBranchAddress("genTrueLepPh_E", &genTrueLepPh_E);
  tree->SetBranchAddress("genTrueLepPh_Px", &genTrueLepPh_Px);
  tree->SetBranchAddress("genTrueLepPh_Py", &genTrueLepPh_Py);
  tree->SetBranchAddress("genTrueLepPh_Pz", &genTrueLepPh_Pz);
  tree->SetBranchAddress("genTrueLepPh_PDGID", &genTrueLepPh_PDGID);
  tree->SetBranchAddress("genTrueLepPh_ANC_PDGID", &genTrueLepPh_ANC_PDGID);
 
  tree->SetBranchAddress("isoPartGenDR10_E", &isoPartGenDR10_E);
  tree->SetBranchAddress("isoPartGenDR10_Px", &isoPartGenDR10_Px);
  tree->SetBranchAddress("isoPartGenDR10_Py", &isoPartGenDR10_Py);
  tree->SetBranchAddress("isoPartGenDR10_Pz", &isoPartGenDR10_Pz);
  tree->SetBranchAddress("isoPartGenDR10_PDGID", &isoPartGenDR10_PDGID);
  tree->SetBranchAddress("isoPartGenDR10_relIso", &isoPartGenDR10_relIso);
  
  tree->SetBranchAddress("isoPartRecoDR10_E", &isoPartRecoDR10_E);
  tree->SetBranchAddress("isoPartRecoDR10_Px", &isoPartRecoDR10_Px);
  tree->SetBranchAddress("isoPartRecoDR10_Py", &isoPartRecoDR10_Py);
  tree->SetBranchAddress("isoPartRecoDR10_Pz", &isoPartRecoDR10_Pz);
  tree->SetBranchAddress("isoPartRecoDR10_PDGID", &isoPartRecoDR10_PDGID);
  tree->SetBranchAddress("isoPartRecoDR10_relIso", &isoPartRecoDR10_relIso);

  tree->SetBranchAddress("genjet_E", &genjet_E);
  tree->SetBranchAddress("genjet_Px", &genjet_Px);
  tree->SetBranchAddress("genjet_Py", &genjet_Py);
  tree->SetBranchAddress("genjet_Pz", &genjet_Pz);

  tree->SetBranchAddress("recojet_E", &recojet_E);
  tree->SetBranchAddress("recojet_Px", &recojet_Px);
  tree->SetBranchAddress("recojet_Py", &recojet_Py);
  tree->SetBranchAddress("recojet_Pz", &recojet_Pz);

  tree->SetBranchAddress("genjet_subjet_E", &genjet_subjet_E);
  tree->SetBranchAddress("genjet_subjet_Px", &genjet_subjet_Px);
  tree->SetBranchAddress("genjet_subjet_Py", &genjet_subjet_Py);
  tree->SetBranchAddress("genjet_subjet_Pz", &genjet_subjet_Pz);
  tree->SetBranchAddress("genjet_subjet_jetindex", &genjet_subjet_jetindex);

  tree->SetBranchAddress("recojet_subjet_E", &recojet_subjet_E);
  tree->SetBranchAddress("recojet_subjet_Px", &recojet_subjet_Px);
  tree->SetBranchAddress("recojet_subjet_Py", &recojet_subjet_Py);
  tree->SetBranchAddress("recojet_subjet_Pz", &recojet_subjet_Pz);
  tree->SetBranchAddress("recojet_subjet_jetindex", &recojet_subjet_jetindex);



  //we should have the following
  //index 0 and 1 are outgoing e+e-
  //index 2 and 3 ISR photons from incoming e+e-
  //index 4 and 5 Z quarks
  //index 6 is Higgs (boosted back to real system)
  //index 7 and 8 are H daughters

  //beam strahlung effects not in history -->take outgoing e+e- for centre of mass system, i.e. use 0 and 1, NOT 2 and 3

  double sqrtS_low=750;
  double sqrtS_high=2750;
  double sqrtS_high_reco=2500;
  double sqrtS_nom=3000;

  float theta_max_Z_q=-100;
  double alpha_min =200;

  for(unsigned int i_entry=0;i_entry<tree->GetEntries();i_entry++){
    //fill jet energy resolution histograms
    tree->GetEntry(i_entry);

    if(i_entry%5000==0){
      std::cout<<"entry "<<i_entry<<std::endl;
    }


    TLorentzVector tempTotEventP4(0,0,0,0);
    TLorentzVector tempTotEventP4HZ(0,0,0,0);
    TLorentzVector tempZP4(0,0,0,0);
    TLorentzVector tempHP4(0,0,0,0);


    bool H_decays_bbar=true;

    if(fabs((*trueME_PDGID)[7])!=5 || fabs((*trueME_PDGID)[8])!=5){
      H_decays_bbar=false;
      std::cout<<"entry "<<i_entry <<" is NO bbar "<<(*trueME_PDGID)[7]<<"/"<<(*trueME_PDGID)[8]<<std::endl;
    }

    if(!H_decays_bbar){
      continue;
    }
    std::cout<<"entry "<<i_entry <<" is bbar "<<(*trueME_PDGID)[7]<<"/"<<(*trueME_PDGID)[8]<<std::endl;

    tempHP4.SetPxPyPzE((*trueME_Px)[6],(*trueME_Py)[6],(*trueME_Pz)[6],(*trueME_E)[6]);
    for(unsigned int i=0;i<trueME_E->size();i++){
      TLorentzVector temp(0,0,0,0);
      temp.SetPxPyPzE((*trueME_Px)[i],(*trueME_Py)[i],(*trueME_Pz)[i],(*trueME_E)[i]);
      if(i<2){
	tempTotEventP4+=temp;
	//std::cout<<i<<" px/py/pz/E "<<tempTotEventP4.Px()<<"/"<<tempTotEventP4.Py()<<"/"<<tempTotEventP4.Pz()<<"/"<<tempTotEventP4.E()<<std::endl;
      }
      //at this point sqrtS should be known
      if(i==1){
	h_hist_vec[0]->Fill(tempTotEventP4.M());
      }
      if(i>3 && i<6){//index 4 and 5 are Z daughters
	tempZP4+=temp;
	theta_max_Z_q=temp.Theta();
	if(i==4){
	  //if(fabs((*trueME_PDGID)[i])==5){
	  TLorentzVector temp2(0,0,0,0);
	  temp2.SetPxPyPzE((*trueME_Px)[i+1],(*trueME_Py)[i+1],(*trueME_Pz)[i+1],(*trueME_E)[i+1]);
	  if(fabs((*trueME_PDGID)[7])==5){
	    TLorentzVector temp2HD(0,0,0,0);
	    temp2HD.SetPxPyPzE((*trueME_Px)[7],(*trueME_Py)[7],(*trueME_Pz)[7],(*trueME_E)[7]);
	    alpha_min=temp.Angle(temp2HD.Vect());
	    if(temp2.Angle(temp2HD.Vect())<alpha_min){
	      alpha_min=temp2.Angle(temp2HD.Vect());  
	    }
	    //both Z quarks checked against first H bottom
	    temp2HD.SetPxPyPzE((*trueME_Px)[8],(*trueME_Py)[8],(*trueME_Pz)[8],(*trueME_E)[8]);
	    if(temp.Angle(temp2HD.Vect())<alpha_min){
	      alpha_min=temp.Angle(temp2HD.Vect());  
	    }
	    if(temp2.Angle(temp2HD.Vect())<alpha_min){
	      alpha_min=temp2.Angle(temp2HD.Vect());  
	    }
	  if(tempTotEventP4.M()<sqrtS_low){
	    h_hist_vec[19]->Fill(alpha_min*TMath::RadToDeg());
	  }else if(tempTotEventP4.M()<sqrtS_high){
	    h_hist_vec[20]->Fill(alpha_min*TMath::RadToDeg());
	  }else{
	    h_hist_vec[21]->Fill(alpha_min*TMath::RadToDeg());
	  }
	  }
	  if(fabs(temp2.Theta())>fabs(temp.Theta())){
	    theta_max_Z_q=temp2.Theta();
	  }
	  if(tempTotEventP4.M()<sqrtS_low){
	    h_hist_vec[10]->Fill(temp.Angle(temp2.Vect())*TMath::RadToDeg());
	    h_hist_vec[16]->Fill(theta_max_Z_q*TMath::RadToDeg());
	  }else if(tempTotEventP4.M()<sqrtS_high){
	    h_hist_vec[11]->Fill(temp.Angle(temp2.Vect())*TMath::RadToDeg());
	    h_hist_vec[17]->Fill(theta_max_Z_q*TMath::RadToDeg());
	  }else{
	    h_hist_vec[12]->Fill(temp.Angle(temp2.Vect())*TMath::RadToDeg());
	    h_hist_vec[18]->Fill(theta_max_Z_q*TMath::RadToDeg());
	  }
	  //}
	}
      }
      if(i==6){//fill higgs histograms
	if(tempTotEventP4.M()<sqrtS_low){
	  h_hist_vec[1]->Fill(temp.Pt());
	  h_hist_vec[22]->Fill(temp.Theta()*TMath::RadToDeg());
	}else if(tempTotEventP4.M()<sqrtS_high){
	  h_hist_vec[2]->Fill(temp.Pt());
	  h_hist_vec[23]->Fill(temp.Theta()*TMath::RadToDeg());
	}else{
	  h_hist_vec[3]->Fill(temp.Pt());
	  h_hist_vec[24]->Fill(temp.Theta()*TMath::RadToDeg());
	}
      }
      if(i==7){//first H daughter, consider for now only b or bbar
	if(fabs((*trueME_PDGID)[i])==5){
	  TLorentzVector temp2(0,0,0,0);
	  temp2.SetPxPyPzE((*trueME_Px)[i+1],(*trueME_Py)[i+1],(*trueME_Pz)[i+1],(*trueME_E)[i+1]);
	  float theta_max_H_q=temp2.Theta();
	  if(fabs(temp.Theta())>fabs(temp2.Theta())){
	    theta_max_H_q=temp.Theta();
	  }
	  if(tempTotEventP4.M()<sqrtS_low){
	    h_hist_vec[7]->Fill(temp.Angle(temp2.Vect())*TMath::RadToDeg());
	    h_hist_vec[13]->Fill(theta_max_H_q*TMath::RadToDeg());
	  }else if(tempTotEventP4.M()<sqrtS_high){
	    h_hist_vec[8]->Fill(temp.Angle(temp2.Vect())*TMath::RadToDeg());
	    h_hist_vec[14]->Fill(theta_max_H_q*TMath::RadToDeg());
	  }else{
	    h_hist_vec[9]->Fill(temp.Angle(temp2.Vect())*TMath::RadToDeg());
	    h_hist_vec[15]->Fill(theta_max_H_q*TMath::RadToDeg());
	  }
	}
      }
      //if(trueME_E->size()==9){
      //std::cout<<"PDGID's "<<i<<"/"<<(*trueME_PDGID)[i]<<"/"<<tempTotEventP4.M()<<std::endl;
      //}
    }
    tempTotEventP4HZ=tempZP4+tempHP4;

    if(tempTotEventP4.M()>sqrtS_high_reco){
      //4 and 5 are Z daughters
      if((*trueME_E)[4]>(*trueME_E)[5]){
	h_hist_vec[96]->Fill((*trueME_E)[4]/((*trueME_E)[4]+(*trueME_E)[5]));
      h_hist_vec[97]->Fill((*trueME_E)[5]/((*trueME_E)[4]+(*trueME_E)[5]));
      }else{
	h_hist_vec[96]->Fill((*trueME_E)[5]/((*trueME_E)[4]+(*trueME_E)[5]));
	h_hist_vec[97]->Fill((*trueME_E)[4]/((*trueME_E)[4]+(*trueME_E)[5]));
      }
      //7 and 8 are H daughters
      if((*trueME_E)[7]>(*trueME_E)[8]){
	h_hist_vec[98]->Fill((*trueME_E)[7]/((*trueME_E)[7]+(*trueME_E)[8]));
	h_hist_vec[99]->Fill((*trueME_E)[8]/((*trueME_E)[7]+(*trueME_E)[8]));
      }else{
	h_hist_vec[98]->Fill((*trueME_E)[8]/((*trueME_E)[7]+(*trueME_E)[8]));
	h_hist_vec[99]->Fill((*trueME_E)[7]/((*trueME_E)[7]+(*trueME_E)[8]));
      }
    }


    TLorentzVector tempTotGenP4(0,0,0,0);
    tempTotGenP4.SetPxPyPzE(true_Px,true_Py,true_Pz,true_E);

    TLorentzVector tempTotInvGenP4(0,0,0,0);
    tempTotInvGenP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);

    TLorentzVector tempRecoMETP4(0,0,0,0);
    tempRecoMETP4.SetPxPyPzE(-totPFO_Px,-totPFO_Py,-totPFO_Pz,sqrtS_nom-totPFO_E);

    TLorentzVector tempTotRecoP4(0,0,0,0);
    tempTotRecoP4.SetPxPyPzE(totPFO_Px,totPFO_Py,totPFO_Pz,totPFO_E);

    h_hist_vec[82]->Fill(tempTotInvGenP4.Pt());
    h_hist_vec[83]->Fill(tempRecoMETP4.Pt());

    TLorentzVector tempGenTruePhP4(0,0,0,0);
    for(unsigned int i=0;i<genTrueLepPh_E->size();i++){
      if((*genTrueLepPh_PDGID)[i]==22){
	//photons from beam
	if((*genTrueLepPh_ANC_PDGID)[i]==0){
	  TLorentzVector temp(0,0,0,0);
	  temp.SetPxPyPzE((*genTrueLepPh_Px)[i],(*genTrueLepPh_Py)[i],(*genTrueLepPh_Pz)[i],(*genTrueLepPh_E)[i]);
	  tempGenTruePhP4+=temp;
	}  
      }
    }
    
    TLorentzVector tempGenIsoPhP4(0,0,0,0);
    TLorentzVector tempGenIsoLepP4(0,0,0,0);
    for(unsigned int i=0;i<isoPartGenDR10_E->size();i++){
      if((*isoPartGenDR10_PDGID)[i]==22){
	//photons are isolated
	if((*isoPartGenDR10_relIso)[i]<0.10){
	  TLorentzVector temp(0,0,0,0);
	  temp.SetPxPyPzE((*isoPartGenDR10_Px)[i],(*isoPartGenDR10_Py)[i],(*isoPartGenDR10_Pz)[i],(*isoPartGenDR10_E)[i]);
	  tempGenIsoPhP4+=temp;
	}  
      }else if(abs((*isoPartGenDR10_PDGID)[i])==11 || abs((*isoPartGenDR10_PDGID)[i])==13){
	if((*isoPartGenDR10_relIso)[i]<0.10){
	  TLorentzVector temp(0,0,0,0);
	  temp.SetPxPyPzE((*isoPartGenDR10_Px)[i],(*isoPartGenDR10_Py)[i],(*isoPartGenDR10_Pz)[i],(*isoPartGenDR10_E)[i]);
	  tempGenIsoLepP4+=temp;
	}
      }
    }
    
    TLorentzVector tempRecoIsoPhP4(0,0,0,0);
    TLorentzVector tempRecoIsoLepP4(0,0,0,0);
    for(unsigned int i=0;i<isoPartRecoDR10_E->size();i++){
      if((*isoPartRecoDR10_PDGID)[i]==22){
	//photons are isolated
	if((*isoPartRecoDR10_relIso)[i]<0.10){
	  TLorentzVector temp(0,0,0,0);
	  temp.SetPxPyPzE(-(*isoPartRecoDR10_Px)[i],-(*isoPartRecoDR10_Py)[i],-(*isoPartRecoDR10_Pz)[i],-(*isoPartRecoDR10_E)[i]);
	  tempRecoIsoPhP4+=temp;
	}else if(abs((*isoPartRecoDR10_PDGID)[i])==11 || abs((*isoPartRecoDR10_PDGID)[i])==13){
	  if((*isoPartRecoDR10_relIso)[i]<0.10){
	    TLorentzVector temp(0,0,0,0);
	    temp.SetPxPyPzE((*isoPartRecoDR10_Px)[i],(*isoPartRecoDR10_Py)[i],(*isoPartRecoDR10_Pz)[i],(*isoPartRecoDR10_E)[i]);
	    tempRecoIsoLepP4+=temp;
	  }
	}  
      }
    }

    //h_hist_vec[37]->Fill((tempTotGenP4+tempGenIsoPhP4).M());
    //h_hist_vec[38]->Fill((tempTotRecoP4+tempRecoIsoPhP4).M());

    h_hist_vec[37]->Fill((tempTotGenP4-tempGenIsoPhP4).M());
    h_hist_vec[38]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M());
    h_hist_vec[39]->Fill(tempTotEventP4HZ.M());
    h_hist_vec[40]->Fill((tempTotGenP4-tempGenIsoPhP4+tempTotInvGenP4).M());
    h_hist_vec[41]->Fill((tempTotRecoP4).M());
    h_hist_vec[42]->Fill((tempTotGenP4).M());
    h_hist_vec[43]->Fill((tempTotGenP4+tempTotInvGenP4).M());
    h_hist_vec[44]->Fill((tempTotGenP4-tempGenTruePhP4+tempTotInvGenP4).M());
    h_hist_vec[45]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempTotInvGenP4).M());
  

  }
}

  


void HZAnalyzerFull(){

  CLICdpStyle();

  gROOT->ProcessLine("#include <vector>");

  //const char* final_histo_name="/afs/cern.ch/user/w/weberma2/performanceHistoFiles/HZqq_testanalyzer_10906_10905_antibbarOnly_DR10_DR12.root";
  const char* final_histo_name="/afs/cern.ch/user/w/weberma2/performanceHistoFiles/HZqq_testanalyzer_10906_10905_SignalOnly_DR10_DR15.root";
  
  TFile* file_CLIC_HZqq = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/180810_gcc62_CT/HZStudy_HZqq_10906_3TeVOverlay_CLIC_o3_v14_PartonInfo_DR10_DR15.root");
  TFile* file_CLIC_HZqq_noBGG = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/180810_gcc62_CT/HZStudy_HZqq_10905_noOverlay_CLIC_o3_v14_PartonInfo_DR10_DR15.root");


  TFile* file_histogram=new TFile(final_histo_name,"recreate");

  int n_bins_high=200;

  double lim_energy_low=0;
  double lim_energy_high=3050.;

  TH1F* h_sqrtS_e1_e2_effective = new TH1F("h_sqrtS_e1_e2_effective","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_H1_H2_E_sum = new TH1F("h_H_H2_E_sum","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
  TH1F* h_H1_E = new TH1F("h_H1_E","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
  TH1F* h_H2_E = new TH1F("h_H1_E","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
  TH1F* h_Z_E = new TH1F("h_Z_E_sum","", n_bins_high, lim_energy_low,0.5*lim_energy_high);

  double lim_mass_low=0;
  double lim_mass_high=500;

  TH1F* h_H1_H2_mass = new TH1F("h_H_H2_E_sum","", n_bins_high, lim_energy_low,0.5*lim_energy_high);

  double lim_dalpha_low=0.;
  double lim_dalpha_high=180.;

  TH1F* h_dalpha_H1_H2_comb_vs_Z = new TH1F("h_dalpha_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dphi_H1_H2_comb_vs_Z = new TH1F("h_dphi_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dtheta_H1_H2_comb_vs_Z = new TH1F("h_dtheta_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_dalpha_H1_H2 = new TH1F("h_dalpha_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dphi_H1_H2 = new TH1F("h_dphi_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dtheta_H1_H2 = new TH1F("h_dtheta_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  double lim_dalpha_qqbar_low=0.;
  double lim_dalpha_qqbar_high=30.;

  TH1F* h_dalpha_H1_bbar = new TH1F("h_dalpha_H1_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dphi_H1_bbar = new TH1F("h_dphi_H1_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dtheta_H1_bbar = new TH1F("h_dtheta_H1_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);

  TH1F* h_dalpha_H2_bbar = new TH1F("h_dalpha_H2_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dphi_H2_bbar = new TH1F("h_dphi_H2_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dtheta_H2_bbar = new TH1F("h_dtheta_H2_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);

  TH1F* h_dalpha_H1_qqbar = new TH1F("h_dalpha_H1_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dphi_H1_qqbar = new TH1F("h_dphi_H1_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dtheta_H1_qqbar = new TH1F("h_dtheta_H1_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);

  TH1F* h_dalpha_H2_qqbar = new TH1F("h_dalpha_H2_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dphi_H2_qqbar = new TH1F("h_dphi_H2_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dtheta_H2_qqbar = new TH1F("h_dtheta_H2_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);

  TH1F* h_dalpha_Z_qqbar = new TH1F("h_dalpha_Z_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dphi_Z_qqbar = new TH1F("h_dphi_Z_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dtheta_Z_qqbar = new TH1F("h_dtheta_Z_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);


  double lim_dalpha_qqbar_allH_low=0.;
  double lim_dalpha_qqbar_allH_high=100.;

  TH1F* h_dalpha_max_allH_qqbar = new TH1F("h_dalpha_max_allH_qqbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
  TH1F* h_dphi_max_allH_qqbar = new TH1F("h_dphi_max_allH_qqbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
  TH1F* h_dtheta_max_allH_qqbar = new TH1F("h_dtheta_max_allH_qqbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);

  TH1F* h_dalpha_max_allH_bbar = new TH1F("h_dalpha_max_allH_bbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
  TH1F* h_dphi_max_allH_bbar = new TH1F("h_dphi_max_allH_bbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
  TH1F* h_dtheta_max_allH_bbar = new TH1F("h_dtheta_max_allH_bbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);

  TH1F* h_dalpha_allH_diboson_qqbar = new TH1F("h_dalpha_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dphi_allH_diboson_qqbar = new TH1F("h_dphi_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
  TH1F* h_dtheta_allH_diboson_qqbar = new TH1F("h_dtheta_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);

  double lim_H1_E_over_allH_E_low=0.5;
  double lim_H1_E_over_allH_E_high=1.0;

  TH1F* h_E_H1_over_E_allH = new TH1F("h_E_H1_over_E_allH","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
  TH1F* h_E_H1_over_E_allH_bbar = new TH1F("h_E_H1_over_E_allH_bbar","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
  TH1F* h_E_H1_over_E_allH_qqbar = new TH1F("h_E_H1_over_E_allH_qqbar","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);

  std::vector<TH1F*> hist_vec_HZ_parton;

  hist_vec_HZ_parton.push_back(h_sqrtS_e1_e2_effective);
  hist_vec_HZ_parton.push_back(h_H1_H2_E_sum);
  hist_vec_HZ_parton.push_back(h_H1_E);
  hist_vec_HZ_parton.push_back(h_H2_E);
  hist_vec_HZ_parton.push_back(h_Z_E);
  hist_vec_HZ_parton.push_back(h_H1_H2_mass);
  hist_vec_HZ_parton.push_back(h_dalpha_H1_H2_comb_vs_Z);
  hist_vec_HZ_parton.push_back(h_dphi_H1_H2_comb_vs_Z);
  hist_vec_HZ_parton.push_back(h_dtheta_H1_H2_comb_vs_Z);
  hist_vec_HZ_parton.push_back(h_dalpha_H1_H2);
  hist_vec_HZ_parton.push_back(h_dphi_H1_H2);
  hist_vec_HZ_parton.push_back(h_dtheta_H1_H2);
  hist_vec_HZ_parton.push_back(h_dalpha_H1_bbar);
  hist_vec_HZ_parton.push_back(h_dphi_H1_bbar);
  hist_vec_HZ_parton.push_back(h_dtheta_H1_bbar);
  hist_vec_HZ_parton.push_back(h_dalpha_H2_bbar);
  hist_vec_HZ_parton.push_back(h_dphi_H2_bbar);
  hist_vec_HZ_parton.push_back(h_dtheta_H2_bbar);
  hist_vec_HZ_parton.push_back(h_dalpha_H1_qqbar);
  hist_vec_HZ_parton.push_back(h_dphi_H1_qqbar);
  hist_vec_HZ_parton.push_back(h_dtheta_H1_qqbar);
  hist_vec_HZ_parton.push_back(h_dalpha_H2_qqbar);
  hist_vec_HZ_parton.push_back(h_dphi_H2_qqbar);
  hist_vec_HZ_parton.push_back(h_dtheta_H2_qqbar);
  hist_vec_HZ_parton.push_back(h_dalpha_Z_qqbar);
  hist_vec_HZ_parton.push_back(h_dphi_Z_qqbar);
  hist_vec_HZ_parton.push_back(h_dtheta_Z_qqbar);
  hist_vec_HZ_parton.push_back(h_dalpha_max_allH_qqbar);
  hist_vec_HZ_parton.push_back(h_dphi_max_allH_qqbar);
  hist_vec_HZ_parton.push_back(h_dtheta_max_allH_qqbar);
  hist_vec_HZ_parton.push_back(h_dalpha_max_allH_bbar);
  hist_vec_HZ_parton.push_back(h_dphi_max_allH_bbar);
  hist_vec_HZ_parton.push_back(h_dtheta_max_allH_bbar);
  hist_vec_HZ_parton.push_back(h_dalpha_allH_diboson_qqbar);
  hist_vec_HZ_parton.push_back(h_dphi_allH_diboson_qqbar);
  hist_vec_HZ_parton.push_back(h_dtheta_allH_diboson_qqbar);
  hist_vec_HZ_parton.push_back(h_E_H1_over_E_allH);
  hist_vec_HZ_parton.push_back(h_E_H1_over_E_allH_bbar);
  hist_vec_HZ_parton.push_back(h_E_H1_over_E_allH_qqbar);

  for(unsigned int i=0;i<hist_vec_HZ_parton.size();i++){
    hist_vec_HZ_parton[i]->Sumw2();
  }
  std::vector<TH2F*> hist_vec_HZ_2DHist;

 fill_HZ_histograms(file_CLIC_HZqq,hist_vec_HZ_parton,hist_vec_HZ_2DHist);



 file_histogram->Write();
 file_histogram->Close();

}
