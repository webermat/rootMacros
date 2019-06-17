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
   gStyle->SetTitleFontSize(0.06);
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


TCanvas* setRatioCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,0,50,800,300);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}


void fill_HZ_histograms(TFile* file, std::vector<TH1F*> h_hist_parton, std::vector<TH1F*> h_hist_vec_gen, std::vector<TH1F*> h_hist_vec_reco, std::vector<TH2F*> h_hist_vec_2D, bool usePartonInfo ,double x_sec, bool fill_partonInfo, bool fill_genInfo){

  std::cout<<"size of histogram vectors "<<h_hist_parton.size()<<"/"<<h_hist_vec_gen.size()<<"/"<<h_hist_vec_reco.size()<<"/"<<h_hist_vec_2D.size()<<std::endl;

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
  vector<float> *recojet_CHFraction=0;
  vector<float> *recojet_PhFraction=0;
  vector<int> *recojet_NCH=0;
  vector<int> *recojet_NPh=0;
  vector<int> *recojet_Mult=0;

  vector<float> *genjet_Px=0;
  vector<float> *genjet_Py=0;
  vector<float> *genjet_Pz=0;
  vector<float> *genjet_E=0;
  vector<float> *genjet_CHFraction=0;
  vector<float> *genjet_PhFraction=0;
  vector<int> *genjet_NCH=0;
  vector<int> *genjet_NPh=0;
  vector<int> *genjet_Mult=0;

  vector<float> *gen_tau_Px=0;
  vector<float> *gen_tau_Py=0;
  vector<float> *gen_tau_Pz=0;
  vector<float> *gen_tau_E=0;
  vector<int> *gen_tau_Charge=0;
  vector<int> *gen_tau_NCH=0;
  vector<int> *gen_tau_NPh=0;
  vector<float> *gen_tau_CHFraction=0;
  vector<float> *gen_tau_PhFraction=0;

  vector<float> *reco_tau_Px=0;
  vector<float> *reco_tau_Py=0;
  vector<float> *reco_tau_Pz=0;
  vector<float> *reco_tau_E=0;
  vector<int> *reco_tau_Charge=0;
  vector<int> *reco_tau_NCH=0;
  vector<int> *reco_tau_NPh=0;
  vector<float> *reco_tau_CHFraction=0;
  vector<float> *reco_tau_PhFraction=0;

  vector<float> *recojet_subjet_Px=0;
  vector<float> *recojet_subjet_Py=0;
  vector<float> *recojet_subjet_Pz=0;
  vector<float> *recojet_subjet_E=0;
  vector<int> *recojet_subjet_jetindex=0;

  vector<float> *recojet_nsubjettiness1=0;
  vector<float> *recojet_nsubjettiness2=0;
  vector<float> *recojet_nsubjettiness3=0;
  vector<float> *recojet_nsubjettiness1_lrz=0;
  vector<float> *recojet_nsubjettiness2_lrz=0;
  vector<float> *recojet_nsubjettiness3_lrz=0;

  vector<float> *recojet_beta1_ECorr2=0;
  vector<float> *recojet_beta1_ECorr3=0;
  vector<float> *recojet_beta1_ECorr2_E_theta=0;
  vector<float> *recojet_beta1_ECorr3_E_theta=0;

  vector<float> *recojet_beta1_N2=0;
  vector<float> *recojet_beta1_N3=0;
  vector<float> *recojet_beta1_C2=0;
  vector<float> *recojet_beta1_C3=0;
  vector<float> *recojet_beta1_D2=0;
  vector<float> *recojet_beta1_N2_E_theta=0;
  vector<float> *recojet_beta1_N3_E_theta=0;
  vector<float> *recojet_beta1_C2_E_theta=0;
  vector<float> *recojet_beta1_C3_E_theta=0;
  vector<float> *recojet_beta1_D2_E_theta=0;

  vector<float> *recojet_beta2_ECorr2=0;
  vector<float> *recojet_beta2_ECorr3=0;
  vector<float> *recojet_beta2_ECorr2_E_theta=0;
  vector<float> *recojet_beta2_ECorr3_E_theta=0;

  vector<float> *recojet_beta2_N2=0;
  vector<float> *recojet_beta2_N3=0;
  vector<float> *recojet_beta2_C2=0;
  vector<float> *recojet_beta2_C3=0;
  vector<float> *recojet_beta2_D2=0;
  vector<float> *recojet_beta2_N2_E_theta=0;
  vector<float> *recojet_beta2_N3_E_theta=0;
  vector<float> *recojet_beta2_C2_E_theta=0;
  vector<float> *recojet_beta2_C3_E_theta=0;
  vector<float> *recojet_beta2_D2_E_theta=0;

  vector<float> *recojet_beta0_5_ECorr2=0;
  vector<float> *recojet_beta0_5_ECorr3=0;
  vector<float> *recojet_beta0_5_ECorr2_E_theta=0;
  vector<float> *recojet_beta0_5_ECorr3_E_theta=0;

  vector<float> *recojet_beta0_5_N2=0;
  vector<float> *recojet_beta0_5_N3=0;
  vector<float> *recojet_beta0_5_C2=0;
  vector<float> *recojet_beta0_5_C3=0;
  vector<float> *recojet_beta0_5_D2=0;
  vector<float> *recojet_beta0_5_N2_E_theta=0;
  vector<float> *recojet_beta0_5_N3_E_theta=0;
  vector<float> *recojet_beta0_5_C2_E_theta=0;
  vector<float> *recojet_beta0_5_C3_E_theta=0;
  vector<float> *recojet_beta0_5_D2_E_theta=0;

  vector<float> *genjet_subjet_Px=0;
  vector<float> *genjet_subjet_Py=0;
  vector<float> *genjet_subjet_Pz=0;
  vector<float> *genjet_subjet_E=0;
  vector<int> *genjet_subjet_jetindex=0;

  vector<float> *genjet_nsubjettiness1=0;
  vector<float> *genjet_nsubjettiness2=0;
  vector<float> *genjet_nsubjettiness3=0;
  vector<float> *genjet_nsubjettiness1_lrz=0;
  vector<float> *genjet_nsubjettiness2_lrz=0;
  vector<float> *genjet_nsubjettiness3_lrz=0;

  vector<float> *genjet_beta1_ECorr2=0;
  vector<float> *genjet_beta1_ECorr3=0;
  vector<float> *genjet_beta1_ECorr2_E_theta=0;
  vector<float> *genjet_beta1_ECorr3_E_theta=0; 
  vector<float> *genjet_beta1_N2=0;
  vector<float> *genjet_beta1_N3=0;
  vector<float> *genjet_beta1_C2=0;
  vector<float> *genjet_beta1_C3=0;
  vector<float> *genjet_beta1_D2=0;
  vector<float> *genjet_beta1_N2_E_theta=0;
  vector<float> *genjet_beta1_N3_E_theta=0;
  vector<float> *genjet_beta1_C2_E_theta=0;
  vector<float> *genjet_beta1_C3_E_theta=0;
  vector<float> *genjet_beta1_D2_E_theta=0;

  vector<float> *genjet_beta2_ECorr2=0;
  vector<float> *genjet_beta2_ECorr3=0;
  vector<float> *genjet_beta2_ECorr2_E_theta=0;
  vector<float> *genjet_beta2_ECorr3_E_theta=0; 
  vector<float> *genjet_beta2_N2=0;
  vector<float> *genjet_beta2_N3=0;
  vector<float> *genjet_beta2_C2=0;
  vector<float> *genjet_beta2_C3=0;
  vector<float> *genjet_beta2_D2=0;
  vector<float> *genjet_beta2_N2_E_theta=0;
  vector<float> *genjet_beta2_N3_E_theta=0;
  vector<float> *genjet_beta2_C2_E_theta=0;
  vector<float> *genjet_beta2_C3_E_theta=0;
  vector<float> *genjet_beta2_D2_E_theta=0;

  vector<float> *genjet_beta0_5_ECorr2=0;
  vector<float> *genjet_beta0_5_ECorr3=0;
  vector<float> *genjet_beta0_5_ECorr2_E_theta=0;
  vector<float> *genjet_beta0_5_ECorr3_E_theta=0; 
  vector<float> *genjet_beta0_5_N2=0;
  vector<float> *genjet_beta0_5_N3=0;
  vector<float> *genjet_beta0_5_C2=0;
  vector<float> *genjet_beta0_5_C3=0;
  vector<float> *genjet_beta0_5_D2=0;
  vector<float> *genjet_beta0_5_N2_E_theta=0;
  vector<float> *genjet_beta0_5_N3_E_theta=0;
  vector<float> *genjet_beta0_5_C2_E_theta=0;
  vector<float> *genjet_beta0_5_C3_E_theta=0;
  vector<float> *genjet_beta0_5_D2_E_theta=0;

  vector<float> *trueME_Px=0;
  vector<float> *trueME_Py=0;
  vector<float> *trueME_Pz=0;
  vector<float> *trueME_E=0;
  vector<int> *trueME_PDGID=0;
 
  if(usePartonInfo){
    tree->SetBranchAddress("trueME_E", &trueME_E);
    tree->SetBranchAddress("trueME_Px", &trueME_Px);
    tree->SetBranchAddress("trueME_Py", &trueME_Py);
    tree->SetBranchAddress("trueME_Pz", &trueME_Pz);
    tree->SetBranchAddress("trueME_PDGID", &trueME_PDGID);
  }
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
   //cuts are set to 10 GeV in creation here
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
  tree->SetBranchAddress("genjet_CHFraction", &genjet_CHFraction);
  tree->SetBranchAddress("genjet_PhFraction", &genjet_PhFraction);
  tree->SetBranchAddress("genjet_NCH", &genjet_NCH);
  tree->SetBranchAddress("genjet_NPh", &genjet_NPh);
  tree->SetBranchAddress("genjet_Mult", &genjet_Mult);

  tree->SetBranchAddress("genjet_nsubjettiness1", &genjet_nsubjettiness1);
  tree->SetBranchAddress("genjet_nsubjettiness2", &genjet_nsubjettiness2);
  tree->SetBranchAddress("genjet_nsubjettiness3", &genjet_nsubjettiness3);
  tree->SetBranchAddress("genjet_nsubjettiness1_lrz", &genjet_nsubjettiness1_lrz);
  tree->SetBranchAddress("genjet_nsubjettiness2_lrz", &genjet_nsubjettiness2_lrz);
  tree->SetBranchAddress("genjet_nsubjettiness3_lrz", &genjet_nsubjettiness3_lrz);

  tree->SetBranchAddress("genjet_beta1_ECorr2", &genjet_beta1_ECorr2);
  tree->SetBranchAddress("genjet_beta1_ECorr3", &genjet_beta1_ECorr3);
  tree->SetBranchAddress("genjet_beta1_ECorr2_E_theta", &genjet_beta1_ECorr2_E_theta);
  tree->SetBranchAddress("genjet_beta1_ECorr3_E_theta", &genjet_beta1_ECorr3_E_theta);

  tree->SetBranchAddress("genjet_beta1_N2", &genjet_beta1_N2);
  tree->SetBranchAddress("genjet_beta1_N3", &genjet_beta1_N3);
  tree->SetBranchAddress("genjet_beta1_C2", &genjet_beta1_C2);
  tree->SetBranchAddress("genjet_beta1_C3", &genjet_beta1_C3);
  tree->SetBranchAddress("genjet_beta1_D2", &genjet_beta1_D2);
  tree->SetBranchAddress("genjet_beta1_N2_E_theta", &genjet_beta1_N2_E_theta);
  tree->SetBranchAddress("genjet_beta1_N3_E_theta", &genjet_beta1_N3_E_theta);
  tree->SetBranchAddress("genjet_beta1_C2_E_theta", &genjet_beta1_C2_E_theta);
  tree->SetBranchAddress("genjet_beta1_C3_E_theta", &genjet_beta1_C3_E_theta);
  tree->SetBranchAddress("genjet_beta1_D2_E_theta", &genjet_beta1_D2_E_theta);

  tree->SetBranchAddress("genjet_beta2_ECorr2", &genjet_beta2_ECorr2);
  tree->SetBranchAddress("genjet_beta2_ECorr3", &genjet_beta2_ECorr3);
  tree->SetBranchAddress("genjet_beta2_ECorr2_E_theta", &genjet_beta2_ECorr2_E_theta);
  tree->SetBranchAddress("genjet_beta2_ECorr3_E_theta", &genjet_beta2_ECorr3_E_theta);

  tree->SetBranchAddress("genjet_beta2_N2", &genjet_beta2_N2);
  tree->SetBranchAddress("genjet_beta2_N3", &genjet_beta2_N3);
  tree->SetBranchAddress("genjet_beta2_C2", &genjet_beta2_C2);
  tree->SetBranchAddress("genjet_beta2_C3", &genjet_beta2_C3);
  tree->SetBranchAddress("genjet_beta2_D2", &genjet_beta2_D2);
  tree->SetBranchAddress("genjet_beta2_N2_E_theta", &genjet_beta2_N2_E_theta);
  tree->SetBranchAddress("genjet_beta2_N3_E_theta", &genjet_beta2_N3_E_theta);
  tree->SetBranchAddress("genjet_beta2_C2_E_theta", &genjet_beta2_C2_E_theta);
  tree->SetBranchAddress("genjet_beta2_C3_E_theta", &genjet_beta2_C3_E_theta);
  tree->SetBranchAddress("genjet_beta2_D2_E_theta", &genjet_beta2_D2_E_theta);

  tree->SetBranchAddress("genjet_beta0_5_ECorr2", &genjet_beta0_5_ECorr2);
  tree->SetBranchAddress("genjet_beta0_5_ECorr3", &genjet_beta0_5_ECorr3);
  tree->SetBranchAddress("genjet_beta0_5_ECorr2_E_theta", &genjet_beta0_5_ECorr2_E_theta);
  tree->SetBranchAddress("genjet_beta0_5_ECorr3_E_theta", &genjet_beta0_5_ECorr3_E_theta);

  tree->SetBranchAddress("genjet_beta0_5_N2", &genjet_beta0_5_N2);
  tree->SetBranchAddress("genjet_beta0_5_N3", &genjet_beta0_5_N3);
  tree->SetBranchAddress("genjet_beta0_5_C2", &genjet_beta0_5_C2);
  tree->SetBranchAddress("genjet_beta0_5_C3", &genjet_beta0_5_C3);
  tree->SetBranchAddress("genjet_beta0_5_D2", &genjet_beta0_5_D2);
  tree->SetBranchAddress("genjet_beta0_5_N2_E_theta", &genjet_beta0_5_N2_E_theta);
  tree->SetBranchAddress("genjet_beta0_5_N3_E_theta", &genjet_beta0_5_N3_E_theta);
  tree->SetBranchAddress("genjet_beta0_5_C2_E_theta", &genjet_beta0_5_C2_E_theta);
  tree->SetBranchAddress("genjet_beta0_5_C3_E_theta", &genjet_beta0_5_C3_E_theta);
  tree->SetBranchAddress("genjet_beta0_5_D2_E_theta", &genjet_beta0_5_D2_E_theta);

  tree->SetBranchAddress("genjet_subjet_E", &genjet_subjet_E);
  tree->SetBranchAddress("genjet_subjet_Px", &genjet_subjet_Px);
  tree->SetBranchAddress("genjet_subjet_Py", &genjet_subjet_Py);
  tree->SetBranchAddress("genjet_subjet_Pz", &genjet_subjet_Pz);
  tree->SetBranchAddress("genjet_subjet_jetindex", &genjet_subjet_jetindex);

  tree->SetBranchAddress("recojet_E", &recojet_E);
  tree->SetBranchAddress("recojet_Px", &recojet_Px);
  tree->SetBranchAddress("recojet_Py", &recojet_Py);
  tree->SetBranchAddress("recojet_Pz", &recojet_Pz);
  tree->SetBranchAddress("recojet_CHFraction", &recojet_CHFraction);
  tree->SetBranchAddress("recojet_PhFraction", &recojet_PhFraction);
  tree->SetBranchAddress("recojet_NCH", &recojet_NCH);
  tree->SetBranchAddress("recojet_NPh", &recojet_NPh);
  tree->SetBranchAddress("recojet_Mult", &recojet_Mult);

  tree->SetBranchAddress("recojet_nsubjettiness1", &recojet_nsubjettiness1);
  tree->SetBranchAddress("recojet_nsubjettiness2", &recojet_nsubjettiness2);
  tree->SetBranchAddress("recojet_nsubjettiness3", &recojet_nsubjettiness3);
  tree->SetBranchAddress("recojet_nsubjettiness1_lrz", &recojet_nsubjettiness1_lrz);
  tree->SetBranchAddress("recojet_nsubjettiness2_lrz", &recojet_nsubjettiness2_lrz);
  tree->SetBranchAddress("recojet_nsubjettiness3_lrz", &recojet_nsubjettiness3_lrz);

  tree->SetBranchAddress("recojet_beta1_ECorr2", &recojet_beta1_ECorr2);
  tree->SetBranchAddress("recojet_beta1_ECorr3", &recojet_beta1_ECorr3);
  tree->SetBranchAddress("recojet_beta1_ECorr2_E_theta", &recojet_beta1_ECorr2_E_theta);
  tree->SetBranchAddress("recojet_beta1_ECorr3_E_theta", &recojet_beta1_ECorr3_E_theta);

  tree->SetBranchAddress("recojet_beta1_N2", &recojet_beta1_N2);
  tree->SetBranchAddress("recojet_beta1_N3", &recojet_beta1_N3);
  tree->SetBranchAddress("recojet_beta1_C2", &recojet_beta1_C2);
  tree->SetBranchAddress("recojet_beta1_C3", &recojet_beta1_C3);
  tree->SetBranchAddress("recojet_beta1_D2", &recojet_beta1_D2);
  tree->SetBranchAddress("recojet_beta1_N2_E_theta", &recojet_beta1_N2_E_theta);
  tree->SetBranchAddress("recojet_beta1_N3_E_theta", &recojet_beta1_N3_E_theta);
  tree->SetBranchAddress("recojet_beta1_C2_E_theta", &recojet_beta1_C2_E_theta);
  tree->SetBranchAddress("recojet_beta1_C3_E_theta", &recojet_beta1_C3_E_theta);
  tree->SetBranchAddress("recojet_beta1_D2_E_theta", &recojet_beta1_D2_E_theta);

  tree->SetBranchAddress("recojet_beta2_ECorr2", &recojet_beta2_ECorr2);
  tree->SetBranchAddress("recojet_beta2_ECorr3", &recojet_beta2_ECorr3);
  tree->SetBranchAddress("recojet_beta2_ECorr2_E_theta", &recojet_beta2_ECorr2_E_theta);
  tree->SetBranchAddress("recojet_beta2_ECorr3_E_theta", &recojet_beta2_ECorr3_E_theta);

  tree->SetBranchAddress("recojet_beta2_N2", &recojet_beta2_N2);
  tree->SetBranchAddress("recojet_beta2_N3", &recojet_beta2_N3);
  tree->SetBranchAddress("recojet_beta2_C2", &recojet_beta2_C2);
  tree->SetBranchAddress("recojet_beta2_C3", &recojet_beta2_C3);
  tree->SetBranchAddress("recojet_beta2_D2", &recojet_beta2_D2);
  tree->SetBranchAddress("recojet_beta2_N2_E_theta", &recojet_beta2_N2_E_theta);
  tree->SetBranchAddress("recojet_beta2_N3_E_theta", &recojet_beta2_N3_E_theta);
  tree->SetBranchAddress("recojet_beta2_C2_E_theta", &recojet_beta2_C2_E_theta);
  tree->SetBranchAddress("recojet_beta2_C3_E_theta", &recojet_beta2_C3_E_theta);
  tree->SetBranchAddress("recojet_beta2_D2_E_theta", &recojet_beta2_D2_E_theta);

  tree->SetBranchAddress("recojet_beta0_5_ECorr2", &recojet_beta0_5_ECorr2);
  tree->SetBranchAddress("recojet_beta0_5_ECorr3", &recojet_beta0_5_ECorr3);
  tree->SetBranchAddress("recojet_beta0_5_ECorr2_E_theta", &recojet_beta0_5_ECorr2_E_theta);
  tree->SetBranchAddress("recojet_beta0_5_ECorr3_E_theta", &recojet_beta0_5_ECorr3_E_theta);

  tree->SetBranchAddress("recojet_beta0_5_N2", &recojet_beta0_5_N2);
  tree->SetBranchAddress("recojet_beta0_5_N3", &recojet_beta0_5_N3);
  tree->SetBranchAddress("recojet_beta0_5_C2", &recojet_beta0_5_C2);
  tree->SetBranchAddress("recojet_beta0_5_C3", &recojet_beta0_5_C3);
  tree->SetBranchAddress("recojet_beta0_5_D2", &recojet_beta0_5_D2);
  tree->SetBranchAddress("recojet_beta0_5_N2_E_theta", &recojet_beta0_5_N2_E_theta);
  tree->SetBranchAddress("recojet_beta0_5_N3_E_theta", &recojet_beta0_5_N3_E_theta);
  tree->SetBranchAddress("recojet_beta0_5_C2_E_theta", &recojet_beta0_5_C2_E_theta);
  tree->SetBranchAddress("recojet_beta0_5_C3_E_theta", &recojet_beta0_5_C3_E_theta);
  tree->SetBranchAddress("recojet_beta0_5_D2_E_theta", &recojet_beta0_5_D2_E_theta);

  tree->SetBranchAddress("recojet_subjet_E", &recojet_subjet_E);
  tree->SetBranchAddress("recojet_subjet_Px", &recojet_subjet_Px);
  tree->SetBranchAddress("recojet_subjet_Py", &recojet_subjet_Py);
  tree->SetBranchAddress("recojet_subjet_Pz", &recojet_subjet_Pz);
  tree->SetBranchAddress("recojet_subjet_jetindex", &recojet_subjet_jetindex);

  tree->SetBranchAddress("gen_tau_E", &gen_tau_E);
  tree->SetBranchAddress("gen_tau_Px", &gen_tau_Px);
  tree->SetBranchAddress("gen_tau_Py", &gen_tau_Py);
  tree->SetBranchAddress("gen_tau_Pz", &gen_tau_Pz);
  tree->SetBranchAddress("gen_tau_Charge", &gen_tau_Charge);
  tree->SetBranchAddress("gen_tau_NCH", &gen_tau_NCH);
  tree->SetBranchAddress("gen_tau_NPh", &gen_tau_NPh);
  tree->SetBranchAddress("gen_tau_CHFraction", &gen_tau_CHFraction);
  tree->SetBranchAddress("gen_tau_PhFraction", &gen_tau_PhFraction);

  tree->SetBranchAddress("reco_tau_E", &reco_tau_E);
  tree->SetBranchAddress("reco_tau_Px", &reco_tau_Px);
  tree->SetBranchAddress("reco_tau_Py", &reco_tau_Py);
  tree->SetBranchAddress("reco_tau_Pz", &reco_tau_Pz);
  tree->SetBranchAddress("reco_tau_Charge", &reco_tau_Charge);
  tree->SetBranchAddress("reco_tau_NCH", &reco_tau_NCH);
  tree->SetBranchAddress("reco_tau_NPh", &reco_tau_NPh);
  tree->SetBranchAddress("reco_tau_CHFraction", &reco_tau_CHFraction);
  tree->SetBranchAddress("reco_tau_PhFraction", &reco_tau_PhFraction);

  //we should have the following
  //index 0 and 1 are outgoing e+e-
  //index 2 and 3 ISR photons from incoming e+e-
  //index 4 and 5 Z quarks
  //index 6 is Higgs (boosted back to real system)
  //index 7 and 8 are H daughters

  //beam strahlung effects not in history -->take outgoing e+e- for centre of mass system, i.e. use 0 and 1, NOT 2 and 3

  double sqrtS_low=750;
  double sqrtS_high=2500;
  double sqrtS_high_reco=2500;
  double sqrtS_nom=3000;

  //number of alloweed isolated leptons (+1), isolation criteria relative 10% in cone of 10 degrees, lower energy limit on leptons and photons set to 10 GeV
  unsigned int m_cut_nLeptons = 2;

  double alpha_min =200;

  double lumi = 5000.;//cross sections are in femtobarn, at 3 TeV we expect 4 inverse attobarn 

  double weight=x_sec*lumi/(double)tree->GetEntries();
  std::cout<<"weight for sample "<<weight<<std::endl;

  for(unsigned int i_entry=0;i_entry<tree->GetEntries();i_entry++){
  //for(unsigned int i_entry=0;i_entry<0;i_entry++){
    //fill jet energy resolution histograms
    tree->GetEntry(i_entry);

    if(i_entry%5000==0){
      std::cout<<"entry "<<i_entry<<std::endl;
    }


    TLorentzVector tempTotEventP4(0,0,0,0);
    TLorentzVector tempTotEventP4HZ(0,0,0,0);
    TLorentzVector tempZP4(0,0,0,0);
    TLorentzVector tempHP4(0,0,0,0);

    TLorentzVector tempZ_q1(0,0,0,0);
    TLorentzVector tempZ_q2(0,0,0,0);

    TLorentzVector tempH_q1(0,0,0,0);
    TLorentzVector tempH_q2(0,0,0,0);

    bool H_decays_bbar=true;
    if(usePartonInfo){
      if(fabs((*trueME_PDGID)[7])!=5 || fabs((*trueME_PDGID)[8])!=5){
	H_decays_bbar=false;
	std::cout<<"entry "<<i_entry <<" is NO bbar "<<(*trueME_PDGID)[7]<<"/"<<(*trueME_PDGID)[8]<<std::endl;
      }
    }
    if(!H_decays_bbar){
      continue;
    }
    
    if(usePartonInfo){
      //4 and 5 are Z daugthers
      if((*trueME_E)[4]>(*trueME_E)[5]){
	tempZ_q1.SetPxPyPzE((*trueME_Px)[4],(*trueME_Py)[4],(*trueME_Pz)[4],(*trueME_E)[4]);
	tempZ_q2.SetPxPyPzE((*trueME_Px)[5],(*trueME_Py)[5],(*trueME_Pz)[5],(*trueME_E)[5]);
      }else{
	tempZ_q1.SetPxPyPzE((*trueME_Px)[5],(*trueME_Py)[5],(*trueME_Pz)[5],(*trueME_E)[5]);
	tempZ_q2.SetPxPyPzE((*trueME_Px)[4],(*trueME_Py)[4],(*trueME_Pz)[4],(*trueME_E)[4]);
      }
      //7 and 8 are H daugthers
      if((*trueME_E)[7]>(*trueME_E)[8]){
	tempH_q1.SetPxPyPzE((*trueME_Px)[7],(*trueME_Py)[7],(*trueME_Pz)[7],(*trueME_E)[7]);
	tempH_q2.SetPxPyPzE((*trueME_Px)[8],(*trueME_Py)[8],(*trueME_Pz)[8],(*trueME_E)[8]);
      }else{
	tempH_q1.SetPxPyPzE((*trueME_Px)[8],(*trueME_Py)[8],(*trueME_Pz)[8],(*trueME_E)[8]);
	tempH_q2.SetPxPyPzE((*trueME_Px)[7],(*trueME_Py)[7],(*trueME_Pz)[7],(*trueME_E)[7]);
      }
      tempHP4.SetPxPyPzE((*trueME_Px)[6],(*trueME_Py)[6],(*trueME_Pz)[6],(*trueME_E)[6]);
      for(unsigned int i=0;i<trueME_E->size();i++){
	TLorentzVector temp(0,0,0,0);
	temp.SetPxPyPzE((*trueME_Px)[i],(*trueME_Py)[i],(*trueME_Pz)[i],(*trueME_E)[i]);
	if(i<2){
	  //determined by outgoing e+e- pair
	  tempTotEventP4+=temp;
	}
	if(i>3 && i<6){//index 4 and 5 are Z daughters
	  tempZP4+=temp;
	}
      }
      tempTotEventP4HZ=tempZP4+tempHP4;
      float H_b_theta_FW=-1;
      float Z_q_theta_FW=-1;
      //calculation
      if(fabs(tempH_q1.Theta()-TMath::Pi()*0.5)>fabs(tempH_q2.Theta()-TMath::Pi()*0.5)){
	H_b_theta_FW=tempH_q1.Theta()*TMath::RadToDeg();
      }else{
	H_b_theta_FW=tempH_q2.Theta()*TMath::RadToDeg();
      }
     //calculation
      if(fabs(tempZ_q1.Theta()-TMath::Pi()*0.5)>fabs(tempZ_q2.Theta()-TMath::Pi()*0.5)){
	Z_q_theta_FW=tempZ_q1.Theta()*TMath::RadToDeg();
      }else{
	Z_q_theta_FW=tempZ_q2.Theta()*TMath::RadToDeg();
      }
      float Z_q_H_b_angle_min=tempH_q1.Angle(tempZ_q1.Vect());
      if(Z_q_H_b_angle_min>tempH_q1.Angle(tempZ_q2.Vect())){
	Z_q_H_b_angle_min=tempH_q1.Angle(tempZ_q2.Vect());
      }
      if(Z_q_H_b_angle_min>tempH_q2.Angle(tempZ_q1.Vect())){
	Z_q_H_b_angle_min=tempH_q2.Angle(tempZ_q1.Vect());
      }
     if(Z_q_H_b_angle_min>tempH_q2.Angle(tempZ_q2.Vect())){
	Z_q_H_b_angle_min=tempH_q2.Angle(tempZ_q2.Vect());
      }
     if(fill_partonInfo){
       h_hist_parton[0]->Fill(tempTotEventP4.M(),weight);  
       h_hist_parton[37]->Fill((tempHP4+tempZP4).M(),weight);
       if(tempZ_q1.E()>tempZ_q2.E()){
	 h_hist_parton[54]->Fill(tempZ_q1.E()/tempZP4.E(),weight);
	 h_hist_parton[55]->Fill(tempZ_q2.E()/tempZP4.E(),weight);
       }else{
	 h_hist_parton[54]->Fill(tempZ_q2.E()/tempZP4.E(),weight);
	 h_hist_parton[55]->Fill(tempZ_q1.E()/tempZP4.E(),weight);
       }
       if(tempH_q1.E()>tempH_q2.E()){
	 h_hist_parton[56]->Fill(tempH_q1.E()/tempHP4.E(),weight);
	 h_hist_parton[57]->Fill(tempH_q2.E()/tempHP4.E(),weight);
       }else{
	 h_hist_parton[56]->Fill(tempH_q2.E()/tempHP4.E(),weight);
	 h_hist_parton[57]->Fill(tempH_q1.E()/tempHP4.E(),weight);
       }
       if(tempTotEventP4.M()<sqrtS_low){
	 h_hist_parton[1]->Fill(tempHP4.Pt(),weight);
	 h_hist_parton[4]->Fill(tempZP4.Pt(),weight);
	 h_hist_parton[7]->Fill(tempH_q1.Angle(tempH_q2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[10]->Fill(tempZ_q1.Angle(tempZ_q2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[13]->Fill(H_b_theta_FW,weight);//angle calculated in degrees previously
	 h_hist_parton[16]->Fill(Z_q_theta_FW,weight);//angle calculated in degrees previously
	 h_hist_parton[19]->Fill(Z_q_H_b_angle_min*TMath::RadToDeg(),weight);
	 h_hist_parton[22]->Fill(tempHP4.Theta()*TMath::RadToDeg(),weight);
	 h_hist_parton[25]->Fill(tempZP4.Theta()*TMath::RadToDeg(),weight);
	 h_hist_parton[28]->Fill(fabs(tempZP4.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	 h_hist_parton[31]->Fill(DeltaPhi(tempZP4.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	 h_hist_parton[34]->Fill(tempZP4.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	 
	 
       }else if(tempTotEventP4.M()<sqrtS_high){
	 h_hist_parton[2]->Fill(tempHP4.Pt(),weight);
	 h_hist_parton[5]->Fill(tempZP4.Pt(),weight);
	 h_hist_parton[8]->Fill(tempH_q1.Angle(tempH_q2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[11]->Fill(tempZ_q1.Angle(tempZ_q2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[14]->Fill(H_b_theta_FW,weight);//angle calculated in degrees previously
	 h_hist_parton[17]->Fill(Z_q_theta_FW,weight);//angle calculated in degrees previously
	 h_hist_parton[20]->Fill(Z_q_H_b_angle_min*TMath::RadToDeg(),weight);
	 h_hist_parton[23]->Fill(tempHP4.Theta()*TMath::RadToDeg(),weight);
	 h_hist_parton[26]->Fill(tempZP4.Theta()*TMath::RadToDeg(),weight);
	 h_hist_parton[29]->Fill(fabs(tempZP4.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	 h_hist_parton[32]->Fill(DeltaPhi(tempZP4.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	 h_hist_parton[35]->Fill(tempZP4.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	 
	 
       }else{
	 h_hist_parton[3]->Fill(tempHP4.Pt(),weight);
	 h_hist_parton[6]->Fill(tempZP4.Pt(),weight);
	 h_hist_parton[9]->Fill(tempH_q1.Angle(tempH_q2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[12]->Fill(tempZ_q1.Angle(tempZ_q2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[15]->Fill(H_b_theta_FW,weight);//angle calculated in degrees previously
	 h_hist_parton[18]->Fill(Z_q_theta_FW,weight);//angle calculated in degrees previously
	 h_hist_parton[21]->Fill(Z_q_H_b_angle_min*TMath::RadToDeg(),weight);
	 h_hist_parton[24]->Fill(tempHP4.Theta()*TMath::RadToDeg(),weight);
	 h_hist_parton[27]->Fill(tempZP4.Theta()*TMath::RadToDeg(),weight);
	 h_hist_parton[30]->Fill(fabs(tempZP4.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	 h_hist_parton[33]->Fill(DeltaPhi(tempZP4.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	 h_hist_parton[36]->Fill(tempZP4.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
       }
     }//if parton info is to be filled      
    }


    TLorentzVector tempTotGenP4(0,0,0,0);
    tempTotGenP4.SetPxPyPzE(true_Px,true_Py,true_Pz,true_E);

    TLorentzVector tempTotInvGenP4(0,0,0,0);
    tempTotInvGenP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);

    TLorentzVector tempRecoMETP4(0,0,0,0);
    tempRecoMETP4.SetPxPyPzE(-totPFO_Px,-totPFO_Py,-totPFO_Pz,sqrtS_nom-totPFO_E);

    TLorentzVector tempTotRecoP4(0,0,0,0);
    tempTotRecoP4.SetPxPyPzE(totPFO_Px,totPFO_Py,totPFO_Pz,totPFO_E);

    h_hist_vec_reco[83]->Fill(tempRecoMETP4.Pt(),weight);

    unsigned int n_IsoLep_gen=0;
    unsigned int n_IsoLep_reco=0;

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
    
    //here the threshold are set at filling the ntuple, to 10 GeV
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
	  n_IsoLep_gen+=1;
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
	    n_IsoLep_reco+=1;
	  }
	}  
      }
    }
    TLorentzVector tempGenJetSum(0,0,0,0);
    for(unsigned int i=0;i<genjet_E->size();i++){
     TLorentzVector temp(0,0,0,0);
     temp.SetPxPyPzE((*genjet_Px)[i],(*genjet_Py)[i],(*genjet_Pz)[i],(*genjet_E)[i]);
     tempGenJetSum+=temp;
    }
    TLorentzVector tempRecoJetSum(0,0,0,0);
    for(unsigned int i=0;i<recojet_E->size();i++){
     TLorentzVector temp(0,0,0,0);
     temp.SetPxPyPzE((*recojet_Px)[i],(*recojet_Py)[i],(*recojet_Pz)[i],(*recojet_E)[i]);
     tempRecoJetSum+=temp;
    }
    TLorentzVector gj_m1(0,0,0,0);
    TLorentzVector gj_m2(0,0,0,0);
    unsigned int ind_gj1=0;
    unsigned int ind_gj2=1;
    if(genjet_E->size()==2){
     gj_m1.SetPxPyPzE((*genjet_Px)[0],(*genjet_Py)[0],(*genjet_Pz)[0],(*genjet_E)[0]);
     gj_m2.SetPxPyPzE((*genjet_Px)[1],(*genjet_Py)[1],(*genjet_Pz)[1],(*genjet_E)[1]);
     if(gj_m2.M()>gj_m1.M()){
       ind_gj1=1;
       ind_gj2=0;
       TLorentzVector temp=gj_m1;
       gj_m1=gj_m2;
       gj_m2=temp;
     }
     if(gj_m2.M()>=gj_m1.M()){
       std::cout<<"gj mass order should have not been the case nowadays "<<gj_m2.M()<<"/"<<gj_m1.M()<<std::endl;
     }
    if(usePartonInfo && fill_partonInfo && fill_genInfo){
      if(tempTotEventP4.M()>sqrtS_high){
	h_hist_parton[38]->Fill(gj_m1.M(),weight);
	h_hist_parton[39]->Fill(gj_m1.M(),weight);
	h_hist_parton[42]->Fill(DeltaPhi(gj_m1.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	h_hist_parton[43]->Fill(fabs(gj_m1.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	h_hist_parton[44]->Fill(DeltaPhi(gj_m1.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	h_hist_parton[45]->Fill(fabs(gj_m1.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	h_hist_parton[50]->Fill(gj_m1.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	h_hist_parton[51]->Fill(gj_m2.Angle(tempZP4.Vect())*TMath::RadToDeg(),weight);
      }
    }
    int ind_sj1_gj1=-1;
    float E_sj1_gj1=-1;//jets forced into two subjets, thus only one check necessary
    int ind_sj2_gj1=-1;
    int ind_sj1_gj2=-1;
    float E_sj1_gj2=-1;
    int ind_sj2_gj2=-1;

    TLorentzVector temp_sj1_gj1(0,0,0,0);
    TLorentzVector temp_sj2_gj1(0,0,0,0);
    TLorentzVector temp_sj1_gj2(0,0,0,0);
    TLorentzVector temp_sj2_gj2(0,0,0,0);
    
    //if too few components in jet, then jet index NOT found -->i.e. remains at -1
    for(unsigned int i=0;i<genjet_subjet_E->size();i++){
      if((*genjet_subjet_jetindex)[i]==ind_gj1){
	if((*genjet_subjet_E)[i]>E_sj1_gj1){
	  ind_sj2_gj1=ind_sj1_gj1;
	  ind_sj1_gj1=i;
	}else{
	  ind_sj2_gj1=i;
	}
      }
      if((*genjet_subjet_jetindex)[i]==ind_gj2){
	if((*genjet_subjet_E)[i]>E_sj1_gj2){
	  ind_sj2_gj2=ind_sj1_gj2;
	  ind_sj1_gj2=i;
	}else{
	  ind_sj2_gj2=i;
	}
      }
    }
    if(ind_sj1_gj1!=-1){
      temp_sj1_gj1.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj1],(*genjet_subjet_Py)[ind_sj1_gj1],(*genjet_subjet_Pz)[ind_sj1_gj1],(*genjet_subjet_E)[ind_sj1_gj1]);
    }
    if(ind_sj2_gj1!=-1){
      temp_sj2_gj1.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj1],(*genjet_subjet_Py)[ind_sj2_gj1],(*genjet_subjet_Pz)[ind_sj2_gj1],(*genjet_subjet_E)[ind_sj2_gj1]);
    }
    if(ind_sj1_gj2!=-1){
      temp_sj1_gj2.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj2],(*genjet_subjet_Py)[ind_sj1_gj2],(*genjet_subjet_Pz)[ind_sj1_gj2],(*genjet_subjet_E)[ind_sj1_gj2]);
    }
    if(ind_sj2_gj2!=-1){
      temp_sj2_gj2.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj2],(*genjet_subjet_Py)[ind_sj2_gj2],(*genjet_subjet_Pz)[ind_sj2_gj2],(*genjet_subjet_E)[ind_sj2_gj2]);
    }

    if(fill_genInfo && genjet_E->size()==2 && n_IsoLep_gen<m_cut_nLeptons){//effectively no cut here, but exactly two genjets
      if((tempTotGenP4-tempGenIsoPhP4).M()<sqrtS_low){
	h_hist_vec_gen[0]->Fill(gj_m1.Angle(gj_m2.Vect())*TMath::RadToDeg(),weight);
	h_hist_vec_gen[3]->Fill(DeltaPhi(gj_m1.Phi(),gj_m2.Phi())*TMath::RadToDeg(),weight);
	h_hist_vec_gen[6]->Fill(fabs(gj_m1.Theta()-gj_m2.Theta())*TMath::RadToDeg(),weight);
	h_hist_vec_gen[9]->Fill(gj_m1.Theta()*TMath::RadToDeg(),weight);
	h_hist_vec_gen[12]->Fill(gj_m2.Theta()*TMath::RadToDeg(),weight);
	h_hist_vec_gen[15]->Fill(gj_m1.M(),weight);
	h_hist_vec_gen[18]->Fill(gj_m2.M(),weight);
	if((*genjet_nsubjettiness2)[ind_gj1]!=-1){//set to -1 if too few constitutents
	  h_hist_vec_gen[21]->Fill((*genjet_nsubjettiness2)[ind_gj1]/(*genjet_nsubjettiness1)[ind_gj1],weight);
	}
	if((*genjet_nsubjettiness2)[ind_gj1]!=-1){//set to -1 if too few constitutents
	  h_hist_vec_gen[24]->Fill((*genjet_nsubjettiness2_lrz)[ind_gj1]/(*genjet_nsubjettiness1_lrz)[ind_gj1],weight);
	}
	if((*genjet_nsubjettiness3)[ind_gj1]!=-1){//set to -1 if too few constitutents
	  h_hist_vec_gen[27]->Fill((*genjet_nsubjettiness3)[ind_gj1]/(*genjet_nsubjettiness2)[ind_gj1],weight);
	}
	if((*genjet_nsubjettiness3)[ind_gj1]!=-1){//set to -1 if too few constitutents
	  h_hist_vec_gen[30]->Fill((*genjet_nsubjettiness3_lrz)[ind_gj1]/(*genjet_nsubjettiness2_lrz)[ind_gj1],weight);
	}
	if((*genjet_nsubjettiness2)[ind_gj2]!=-1){//set to -1 if too few constitutents
	  h_hist_vec_gen[33]->Fill((*genjet_nsubjettiness2)[ind_gj2]/(*genjet_nsubjettiness1)[ind_gj2],weight);
	}
	if((*genjet_nsubjettiness2)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[36]->Fill((*genjet_nsubjettiness2_lrz)[ind_gj2]/(*genjet_nsubjettiness1_lrz)[ind_gj2],weight);
	 }
	 if((*genjet_nsubjettiness3)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[39]->Fill((*genjet_nsubjettiness3)[ind_gj2]/(*genjet_nsubjettiness2)[ind_gj2],weight);
	 }
	 if((*genjet_nsubjettiness3)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[42]->Fill((*genjet_nsubjettiness3_lrz)[ind_gj2]/(*genjet_nsubjettiness2_lrz)[ind_gj2],weight);
	 }
	 h_hist_vec_gen[45]->Fill((*genjet_beta1_N2)[ind_gj1],weight);
	 h_hist_vec_gen[48]->Fill((*genjet_beta1_N2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[51]->Fill((*genjet_beta2_N2)[ind_gj1],weight);
	 h_hist_vec_gen[54]->Fill((*genjet_beta2_N2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[57]->Fill((*genjet_beta0_5_N2)[ind_gj1],weight);
	 h_hist_vec_gen[60]->Fill((*genjet_beta0_5_N2_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[63]->Fill((*genjet_beta1_N2)[ind_gj2],weight);
	 h_hist_vec_gen[66]->Fill((*genjet_beta1_N2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[69]->Fill((*genjet_beta2_N2)[ind_gj2],weight);
	 h_hist_vec_gen[72]->Fill((*genjet_beta2_N2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[75]->Fill((*genjet_beta0_5_N2)[ind_gj2],weight);
	 h_hist_vec_gen[78]->Fill((*genjet_beta0_5_N2_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[81]->Fill((*genjet_beta1_N3)[ind_gj1],weight);
	 h_hist_vec_gen[84]->Fill((*genjet_beta1_N3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[87]->Fill((*genjet_beta2_N3)[ind_gj1],weight);
	 h_hist_vec_gen[90]->Fill((*genjet_beta2_N3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[93]->Fill((*genjet_beta0_5_N3)[ind_gj1],weight);
	 h_hist_vec_gen[96]->Fill((*genjet_beta0_5_N3_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[99]->Fill((*genjet_beta1_N3)[ind_gj2],weight);
	 h_hist_vec_gen[102]->Fill((*genjet_beta1_N3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[105]->Fill((*genjet_beta2_N3)[ind_gj2],weight);
	 h_hist_vec_gen[108]->Fill((*genjet_beta2_N3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[111]->Fill((*genjet_beta0_5_N3)[ind_gj2],weight);
	 h_hist_vec_gen[114]->Fill((*genjet_beta0_5_N3_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[117]->Fill((*genjet_beta1_C2)[ind_gj1],weight);
	 h_hist_vec_gen[120]->Fill((*genjet_beta1_C2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[123]->Fill((*genjet_beta2_C2)[ind_gj1],weight);
	 h_hist_vec_gen[126]->Fill((*genjet_beta2_C2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[129]->Fill((*genjet_beta0_5_C2)[ind_gj1],weight);
	 h_hist_vec_gen[132]->Fill((*genjet_beta0_5_C2_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[135]->Fill((*genjet_beta1_C2)[ind_gj2],weight);
	 h_hist_vec_gen[138]->Fill((*genjet_beta1_C2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[141]->Fill((*genjet_beta2_C2)[ind_gj2],weight);
	 h_hist_vec_gen[144]->Fill((*genjet_beta2_C2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[147]->Fill((*genjet_beta0_5_C2)[ind_gj2],weight);
	 h_hist_vec_gen[150]->Fill((*genjet_beta0_5_C2_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[153]->Fill((*genjet_beta1_C3)[ind_gj1],weight);
	 h_hist_vec_gen[156]->Fill((*genjet_beta1_C3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[159]->Fill((*genjet_beta2_C3)[ind_gj1],weight);
	 h_hist_vec_gen[162]->Fill((*genjet_beta2_C3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[165]->Fill((*genjet_beta0_5_C3)[ind_gj1],weight);
	 h_hist_vec_gen[168]->Fill((*genjet_beta0_5_C3_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[171]->Fill((*genjet_beta1_C3)[ind_gj2],weight);
	 h_hist_vec_gen[174]->Fill((*genjet_beta1_C3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[177]->Fill((*genjet_beta2_C3)[ind_gj2],weight);
	 h_hist_vec_gen[180]->Fill((*genjet_beta2_C3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[183]->Fill((*genjet_beta0_5_C3)[ind_gj2],weight);
	 h_hist_vec_gen[186]->Fill((*genjet_beta0_5_C3_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[189]->Fill((*genjet_beta1_D2)[ind_gj1],weight);
	 h_hist_vec_gen[192]->Fill((*genjet_beta1_D2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[195]->Fill((*genjet_beta2_D2)[ind_gj1],weight);
	 h_hist_vec_gen[198]->Fill((*genjet_beta2_D2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[201]->Fill((*genjet_beta0_5_D2)[ind_gj1],weight);
	 h_hist_vec_gen[204]->Fill((*genjet_beta0_5_D2_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[207]->Fill((*genjet_beta1_D2)[ind_gj2],weight);
	 h_hist_vec_gen[210]->Fill((*genjet_beta1_D2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[213]->Fill((*genjet_beta2_D2)[ind_gj2],weight);
	 h_hist_vec_gen[216]->Fill((*genjet_beta2_D2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[219]->Fill((*genjet_beta0_5_D2)[ind_gj2],weight);
	 h_hist_vec_gen[222]->Fill((*genjet_beta0_5_D2_E_theta)[ind_gj2],weight);
	 //D_{2}^{(1,2)} defined as e_3^{(1)}/(e_2^{(2)})^3/2
	 if((*genjet_beta2_ECorr2)[ind_gj1]!=0){
	   h_hist_vec_gen[225]->Fill((*genjet_beta1_ECorr3)[ind_gj1]/pow((*genjet_beta2_ECorr2)[ind_gj1],3./2.),weight);
	 }else{
	   h_hist_vec_gen[225]->Fill(1.e-5);
	 }
	 if((*genjet_beta2_ECorr2_E_theta)[ind_gj1]!=0){
	   h_hist_vec_gen[228]->Fill((*genjet_beta1_ECorr3_E_theta)[ind_gj1]/pow((*genjet_beta2_ECorr2_E_theta)[ind_gj1],3./2),weight);
	 }else{
	   h_hist_vec_gen[228]->Fill(1.e-5);
	 }
	 if((*genjet_beta2_ECorr2)[ind_gj2]!=0){
	   h_hist_vec_gen[231]->Fill((*genjet_beta1_ECorr3)[ind_gj2]/pow((*genjet_beta2_ECorr2)[ind_gj2],3./2.),weight);
	 }else{
	   h_hist_vec_gen[231]->Fill(1.e-5);
	 }
	 if((*genjet_beta2_ECorr2_E_theta)[ind_gj2]!=0){
	   h_hist_vec_gen[234]->Fill((*genjet_beta1_ECorr3_E_theta)[ind_gj2]/pow((*genjet_beta2_ECorr2_E_theta)[ind_gj2],3./2),weight);
	 }else{
	   h_hist_vec_gen[234]->Fill(1.e-5);
	 }
	 if(ind_sj1_gj1!=-1){
	   h_hist_vec_gen[237]->Fill((*genjet_subjet_E)[ind_sj1_gj1],weight);
	   h_hist_vec_gen[249]->Fill((*genjet_subjet_E)[ind_sj1_gj1]/gj_m1.E(),weight);
	   if(ind_sj2_gj1!=-1){
	     h_hist_vec_gen[255]->Fill(temp_sj1_gj1.Angle(temp_sj2_gj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj1!=-1){
	   h_hist_vec_gen[240]->Fill((*genjet_subjet_E)[ind_sj2_gj1],weight);
	 }
	 if(ind_sj1_gj2!=-1){
	   h_hist_vec_gen[243]->Fill((*genjet_subjet_E)[ind_sj1_gj2],weight);
	   h_hist_vec_gen[252]->Fill((*genjet_subjet_E)[ind_sj1_gj2]/gj_m2.E(),weight);
	   if(ind_sj2_gj2!=-1){
	     h_hist_vec_gen[258]->Fill(temp_sj1_gj2.Angle(temp_sj2_gj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj2!=-1){
	   h_hist_vec_gen[246]->Fill((*genjet_subjet_E)[ind_sj2_gj2],weight);
	 }

       }else if((tempTotGenP4-tempGenIsoPhP4).M()<sqrtS_high_reco){
	 h_hist_vec_gen[1]->Fill(gj_m1.Angle(gj_m2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[4]->Fill(DeltaPhi(gj_m1.Phi(),gj_m2.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[7]->Fill(fabs(gj_m1.Theta()-gj_m2.Theta())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[10]->Fill(gj_m1.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[13]->Fill(gj_m2.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[16]->Fill(gj_m1.M(),weight);
	 h_hist_vec_gen[19]->Fill(gj_m2.M(),weight);
	 if((*genjet_nsubjettiness2)[ind_gj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[22]->Fill((*genjet_nsubjettiness2)[ind_gj1]/(*genjet_nsubjettiness1)[ind_gj1],weight);
	 }
	 if((*genjet_nsubjettiness2)[ind_gj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[25]->Fill((*genjet_nsubjettiness2_lrz)[ind_gj1]/(*genjet_nsubjettiness1_lrz)[ind_gj1],weight);
	 }
	 if((*genjet_nsubjettiness3)[ind_gj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[28]->Fill((*genjet_nsubjettiness3)[ind_gj1]/(*genjet_nsubjettiness2)[ind_gj1],weight);
	 }
	 if((*genjet_nsubjettiness3)[ind_gj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[31]->Fill((*genjet_nsubjettiness3_lrz)[ind_gj1]/(*genjet_nsubjettiness2_lrz)[ind_gj1],weight);
	 }
	 if((*genjet_nsubjettiness2)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[34]->Fill((*genjet_nsubjettiness2)[ind_gj2]/(*genjet_nsubjettiness1)[ind_gj2],weight);
	 }
	 if((*genjet_nsubjettiness2)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[37]->Fill((*genjet_nsubjettiness2_lrz)[ind_gj2]/(*genjet_nsubjettiness1_lrz)[ind_gj2],weight);
	 }
	 if((*genjet_nsubjettiness3)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[40]->Fill((*genjet_nsubjettiness3)[ind_gj2]/(*genjet_nsubjettiness2)[ind_gj2],weight);
	 }
	 if((*genjet_nsubjettiness3)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[43]->Fill((*genjet_nsubjettiness3_lrz)[ind_gj2]/(*genjet_nsubjettiness2_lrz)[ind_gj2],weight);
	 }
	 h_hist_vec_gen[46]->Fill((*genjet_beta1_N2)[ind_gj1],weight);
	 h_hist_vec_gen[49]->Fill((*genjet_beta1_N2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[52]->Fill((*genjet_beta2_N2)[ind_gj1],weight);
	 h_hist_vec_gen[55]->Fill((*genjet_beta2_N2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[58]->Fill((*genjet_beta0_5_N2)[ind_gj1],weight);
	 h_hist_vec_gen[61]->Fill((*genjet_beta0_5_N2_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[64]->Fill((*genjet_beta1_N2)[ind_gj2],weight);
	 h_hist_vec_gen[67]->Fill((*genjet_beta1_N2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[70]->Fill((*genjet_beta2_N2)[ind_gj2],weight);
	 h_hist_vec_gen[73]->Fill((*genjet_beta2_N2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[76]->Fill((*genjet_beta0_5_N2)[ind_gj2],weight);
	 h_hist_vec_gen[79]->Fill((*genjet_beta0_5_N2_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[82]->Fill((*genjet_beta1_N3)[ind_gj1],weight);
	 h_hist_vec_gen[85]->Fill((*genjet_beta1_N3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[88]->Fill((*genjet_beta2_N3)[ind_gj1],weight);
	 h_hist_vec_gen[91]->Fill((*genjet_beta2_N3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[94]->Fill((*genjet_beta0_5_N3)[ind_gj1],weight);
	 h_hist_vec_gen[97]->Fill((*genjet_beta0_5_N3_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[100]->Fill((*genjet_beta1_N3)[ind_gj2],weight);
	 h_hist_vec_gen[103]->Fill((*genjet_beta1_N3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[106]->Fill((*genjet_beta2_N3)[ind_gj2],weight);
	 h_hist_vec_gen[109]->Fill((*genjet_beta2_N3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[112]->Fill((*genjet_beta0_5_N3)[ind_gj2],weight);
	 h_hist_vec_gen[115]->Fill((*genjet_beta0_5_N3_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[118]->Fill((*genjet_beta1_C2)[ind_gj1],weight);
	 h_hist_vec_gen[121]->Fill((*genjet_beta1_C2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[124]->Fill((*genjet_beta2_C2)[ind_gj1],weight);
	 h_hist_vec_gen[127]->Fill((*genjet_beta2_C2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[130]->Fill((*genjet_beta0_5_C2)[ind_gj1],weight);
	 h_hist_vec_gen[133]->Fill((*genjet_beta0_5_C2_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[136]->Fill((*genjet_beta1_C2)[ind_gj2],weight);
	 h_hist_vec_gen[139]->Fill((*genjet_beta1_C2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[142]->Fill((*genjet_beta2_C2)[ind_gj2],weight);
	 h_hist_vec_gen[145]->Fill((*genjet_beta2_C2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[148]->Fill((*genjet_beta0_5_C2)[ind_gj2],weight);
	 h_hist_vec_gen[151]->Fill((*genjet_beta0_5_C2_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[154]->Fill((*genjet_beta1_C3)[ind_gj1],weight);
	 h_hist_vec_gen[157]->Fill((*genjet_beta1_C3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[160]->Fill((*genjet_beta2_C3)[ind_gj1],weight);
	 h_hist_vec_gen[163]->Fill((*genjet_beta2_C3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[166]->Fill((*genjet_beta0_5_C3)[ind_gj1],weight);
	 h_hist_vec_gen[169]->Fill((*genjet_beta0_5_C3_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[172]->Fill((*genjet_beta1_C3)[ind_gj2],weight);
	 h_hist_vec_gen[175]->Fill((*genjet_beta1_C3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[178]->Fill((*genjet_beta2_C3)[ind_gj2],weight);
	 h_hist_vec_gen[181]->Fill((*genjet_beta2_C3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[184]->Fill((*genjet_beta0_5_C3)[ind_gj2],weight);
	 h_hist_vec_gen[187]->Fill((*genjet_beta0_5_C3_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[190]->Fill((*genjet_beta1_D2)[ind_gj1],weight);
	 h_hist_vec_gen[193]->Fill((*genjet_beta1_D2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[196]->Fill((*genjet_beta2_D2)[ind_gj1],weight);
	 h_hist_vec_gen[199]->Fill((*genjet_beta2_D2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[202]->Fill((*genjet_beta0_5_D2)[ind_gj1],weight);
	 h_hist_vec_gen[205]->Fill((*genjet_beta0_5_D2_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[208]->Fill((*genjet_beta1_D2)[ind_gj2],weight);
	 h_hist_vec_gen[211]->Fill((*genjet_beta1_D2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[214]->Fill((*genjet_beta2_D2)[ind_gj2],weight);
	 h_hist_vec_gen[217]->Fill((*genjet_beta2_D2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[220]->Fill((*genjet_beta0_5_D2)[ind_gj2],weight);
	 h_hist_vec_gen[223]->Fill((*genjet_beta0_5_D2_E_theta)[ind_gj2],weight);

	 //D_{2}^{(1,2)} defined as e_3^{(1)}/(e_2^{(2)})^3/2
	 if((*genjet_beta2_ECorr2)[ind_gj1]!=0){
	   h_hist_vec_gen[226]->Fill((*genjet_beta1_ECorr3)[ind_gj1]/pow((*genjet_beta2_ECorr2)[ind_gj1],3./2.),weight);
	 }else{
	   h_hist_vec_gen[226]->Fill(1.e-5);
	 }
	 if((*genjet_beta2_ECorr2_E_theta)[ind_gj1]!=0){
	   h_hist_vec_gen[229]->Fill((*genjet_beta1_ECorr3_E_theta)[ind_gj1]/pow((*genjet_beta2_ECorr2_E_theta)[ind_gj1],3./2),weight);
	 }else{
	   h_hist_vec_gen[229]->Fill(1.e-5);
	 }
	 if((*genjet_beta2_ECorr2)[ind_gj2]!=0){
	   h_hist_vec_gen[232]->Fill((*genjet_beta1_ECorr3)[ind_gj2]/pow((*genjet_beta2_ECorr2)[ind_gj2],3./2.),weight);
	 }else{
	   h_hist_vec_gen[232]->Fill(1.e-5);
	 }
	 if((*genjet_beta2_ECorr2_E_theta)[ind_gj2]!=0){
	   h_hist_vec_gen[235]->Fill((*genjet_beta1_ECorr3_E_theta)[ind_gj2]/pow((*genjet_beta2_ECorr2_E_theta)[ind_gj2],3./2),weight);
	 }else{
	   h_hist_vec_gen[235]->Fill(1.e-5);
	 }
	 if(ind_sj1_gj1!=-1){
	   h_hist_vec_gen[238]->Fill((*genjet_subjet_E)[ind_sj1_gj1],weight);
	   h_hist_vec_gen[250]->Fill((*genjet_subjet_E)[ind_sj1_gj1]/gj_m1.E(),weight);
	   if(ind_sj2_gj1!=-1){
	     h_hist_vec_gen[256]->Fill(temp_sj1_gj1.Angle(temp_sj2_gj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj1!=-1){
	   h_hist_vec_gen[241]->Fill((*genjet_subjet_E)[ind_sj2_gj1],weight);
	 }
	 if(ind_sj1_gj2!=-1){
	   h_hist_vec_gen[244]->Fill((*genjet_subjet_E)[ind_sj1_gj2],weight);
	   h_hist_vec_gen[253]->Fill((*genjet_subjet_E)[ind_sj1_gj2]/gj_m2.E(),weight);
	   if(ind_sj2_gj2!=-1){
	     h_hist_vec_gen[259]->Fill(temp_sj1_gj2.Angle(temp_sj2_gj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj2!=-1){
	   h_hist_vec_gen[247]->Fill((*genjet_subjet_E)[ind_sj2_gj2],weight);
	 }

       }else{
	 h_hist_vec_gen[2]->Fill(gj_m1.Angle(gj_m2.Vect())*TMath::RadToDeg(),weight);    
	 h_hist_vec_gen[5]->Fill(DeltaPhi(gj_m1.Phi(),gj_m2.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[8]->Fill(fabs(gj_m1.Theta()-gj_m2.Theta())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[11]->Fill(gj_m1.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[14]->Fill(gj_m2.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[17]->Fill(gj_m1.M(),weight);
	 h_hist_vec_gen[20]->Fill(gj_m2.M(),weight);
	 if((*genjet_nsubjettiness2)[ind_gj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[23]->Fill((*genjet_nsubjettiness2)[ind_gj1]/(*genjet_nsubjettiness1)[ind_gj1],weight);
	 }
	 if((*genjet_nsubjettiness2)[ind_gj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[26]->Fill((*genjet_nsubjettiness2_lrz)[ind_gj1]/(*genjet_nsubjettiness1_lrz)[ind_gj1],weight);
	 }
	 if((*genjet_nsubjettiness3)[ind_gj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[29]->Fill((*genjet_nsubjettiness3)[ind_gj1]/(*genjet_nsubjettiness2)[ind_gj1],weight);
	 }
	 if((*genjet_nsubjettiness3)[ind_gj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[32]->Fill((*genjet_nsubjettiness3_lrz)[ind_gj1]/(*genjet_nsubjettiness2_lrz)[ind_gj1],weight);
	 }
	 if((*genjet_nsubjettiness2)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[35]->Fill((*genjet_nsubjettiness2)[ind_gj2]/(*genjet_nsubjettiness1)[ind_gj2],weight);
	 }
	 if((*genjet_nsubjettiness2)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[38]->Fill((*genjet_nsubjettiness2_lrz)[ind_gj2]/(*genjet_nsubjettiness1_lrz)[ind_gj2],weight);
	 }
	 if((*genjet_nsubjettiness3)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[41]->Fill((*genjet_nsubjettiness3)[ind_gj2]/(*genjet_nsubjettiness2)[ind_gj2],weight);
	 }
	 if((*genjet_nsubjettiness3)[ind_gj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_gen[44]->Fill((*genjet_nsubjettiness3_lrz)[ind_gj2]/(*genjet_nsubjettiness2_lrz)[ind_gj2],weight);
	 }
	 h_hist_vec_gen[47]->Fill((*genjet_beta1_N2)[ind_gj1],weight);
	 h_hist_vec_gen[50]->Fill((*genjet_beta1_N2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[53]->Fill((*genjet_beta2_N2)[ind_gj1],weight);
	 h_hist_vec_gen[56]->Fill((*genjet_beta2_N2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[59]->Fill((*genjet_beta0_5_N2)[ind_gj1],weight);
	 h_hist_vec_gen[62]->Fill((*genjet_beta0_5_N2_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[65]->Fill((*genjet_beta1_N2)[ind_gj2],weight);
	 h_hist_vec_gen[68]->Fill((*genjet_beta1_N2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[71]->Fill((*genjet_beta2_N2)[ind_gj2],weight);
	 h_hist_vec_gen[74]->Fill((*genjet_beta2_N2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[77]->Fill((*genjet_beta0_5_N2)[ind_gj2],weight);
	 h_hist_vec_gen[80]->Fill((*genjet_beta0_5_N2_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[83]->Fill((*genjet_beta1_N3)[ind_gj1],weight);
	 h_hist_vec_gen[86]->Fill((*genjet_beta1_N3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[89]->Fill((*genjet_beta2_N3)[ind_gj1],weight);
	 h_hist_vec_gen[92]->Fill((*genjet_beta2_N3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[95]->Fill((*genjet_beta0_5_N3)[ind_gj1],weight);
	 h_hist_vec_gen[98]->Fill((*genjet_beta0_5_N3_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[101]->Fill((*genjet_beta1_N3)[ind_gj2],weight);
	 h_hist_vec_gen[104]->Fill((*genjet_beta1_N3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[107]->Fill((*genjet_beta2_N3)[ind_gj2],weight);
	 h_hist_vec_gen[110]->Fill((*genjet_beta2_N3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[113]->Fill((*genjet_beta0_5_N3)[ind_gj2],weight);
	 h_hist_vec_gen[116]->Fill((*genjet_beta0_5_N3_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[119]->Fill((*genjet_beta1_C2)[ind_gj1],weight);
	 h_hist_vec_gen[122]->Fill((*genjet_beta1_C2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[125]->Fill((*genjet_beta2_C2)[ind_gj1],weight);
	 h_hist_vec_gen[128]->Fill((*genjet_beta2_C2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[131]->Fill((*genjet_beta0_5_C2)[ind_gj1],weight);
	 h_hist_vec_gen[134]->Fill((*genjet_beta0_5_C2_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[137]->Fill((*genjet_beta1_C2)[ind_gj2],weight);
	 h_hist_vec_gen[140]->Fill((*genjet_beta1_C2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[143]->Fill((*genjet_beta2_C2)[ind_gj2],weight);
	 h_hist_vec_gen[146]->Fill((*genjet_beta2_C2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[149]->Fill((*genjet_beta0_5_C2)[ind_gj2],weight);
	 h_hist_vec_gen[152]->Fill((*genjet_beta0_5_C2_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[155]->Fill((*genjet_beta1_C3)[ind_gj1],weight);
	 h_hist_vec_gen[158]->Fill((*genjet_beta1_C3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[161]->Fill((*genjet_beta2_C3)[ind_gj1],weight);
	 h_hist_vec_gen[164]->Fill((*genjet_beta2_C3_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[167]->Fill((*genjet_beta0_5_C3)[ind_gj1],weight);
	 h_hist_vec_gen[170]->Fill((*genjet_beta0_5_C3_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[173]->Fill((*genjet_beta1_C3)[ind_gj2],weight);
	 h_hist_vec_gen[176]->Fill((*genjet_beta1_C3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[179]->Fill((*genjet_beta2_C3)[ind_gj2],weight);
	 h_hist_vec_gen[182]->Fill((*genjet_beta2_C3_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[185]->Fill((*genjet_beta0_5_C3)[ind_gj2],weight);
	 h_hist_vec_gen[188]->Fill((*genjet_beta0_5_C3_E_theta)[ind_gj2],weight);

	 h_hist_vec_gen[191]->Fill((*genjet_beta1_D2)[ind_gj1],weight);
	 h_hist_vec_gen[194]->Fill((*genjet_beta1_D2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[197]->Fill((*genjet_beta2_D2)[ind_gj1],weight);
	 h_hist_vec_gen[200]->Fill((*genjet_beta2_D2_E_theta)[ind_gj1],weight);
	 h_hist_vec_gen[203]->Fill((*genjet_beta0_5_D2)[ind_gj1],weight);
	 h_hist_vec_gen[206]->Fill((*genjet_beta0_5_D2_E_theta)[ind_gj1],weight);

	 h_hist_vec_gen[209]->Fill((*genjet_beta1_D2)[ind_gj2],weight);
	 h_hist_vec_gen[212]->Fill((*genjet_beta1_D2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[215]->Fill((*genjet_beta2_D2)[ind_gj2],weight);
	 h_hist_vec_gen[218]->Fill((*genjet_beta2_D2_E_theta)[ind_gj2],weight);
	 h_hist_vec_gen[221]->Fill((*genjet_beta0_5_D2)[ind_gj2],weight);
	 h_hist_vec_gen[224]->Fill((*genjet_beta0_5_D2_E_theta)[ind_gj2],weight);

	 //D_{2}^{(1,2)} defined as e_3^{(1)}/(e_2^{(2)})^3/2
	 if((*genjet_beta2_ECorr2)[ind_gj1]!=0){
	   h_hist_vec_gen[227]->Fill((*genjet_beta1_ECorr3)[ind_gj1]/pow((*genjet_beta2_ECorr2)[ind_gj1],3./2.),weight);
	 }else{
	   h_hist_vec_gen[227]->Fill(1.e-5);
	 }
	 if((*genjet_beta2_ECorr2_E_theta)[ind_gj1]!=0){
	   h_hist_vec_gen[230]->Fill((*genjet_beta1_ECorr3_E_theta)[ind_gj1]/pow((*genjet_beta2_ECorr2_E_theta)[ind_gj1],3./2),weight);
	 }else{
	   h_hist_vec_gen[230]->Fill(1.e-5);
	 }
	 if((*genjet_beta2_ECorr2)[ind_gj2]!=0){
	   h_hist_vec_gen[233]->Fill((*genjet_beta1_ECorr3)[ind_gj2]/pow((*genjet_beta2_ECorr2)[ind_gj2],3./2.),weight);
	 }else{
	   h_hist_vec_gen[233]->Fill(1.e-5);
	 }
	 if((*genjet_beta2_ECorr2_E_theta)[ind_gj2]!=0){
	   h_hist_vec_gen[236]->Fill((*genjet_beta1_ECorr3_E_theta)[ind_gj2]/pow((*genjet_beta2_ECorr2_E_theta)[ind_gj2],3./2),weight);
	 }else{
	   h_hist_vec_gen[236]->Fill(1.e-5);
	 }
	 if(ind_sj1_gj1!=-1){
	   h_hist_vec_gen[239]->Fill((*genjet_subjet_E)[ind_sj1_gj1],weight);
	   h_hist_vec_gen[251]->Fill((*genjet_subjet_E)[ind_sj1_gj1]/gj_m1.E(),weight);
	   if(ind_sj2_gj1!=-1){
	     h_hist_vec_gen[257]->Fill(temp_sj1_gj1.Angle(temp_sj2_gj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj1!=-1){
	   h_hist_vec_gen[242]->Fill((*genjet_subjet_E)[ind_sj2_gj1],weight);
	 }
	 if(ind_sj1_gj2!=-1){
	   h_hist_vec_gen[245]->Fill((*genjet_subjet_E)[ind_sj1_gj2],weight);
	   h_hist_vec_gen[254]->Fill((*genjet_subjet_E)[ind_sj1_gj2]/gj_m2.E(),weight);
	   if(ind_sj2_gj2!=-1){
	     h_hist_vec_gen[260]->Fill(temp_sj1_gj2.Angle(temp_sj2_gj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj2!=-1){
	   h_hist_vec_gen[248]->Fill((*genjet_subjet_E)[ind_sj2_gj2],weight);
	 }
       }
       h_hist_vec_gen[261]->Fill((tempTotGenP4-tempGenIsoPhP4).M());
       h_hist_vec_gen[262]->Fill((tempTotGenP4-tempGenIsoPhP4+tempTotInvGenP4).M());
       h_hist_vec_gen[263]->Fill(tempTotGenP4.M());
       h_hist_vec_gen[264]->Fill((tempTotGenP4+tempTotInvGenP4).M());
       h_hist_vec_gen[265]->Fill((tempTotGenP4-tempGenTruePhP4+tempTotInvGenP4).M());
     }//bracket closed on two genjets plus isolated lepton cut
    }//two genjet loop AGAIN 
    
    TLorentzVector rj_m1(0,0,0,0);
    TLorentzVector rj_m2(0,0,0,0);
    unsigned int ind_rj1=0;
    unsigned int ind_rj2=1;
    if(recojet_E->size()==2){
      rj_m1.SetPxPyPzE((*recojet_Px)[0],(*recojet_Py)[0],(*recojet_Pz)[0],(*recojet_E)[0]);
      rj_m2.SetPxPyPzE((*recojet_Px)[1],(*recojet_Py)[1],(*recojet_Pz)[1],(*recojet_E)[1]);
      if(rj_m2.M()>rj_m1.M()){
	ind_rj1=1;
	ind_rj2=0;
	TLorentzVector temp=rj_m1;
	rj_m1=rj_m2;
	rj_m2=temp;
      }
      if(rj_m2.M()>=rj_m1.M()){
	std::cout<<"rj mass order should have not been the case nowadays "<<rj_m2.M()<<"/"<<rj_m1.M()<<std::endl;
      }
      if(usePartonInfo && fill_partonInfo){
	if(tempTotEventP4.M()>sqrtS_high){
	  h_hist_parton[40]->Fill(rj_m1.M(),weight);
	  h_hist_parton[41]->Fill(rj_m1.M(),weight);
	  h_hist_parton[46]->Fill(DeltaPhi(rj_m1.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	  h_hist_parton[47]->Fill(fabs(rj_m1.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	  h_hist_parton[48]->Fill(DeltaPhi(rj_m1.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	  h_hist_parton[49]->Fill(fabs(rj_m1.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	  h_hist_parton[52]->Fill(rj_m1.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	  h_hist_parton[53]->Fill(rj_m2.Angle(tempZP4.Vect())*TMath::RadToDeg(),weight);
	}
      }
      int ind_sj1_rj1=-1;
      float E_sj1_rj1=-1;//jets forced into two subjets, thus only one check necessary
      int ind_sj2_rj1=-1;
      int ind_sj1_rj2=-1;
      float E_sj1_rj2=-1;
      int ind_sj2_rj2=-1;
      
      TLorentzVector temp_sj1_rj1(0,0,0,0);
      TLorentzVector temp_sj2_rj1(0,0,0,0);
      TLorentzVector temp_sj1_rj2(0,0,0,0);
      TLorentzVector temp_sj2_rj2(0,0,0,0);
      
      //if too few components in jet, then jet index NOT found -->i.e. remains at -1
      for(unsigned int i=0;i<recojet_subjet_E->size();i++){
	if((*recojet_subjet_jetindex)[i]==ind_rj1){
	  if((*recojet_subjet_E)[i]>E_sj1_rj1){
	    ind_sj2_rj1=ind_sj1_rj1;
	    ind_sj1_rj1=i;
	  }else{
	    ind_sj2_rj1=i;
	}
	}
	if((*recojet_subjet_jetindex)[i]==ind_rj2){
	  if((*recojet_subjet_E)[i]>E_sj1_rj2){
	    ind_sj2_rj2=ind_sj1_rj2;
	    ind_sj1_rj2=i;
	  }else{
	    ind_sj2_rj2=i;
	  }
	}
      }
      if(ind_sj1_rj1!=-1){
	temp_sj1_rj1.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj1],(*recojet_subjet_Py)[ind_sj1_rj1],(*recojet_subjet_Pz)[ind_sj1_rj1],(*recojet_subjet_E)[ind_sj1_rj1]);
      }
      if(ind_sj2_rj1!=-1){
	temp_sj2_rj1.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj1],(*recojet_subjet_Py)[ind_sj2_rj1],(*recojet_subjet_Pz)[ind_sj2_rj1],(*recojet_subjet_E)[ind_sj2_rj1]);
      }
      if(ind_sj1_rj2!=-1){
	temp_sj1_rj2.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj2],(*recojet_subjet_Py)[ind_sj1_rj2],(*recojet_subjet_Pz)[ind_sj1_rj2],(*recojet_subjet_E)[ind_sj1_rj2]);
      }
      if(ind_sj2_rj2!=-1){
	temp_sj2_rj2.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj2],(*recojet_subjet_Py)[ind_sj2_rj2],(*recojet_subjet_Pz)[ind_sj2_rj2],(*recojet_subjet_E)[ind_sj2_rj2]);
      }
      if(recojet_E->size()==2 && n_IsoLep_reco<m_cut_nLeptons){//effectively no cut here, but exactly two recojets
       if((tempTotRecoP4-tempRecoIsoPhP4).M()<sqrtS_low){
	 h_hist_vec_reco[0]->Fill(rj_m1.Angle(rj_m2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[3]->Fill(DeltaPhi(rj_m1.Phi(),rj_m2.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[6]->Fill(fabs(rj_m1.Theta()-rj_m2.Theta())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[9]->Fill(rj_m1.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[12]->Fill(rj_m2.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[15]->Fill(rj_m1.M(),weight);
	 h_hist_vec_reco[18]->Fill(rj_m2.M(),weight);
	 if((*recojet_nsubjettiness2)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[21]->Fill((*recojet_nsubjettiness2)[ind_rj1]/(*recojet_nsubjettiness1)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness2)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[24]->Fill((*recojet_nsubjettiness2_lrz)[ind_rj1]/(*recojet_nsubjettiness1_lrz)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[27]->Fill((*recojet_nsubjettiness3)[ind_rj1]/(*recojet_nsubjettiness2)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[30]->Fill((*recojet_nsubjettiness3_lrz)[ind_rj1]/(*recojet_nsubjettiness2_lrz)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness2)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[33]->Fill((*recojet_nsubjettiness2)[ind_rj2]/(*recojet_nsubjettiness1)[ind_rj2],weight);
	 }
	 if((*recojet_nsubjettiness2)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[36]->Fill((*recojet_nsubjettiness2_lrz)[ind_rj2]/(*recojet_nsubjettiness1_lrz)[ind_rj2],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[39]->Fill((*recojet_nsubjettiness3)[ind_rj2]/(*recojet_nsubjettiness2)[ind_rj2],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[42]->Fill((*recojet_nsubjettiness3_lrz)[ind_rj2]/(*recojet_nsubjettiness2_lrz)[ind_rj2],weight);
	 }
	 h_hist_vec_reco[45]->Fill((*recojet_beta1_N2)[ind_rj1],weight);
	 h_hist_vec_reco[48]->Fill((*recojet_beta1_N2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[51]->Fill((*recojet_beta2_N2)[ind_rj1],weight);
	 h_hist_vec_reco[54]->Fill((*recojet_beta2_N2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[57]->Fill((*recojet_beta0_5_N2)[ind_rj1],weight);
	 h_hist_vec_reco[60]->Fill((*recojet_beta0_5_N2_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[63]->Fill((*recojet_beta1_N2)[ind_rj2],weight);
	 h_hist_vec_reco[66]->Fill((*recojet_beta1_N2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[69]->Fill((*recojet_beta2_N2)[ind_rj2],weight);
	 h_hist_vec_reco[72]->Fill((*recojet_beta2_N2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[75]->Fill((*recojet_beta0_5_N2)[ind_rj2],weight);
	 h_hist_vec_reco[78]->Fill((*recojet_beta0_5_N2_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[81]->Fill((*recojet_beta1_N3)[ind_rj1],weight);
	 h_hist_vec_reco[84]->Fill((*recojet_beta1_N3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[87]->Fill((*recojet_beta2_N3)[ind_rj1],weight);
	 h_hist_vec_reco[90]->Fill((*recojet_beta2_N3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[93]->Fill((*recojet_beta0_5_N3)[ind_rj1],weight);
	 h_hist_vec_reco[96]->Fill((*recojet_beta0_5_N3_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[99]->Fill((*recojet_beta1_N3)[ind_rj2],weight);
	 h_hist_vec_reco[102]->Fill((*recojet_beta1_N3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[105]->Fill((*recojet_beta2_N3)[ind_rj2],weight);
	 h_hist_vec_reco[108]->Fill((*recojet_beta2_N3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[111]->Fill((*recojet_beta0_5_N3)[ind_rj2],weight);
	 h_hist_vec_reco[114]->Fill((*recojet_beta0_5_N3_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[117]->Fill((*recojet_beta1_C2)[ind_rj1],weight);
	 h_hist_vec_reco[120]->Fill((*recojet_beta1_C2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[123]->Fill((*recojet_beta2_C2)[ind_rj1],weight);
	 h_hist_vec_reco[126]->Fill((*recojet_beta2_C2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[129]->Fill((*recojet_beta0_5_C2)[ind_rj1],weight);
	 h_hist_vec_reco[132]->Fill((*recojet_beta0_5_C2_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[135]->Fill((*recojet_beta1_C2)[ind_rj2],weight);
	 h_hist_vec_reco[138]->Fill((*recojet_beta1_C2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[141]->Fill((*recojet_beta2_C2)[ind_rj2],weight);
	 h_hist_vec_reco[144]->Fill((*recojet_beta2_C2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[147]->Fill((*recojet_beta0_5_C2)[ind_rj2],weight);
	 h_hist_vec_reco[150]->Fill((*recojet_beta0_5_C2_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[153]->Fill((*recojet_beta1_C3)[ind_rj1],weight);
	 h_hist_vec_reco[156]->Fill((*recojet_beta1_C3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[159]->Fill((*recojet_beta2_C3)[ind_rj1],weight);
	 h_hist_vec_reco[162]->Fill((*recojet_beta2_C3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[165]->Fill((*recojet_beta0_5_C3)[ind_rj1],weight);
	 h_hist_vec_reco[168]->Fill((*recojet_beta0_5_C3_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[171]->Fill((*recojet_beta1_C3)[ind_rj2],weight);
	 h_hist_vec_reco[174]->Fill((*recojet_beta1_C3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[177]->Fill((*recojet_beta2_C3)[ind_rj2],weight);
	 h_hist_vec_reco[180]->Fill((*recojet_beta2_C3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[183]->Fill((*recojet_beta0_5_C3)[ind_rj2],weight);
	 h_hist_vec_reco[186]->Fill((*recojet_beta0_5_C3_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[189]->Fill((*recojet_beta1_D2)[ind_rj1],weight);
	 h_hist_vec_reco[192]->Fill((*recojet_beta1_D2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[195]->Fill((*recojet_beta2_D2)[ind_rj1],weight);
	 h_hist_vec_reco[198]->Fill((*recojet_beta2_D2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[201]->Fill((*recojet_beta0_5_D2)[ind_rj1],weight);
	 h_hist_vec_reco[204]->Fill((*recojet_beta0_5_D2_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[207]->Fill((*recojet_beta1_D2)[ind_rj2],weight);
	 h_hist_vec_reco[210]->Fill((*recojet_beta1_D2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[213]->Fill((*recojet_beta2_D2)[ind_rj2],weight);
	 h_hist_vec_reco[216]->Fill((*recojet_beta2_D2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[219]->Fill((*recojet_beta0_5_D2)[ind_rj2],weight);
	 h_hist_vec_reco[222]->Fill((*recojet_beta0_5_D2_E_theta)[ind_rj2],weight);
	 //D_{2}^{(1,2)} defined as e_3^{(1)}/(e_2^{(2)})^3/2
	 if((*recojet_beta2_ECorr2)[ind_rj1]!=0){
	   h_hist_vec_reco[225]->Fill((*recojet_beta1_ECorr3)[ind_rj1]/pow((*recojet_beta2_ECorr2)[ind_rj1],3./2.),weight);
	 }else{
	   h_hist_vec_reco[225]->Fill(1.e-5);
	 }
	 if((*recojet_beta2_ECorr2_E_theta)[ind_rj1]!=0){
	   h_hist_vec_reco[228]->Fill((*recojet_beta1_ECorr3_E_theta)[ind_rj1]/pow((*recojet_beta2_ECorr2_E_theta)[ind_rj1],3./2),weight);
	 }else{
	   h_hist_vec_reco[228]->Fill(1.e-5);
	 }
	 if((*recojet_beta2_ECorr2)[ind_rj2]!=0){
	   h_hist_vec_reco[231]->Fill((*recojet_beta1_ECorr3)[ind_rj2]/pow((*recojet_beta2_ECorr2)[ind_rj2],3./2.),weight);
	 }else{
	   h_hist_vec_reco[231]->Fill(1.e-5);
	 }
	 if((*recojet_beta2_ECorr2_E_theta)[ind_rj2]!=0){
	   h_hist_vec_reco[234]->Fill((*recojet_beta1_ECorr3_E_theta)[ind_rj2]/pow((*recojet_beta2_ECorr2_E_theta)[ind_rj2],3./2),weight);
	 }else{
	   h_hist_vec_reco[234]->Fill(1.e-5);
	 }
	 if(ind_sj1_rj1!=-1){
	   h_hist_vec_reco[237]->Fill((*recojet_subjet_E)[ind_sj1_rj1],weight);
	   h_hist_vec_reco[249]->Fill((*recojet_subjet_E)[ind_sj1_rj1]/rj_m1.E(),weight);
	   if(ind_sj2_rj1!=-1){
	     h_hist_vec_reco[255]->Fill(temp_sj1_rj1.Angle(temp_sj2_rj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj1!=-1){
	   h_hist_vec_reco[240]->Fill((*recojet_subjet_E)[ind_sj2_rj1],weight);
	 }
	 if(ind_sj1_rj2!=-1){
	   h_hist_vec_reco[243]->Fill((*recojet_subjet_E)[ind_sj1_rj2],weight);
	   h_hist_vec_reco[252]->Fill((*recojet_subjet_E)[ind_sj1_rj2]/rj_m2.E(),weight);
	   if(ind_sj2_rj2!=-1){
	     h_hist_vec_reco[258]->Fill(temp_sj1_rj2.Angle(temp_sj2_rj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj2!=-1){
	   h_hist_vec_reco[246]->Fill((*recojet_subjet_E)[ind_sj2_rj2],weight);
	 }

       }else if((tempTotRecoP4-tempRecoIsoPhP4).M()<sqrtS_high_reco){
	 h_hist_vec_reco[1]->Fill(rj_m1.Angle(rj_m2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[4]->Fill(DeltaPhi(rj_m1.Phi(),rj_m2.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[7]->Fill(fabs(rj_m1.Theta()-rj_m2.Theta())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[10]->Fill(rj_m1.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[13]->Fill(rj_m2.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[16]->Fill(rj_m1.M(),weight);
	 h_hist_vec_reco[19]->Fill(rj_m2.M(),weight);
	 if((*recojet_nsubjettiness2)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[22]->Fill((*recojet_nsubjettiness2)[ind_rj1]/(*recojet_nsubjettiness1)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness2)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[25]->Fill((*recojet_nsubjettiness2_lrz)[ind_rj1]/(*recojet_nsubjettiness1_lrz)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[28]->Fill((*recojet_nsubjettiness3)[ind_rj1]/(*recojet_nsubjettiness2)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[31]->Fill((*recojet_nsubjettiness3_lrz)[ind_rj1]/(*recojet_nsubjettiness2_lrz)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness2)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[34]->Fill((*recojet_nsubjettiness2)[ind_rj2]/(*recojet_nsubjettiness1)[ind_rj2],weight);
	 }
	 if((*recojet_nsubjettiness2)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[37]->Fill((*recojet_nsubjettiness2_lrz)[ind_rj2]/(*recojet_nsubjettiness1_lrz)[ind_rj2],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[40]->Fill((*recojet_nsubjettiness3)[ind_rj2]/(*recojet_nsubjettiness2)[ind_rj2],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[43]->Fill((*recojet_nsubjettiness3_lrz)[ind_rj2]/(*recojet_nsubjettiness2_lrz)[ind_rj2],weight);
	 }
	 h_hist_vec_reco[46]->Fill((*recojet_beta1_N2)[ind_rj1],weight);
	 h_hist_vec_reco[49]->Fill((*recojet_beta1_N2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[52]->Fill((*recojet_beta2_N2)[ind_rj1],weight);
	 h_hist_vec_reco[55]->Fill((*recojet_beta2_N2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[58]->Fill((*recojet_beta0_5_N2)[ind_rj1],weight);
	 h_hist_vec_reco[61]->Fill((*recojet_beta0_5_N2_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[64]->Fill((*recojet_beta1_N2)[ind_rj2],weight);
	 h_hist_vec_reco[67]->Fill((*recojet_beta1_N2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[70]->Fill((*recojet_beta2_N2)[ind_rj2],weight);
	 h_hist_vec_reco[73]->Fill((*recojet_beta2_N2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[76]->Fill((*recojet_beta0_5_N2)[ind_rj2],weight);
	 h_hist_vec_reco[79]->Fill((*recojet_beta0_5_N2_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[82]->Fill((*recojet_beta1_N3)[ind_rj1],weight);
	 h_hist_vec_reco[85]->Fill((*recojet_beta1_N3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[88]->Fill((*recojet_beta2_N3)[ind_rj1],weight);
	 h_hist_vec_reco[91]->Fill((*recojet_beta2_N3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[94]->Fill((*recojet_beta0_5_N3)[ind_rj1],weight);
	 h_hist_vec_reco[97]->Fill((*recojet_beta0_5_N3_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[100]->Fill((*recojet_beta1_N3)[ind_rj2],weight);
	 h_hist_vec_reco[103]->Fill((*recojet_beta1_N3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[106]->Fill((*recojet_beta2_N3)[ind_rj2],weight);
	 h_hist_vec_reco[109]->Fill((*recojet_beta2_N3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[112]->Fill((*recojet_beta0_5_N3)[ind_rj2],weight);
	 h_hist_vec_reco[115]->Fill((*recojet_beta0_5_N3_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[118]->Fill((*recojet_beta1_C2)[ind_rj1],weight);
	 h_hist_vec_reco[121]->Fill((*recojet_beta1_C2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[124]->Fill((*recojet_beta2_C2)[ind_rj1],weight);
	 h_hist_vec_reco[127]->Fill((*recojet_beta2_C2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[130]->Fill((*recojet_beta0_5_C2)[ind_rj1],weight);
	 h_hist_vec_reco[133]->Fill((*recojet_beta0_5_C2_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[136]->Fill((*recojet_beta1_C2)[ind_rj2],weight);
	 h_hist_vec_reco[139]->Fill((*recojet_beta1_C2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[142]->Fill((*recojet_beta2_C2)[ind_rj2],weight);
	 h_hist_vec_reco[145]->Fill((*recojet_beta2_C2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[148]->Fill((*recojet_beta0_5_C2)[ind_rj2],weight);
	 h_hist_vec_reco[151]->Fill((*recojet_beta0_5_C2_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[154]->Fill((*recojet_beta1_C3)[ind_rj1],weight);
	 h_hist_vec_reco[157]->Fill((*recojet_beta1_C3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[160]->Fill((*recojet_beta2_C3)[ind_rj1],weight);
	 h_hist_vec_reco[163]->Fill((*recojet_beta2_C3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[166]->Fill((*recojet_beta0_5_C3)[ind_rj1],weight);
	 h_hist_vec_reco[169]->Fill((*recojet_beta0_5_C3_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[172]->Fill((*recojet_beta1_C3)[ind_rj2],weight);
	 h_hist_vec_reco[175]->Fill((*recojet_beta1_C3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[178]->Fill((*recojet_beta2_C3)[ind_rj2],weight);
	 h_hist_vec_reco[181]->Fill((*recojet_beta2_C3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[184]->Fill((*recojet_beta0_5_C3)[ind_rj2],weight);
	 h_hist_vec_reco[187]->Fill((*recojet_beta0_5_C3_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[190]->Fill((*recojet_beta1_D2)[ind_rj1],weight);
	 h_hist_vec_reco[193]->Fill((*recojet_beta1_D2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[196]->Fill((*recojet_beta2_D2)[ind_rj1],weight);
	 h_hist_vec_reco[199]->Fill((*recojet_beta2_D2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[202]->Fill((*recojet_beta0_5_D2)[ind_rj1],weight);
	 h_hist_vec_reco[205]->Fill((*recojet_beta0_5_D2_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[208]->Fill((*recojet_beta1_D2)[ind_rj2],weight);
	 h_hist_vec_reco[211]->Fill((*recojet_beta1_D2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[214]->Fill((*recojet_beta2_D2)[ind_rj2],weight);
	 h_hist_vec_reco[217]->Fill((*recojet_beta2_D2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[220]->Fill((*recojet_beta0_5_D2)[ind_rj2],weight);
	 h_hist_vec_reco[223]->Fill((*recojet_beta0_5_D2_E_theta)[ind_rj2],weight);

	 //D_{2}^{(1,2)} defined as e_3^{(1)}/(e_2^{(2)})^3/2
	 if((*recojet_beta2_ECorr2)[ind_rj1]!=0){
	   h_hist_vec_reco[226]->Fill((*recojet_beta1_ECorr3)[ind_rj1]/pow((*recojet_beta2_ECorr2)[ind_rj1],3./2.),weight);
	 }else{
	   h_hist_vec_reco[226]->Fill(1.e-5);
	 }
	 if((*recojet_beta2_ECorr2_E_theta)[ind_rj1]!=0){
	   h_hist_vec_reco[229]->Fill((*recojet_beta1_ECorr3_E_theta)[ind_rj1]/pow((*recojet_beta2_ECorr2_E_theta)[ind_rj1],3./2),weight);
	 }else{
	   h_hist_vec_reco[229]->Fill(1.e-5);
	 }
	 if((*recojet_beta2_ECorr2)[ind_rj2]!=0){
	   h_hist_vec_reco[232]->Fill((*recojet_beta1_ECorr3)[ind_rj2]/pow((*recojet_beta2_ECorr2)[ind_rj2],3./2.),weight);
	 }else{
	   h_hist_vec_reco[232]->Fill(1.e-5);
	 }
	 if((*recojet_beta2_ECorr2_E_theta)[ind_rj2]!=0){
	   h_hist_vec_reco[235]->Fill((*recojet_beta1_ECorr3_E_theta)[ind_rj2]/pow((*recojet_beta2_ECorr2_E_theta)[ind_rj2],3./2),weight);
	 }else{
	   h_hist_vec_reco[235]->Fill(1.e-5);
	 }
	 if(ind_sj1_rj1!=-1){
	   h_hist_vec_reco[238]->Fill((*recojet_subjet_E)[ind_sj1_rj1],weight);
	   h_hist_vec_reco[250]->Fill((*recojet_subjet_E)[ind_sj1_rj1]/rj_m1.E(),weight);
	   if(ind_sj2_rj1!=-1){
	     h_hist_vec_reco[256]->Fill(temp_sj1_rj1.Angle(temp_sj2_rj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj1!=-1){
	   h_hist_vec_reco[241]->Fill((*recojet_subjet_E)[ind_sj2_rj1],weight);
	 }
	 if(ind_sj1_rj2!=-1){
	   h_hist_vec_reco[244]->Fill((*recojet_subjet_E)[ind_sj1_rj2],weight);
	   h_hist_vec_reco[253]->Fill((*recojet_subjet_E)[ind_sj1_rj2]/rj_m2.E(),weight);
	   if(ind_sj2_rj2!=-1){
	     h_hist_vec_reco[259]->Fill(temp_sj1_rj2.Angle(temp_sj2_rj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj2!=-1){
	   h_hist_vec_reco[247]->Fill((*recojet_subjet_E)[ind_sj2_rj2],weight);
	 }

       }else{
	 h_hist_vec_reco[2]->Fill(rj_m1.Angle(rj_m2.Vect())*TMath::RadToDeg(),weight);    
	 h_hist_vec_reco[5]->Fill(DeltaPhi(rj_m1.Phi(),rj_m2.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[8]->Fill(fabs(rj_m1.Theta()-rj_m2.Theta())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[11]->Fill(rj_m1.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[14]->Fill(rj_m2.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[17]->Fill(rj_m1.M(),weight);
	 h_hist_vec_reco[20]->Fill(rj_m2.M(),weight);
	 if((*recojet_nsubjettiness2)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[23]->Fill((*recojet_nsubjettiness2)[ind_rj1]/(*recojet_nsubjettiness1)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness2)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[26]->Fill((*recojet_nsubjettiness2_lrz)[ind_rj1]/(*recojet_nsubjettiness1_lrz)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[29]->Fill((*recojet_nsubjettiness3)[ind_rj1]/(*recojet_nsubjettiness2)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj1]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[32]->Fill((*recojet_nsubjettiness3_lrz)[ind_rj1]/(*recojet_nsubjettiness2_lrz)[ind_rj1],weight);
	 }
	 if((*recojet_nsubjettiness2)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[35]->Fill((*recojet_nsubjettiness2)[ind_rj2]/(*recojet_nsubjettiness1)[ind_rj2],weight);
	 }
	 if((*recojet_nsubjettiness2)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[38]->Fill((*recojet_nsubjettiness2_lrz)[ind_rj2]/(*recojet_nsubjettiness1_lrz)[ind_rj2],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[41]->Fill((*recojet_nsubjettiness3)[ind_rj2]/(*recojet_nsubjettiness2)[ind_rj2],weight);
	 }
	 if((*recojet_nsubjettiness3)[ind_rj2]!=-1){//set to -1 if too few constitutents
	   h_hist_vec_reco[44]->Fill((*recojet_nsubjettiness3_lrz)[ind_rj2]/(*recojet_nsubjettiness2_lrz)[ind_rj2],weight);
	 }
	 h_hist_vec_reco[47]->Fill((*recojet_beta1_N2)[ind_rj1],weight);
	 h_hist_vec_reco[50]->Fill((*recojet_beta1_N2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[53]->Fill((*recojet_beta2_N2)[ind_rj1],weight);
	 h_hist_vec_reco[56]->Fill((*recojet_beta2_N2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[59]->Fill((*recojet_beta0_5_N2)[ind_rj1],weight);
	 h_hist_vec_reco[62]->Fill((*recojet_beta0_5_N2_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[65]->Fill((*recojet_beta1_N2)[ind_rj2],weight);
	 h_hist_vec_reco[68]->Fill((*recojet_beta1_N2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[71]->Fill((*recojet_beta2_N2)[ind_rj2],weight);
	 h_hist_vec_reco[74]->Fill((*recojet_beta2_N2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[77]->Fill((*recojet_beta0_5_N2)[ind_rj2],weight);
	 h_hist_vec_reco[80]->Fill((*recojet_beta0_5_N2_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[83]->Fill((*recojet_beta1_N3)[ind_rj1],weight);
	 h_hist_vec_reco[86]->Fill((*recojet_beta1_N3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[89]->Fill((*recojet_beta2_N3)[ind_rj1],weight);
	 h_hist_vec_reco[92]->Fill((*recojet_beta2_N3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[95]->Fill((*recojet_beta0_5_N3)[ind_rj1],weight);
	 h_hist_vec_reco[98]->Fill((*recojet_beta0_5_N3_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[101]->Fill((*recojet_beta1_N3)[ind_rj2],weight);
	 h_hist_vec_reco[104]->Fill((*recojet_beta1_N3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[107]->Fill((*recojet_beta2_N3)[ind_rj2],weight);
	 h_hist_vec_reco[110]->Fill((*recojet_beta2_N3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[113]->Fill((*recojet_beta0_5_N3)[ind_rj2],weight);
	 h_hist_vec_reco[116]->Fill((*recojet_beta0_5_N3_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[119]->Fill((*recojet_beta1_C2)[ind_rj1],weight);
	 h_hist_vec_reco[122]->Fill((*recojet_beta1_C2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[125]->Fill((*recojet_beta2_C2)[ind_rj1],weight);
	 h_hist_vec_reco[128]->Fill((*recojet_beta2_C2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[131]->Fill((*recojet_beta0_5_C2)[ind_rj1],weight);
	 h_hist_vec_reco[134]->Fill((*recojet_beta0_5_C2_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[137]->Fill((*recojet_beta1_C2)[ind_rj2],weight);
	 h_hist_vec_reco[140]->Fill((*recojet_beta1_C2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[143]->Fill((*recojet_beta2_C2)[ind_rj2],weight);
	 h_hist_vec_reco[146]->Fill((*recojet_beta2_C2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[149]->Fill((*recojet_beta0_5_C2)[ind_rj2],weight);
	 h_hist_vec_reco[152]->Fill((*recojet_beta0_5_C2_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[155]->Fill((*recojet_beta1_C3)[ind_rj1],weight);
	 h_hist_vec_reco[158]->Fill((*recojet_beta1_C3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[161]->Fill((*recojet_beta2_C3)[ind_rj1],weight);
	 h_hist_vec_reco[164]->Fill((*recojet_beta2_C3_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[167]->Fill((*recojet_beta0_5_C3)[ind_rj1],weight);
	 h_hist_vec_reco[170]->Fill((*recojet_beta0_5_C3_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[173]->Fill((*recojet_beta1_C3)[ind_rj2],weight);
	 h_hist_vec_reco[176]->Fill((*recojet_beta1_C3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[179]->Fill((*recojet_beta2_C3)[ind_rj2],weight);
	 h_hist_vec_reco[182]->Fill((*recojet_beta2_C3_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[185]->Fill((*recojet_beta0_5_C3)[ind_rj2],weight);
	 h_hist_vec_reco[188]->Fill((*recojet_beta0_5_C3_E_theta)[ind_rj2],weight);

	 h_hist_vec_reco[191]->Fill((*recojet_beta1_D2)[ind_rj1],weight);
	 h_hist_vec_reco[194]->Fill((*recojet_beta1_D2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[197]->Fill((*recojet_beta2_D2)[ind_rj1],weight);
	 h_hist_vec_reco[200]->Fill((*recojet_beta2_D2_E_theta)[ind_rj1],weight);
	 h_hist_vec_reco[203]->Fill((*recojet_beta0_5_D2)[ind_rj1],weight);
	 h_hist_vec_reco[206]->Fill((*recojet_beta0_5_D2_E_theta)[ind_rj1],weight);

	 h_hist_vec_reco[209]->Fill((*recojet_beta1_D2)[ind_rj2],weight);
	 h_hist_vec_reco[212]->Fill((*recojet_beta1_D2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[215]->Fill((*recojet_beta2_D2)[ind_rj2],weight);
	 h_hist_vec_reco[218]->Fill((*recojet_beta2_D2_E_theta)[ind_rj2],weight);
	 h_hist_vec_reco[221]->Fill((*recojet_beta0_5_D2)[ind_rj2],weight);
	 h_hist_vec_reco[224]->Fill((*recojet_beta0_5_D2_E_theta)[ind_rj2],weight);

	 //D_{2}^{(1,2)} defined as e_3^{(1)}/(e_2^{(2)})^3/2
	 if((*recojet_beta2_ECorr2)[ind_rj1]!=0){
	   h_hist_vec_reco[227]->Fill((*recojet_beta1_ECorr3)[ind_rj1]/pow((*recojet_beta2_ECorr2)[ind_rj1],3./2.),weight);
	 }else{
	   h_hist_vec_reco[227]->Fill(1.e-5);
	 }
	 if((*recojet_beta2_ECorr2_E_theta)[ind_rj1]!=0){
	   h_hist_vec_reco[230]->Fill((*recojet_beta1_ECorr3_E_theta)[ind_rj1]/pow((*recojet_beta2_ECorr2_E_theta)[ind_rj1],3./2),weight);
	 }else{
	   h_hist_vec_reco[230]->Fill(1.e-5);
	 }
	 if((*recojet_beta2_ECorr2)[ind_rj2]!=0){
	   h_hist_vec_reco[233]->Fill((*recojet_beta1_ECorr3)[ind_rj2]/pow((*recojet_beta2_ECorr2)[ind_rj2],3./2.),weight);
	 }else{
	   h_hist_vec_reco[233]->Fill(1.e-5);
	 }
	 if((*recojet_beta2_ECorr2_E_theta)[ind_rj2]!=0){
	   h_hist_vec_reco[236]->Fill((*recojet_beta1_ECorr3_E_theta)[ind_rj2]/pow((*recojet_beta2_ECorr2_E_theta)[ind_rj2],3./2),weight);
	 }else{
	   h_hist_vec_reco[236]->Fill(1.e-5);
	 }
	 if(ind_sj1_rj1!=-1){
	   h_hist_vec_reco[239]->Fill((*recojet_subjet_E)[ind_sj1_rj1],weight);
	   h_hist_vec_reco[251]->Fill((*recojet_subjet_E)[ind_sj1_rj1]/rj_m1.E(),weight);
	   if(ind_sj2_rj1!=-1){
	     h_hist_vec_reco[257]->Fill(temp_sj1_rj1.Angle(temp_sj2_rj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj1!=-1){
	   h_hist_vec_reco[242]->Fill((*recojet_subjet_E)[ind_sj2_rj1],weight);
	 }
	 if(ind_sj1_rj2!=-1){
	   h_hist_vec_reco[245]->Fill((*recojet_subjet_E)[ind_sj1_rj2],weight);
	   h_hist_vec_reco[254]->Fill((*recojet_subjet_E)[ind_sj1_rj2]/rj_m2.E(),weight);
	   if(ind_sj2_rj2!=-1){
	     h_hist_vec_reco[260]->Fill(temp_sj1_rj2.Angle(temp_sj2_rj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj2!=-1){
	   h_hist_vec_reco[248]->Fill((*recojet_subjet_E)[ind_sj2_rj2],weight);
	 }
       }
      }//bracket closed on two recojets plus isolated lepton cut
    }//recojets are closed

    if((tempTotRecoP4-tempRecoIsoPhP4).M()<sqrtS_low){
      if(fill_genInfo){
	if((tempTotGenP4-tempGenIsoPhP4).M()<sqrtS_low){
	  h_hist_vec_reco[261]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else if((tempTotGenP4-tempGenIsoPhP4).M()<sqrtS_high_reco) {
	  h_hist_vec_reco[262]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else{
	  h_hist_vec_reco[263]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}
      }
      if(fill_partonInfo){
	if(tempTotEventP4.M()<sqrtS_low){
	  h_hist_vec_reco[264]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else if(tempTotEventP4.M()<sqrtS_high_reco) {
	  h_hist_vec_reco[265]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else{
	  h_hist_vec_reco[266]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}
      }
    }else if((tempTotRecoP4-tempRecoIsoPhP4).M()<sqrtS_high_reco){
      if(fill_genInfo){
	if((tempTotGenP4-tempGenIsoPhP4).M()<sqrtS_low){
	  h_hist_vec_reco[267]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else if((tempTotGenP4-tempGenIsoPhP4).M()<sqrtS_high_reco) {
	  h_hist_vec_reco[268]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else{
	  h_hist_vec_reco[269]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}
      }
      if(fill_partonInfo){
	if(tempTotEventP4.M()<sqrtS_low){
	  h_hist_vec_reco[270]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else if(tempTotEventP4.M()<sqrtS_high_reco) {
	  h_hist_vec_reco[271]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else{
	  h_hist_vec_reco[272]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}
      }
    }else{
      if(fill_genInfo){
	if((tempTotGenP4-tempGenIsoPhP4).M()<sqrtS_low){
	  h_hist_vec_reco[273]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else if((tempTotGenP4-tempGenIsoPhP4).M()<sqrtS_high_reco) {
	  h_hist_vec_reco[274]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else{
	  h_hist_vec_reco[275]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}
      }
      if(fill_partonInfo){
	if(tempTotEventP4.M()<sqrtS_low){
	  h_hist_vec_reco[276]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else if(tempTotEventP4.M()<sqrtS_high_reco) {
	  h_hist_vec_reco[277]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else{
	  h_hist_vec_reco[278]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}
      }
    }
    if(fill_genInfo && fill_partonInfo){
      h_hist_vec_2D[0]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenTruePhP4).M(),weight);
      h_hist_vec_2D[1]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4).M(),weight);
      h_hist_vec_2D[2]->Fill(tempTotEventP4.M(),tempGenJetSum.M(),weight);    
      h_hist_vec_2D[3]->Fill(tempTotEventP4.M(),(tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      h_hist_vec_2D[4]->Fill(tempTotEventP4.M(),tempRecoJetSum.M(),weight);
      h_hist_vec_2D[5]->Fill((tempTotGenP4-tempGenIsoPhP4).M(),(tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      h_hist_vec_2D[6]->Fill(tempGenJetSum.M(),tempRecoJetSum.M(),weight);
      h_hist_vec_2D[7]->Fill(tempTotEventP4.M(),tempTotGenP4.M(),weight);
      h_hist_vec_2D[8]->Fill(tempTotEventP4.M(),tempTotRecoP4.M(),weight);
      h_hist_vec_2D[9]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4+tempTotInvGenP4).M(),weight);
      h_hist_vec_2D[10]->Fill(tempTotEventP4.M(),(tempTotGenP4+tempTotInvGenP4).M(),weight);
    }
    }//loop over tree
  }

  


void HZAnalyzerFull(){

  CLICdpStyle();

  gROOT->ProcessLine("#include <vector>");

  //const char* final_histo_name="/afs/cern.ch/user/w/weberma2/performanceHistoFiles/HZqq_testanalyzer_10906_10905_antibbarOnly_DR12_DR12.root";
  const char* final_histo_name="/afs/cern.ch/user/w/weberma2/performanceHistoFiles/HZqq_testanalyzer_12075_12074_SignalOnly_DR7.root";
  
  TFile* file_CLIC_HZqq = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/181011_gcc62_CT/3TeVBackground/DR7/HZStudy_hz_qq_12075_3TeVOverlay_CLIC_o3_v14_DR7_betaVar.root");
  TFile* file_CLIC_HZqq_noBGG = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/181011_gcc62_CT/noBackground/DR7/HZStudy_hz_qq_12074_noOverlay_CLIC_o3_v14_DR7_betaVar.root");

  double xsec_qqqq=548.5;//+-2.02 fb, --lumi 0.1823fb-1)
  double xsec_ee_qq=2849.5;//+-48.4 --lumi 0.03509 fb-1)
  double xsec_hz_qq=3.670;//+-0.0035 fb, --lumi 27.25 fb-1
  double xsec_tt=52.569;//+-0.031, lumi 0.9511 fb-1)

  bool usePartonInfo = true;
  bool fillPartonInfo = true;
  bool fillGenInfo = true;


  TFile* file_histogram=new TFile(final_histo_name,"recreate");

  int n_bins_high=200;

  double lim_energy_low=-200;
  double lim_energy_high=3500.;

  double lim_mass_low=0;
  double lim_mass_high=500;

  TH1F* h_sqrtS_e1_e2_effective = new TH1F("h_sqrtS_e1_e2_effective","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_H_pt_sqrt_s_0_750 = new TH1F("h_H_pt_sqrt_s_0_750","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_H_pt_sqrt_s_750_2500 = new TH1F("h_H_pt_sqrt_s_750_2500","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_H_pt_sqrt_s_2500_3000 = new TH1F("h_H_pt_sqrt_s_2500_3000","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_Z_pt_sqrt_s_0_750 = new TH1F("h_Z_pt_sqrt_s_0_750","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_Z_pt_sqrt_s_750_2500 = new TH1F("h_Z_pt_sqrt_s_750_2500","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_Z_pt_sqrt_s_2500_3000 = new TH1F("h_Z_pt_sqrt_s_2500_3000","", n_bins_high, lim_energy_low,lim_energy_high);

  double lim_dalpha_low=0.;
  double lim_dalpha_high=180.;

  TH1F* h_H_dalpha_bb_sqrt_s_0_750 = new TH1F("h_H_dalpha_bb_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_H_dalpha_bb_sqrt_s_750_2500 = new TH1F("h_H_dalpha_bb_sqrt_s_750_2500","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_H_dalpha_bb_sqrt_s_2500_3000 = new TH1F("h_H_dalpha_bb_sqrt_s_0_2500_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_Z_dalpha_qq_sqrt_s_0_750 = new TH1F("h_Z_dalpha_qq_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_Z_dalpha_qq_sqrt_s_750_2500 = new TH1F("h_Z_dalpha_qq_sqrt_s_750_2500","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_Z_dalpha_qq_sqrt_s_2500_3000 = new TH1F("h_Z_dalpha_qq_sqrt_s_0_2500_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  double lim_theta_low=0.;
  double lim_theta_high=180.;

  TH1F* h_H_theta_b_fw_sqrt_s_0_750 = new TH1F("h_H_theta_b_fw_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_H_theta_b_fw_sqrt_s_750_2500 = new TH1F("h_H_theta_b_fw_sqrt_s_750_2500","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_H_theta_b_fw_sqrt_s_2500_3000 = new TH1F("h_H_theta_b_fw_sqrt_s_0_2500_3000","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_q_fw_sqrt_s_0_750 = new TH1F("h_Z_theta_q_fw_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_q_fw_sqrt_s_750_2500 = new TH1F("h_Z_theta_q_fw_sqrt_s_750_2500","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_q_fw_sqrt_s_2500_3000 = new TH1F("h_Z_theta_q_fw_sqrt_s_0_2500_3000","", n_bins_high, lim_theta_low,lim_theta_high);

  TH1F* h_dalpha_min_Z_q_H_b_sqrt_s_0_750 = new TH1F("h_dalpha_min_Z_q_H_b_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dalpha_min_Z_q_H_b_sqrt_s_750_2500 = new TH1F("h_dalpha_min_Z_q_H_b_sqrt_s_750_2500","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dalpha_min_Z_q_H_b_sqrt_s_2500_3000 = new TH1F("h_dalpha_min_Z_q_H_b_sqrt_s_0_2500_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_H_theta_sqrt_s_0_750 = new TH1F("h_H_theta_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_H_theta_sqrt_s_750_2500 = new TH1F("h_H_theta_sqrt_s_750_2500","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_H_theta_sqrt_s_2500_3000 = new TH1F("h_H_theta_sqrt_s_0_2500_3000","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_sqrt_s_0_750 = new TH1F("h_Z_theta_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_sqrt_s_750_2500 = new TH1F("h_Z_theta_sqrt_s_750_2500","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_sqrt_s_2500_3000 = new TH1F("h_Z_theta_sqrt_s_0_2500_3000","", n_bins_high, lim_theta_low,lim_theta_high);

  TH1F* h_dtheta_H_Z_sqrt_s_0_750 = new TH1F("h_dtheta_H_Z_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dtheta_H_Z_sqrt_s_750_2500 = new TH1F("h_dtheta_H_Z_sqrt_s_750_2500","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dtheta_H_Z_sqrt_s_2500_3000 = new TH1F("h_dtheta_H_Z_sqrt_s_0_2500_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_dphi_H_Z_sqrt_s_0_750 = new TH1F("h_dphi_H_Z_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dphi_H_Z_sqrt_s_750_2500 = new TH1F("h_dphi_H_Z_sqrt_s_750_2500","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dphi_H_Z_sqrt_s_2500_3000 = new TH1F("h_dphi_H_Z_sqrt_s_0_2500_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_dalpha_H_Z_sqrt_s_0_750 = new TH1F("h_dalpha_H_Z_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dalpha_H_Z_sqrt_s_750_2500 = new TH1F("h_dalpha_H_Z_sqrt_s_750_2500","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dalpha_H_Z_sqrt_s_2500_3000 = new TH1F("h_dalpha_H_Z_sqrt_s_0_2500_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
 

  TH1F* h_sqrtS_HZ_effective = new TH1F("h_sqrtS_HZ_effective","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_mass1_gj_sqrtS_2500 = new TH1F("h_mass1_gj_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_mass2_gj_sqrtS_2500 = new TH1F("h_mass2_gj_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);

  TH1F* h_mass1_rj_sqrtS_2500 = new TH1F("h_mass1_rj_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_mass2_rj_sqrtS_2500 = new TH1F("h_mass2_rj_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);

  double lim_phi_low=-5.;
  double lim_phi_high=185.;

  TH1F* h_dPhi_mass1_gj_H_sqrtS_2500 = new TH1F("h_dPhi_mass1_gj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass1_gj_H_sqrtS_2500 = new TH1F("h_dTheta_mass1_gj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_mass2_gj_H_sqrtS_2500 = new TH1F("h_dPhi_mass2_gj_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass2_gj_H_sqrtS_2500 = new TH1F("h_dTheta_mass2_gj_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);


  TH1F* h_dPhi_mass1_rj_H_sqrtS_2500 = new TH1F("h_dPhi_mass1_rj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass1_rj_H_sqrtS_2500 = new TH1F("h_dTheta_mass1_rj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_mass2_rj_H_sqrtS_2500 = new TH1F("h_dPhi_mass2_rj_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass2_rj_H_sqrtS_2500 = new TH1F("h_dTheta_mass2_rj_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);


  TH1F* h_dAlpha_mass1_gj_H_sqrtS_2500 = new TH1F("h_dAlpha_mass1_gj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dAlpha_mass2_gj_Z_sqrtS_2500 = new TH1F("h_dAlpha_mass2_gj_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_dAlpha_mass1_rj_H_sqrtS_2500 = new TH1F("h_dAlpha_mass1_rj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dAlpha_mass2_rj_Z_sqrtS_2500 = new TH1F("h_dAlpha_mass2_rj_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_dPhi_MET_mass1_gj_H_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_gj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass1_rj_H_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_rj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  double lim_MET_low=0.;
  double lim_MET_high=250.;
  TH1F* h_genMET_sqrtS_2500 = new TH1F("h_genMET_sqrtS_2500","", n_bins_high, lim_MET_low,lim_MET_high);
  TH1F* h_recoMET_sqrtS_2500 = new TH1F("h_recoMET_sqrtS_2500","", n_bins_high, lim_MET_low,lim_MET_high);

  TH1F* h_dPhi_MET50_mass1_gj_H_sqrtS_2500 = new TH1F("h_dPhi_MET50_mass1_gj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET50_mass1_rj_H_sqrtS_2500 = new TH1F("h_dPhi_MET50_mass1_rj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_dPhi_MET100_mass1_gj_H_sqrtS_2500 = new TH1F("h_dPhi_MET100_mass1_gj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET100_mass1_rj_H_sqrtS_2500 = new TH1F("h_dPhi_MET100_mass1_rj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_dPhi_MET150_mass1_gj_H_sqrtS_2500 = new TH1F("h_dPhi_MET150_mass1_gj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET150_mass1_rj_H_sqrtS_2500 = new TH1F("h_dPhi_MET150_mass1_rj_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  double lim_sj1ratio_low=0.5;
  double lim_sj1ratio_high=1.;

  double lim_sj2ratio_low=0;
  double lim_sj2ratio_high=0.5;

  TH1F* h_E_q1_over_Z = new TH1F("h_E_q1_over_Z","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_q2_over_Z = new TH1F("h_E_q2_over_Z","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_E_b1_over_H = new TH1F("h_E_b1_over_H","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_b2_over_H = new TH1F("h_E_b2_over_H","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);

  
  std::vector<TH1F*> hist_vec_HZ_parton;
  hist_vec_HZ_parton.push_back(h_sqrtS_e1_e2_effective);
  hist_vec_HZ_parton.push_back(h_H_pt_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_pt_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_H_pt_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_Z_pt_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_Z_pt_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_Z_pt_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_H_dalpha_bb_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_dalpha_bb_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_H_dalpha_bb_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_Z_dalpha_qq_sqrt_s_0_750);//10
  hist_vec_HZ_parton.push_back(h_Z_dalpha_qq_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_Z_dalpha_qq_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_H_theta_b_fw_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_theta_b_fw_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_H_theta_b_fw_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_Z_theta_q_fw_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_Z_theta_q_fw_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_Z_theta_q_fw_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_dalpha_min_Z_q_H_b_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_dalpha_min_Z_q_H_b_sqrt_s_750_2500);//20
  hist_vec_HZ_parton.push_back(h_dalpha_min_Z_q_H_b_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_H_theta_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_theta_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_H_theta_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_Z_theta_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_Z_theta_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_Z_theta_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_dtheta_H_Z_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_dtheta_H_Z_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_dtheta_H_Z_sqrt_s_2500_3000);//30
  hist_vec_HZ_parton.push_back(h_dphi_H_Z_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_dphi_H_Z_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_dphi_H_Z_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_dalpha_H_Z_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_dalpha_H_Z_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_dalpha_H_Z_sqrt_s_2500_3000);
  hist_vec_HZ_parton.push_back(h_sqrtS_HZ_effective);
  hist_vec_HZ_parton.push_back(h_mass1_gj_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_mass2_gj_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_mass1_rj_sqrtS_2500);//40
  hist_vec_HZ_parton.push_back(h_mass2_rj_sqrtS_2500);  
  hist_vec_HZ_parton.push_back(h_dPhi_mass1_gj_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass1_gj_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass2_gj_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass2_gj_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass1_rj_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass1_rj_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass2_rj_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass2_rj_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass1_gj_H_sqrtS_2500);//50
  hist_vec_HZ_parton.push_back(h_dAlpha_mass2_gj_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass1_rj_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass2_rj_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_E_q1_over_Z);
  hist_vec_HZ_parton.push_back(h_E_q2_over_Z);
  hist_vec_HZ_parton.push_back(h_E_b1_over_H);
  hist_vec_HZ_parton.push_back(h_E_b2_over_H);  //57


  for(unsigned int i=0;i<hist_vec_HZ_parton.size();i++){
    hist_vec_HZ_parton[i]->Sumw2();
  }



 float n_bins_high_gen = 100;
 float n_bins_high_gen_mass = 500;

 float n_bins_low_dangle = 00;
 float n_bins_high_dangle = 180;

 TH1F* h_HZ_dAlpha_j1_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_dAlpha_j1_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dAlpha_j1_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_dAlpha_j1_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dAlpha_j1_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_dAlpha_j1_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);

 TH1F* h_HZ_dPhi_j1_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_dPhi_j1_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dPhi_j1_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_dPhi_j1_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dPhi_j1_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_dPhi_j1_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);

 TH1F* h_HZ_dTheta_j1_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_dTheta_j1_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dTheta_j1_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_dTheta_j1_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dTheta_j1_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_dTheta_j1_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);

 TH1F* h_HZ_Theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_Theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_Theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_Theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_Theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_Theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);

 TH1F* h_HZ_Theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_Theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_Theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_Theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_Theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_Theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_dangle, n_bins_high_dangle);

 float n_bins_low_jetmass = 00;
 float n_bins_high_jetmass = 500;

 
 TH1F* h_HZ_mass_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_mass_j1_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j1_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_mass_j1_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 
 TH1F* h_HZ_mass_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_mass_j2_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j2_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_mass_j2_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);

 float n_bins_low_tau_rat = 0.0;
 float n_bins_high_tau_rat = 2.0;

 TH1F* h_HZ_tau21_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_tau21_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_tau21_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_tau21_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);

 TH1F* h_HZ_tau21_lrz_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_tau21_lrz_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_lrz_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_tau21_lrz_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_lrz_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_tau21_lrz_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);

 TH1F* h_HZ_tau32_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_tau32_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_tau32_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_tau32_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);

 TH1F* h_HZ_tau32_lrz_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_tau32_lrz_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_lrz_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_tau32_lrz_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_lrz_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_tau32_lrz_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);

 TH1F* h_HZ_tau21_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_tau21_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_tau21_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_tau21_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);

 TH1F* h_HZ_tau21_lrz_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_tau21_lrz_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_lrz_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_tau21_lrz_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_lrz_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_tau21_lrz_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);

 TH1F* h_HZ_tau32_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_tau32_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_tau32_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_tau32_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);

 TH1F* h_HZ_tau32_lrz_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_tau32_lrz_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_lrz_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_tau32_lrz_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_lrz_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_tau32_lrz_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_tau_rat, n_bins_high_tau_rat);

 float n_bins_low_N2 = 00;
 float n_bins_high_N2 = 0.6;

 TH1F* h_HZ_N2_beta1_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta1_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta1_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta1_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);

 TH1F* h_HZ_N2_beta1_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta1_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta1_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta1_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);

 TH1F* h_HZ_N2_beta2_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta2_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta2_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta2_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);

 TH1F* h_HZ_N2_beta2_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta2_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta2_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta2_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);

 TH1F* h_HZ_N2_beta0_5_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta0_5_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta0_5_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta0_5_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);

 TH1F* h_HZ_N2_beta0_5_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta0_5_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta0_5_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta0_5_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);




 TH1F* h_HZ_N2_beta1_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta1_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta1_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta1_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);

 TH1F* h_HZ_N2_beta1_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta1_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta1_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta1_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);

 TH1F* h_HZ_N2_beta2_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta2_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta2_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta2_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);

 TH1F* h_HZ_N2_beta2_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta2_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta2_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta2_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);

 TH1F* h_HZ_N2_beta0_5_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta0_5_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta0_5_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta0_5_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);

 TH1F* h_HZ_N2_beta0_5_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta0_5_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta0_5_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N2_beta0_5_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N2, n_bins_high_N2);


 float n_bins_low_N3 = 0.0;
 float n_bins_high_N3 = 4.5;

 TH1F* h_HZ_N3_beta1_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta1_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta1_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta1_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);

 TH1F* h_HZ_N3_beta1_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta1_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta1_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta1_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);

 TH1F* h_HZ_N3_beta2_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta2_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta2_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta2_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);

 TH1F* h_HZ_N3_beta2_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta2_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta2_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta2_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);

 TH1F* h_HZ_N3_beta0_5_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta0_5_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta0_5_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta0_5_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);

 TH1F* h_HZ_N3_beta0_5_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta0_5_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta0_5_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta0_5_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);


 TH1F* h_HZ_N3_beta1_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta1_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta1_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta1_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);

 TH1F* h_HZ_N3_beta1_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta1_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta1_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta1_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);

 TH1F* h_HZ_N3_beta2_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta2_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta2_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta2_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);

 TH1F* h_HZ_N3_beta2_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta2_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta2_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta2_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);

 TH1F* h_HZ_N3_beta0_5_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta0_5_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta0_5_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta0_5_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);

 TH1F* h_HZ_N3_beta0_5_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta0_5_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta0_5_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_N3_beta0_5_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_N3, n_bins_high_N3);




 float n_bins_low_C2 = 0.0;
 float n_bins_high_C2 = 0.6;

 TH1F* h_HZ_C2_beta1_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta1_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta1_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta1_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 TH1F* h_HZ_C2_beta1_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta1_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta1_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta1_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 TH1F* h_HZ_C2_beta2_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta2_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta2_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta2_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 TH1F* h_HZ_C2_beta2_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta2_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta2_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta2_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 TH1F* h_HZ_C2_beta0_5_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta0_5_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta0_5_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta0_5_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 TH1F* h_HZ_C2_beta0_5_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta0_5_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta0_5_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta0_5_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);



 TH1F* h_HZ_C2_beta1_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta1_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta1_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta1_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 TH1F* h_HZ_C2_beta1_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta1_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta1_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta1_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 TH1F* h_HZ_C2_beta2_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta2_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta2_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta2_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 TH1F* h_HZ_C2_beta2_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta2_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta2_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta2_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 TH1F* h_HZ_C2_beta0_5_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta0_5_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta0_5_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta0_5_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 TH1F* h_HZ_C2_beta0_5_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta0_5_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta0_5_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C2_beta0_5_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C2, n_bins_high_C2);

 float n_bins_low_C3 = 0.0;
 float n_bins_high_C3 = 1.1;
 
 TH1F* h_HZ_C3_beta1_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta1_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta1_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta1_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 TH1F* h_HZ_C3_beta1_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta1_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta1_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta1_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 TH1F* h_HZ_C3_beta2_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta2_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta2_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta2_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 TH1F* h_HZ_C3_beta2_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta2_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta2_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta2_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 TH1F* h_HZ_C3_beta0_5_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta0_5_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta0_5_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta0_5_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 TH1F* h_HZ_C3_beta0_5_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta0_5_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta0_5_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta0_5_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);


 TH1F* h_HZ_C3_beta1_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta1_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta1_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta1_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 TH1F* h_HZ_C3_beta1_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta1_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta1_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta1_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 TH1F* h_HZ_C3_beta2_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta2_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta2_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta2_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 TH1F* h_HZ_C3_beta2_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta2_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta2_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta2_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 TH1F* h_HZ_C3_beta0_5_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta0_5_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta0_5_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta0_5_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 TH1F* h_HZ_C3_beta0_5_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta0_5_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta0_5_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_C3_beta0_5_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_C3, n_bins_high_C3);

 float n_bins_low_D2 = 0.0;
 float n_bins_high_D2 = 5.5;

 TH1F* h_HZ_D2_beta1_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta1_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta1_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta1_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 TH1F* h_HZ_D2_beta1_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta1_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta1_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta1_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 TH1F* h_HZ_D2_beta2_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta2_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta2_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta2_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 TH1F* h_HZ_D2_beta2_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta2_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta2_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta2_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 TH1F* h_HZ_D2_beta0_5_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta0_5_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta0_5_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta0_5_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 TH1F* h_HZ_D2_beta0_5_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta0_5_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta0_5_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta0_5_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);


 TH1F* h_HZ_D2_beta1_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta1_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta1_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta1_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 TH1F* h_HZ_D2_beta1_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta1_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta1_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta1_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 TH1F* h_HZ_D2_beta2_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta2_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta2_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta2_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 TH1F* h_HZ_D2_beta2_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta2_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta2_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta2_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 TH1F* h_HZ_D2_beta0_5_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta0_5_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta0_5_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta0_5_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 TH1F* h_HZ_D2_beta0_5_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta0_5_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta0_5_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_beta0_5_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2, n_bins_high_D2);

 float n_bins_low_D2_1_2 = 0.0;
 float n_bins_high_D2_1_2 = 0.4;
 //defined in paper as D2^(1,2)=e_3_beta1/(e_2_beta2)^3/2

 TH1F* h_HZ_D2_1_2_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_1_2_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_1_2_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_1_2_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 TH1F* h_HZ_D2_1_2_E_theta_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_1_2_E_theta_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_E_theta_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_1_2_E_theta_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_E_theta_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_1_2_E_theta_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 TH1F* h_HZ_D2_1_2_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_1_2_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_1_2_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_1_2_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 TH1F* h_HZ_D2_1_2_E_theta_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_D2_1_2_E_theta_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_E_theta_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_D2_1_2_E_theta_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_E_theta_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_D2_1_2_E_theta_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 //reco and gen ee_qq
 //r_2 E_3/E_2
 //C1 default=E_2*E_0/E_1^{2}=E_2*1/(sum pt)^2 = E_2/jet_pt^2 --> 0->1
 //N2 0-0.5
 //N3 0-7/160 -220 HZ
 //C2 0-1.5   _E_theta 0.6
 //C3 0-1.6   _E_thata 0.8 
 //D2 0-15/350
 float n_bins_low_subjetE= 10;
 float n_bins_high_subjetE= 1500;

 TH1F* h_HZ_subjet1_E_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_subjet1_E_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet1_E_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_subjet1_E_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet1_E_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_subjet1_E_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);

 TH1F* h_HZ_subjet2_E_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_subjet2_E_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet2_E_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_subjet2_E_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet2_E_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_subjet2_E_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);

 TH1F* h_HZ_subjet1_E_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_subjet1_E_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet1_E_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_subjet1_E_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet1_E_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_subjet1_E_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);

 TH1F* h_HZ_subjet2_E_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_subjet2_E_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet2_E_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_subjet2_E_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet2_E_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_subjet2_E_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_subjetE, n_bins_high_subjetE);

 //higher energetic subjet
 float n_bins_low_subjetE_rat= 0.5;
 float n_bins_high_subjetE_rat= 1.0;

 TH1F* h_HZ_subjet1_E_over_jetE_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_subjet1_E_over_jetE_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_subjet1_E_over_jetE_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_subjet1_E_over_jetE_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_subjet1_E_over_jetE_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_subjet1_E_over_jetE_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);

 TH1F* h_HZ_subjet1_E_over_jetE_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_subjet1_E_over_jetE_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_subjet1_E_over_jetE_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_subjet1_E_over_jetE_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_subjet1_E_over_jetE_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_subjet1_E_over_jetE_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);

 float n_bins_low_dangle_subjet=0;
 float n_bins_high_dangle_subjet=30;


 TH1F* h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_0_750 = new TH1F("h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_750_2500 = new TH1F("h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1 [#circ]");

 TH1F* h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_0_750 = new TH1F("h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_0_750","", n_bins_high_gen,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_750_2500 = new TH1F("h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_750_2500","", n_bins_high_gen,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_2500","", n_bins_high_gen,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");


 TH1F* h_sqrtS_gen_isoPh = new TH1F("h_sqrtS_gen_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);//default
 TH1F* h_sqrtS_gen_isoPh_inv = new TH1F("h_sqrtS_gen_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
 TH1F* h_sqrtS_gen = new TH1F("h_sqrtS_gen","", n_bins_high, lim_energy_low,lim_energy_high);
 TH1F* h_sqrtS_gen_inv = new TH1F("h_sqrtS_gen_inv","", n_bins_high, lim_energy_low,lim_energy_high);
 TH1F* h_sqrtS_gen_truePh_inv = new TH1F("h_sqrtS_gen_truePh_inv","", n_bins_high, lim_energy_low,lim_energy_high);

 h_sqrtS_gen_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_sqrtS_gen_isoPh_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_sqrtS_gen->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_sqrtS_gen_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_sqrtS_gen_truePh_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 std::vector<TH1F*> hist_vec_gen_HZ_1D;  
 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_j1_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_j1_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_j1_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_dPhi_j1_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dPhi_j1_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dPhi_j1_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_dTheta_j1_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dTheta_j1_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dTheta_j1_j2_gen_sqrt_s_2500);//8
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_Theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_Theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_Theta_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_Theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_Theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_Theta_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_gen_sqrt_s_2500);//17
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_lrz_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_lrz_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_lrz_j1_gen_sqrt_s_2500);//26
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_lrz_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_lrz_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_lrz_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_j2_gen_sqrt_s_2500);//35
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_lrz_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_lrz_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau21_lrz_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_lrz_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_lrz_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_tau32_lrz_j2_gen_sqrt_s_2500);//44
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_j1_gen_sqrt_s_2500);//53
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j1_gen_sqrt_s_2500);//62

 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_j2_gen_sqrt_s_2500);//71
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j2_gen_sqrt_s_2500);//80
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_j1_gen_sqrt_s_2500);//89
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j1_gen_sqrt_s_2500);//98
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_j2_gen_sqrt_s_2500);//107
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j2_gen_sqrt_s_2500);//116
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_j1_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_j1_gen_sqrt_s_2500);//125

 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j1_gen_sqrt_s_2500);//134

 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_j2_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_j2_gen_sqrt_s_2500);//143

 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j2_gen_sqrt_s_2500);//152
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_j1_gen_sqrt_s_2500);//161

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j1_gen_sqrt_s_2500);//170

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_j2_gen_sqrt_s_2500);//179

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j2_gen_sqrt_s_2500);//188

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_j1_gen_sqrt_s_2500);//197

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j1_gen_sqrt_s_2500);//206

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_j2_gen_sqrt_s_2500);//215

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j2_gen_sqrt_s_2500);//224

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j1_gen_sqrt_s_2500);


 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_j2_gen_sqrt_s_2500);//233

 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet2_E_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet2_E_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet2_E_j1_gen_sqrt_s_2500);//242

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet2_E_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet2_E_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet2_E_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j1_gen_sqrt_s_2500);//251

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j2_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j1_gen_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j2_gen_sqrt_s_2500);//260

 hist_vec_gen_HZ_1D.push_back(h_sqrtS_gen_isoPh);
 hist_vec_gen_HZ_1D.push_back(h_sqrtS_gen_isoPh_inv);
 hist_vec_gen_HZ_1D.push_back(h_sqrtS_gen);
 hist_vec_gen_HZ_1D.push_back(h_sqrtS_gen_inv);
 hist_vec_gen_HZ_1D.push_back(h_sqrtS_gen_truePh_inv);//265


 for(unsigned int i=0;i<hist_vec_gen_HZ_1D.size();i++){
   hist_vec_gen_HZ_1D[i]->Sumw2();
   hist_vec_gen_HZ_1D[i]->SetLineColor(kBlack);
 }




 int n_bins_high_reco=100;
 int n_bins_high_reco_mass=500;

 TH1F* h_HZ_dAlpha_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_dAlpha_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dAlpha_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_dAlpha_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dAlpha_j1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_dAlpha_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_HZ_dAlpha_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");
 h_HZ_dAlpha_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");
 h_HZ_dAlpha_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");

 TH1F* h_HZ_dPhi_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_dPhi_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dPhi_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_dPhi_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dPhi_j1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_dPhi_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_HZ_dPhi_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");
 h_HZ_dPhi_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");
 h_HZ_dPhi_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");

 TH1F* h_HZ_dTheta_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_dTheta_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dTheta_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_dTheta_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_dTheta_j1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_dTheta_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_HZ_dTheta_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");
 h_HZ_dTheta_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");
 h_HZ_dTheta_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");

 TH1F* h_HZ_Theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_Theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_Theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_Theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_Theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_Theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_HZ_Theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#theta(j1)[#circ]");
 h_HZ_Theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#theta(j1)[#circ]");
 h_HZ_Theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("#theta(j1)[#circ]");

 TH1F* h_HZ_Theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_Theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_Theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_Theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_Theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_Theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_HZ_Theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 h_HZ_Theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 h_HZ_Theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 
 TH1F* h_HZ_mass_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_mass_j1_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j1_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_mass_j1_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 
 h_HZ_mass_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_HZ_mass_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_HZ_mass_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");

 TH1F* h_HZ_mass_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_mass_j2_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j2_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_mass_j2_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);

 h_HZ_mass_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_HZ_mass_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_HZ_mass_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_HZ_tau21_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_tau21_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_tau21_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_tau21_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_tau21_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_HZ_tau21_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_HZ_tau21_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j1)");

 TH1F* h_HZ_tau21_lrz_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_tau21_lrz_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_lrz_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_tau21_lrz_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_lrz_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_tau21_lrz_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_tau21_lrz_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_HZ_tau21_lrz_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_HZ_tau21_lrz_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j1)");

 TH1F* h_HZ_tau32_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_tau32_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_tau32_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_tau32_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_tau32_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_HZ_tau32_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_HZ_tau32_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j1)");

 TH1F* h_HZ_tau32_lrz_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_tau32_lrz_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_lrz_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_tau32_lrz_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_lrz_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_tau32_lrz_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_tau32_lrz_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_HZ_tau32_lrz_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_HZ_tau32_lrz_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j1)");

 TH1F* h_HZ_tau21_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_tau21_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_tau21_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_tau21_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_tau21_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_HZ_tau21_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_HZ_tau21_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j2)");

 TH1F* h_HZ_tau21_lrz_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_tau21_lrz_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_lrz_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_tau21_lrz_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau21_lrz_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_tau21_lrz_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_tau21_lrz_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_HZ_tau21_lrz_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_HZ_tau21_lrz_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j2)");

 TH1F* h_HZ_tau32_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_tau32_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_tau32_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_tau32_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_tau32_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_HZ_tau32_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_HZ_tau32_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j2)");

 TH1F* h_HZ_tau32_lrz_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_tau32_lrz_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_lrz_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_tau32_lrz_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_tau32_lrz_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_tau32_lrz_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_tau32_lrz_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_HZ_tau32_lrz_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_HZ_tau32_lrz_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 //N2 j1
 TH1F* h_HZ_N2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_HZ_N2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_HZ_N2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");

 TH1F* h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");

 TH1F* h_HZ_N2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_HZ_N2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_HZ_N2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");

 TH1F* h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");

 TH1F* h_HZ_N2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_HZ_N2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_HZ_N2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");

 TH1F* h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 //N2 j2
 TH1F* h_HZ_N2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_HZ_N2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_HZ_N2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");

 TH1F* h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");

 TH1F* h_HZ_N2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_HZ_N2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_HZ_N2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");

 TH1F* h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");

 TH1F* h_HZ_N2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_HZ_N2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_HZ_N2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");

 TH1F* h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");



 //N3 j1
 TH1F* h_HZ_N3_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_HZ_N3_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_HZ_N3_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");

 TH1F* h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");

 TH1F* h_HZ_N3_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_HZ_N3_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_HZ_N3_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");

 TH1F* h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");

 TH1F* h_HZ_N3_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_HZ_N3_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_HZ_N3_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");

 TH1F* h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 //N3 j2
 TH1F* h_HZ_N3_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_HZ_N3_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_HZ_N3_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");

 TH1F* h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");

 TH1F* h_HZ_N3_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_HZ_N3_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_HZ_N3_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");

 TH1F* h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");

 TH1F* h_HZ_N3_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_HZ_N3_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_HZ_N3_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");

 TH1F* h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");


 //now the C2 and C3 series
 //C2 j1
 TH1F* h_HZ_C2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_HZ_C2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_HZ_C2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");

 TH1F* h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");

 TH1F* h_HZ_C2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_HZ_C2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_HZ_C2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");

 TH1F* h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");

 TH1F* h_HZ_C2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_HZ_C2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_HZ_C2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");

 TH1F* h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 //C2 j2
 TH1F* h_HZ_C2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_HZ_C2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_HZ_C2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");

 TH1F* h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");

 TH1F* h_HZ_C2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_HZ_C2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_HZ_C2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");

 TH1F* h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");

 TH1F* h_HZ_C2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_HZ_C2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_HZ_C2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");

 TH1F* h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");



 //C3 j1
 TH1F* h_HZ_C3_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_HZ_C3_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_HZ_C3_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");

 TH1F* h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");

 TH1F* h_HZ_C3_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_HZ_C3_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_HZ_C3_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");

 TH1F* h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");

 TH1F* h_HZ_C3_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_HZ_C3_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_HZ_C3_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");

 TH1F* h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 //C3 j2
 TH1F* h_HZ_C3_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_HZ_C3_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_HZ_C3_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");

 TH1F* h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");

 TH1F* h_HZ_C3_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_HZ_C3_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_HZ_C3_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");

 TH1F* h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");

 TH1F* h_HZ_C3_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_HZ_C3_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_HZ_C3_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");

 TH1F* h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");

 //D2 series
 //D2 j1
 TH1F* h_HZ_D2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_HZ_D2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_HZ_D2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");

 TH1F* h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");

 TH1F* h_HZ_D2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_HZ_D2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_HZ_D2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");

 TH1F* h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");

 TH1F* h_HZ_D2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_HZ_D2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_HZ_D2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");

 TH1F* h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 //D2 j2
 TH1F* h_HZ_D2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_HZ_D2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_HZ_D2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");

 TH1F* h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");

 TH1F* h_HZ_D2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_HZ_D2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_HZ_D2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");

 TH1F* h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");

 TH1F* h_HZ_D2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_HZ_D2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_HZ_D2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");

 TH1F* h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");




 //D2 (1,2) series
 //D2 (1,2) j1
 TH1F* h_HZ_D2_1_2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_1_2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_1_2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_1_2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_HZ_D2_1_2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_HZ_D2_1_2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_HZ_D2_1_2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");

 TH1F* h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");

 //D2_1_2 (1,2) j2
 TH1F* h_HZ_D2_1_2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_1_2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_1_2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_1_2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_HZ_D2_1_2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_HZ_D2_1_2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_HZ_D2_1_2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");

 TH1F* h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");



 TH1F* h_HZ_subjet1_E_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_subjet1_E_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet1_E_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_subjet1_E_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet1_E_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_subjet1_E_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_HZ_subjet1_E_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");
 h_HZ_subjet1_E_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");
 h_HZ_subjet1_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");

 TH1F* h_HZ_subjet2_E_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_subjet2_E_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet2_E_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_subjet2_E_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet2_E_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_subjet2_E_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_HZ_subjet2_E_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");
 h_HZ_subjet2_E_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");
 h_HZ_subjet2_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");

 TH1F* h_HZ_subjet1_E_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_subjet1_E_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet1_E_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_subjet1_E_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet1_E_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_subjet1_E_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_HZ_subjet1_E_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");
 h_HZ_subjet1_E_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");
 h_HZ_subjet1_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");

 TH1F* h_HZ_subjet2_E_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_subjet2_E_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet2_E_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_subjet2_E_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_subjet2_E_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_subjet2_E_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_HZ_subjet2_E_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");
 h_HZ_subjet2_E_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");
 h_HZ_subjet2_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");

 TH1F* h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}/E (j1)");
 h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}/E (j1)");
 h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}/E (j1)");

 TH1F* h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}/E (j2)");
 h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}/E (j2)");
 h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}/E (j2)");


 TH1F* h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1 [#circ]");

 TH1F* h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");

 float n_bins_low_sqrt_0_750=0;
 float n_bins_high_sqrt_0_750=750;

 TH1F* h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 TH1F* h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 float n_bins_low_sqrt_750_2500=750;
 float n_bins_high_sqrt_750_2500=2500;

 TH1F* h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 TH1F* h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 float n_bins_low_sqrt_2500=2500;
 float n_bins_high_sqrt_2500=3500;

 TH1F* h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 TH1F* h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

  TH1F* h_HZ_sqrtS_reco_isoPh = new TH1F("h_HZ_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);//default
  TH1F* h_HZ_sqrtS_reco = new TH1F("h_HZ_sqrtS_reco","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_HZ_sqrtS_reco_isoPh_inv = new TH1F("h_HZ_sqrtS_reco_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_sqrtS_reco_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_isoPh_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 std::vector<TH1F*> hist_vec_reco_HZ_1D;  
 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_j1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_dPhi_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dPhi_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dPhi_j1_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_dTheta_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dTheta_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dTheta_j1_j2_reco_sqrt_s_2500);//8
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_Theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_Theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_Theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_Theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_Theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_Theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j1_reco_sqrt_s_2500);//17
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_lrz_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_lrz_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_lrz_j1_reco_sqrt_s_2500);//26
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_lrz_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_lrz_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_lrz_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_j2_reco_sqrt_s_2500);//35
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_lrz_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_lrz_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau21_lrz_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_lrz_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_lrz_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_tau32_lrz_j2_reco_sqrt_s_2500);//44
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_j1_reco_sqrt_s_2500);//53
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//62

 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta1_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_j2_reco_sqrt_s_2500);//71
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta2_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//80
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_j1_reco_sqrt_s_2500);//89
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500);//98
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta1_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_j2_reco_sqrt_s_2500);//107
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta2_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500);//116
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_j1_reco_sqrt_s_2500);//125

 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//134

 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_j2_reco_sqrt_s_2500);//143

 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//152
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_j1_reco_sqrt_s_2500);//161

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500);//170

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_j2_reco_sqrt_s_2500);//179

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500);//188

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_j1_reco_sqrt_s_2500);//197

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//206

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_j2_reco_sqrt_s_2500);//215

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//224

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j1_reco_sqrt_s_2500);


 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_j2_reco_sqrt_s_2500);//233

 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_D2_1_2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet2_E_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet2_E_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet2_E_j1_reco_sqrt_s_2500);//242

 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet2_E_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet2_E_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet2_E_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j1_reco_sqrt_s_2500);//251

 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet1_E_over_jetE_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500);//260

 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500);//269

 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500);//278

 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_isoPh);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_isoPh_inv);//281

 for(unsigned int i=0;i<hist_vec_reco_HZ_1D.size();i++){
   hist_vec_reco_HZ_1D[i]->Sumw2();
   hist_vec_reco_HZ_1D[i]->SetLineColor(kBlack);
 }

 
 //parton vs gen true --> here we calculat all true minus the gentrue photon
 TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_genTrue = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_genTrue","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_genTrue->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_genTrue->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus iso photons
 TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_gen_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
  //parton vs genJet sum
 TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_genJet = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_genJet","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_genJet->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_genJet->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");

 //parton vs reco true --> here we calculat all true minus iso photons
  TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_isoPh = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_isoPh->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_isoPh->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs recoJet sum
  TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_recoJet = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_recoJet","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_recoJet->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_recoJet->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");

 //parton vs reco true --> here we calculat all true minus iso photons
  TH2F* h_HZ_srqtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh = new TH2F("h_HZ_srqtS_gen_isoPh_vs_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_HZ_srqtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
  h_HZ_srqtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs recoJet sum
  TH2F* h_HZ_srqtS_genJet_eff_vs_sqrtS_recoJet = new TH2F("h_HZ_srqtS_genJet_vs_sqrtS_recoJet","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_HZ_srqtS_genJet_eff_vs_sqrtS_recoJet->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
  h_HZ_srqtS_genJet_eff_vs_sqrtS_recoJet->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus the gentrue photon
  TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_genAll = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_genAll","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_genAll->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_genAll->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus the gentrue photon
  TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_recoAll = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_recoAll","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_recoAll->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_recoAll->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus iso photons
  TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_gen_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus iso photons
  TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_inv = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_gen_inv","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_inv->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_inv->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");

  std::vector<TH2F*> hist_vec_HZ_2DHist;  
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_genTrue);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_genJet);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_isoPh);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_recoJet);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_genJet_eff_vs_sqrtS_recoJet);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_genAll);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_recoAll);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_inv);//10
  for(unsigned int i=0;i<hist_vec_HZ_2DHist.size();i++){
    hist_vec_HZ_2DHist[i]->Sumw2();
  }


  fill_HZ_histograms(file_CLIC_HZqq, hist_vec_HZ_parton, hist_vec_gen_HZ_1D, hist_vec_reco_HZ_1D, hist_vec_HZ_2DHist, usePartonInfo ,xsec_hz_qq,fillPartonInfo,fillGenInfo);


  //NOW HZ no BG

 TH1F* h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");
 h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");
 h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");

 TH1F* h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");
 h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");
 h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");

 TH1F* h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");
 h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");
 h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");

 TH1F* h_HZ_noBG_Theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_Theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_noBG_Theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_Theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_noBG_Theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_Theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_HZ_noBG_Theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#theta(j1)[#circ]");
 h_HZ_noBG_Theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#theta(j1)[#circ]");
 h_HZ_noBG_Theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("#theta(j1)[#circ]");

 TH1F* h_HZ_noBG_Theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_Theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_noBG_Theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_Theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_HZ_noBG_Theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_Theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_HZ_noBG_Theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 h_HZ_noBG_Theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 h_HZ_noBG_Theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 
 TH1F* h_HZ_noBG_mass_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_mass_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_noBG_mass_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_mass_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_noBG_mass_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_mass_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_jetmass, n_bins_high_jetmass);
 
 h_HZ_noBG_mass_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_HZ_noBG_mass_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_HZ_noBG_mass_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");

 TH1F* h_HZ_noBG_mass_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_mass_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_noBG_mass_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_mass_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_noBG_mass_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_mass_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_jetmass, n_bins_high_jetmass);

 h_HZ_noBG_mass_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_HZ_noBG_mass_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_HZ_noBG_mass_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_HZ_noBG_tau21_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_tau21_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau21_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_tau21_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau21_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_tau21_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_noBG_tau21_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_HZ_noBG_tau21_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_HZ_noBG_tau21_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j1)");

 TH1F* h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j1)");

 TH1F* h_HZ_noBG_tau32_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_tau32_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau32_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_tau32_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau32_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_tau32_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_noBG_tau32_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_HZ_noBG_tau32_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_HZ_noBG_tau32_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j1)");

 TH1F* h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j1)");

 TH1F* h_HZ_noBG_tau21_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_tau21_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau21_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_tau21_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau21_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_tau21_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_noBG_tau21_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_HZ_noBG_tau21_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_HZ_noBG_tau21_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j2)");

 TH1F* h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j2)");

 TH1F* h_HZ_noBG_tau32_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_tau32_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau32_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_tau32_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau32_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_tau32_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_noBG_tau32_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_HZ_noBG_tau32_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_HZ_noBG_tau32_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j2)");

 TH1F* h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 //N2 j1
 TH1F* h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");

 TH1F* h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");

 TH1F* h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");

 TH1F* h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");

 TH1F* h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");

 TH1F* h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 //N2 j2
 TH1F* h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");

 TH1F* h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");

 TH1F* h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");

 TH1F* h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");

 TH1F* h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");

 TH1F* h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");



 //N3 j1
 TH1F* h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");

 TH1F* h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");

 TH1F* h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");

 TH1F* h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");

 TH1F* h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");

 TH1F* h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 //N3 j2
 TH1F* h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");

 TH1F* h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");

 TH1F* h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");

 TH1F* h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");

 TH1F* h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");

 TH1F* h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");


 //now the C2 and C3 series
 //C2 j1
 TH1F* h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");

 TH1F* h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");

 TH1F* h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");

 TH1F* h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");

 TH1F* h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");

 TH1F* h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 //C2 j2
 TH1F* h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");

 TH1F* h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");

 TH1F* h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");

 TH1F* h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");

 TH1F* h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");

 TH1F* h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");



 //C3 j1
 TH1F* h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");

 TH1F* h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");

 TH1F* h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");

 TH1F* h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");

 TH1F* h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");

 TH1F* h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 //C3 j2
 TH1F* h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");

 TH1F* h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");

 TH1F* h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");

 TH1F* h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");

 TH1F* h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");

 TH1F* h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");

 //D2 series
 //D2 j1
 TH1F* h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");

 TH1F* h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");

 TH1F* h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");

 TH1F* h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");

 TH1F* h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");

 TH1F* h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 //D2 j2
 TH1F* h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");

 TH1F* h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");

 TH1F* h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");

 TH1F* h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");

 TH1F* h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");

 TH1F* h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");




 //D2 (1,2) series
 //D2 (1,2) j1
 TH1F* h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");

 TH1F* h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");

 //D2_1_2 (1,2) j2
 TH1F* h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");

 TH1F* h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");



 TH1F* h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");
 h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");
 h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");

 TH1F* h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");
 h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");
 h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");

 TH1F* h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");
 h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");
 h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");

 TH1F* h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");
 h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");
 h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");

 TH1F* h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}/E (j1)");
 h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}/E (j1)");
 h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}/E (j1)");

 TH1F* h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}/E (j2)");
 h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}/E (j2)");
 h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}/E (j2)");


 TH1F* h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1 [#circ]");

 TH1F* h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750 = new TH1F("h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");


 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500 = new TH1F("h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

  TH1F* h_HZ_noBG_sqrtS_reco_isoPh = new TH1F("h_HZ_noBG_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);//default
  TH1F* h_HZ_noBG_sqrtS_reco = new TH1F("h_HZ_noBG_sqrtS_reco","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_HZ_noBG_sqrtS_reco_isoPh_inv = new TH1F("h_HZ_noBG_sqrtS_reco_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_noBG_sqrtS_reco_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_noBG_sqrtS_reco_isoPh_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 std::vector<TH1F*> hist_vec_reco_HZ_noBG_1D;  
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dAlpha_j1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dPhi_j1_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dTheta_j1_j2_reco_sqrt_s_2500);//8
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_Theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_Theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_Theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_Theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_Theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_Theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_mass_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_mass_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_mass_j1_reco_sqrt_s_2500);//17
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_mass_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_mass_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_mass_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_lrz_j1_reco_sqrt_s_2500);//26
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_lrz_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_j2_reco_sqrt_s_2500);//35
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau21_lrz_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_tau32_lrz_j2_reco_sqrt_s_2500);//44
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_j1_reco_sqrt_s_2500);//53
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//62

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta1_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_j2_reco_sqrt_s_2500);//71
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta2_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//80
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_j1_reco_sqrt_s_2500);//89
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500);//98
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta1_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_j2_reco_sqrt_s_2500);//107
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta2_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500);//116
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_j1_reco_sqrt_s_2500);//125

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//134

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_j2_reco_sqrt_s_2500);//143

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//152
 
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_j1_reco_sqrt_s_2500);//161

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500);//170

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_j2_reco_sqrt_s_2500);//179

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500);//188

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_j1_reco_sqrt_s_2500);//197

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//206

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_j2_reco_sqrt_s_2500);//215

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//224

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_E_theta_j1_reco_sqrt_s_2500);


 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_j2_reco_sqrt_s_2500);//233

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_D2_1_2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet2_E_j1_reco_sqrt_s_2500);//242

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet2_E_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_over_jetE_j1_reco_sqrt_s_2500);//251

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_subjet1_E_over_jetE_j2_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500);//260

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500);//269

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500);

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500);//278

 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_isoPh);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco);
 hist_vec_reco_HZ_noBG_1D.push_back(h_HZ_noBG_sqrtS_reco_isoPh_inv);//281

 for(unsigned int i=0;i<hist_vec_reco_HZ_noBG_1D.size();i++){
   hist_vec_reco_HZ_noBG_1D[i]->Sumw2();
   hist_vec_reco_HZ_noBG_1D[i]->SetLineColor(kBlack);
 }
 std::vector<TH1F*>hist_vec_HZ_noBG_parton;
 std::vector<TH1F*>hist_vec_gen_HZ_noBG_1D;
 std::vector<TH2F*>hist_vec_HZ_noBG_2DHist;

 fillPartonInfo=false;
 fillGenInfo=false;

  fill_HZ_histograms(file_CLIC_HZqq_noBGG, hist_vec_HZ_noBG_parton, hist_vec_gen_HZ_noBG_1D, hist_vec_reco_HZ_noBG_1D, hist_vec_HZ_noBG_2DHist, usePartonInfo ,xsec_hz_qq,fillPartonInfo,fillGenInfo);

 
 unsigned int ind50=h_HZ_mass_j1_reco_sqrt_s_2500->GetXaxis()->FindBin(50.);
 unsigned int ind200=h_HZ_mass_j1_reco_sqrt_s_2500->GetXaxis()->FindBin(200.);
 unsigned int ind_int=0;
 for(unsigned int i=ind50;i<ind200;i++){
   if(h_HZ_mass_j1_reco_sqrt_s_2500->GetBinContent(i)>h_HZ_mass_j2_reco_sqrt_s_2500->GetBinContent(i)){
     ind_int=i;
     break;
   }
 }
 std::cout<<"bin 50/int/200: "<<ind50<<"/"<<ind_int<<"/"<<ind200<<std::endl;
 std::cout<<"efficiency integral 50-int/int-200/overlap: "<<h_HZ_mass_j2_reco_sqrt_s_2500->Integral(ind50,ind_int)/h_HZ_mass_j2_reco_sqrt_s_2500->Integral(ind50,ind200)<<"/"<<h_HZ_mass_j1_reco_sqrt_s_2500->Integral(ind_int,200)/h_HZ_mass_j2_reco_sqrt_s_2500->Integral(ind50,ind200)<<"/"<<0.5*(h_HZ_mass_j2_reco_sqrt_s_2500->Integral(ind_int,ind200)/h_HZ_mass_j2_reco_sqrt_s_2500->Integral(ind50,ind200)+h_HZ_mass_j1_reco_sqrt_s_2500->Integral(ind50,ind_int)/h_HZ_mass_j1_reco_sqrt_s_2500->Integral(ind50,ind200))<<std::endl;

  file_histogram->Write();
  file_histogram->Close();

}
