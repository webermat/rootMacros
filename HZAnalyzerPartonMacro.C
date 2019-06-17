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
 
  vector<float> *recojet10_Px=0;
  vector<float> *recojet10_Py=0;
  vector<float> *recojet10_Pz=0;
  vector<float> *recojet10_E=0;

  vector<float> *recojet12_Px=0;
  vector<float> *recojet12_Py=0;
  vector<float> *recojet12_Pz=0;
  vector<float> *recojet12_E=0;

  vector<float> *genjet10_Px=0;
  vector<float> *genjet10_Py=0;
  vector<float> *genjet10_Pz=0;
  vector<float> *genjet10_E=0;

  vector<float> *genjet12_Px=0;
  vector<float> *genjet12_Py=0;
  vector<float> *genjet12_Pz=0;
  vector<float> *genjet12_E=0;

  vector<float> *recojet10_subjet_Px=0;
  vector<float> *recojet10_subjet_Py=0;
  vector<float> *recojet10_subjet_Pz=0;
  vector<float> *recojet10_subjet_E=0;
  vector<int> *recojet10_subjet_jetindex=0;

  vector<float> *recojet12_subjet_Px=0;
  vector<float> *recojet12_subjet_Py=0;
  vector<float> *recojet12_subjet_Pz=0;
  vector<float> *recojet12_subjet_E=0;
  vector<int> *recojet12_subjet_jetindex=0;

  vector<float> *genjet10_subjet_Px=0;
  vector<float> *genjet10_subjet_Py=0;
  vector<float> *genjet10_subjet_Pz=0;
  vector<float> *genjet10_subjet_E=0;
  vector<int> *genjet10_subjet_jetindex=0;

  vector<float> *genjet12_subjet_Px=0;
  vector<float> *genjet12_subjet_Py=0;
  vector<float> *genjet12_subjet_Pz=0;
  vector<float> *genjet12_subjet_E=0;
  vector<int> *genjet12_subjet_jetindex=0;


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

  tree->SetBranchAddress("genjet10_E", &genjet10_E);
  tree->SetBranchAddress("genjet10_Px", &genjet10_Px);
  tree->SetBranchAddress("genjet10_Py", &genjet10_Py);
  tree->SetBranchAddress("genjet10_Pz", &genjet10_Pz);

  tree->SetBranchAddress("genjet12_E", &genjet12_E);
  tree->SetBranchAddress("genjet12_Px", &genjet12_Px);
  tree->SetBranchAddress("genjet12_Py", &genjet12_Py);
  tree->SetBranchAddress("genjet12_Pz", &genjet12_Pz);

  tree->SetBranchAddress("recojet10_E", &recojet10_E);
  tree->SetBranchAddress("recojet10_Px", &recojet10_Px);
  tree->SetBranchAddress("recojet10_Py", &recojet10_Py);
  tree->SetBranchAddress("recojet10_Pz", &recojet10_Pz);

  tree->SetBranchAddress("recojet12_E", &recojet12_E);
  tree->SetBranchAddress("recojet12_Px", &recojet12_Px);
  tree->SetBranchAddress("recojet12_Py", &recojet12_Py);
  tree->SetBranchAddress("recojet12_Pz", &recojet12_Pz);

  tree->SetBranchAddress("genjet10_subjet_E", &genjet10_subjet_E);
  tree->SetBranchAddress("genjet10_subjet_Px", &genjet10_subjet_Px);
  tree->SetBranchAddress("genjet10_subjet_Py", &genjet10_subjet_Py);
  tree->SetBranchAddress("genjet10_subjet_Pz", &genjet10_subjet_Pz);
  tree->SetBranchAddress("genjet10_subjet_jetindex", &genjet10_subjet_jetindex);

  tree->SetBranchAddress("genjet12_subjet_E", &genjet12_subjet_E);
  tree->SetBranchAddress("genjet12_subjet_Px", &genjet12_subjet_Px);
  tree->SetBranchAddress("genjet12_subjet_Py", &genjet12_subjet_Py);
  tree->SetBranchAddress("genjet12_subjet_Pz", &genjet12_subjet_Pz);
  tree->SetBranchAddress("genjet12_subjet_jetindex", &genjet12_subjet_jetindex);

  tree->SetBranchAddress("recojet10_subjet_E", &recojet10_subjet_E);
  tree->SetBranchAddress("recojet10_subjet_Px", &recojet10_subjet_Px);
  tree->SetBranchAddress("recojet10_subjet_Py", &recojet10_subjet_Py);
  tree->SetBranchAddress("recojet10_subjet_Pz", &recojet10_subjet_Pz);
  tree->SetBranchAddress("recojet10_subjet_jetindex", &recojet10_subjet_jetindex);

  tree->SetBranchAddress("recojet12_subjet_E", &recojet12_subjet_E);
  tree->SetBranchAddress("recojet12_subjet_Px", &recojet12_subjet_Px);
  tree->SetBranchAddress("recojet12_subjet_Py", &recojet12_subjet_Py);
  tree->SetBranchAddress("recojet12_subjet_Pz", &recojet12_subjet_Pz);
  tree->SetBranchAddress("recojet12_subjet_jetindex", &recojet12_subjet_jetindex);


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
  

    TLorentzVector tempGenDR10JetSum(0,0,0,0);
    for(unsigned int i=0;i<genjet10_E->size();i++){
     TLorentzVector temp(0,0,0,0);
     temp.SetPxPyPzE((*genjet10_Px)[i],(*genjet10_Py)[i],(*genjet10_Pz)[i],(*genjet10_E)[i]);
     tempGenDR10JetSum+=temp;
    }
    TLorentzVector tempGenDR12JetSum(0,0,0,0);
    for(unsigned int i=0;i<genjet12_E->size();i++){
     TLorentzVector temp(0,0,0,0);
     temp.SetPxPyPzE((*genjet12_Px)[i],(*genjet12_Py)[i],(*genjet12_Pz)[i],(*genjet12_E)[i]);
     tempGenDR12JetSum+=temp;
    }
    TLorentzVector tempRecoDR10JetSum(0,0,0,0);
    for(unsigned int i=0;i<recojet10_E->size();i++){
     TLorentzVector temp(0,0,0,0);
     temp.SetPxPyPzE((*recojet10_Px)[i],(*recojet10_Py)[i],(*recojet10_Pz)[i],(*recojet10_E)[i]);
     tempRecoDR10JetSum+=temp;
    }
    TLorentzVector tempRecoDR12JetSum(0,0,0,0);
    for(unsigned int i=0;i<recojet12_E->size();i++){
     TLorentzVector temp(0,0,0,0);
     temp.SetPxPyPzE((*recojet12_Px)[i],(*recojet12_Py)[i],(*recojet12_Pz)[i],(*recojet12_E)[i]);
     tempRecoDR12JetSum+=temp;
    }
    TLorentzVector gj10_m1(0,0,0,0);
    TLorentzVector gj10_m2(0,0,0,0);
    if(genjet10_E->size()==2){
     gj10_m1.SetPxPyPzE((*genjet10_Px)[0],(*genjet10_Py)[0],(*genjet10_Pz)[0],(*genjet10_E)[0]);
     gj10_m2.SetPxPyPzE((*genjet10_Px)[1],(*genjet10_Py)[1],(*genjet10_Pz)[1],(*genjet10_E)[1]);
     if(gj10_m2.M()>gj10_m1.M()){
       TLorentzVector temp=gj10_m1;
       gj10_m1=gj10_m2;
       gj10_m2=temp;
     }
     if(gj10_m2.M()>=gj10_m1.M()){
       std::cout<<"gj10 mass order should have not been the case nowadays "<<gj10_m2.M()<<"/"<<gj10_m1.M()<<std::endl;
     }
     if((tempTotGenP4-tempGenIsoPhP4).M()>sqrtS_high_reco){
       h_hist_vec[46]->Fill(gj10_m1.M());
       h_hist_vec[47]->Fill(gj10_m2.M());
       h_hist_vec[54]->Fill(DeltaPhi(tempHP4.Phi(),gj10_m1.Phi())*TMath::RadToDeg());
       h_hist_vec[55]->Fill(fabs(tempHP4.Theta()-gj10_m1.Theta())*TMath::RadToDeg());
       h_hist_vec[56]->Fill(DeltaPhi(tempZP4.Phi(),gj10_m2.Phi())*TMath::RadToDeg());
       h_hist_vec[57]->Fill(fabs(tempZP4.Theta()-gj10_m2.Theta())*TMath::RadToDeg());
       h_hist_vec[70]->Fill(tempHP4.Angle(gj10_m1.Vect())*TMath::RadToDeg());
       h_hist_vec[71]->Fill(tempZP4.Angle(gj10_m2.Vect())*TMath::RadToDeg());
       h_hist_vec[78]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),gj10_m1.Phi())*TMath::RadToDeg());
       if(tempTotInvGenP4.Pt()>50){
	 h_hist_vec[84]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),gj10_m1.Phi())*TMath::RadToDeg());
	 if(tempTotInvGenP4.Pt()>100){
	   h_hist_vec[88]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),gj10_m1.Phi())*TMath::RadToDeg());
	   if(tempTotInvGenP4.Pt()>150){
	     h_hist_vec[92]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),gj10_m1.Phi())*TMath::RadToDeg());
	   }
	 }
       }

       if(genjet10_subjet_E->size()==4){
	 if((*genjet10_subjet_jetindex)[0]!=0 || (*genjet10_subjet_jetindex)[1]!=0 || (*genjet10_subjet_jetindex)[2]!=1 || (*genjet10_subjet_jetindex)[3]!=1){
	   std::cout<<"unexpected subjet ordering "<<(*genjet10_subjet_jetindex)[0]<<"/"<<(*genjet10_subjet_jetindex)[1]<<"/"<<(*genjet10_subjet_jetindex)[2]<<"/"<<(*genjet10_subjet_jetindex)[3]<<std::endl;
	 }
	 if(gj10_m1.E()==(*genjet10_E)[0]){
	   //jet1 is Z jet
	   //check subjet indices 2 and 3
	   if((*genjet10_subjet_E)[2]>(*genjet10_subjet_E)[3]){
	     h_hist_vec[100]->Fill((*genjet10_subjet_E)[2]/((*genjet10_subjet_E)[2]+(*genjet10_subjet_E)[3]));
	     h_hist_vec[101]->Fill((*genjet10_subjet_E)[3]/((*genjet10_subjet_E)[2]+(*genjet10_subjet_E)[3]));
	   }else{
	     h_hist_vec[100]->Fill((*genjet10_subjet_E)[3]/((*genjet10_subjet_E)[2]+(*genjet10_subjet_E)[3]));
	     h_hist_vec[101]->Fill((*genjet10_subjet_E)[2]/((*genjet10_subjet_E)[2]+(*genjet10_subjet_E)[3]));
	   }
	   //jet 0 is H jet, check subjet indices 0 and 1
	   if((*genjet10_subjet_E)[0]>(*genjet10_subjet_E)[1]){
	     h_hist_vec[102]->Fill((*genjet10_subjet_E)[0]/((*genjet10_subjet_E)[0]+(*genjet10_subjet_E)[1]));
	     h_hist_vec[103]->Fill((*genjet10_subjet_E)[1]/((*genjet10_subjet_E)[0]+(*genjet10_subjet_E)[1]));
	   }else{
	     h_hist_vec[102]->Fill((*genjet10_subjet_E)[1]/((*genjet10_subjet_E)[0]+(*genjet10_subjet_E)[1]));
	     h_hist_vec[103]->Fill((*genjet10_subjet_E)[0]/((*genjet10_subjet_E)[0]+(*genjet10_subjet_E)[1]));
	   }
	 }else{
	   //jet 0 is Z jet, check subjet indices 0 and 1
	   if((*genjet10_subjet_E)[0]>(*genjet10_subjet_E)[1]){
	     h_hist_vec[100]->Fill((*genjet10_subjet_E)[0]/((*genjet10_subjet_E)[0]+(*genjet10_subjet_E)[1]));
	     h_hist_vec[101]->Fill((*genjet10_subjet_E)[1]/((*genjet10_subjet_E)[0]+(*genjet10_subjet_E)[1]));
	   }else{
	     h_hist_vec[100]->Fill((*genjet10_subjet_E)[1]/((*genjet10_subjet_E)[0]+(*genjet10_subjet_E)[1]));
	     h_hist_vec[101]->Fill((*genjet10_subjet_E)[0]/((*genjet10_subjet_E)[0]+(*genjet10_subjet_E)[1]));
	   }
	   //jet1 is H jet
	   //check subjet indices 2 and 3
	   if((*genjet10_subjet_E)[2]>(*genjet10_subjet_E)[3]){
	     h_hist_vec[102]->Fill((*genjet10_subjet_E)[2]/((*genjet10_subjet_E)[2]+(*genjet10_subjet_E)[3]));
	     h_hist_vec[103]->Fill((*genjet10_subjet_E)[3]/((*genjet10_subjet_E)[2]+(*genjet10_subjet_E)[3]));
	   }else{
	     h_hist_vec[102]->Fill((*genjet10_subjet_E)[3]/((*genjet10_subjet_E)[2]+(*genjet10_subjet_E)[3]));
	     h_hist_vec[103]->Fill((*genjet10_subjet_E)[2]/((*genjet10_subjet_E)[2]+(*genjet10_subjet_E)[3]));
	   }
	 }
       }
     }
    }
    TLorentzVector gj12_m1(0,0,0,0);
    TLorentzVector gj12_m2(0,0,0,0);
    if(genjet12_E->size()==2){
     gj12_m1.SetPxPyPzE((*genjet12_Px)[0],(*genjet12_Py)[0],(*genjet12_Pz)[0],(*genjet12_E)[0]);
     gj12_m2.SetPxPyPzE((*genjet12_Px)[1],(*genjet12_Py)[1],(*genjet12_Pz)[1],(*genjet12_E)[1]);
     if(gj12_m2.M()>gj12_m1.M()){
       TLorentzVector temp=gj12_m1;
       gj12_m1=gj12_m2;
       gj12_m2=temp;
     }
     if(gj12_m2.M()>=gj12_m1.M()){
       std::cout<<"gj12 mass order should have not been the case nowadays "<<gj12_m2.M()<<"/"<<gj12_m1.M()<<std::endl;
     }
     if((tempTotGenP4-tempGenIsoPhP4).M()>sqrtS_high_reco){
       h_hist_vec[48]->Fill(gj12_m1.M());
       h_hist_vec[49]->Fill(gj12_m2.M());
       h_hist_vec[58]->Fill(DeltaPhi(tempHP4.Phi(),gj12_m1.Phi())*TMath::RadToDeg());
       h_hist_vec[59]->Fill(fabs(tempHP4.Theta()-gj12_m1.Theta())*TMath::RadToDeg());
       h_hist_vec[60]->Fill(DeltaPhi(tempZP4.Phi(),gj12_m2.Phi())*TMath::RadToDeg());
       h_hist_vec[61]->Fill(fabs(tempZP4.Theta()-gj12_m2.Theta())*TMath::RadToDeg());
       h_hist_vec[72]->Fill(tempHP4.Angle(gj12_m1.Vect())*TMath::RadToDeg());
       h_hist_vec[73]->Fill(tempZP4.Angle(gj12_m2.Vect())*TMath::RadToDeg());
       h_hist_vec[79]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),gj12_m1.Phi())*TMath::RadToDeg());
       if(tempTotInvGenP4.Pt()>50){
	 h_hist_vec[85]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),gj12_m1.Phi())*TMath::RadToDeg());
	 if(tempTotInvGenP4.Pt()>100){
	   h_hist_vec[89]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),gj12_m1.Phi())*TMath::RadToDeg());
	   if(tempTotInvGenP4.Pt()>150){
	     h_hist_vec[93]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),gj12_m1.Phi())*TMath::RadToDeg());
	   }
	 }
       }
       if(genjet12_subjet_E->size()==4){
	 if((*genjet12_subjet_jetindex)[0]!=0 || (*genjet12_subjet_jetindex)[1]!=0 || (*genjet12_subjet_jetindex)[2]!=1 || (*genjet12_subjet_jetindex)[3]!=1){
	   std::cout<<"unexpected subjet ordering "<<(*genjet12_subjet_jetindex)[0]<<"/"<<(*genjet12_subjet_jetindex)[1]<<"/"<<(*genjet12_subjet_jetindex)[2]<<"/"<<(*genjet12_subjet_jetindex)[3]<<std::endl;
	 }
	 if(gj12_m1.E()==(*genjet12_E)[0]){
	   //jet1 is Z jet
	   //check subjet indices 2 and 3
	   if((*genjet12_subjet_E)[2]>(*genjet12_subjet_E)[3]){
	     h_hist_vec[104]->Fill((*genjet12_subjet_E)[2]/((*genjet12_subjet_E)[2]+(*genjet12_subjet_E)[3]));
	     h_hist_vec[105]->Fill((*genjet12_subjet_E)[3]/((*genjet12_subjet_E)[2]+(*genjet12_subjet_E)[3]));
	   }else{
	     h_hist_vec[104]->Fill((*genjet12_subjet_E)[3]/((*genjet12_subjet_E)[2]+(*genjet12_subjet_E)[3]));
	     h_hist_vec[105]->Fill((*genjet12_subjet_E)[2]/((*genjet12_subjet_E)[2]+(*genjet12_subjet_E)[3]));
	   }
	   //jet 0 is H jet, check subjet indices 0 and 1
	   if((*genjet12_subjet_E)[0]>(*genjet12_subjet_E)[1]){
	     h_hist_vec[106]->Fill((*genjet12_subjet_E)[0]/((*genjet12_subjet_E)[0]+(*genjet12_subjet_E)[1]));
	     h_hist_vec[107]->Fill((*genjet12_subjet_E)[1]/((*genjet12_subjet_E)[0]+(*genjet12_subjet_E)[1]));
	   }else{
	     h_hist_vec[106]->Fill((*genjet12_subjet_E)[1]/((*genjet12_subjet_E)[0]+(*genjet12_subjet_E)[1]));
	     h_hist_vec[107]->Fill((*genjet12_subjet_E)[0]/((*genjet12_subjet_E)[0]+(*genjet12_subjet_E)[1]));
	   }
	 }else{
	   //jet 0 is Z jet, check subjet indices 0 and 1
	   if((*genjet12_subjet_E)[0]>(*genjet12_subjet_E)[1]){
	     h_hist_vec[104]->Fill((*genjet12_subjet_E)[0]/((*genjet12_subjet_E)[0]+(*genjet12_subjet_E)[1]));
	     h_hist_vec[105]->Fill((*genjet12_subjet_E)[1]/((*genjet12_subjet_E)[0]+(*genjet12_subjet_E)[1]));
	   }else{
	     h_hist_vec[106]->Fill((*genjet12_subjet_E)[1]/((*genjet12_subjet_E)[0]+(*genjet12_subjet_E)[1]));
	     h_hist_vec[107]->Fill((*genjet12_subjet_E)[0]/((*genjet12_subjet_E)[0]+(*genjet12_subjet_E)[1]));
	   }
	   //jet1 is H jet
	   //check subjet indices 2 and 3
	   if((*genjet12_subjet_E)[2]>(*genjet12_subjet_E)[3]){
	     h_hist_vec[106]->Fill((*genjet12_subjet_E)[2]/((*genjet12_subjet_E)[2]+(*genjet12_subjet_E)[3]));
	     h_hist_vec[107]->Fill((*genjet12_subjet_E)[3]/((*genjet12_subjet_E)[2]+(*genjet12_subjet_E)[3]));
	   }else{
	     h_hist_vec[106]->Fill((*genjet12_subjet_E)[3]/((*genjet12_subjet_E)[2]+(*genjet12_subjet_E)[3]));
	     h_hist_vec[107]->Fill((*genjet12_subjet_E)[2]/((*genjet12_subjet_E)[2]+(*genjet12_subjet_E)[3]));
	   }
	 }
       }
     }
    }
    TLorentzVector rj10_m1(0,0,0,0);
    TLorentzVector rj10_m2(0,0,0,0);
    if(recojet10_E->size()==2){
      rj10_m1.SetPxPyPzE((*recojet10_Px)[0],(*recojet10_Py)[0],(*recojet10_Pz)[0],(*recojet10_E)[0]);
      rj10_m2.SetPxPyPzE((*recojet10_Px)[1],(*recojet10_Py)[1],(*recojet10_Pz)[1],(*recojet10_E)[1]);
      if(rj10_m2.M()>rj10_m1.M()){
	TLorentzVector temp=rj10_m1;
	rj10_m1=rj10_m2;
	rj10_m2=temp;
      }
      if(rj10_m2.M()>=rj10_m1.M()){
	std::cout<<"rj10 mass order should have not been the case nowadays "<<rj10_m2.M()<<"/"<<rj10_m1.M()<<std::endl;
      }
      if((tempTotRecoP4-tempRecoIsoPhP4).M()>sqrtS_high_reco){
	h_hist_vec[50]->Fill(rj10_m1.M());
	h_hist_vec[51]->Fill(rj10_m2.M());
	h_hist_vec[62]->Fill(DeltaPhi(tempHP4.Phi(),rj10_m1.Phi())*TMath::RadToDeg());
	h_hist_vec[63]->Fill(fabs(tempHP4.Theta()-rj10_m1.Theta())*TMath::RadToDeg());
	h_hist_vec[64]->Fill(DeltaPhi(tempZP4.Phi(),rj10_m2.Phi())*TMath::RadToDeg());
	h_hist_vec[65]->Fill(fabs(tempZP4.Theta()-rj10_m2.Theta())*TMath::RadToDeg());
	h_hist_vec[74]->Fill(tempHP4.Angle(rj10_m1.Vect())*TMath::RadToDeg());
	h_hist_vec[75]->Fill(tempZP4.Angle(rj10_m2.Vect())*TMath::RadToDeg());
	h_hist_vec[80]->Fill(DeltaPhi(tempRecoMETP4.Phi(),rj10_m1.Phi())*TMath::RadToDeg());
	if(tempRecoMETP4.Pt()>50){
	 h_hist_vec[86]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),rj10_m1.Phi())*TMath::RadToDeg());
	 if(tempRecoMETP4.Pt()>100){
	   h_hist_vec[90]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),rj10_m1.Phi())*TMath::RadToDeg());
	   if(tempRecoMETP4.Pt()>150){
	     h_hist_vec[94]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),rj10_m1.Phi())*TMath::RadToDeg());
	   }
	 }
	}
	if(recojet10_subjet_E->size()==4){
	  if((*recojet10_subjet_jetindex)[0]!=0 || (*recojet10_subjet_jetindex)[1]!=0 || (*recojet10_subjet_jetindex)[2]!=1 || (*recojet10_subjet_jetindex)[3]!=1){
	    std::cout<<"unexpected subjet ordering "<<(*recojet10_subjet_jetindex)[0]<<"/"<<(*recojet10_subjet_jetindex)[1]<<"/"<<(*recojet10_subjet_jetindex)[2]<<"/"<<(*recojet10_subjet_jetindex)[3]<<std::endl;
	  }
	  if(rj10_m1.E()==(*recojet10_E)[0]){
	    //jet1 is Z jet
	    //check subjet indices 2 and 3
	    if((*recojet10_subjet_E)[2]>(*recojet10_subjet_E)[3]){
	      h_hist_vec[108]->Fill((*recojet10_subjet_E)[2]/((*recojet10_subjet_E)[2]+(*recojet10_subjet_E)[3]));
	      h_hist_vec[109]->Fill((*recojet10_subjet_E)[3]/((*recojet10_subjet_E)[2]+(*recojet10_subjet_E)[3]));
	    }else{
	      h_hist_vec[108]->Fill((*recojet10_subjet_E)[3]/((*recojet10_subjet_E)[2]+(*recojet10_subjet_E)[3]));
	      h_hist_vec[109]->Fill((*recojet10_subjet_E)[2]/((*recojet10_subjet_E)[2]+(*recojet10_subjet_E)[3]));
	    }
	    //jet 0 is H jet, check subjet indices 0 and 1
	    if((*recojet10_subjet_E)[0]>(*recojet10_subjet_E)[1]){
	      h_hist_vec[110]->Fill((*recojet10_subjet_E)[0]/((*recojet10_subjet_E)[0]+(*recojet10_subjet_E)[1]));
	      h_hist_vec[111]->Fill((*recojet10_subjet_E)[1]/((*recojet10_subjet_E)[0]+(*recojet10_subjet_E)[1]));
	    }else{
	      h_hist_vec[110]->Fill((*recojet10_subjet_E)[1]/((*recojet10_subjet_E)[0]+(*recojet10_subjet_E)[1]));
	      h_hist_vec[111]->Fill((*recojet10_subjet_E)[0]/((*recojet10_subjet_E)[0]+(*recojet10_subjet_E)[1]));
	    }
	  }else{
	    //jet 0 is Z jet, check subjet indices 0 and 1
	    if((*recojet10_subjet_E)[0]>(*recojet10_subjet_E)[1]){
	      h_hist_vec[108]->Fill((*recojet10_subjet_E)[0]/((*recojet10_subjet_E)[0]+(*recojet10_subjet_E)[1]));
	      h_hist_vec[109]->Fill((*recojet10_subjet_E)[1]/((*recojet10_subjet_E)[0]+(*recojet10_subjet_E)[1]));
	    }else{
	      h_hist_vec[108]->Fill((*recojet10_subjet_E)[1]/((*recojet10_subjet_E)[0]+(*recojet10_subjet_E)[1]));
	      h_hist_vec[109]->Fill((*recojet10_subjet_E)[0]/((*recojet10_subjet_E)[0]+(*recojet10_subjet_E)[1]));
	    }
	    //jet1 is H jet
	    //check subjet indices 2 and 3
	    if((*recojet10_subjet_E)[2]>(*recojet10_subjet_E)[3]){
	      h_hist_vec[110]->Fill((*recojet10_subjet_E)[2]/((*recojet10_subjet_E)[2]+(*recojet10_subjet_E)[3]));
	      h_hist_vec[111]->Fill((*recojet10_subjet_E)[3]/((*recojet10_subjet_E)[2]+(*recojet10_subjet_E)[3]));
	    }else{
	      h_hist_vec[110]->Fill((*recojet10_subjet_E)[3]/((*recojet10_subjet_E)[2]+(*recojet10_subjet_E)[3]));
	      h_hist_vec[111]->Fill((*recojet10_subjet_E)[2]/((*recojet10_subjet_E)[2]+(*recojet10_subjet_E)[3]));
	    }
	  }
	}
      }
    }
    TLorentzVector rj12_m1(0,0,0,0);
    TLorentzVector rj12_m2(0,0,0,0);
    if(recojet12_E->size()==2){
     rj12_m1.SetPxPyPzE((*recojet12_Px)[0],(*recojet12_Py)[0],(*recojet12_Pz)[0],(*recojet12_E)[0]);
     rj12_m2.SetPxPyPzE((*recojet12_Px)[1],(*recojet12_Py)[1],(*recojet12_Pz)[1],(*recojet12_E)[1]);
     if(rj12_m2.M()>rj12_m1.M()){
       TLorentzVector temp=rj12_m1;
       rj12_m1=rj12_m2;
       rj12_m2=temp;
     }
     if(rj12_m2.M()>=rj12_m1.M()){
       std::cout<<"rj12 mass order should have not been the case nowadays "<<rj12_m2.M()<<"/"<<rj12_m1.M()<<std::endl;
     }
     if((tempTotRecoP4-tempRecoIsoPhP4).M()>sqrtS_high_reco){
       h_hist_vec[52]->Fill(rj12_m1.M());
       h_hist_vec[53]->Fill(rj12_m2.M());
       h_hist_vec[66]->Fill(DeltaPhi(tempHP4.Phi(),rj12_m1.Phi())*TMath::RadToDeg());
       h_hist_vec[67]->Fill(fabs(tempHP4.Theta()-rj12_m1.Theta())*TMath::RadToDeg());
       h_hist_vec[68]->Fill(DeltaPhi(tempZP4.Phi(),rj12_m2.Phi())*TMath::RadToDeg());
       h_hist_vec[69]->Fill(fabs(tempZP4.Theta()-rj12_m2.Theta())*TMath::RadToDeg());
       h_hist_vec[76]->Fill(tempHP4.Angle(rj12_m1.Vect())*TMath::RadToDeg());
       h_hist_vec[77]->Fill(tempZP4.Angle(rj12_m2.Vect())*TMath::RadToDeg());
       h_hist_vec[81]->Fill(DeltaPhi(tempRecoMETP4.Phi(),rj12_m1.Phi())*TMath::RadToDeg());
       if(tempRecoMETP4.Pt()>50){
	 h_hist_vec[87]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),rj12_m1.Phi())*TMath::RadToDeg());
	 if(tempRecoMETP4.Pt()>100){
	   h_hist_vec[91]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),rj12_m1.Phi())*TMath::RadToDeg());
	   if(tempRecoMETP4.Pt()>150){
	     h_hist_vec[95]->Fill(DeltaPhi(tempTotInvGenP4.Phi(),rj12_m1.Phi())*TMath::RadToDeg());
	   }
	 }
       }


       if(recojet12_subjet_E->size()==4){
	 if((*recojet12_subjet_jetindex)[0]!=0 || (*recojet12_subjet_jetindex)[1]!=0 || (*recojet12_subjet_jetindex)[2]!=1 || (*recojet12_subjet_jetindex)[3]!=1){
	   std::cout<<"unexpected subjet ordering "<<(*recojet12_subjet_jetindex)[0]<<"/"<<(*recojet12_subjet_jetindex)[1]<<"/"<<(*recojet12_subjet_jetindex)[2]<<"/"<<(*recojet12_subjet_jetindex)[3]<<std::endl;
	 }
	 if(rj12_m1.E()==(*recojet12_E)[0]){
	   //jet1 is Z jet
	   //check subjet indices 2 and 3
	   if((*recojet12_subjet_E)[2]>(*recojet12_subjet_E)[3]){
	     h_hist_vec[112]->Fill((*recojet12_subjet_E)[2]/((*recojet12_subjet_E)[2]+(*recojet12_subjet_E)[3]));
	     h_hist_vec[113]->Fill((*recojet12_subjet_E)[3]/((*recojet12_subjet_E)[2]+(*recojet12_subjet_E)[3]));
	   }else{
	     h_hist_vec[112]->Fill((*recojet12_subjet_E)[3]/((*recojet12_subjet_E)[2]+(*recojet12_subjet_E)[3]));
	     h_hist_vec[113]->Fill((*recojet12_subjet_E)[2]/((*recojet12_subjet_E)[2]+(*recojet12_subjet_E)[3]));
	   }
	   //jet 0 is H jet, check subjet indices 0 and 1
	   if((*recojet12_subjet_E)[0]>(*recojet12_subjet_E)[1]){
	     h_hist_vec[114]->Fill((*recojet12_subjet_E)[0]/((*recojet12_subjet_E)[0]+(*recojet12_subjet_E)[1]));
	     h_hist_vec[115]->Fill((*recojet12_subjet_E)[1]/((*recojet12_subjet_E)[0]+(*recojet12_subjet_E)[1]));
	   }else{
	     h_hist_vec[114]->Fill((*recojet12_subjet_E)[1]/((*recojet12_subjet_E)[0]+(*recojet12_subjet_E)[1]));
	     h_hist_vec[115]->Fill((*recojet12_subjet_E)[0]/((*recojet12_subjet_E)[0]+(*recojet12_subjet_E)[1]));
	   }
	 }else{
	   //jet 0 is Z jet, check subjet indices 0 and 1
	   if((*recojet12_subjet_E)[0]>(*recojet12_subjet_E)[1]){
	     h_hist_vec[112]->Fill((*recojet12_subjet_E)[0]/((*recojet12_subjet_E)[0]+(*recojet12_subjet_E)[1]));
	     h_hist_vec[113]->Fill((*recojet12_subjet_E)[1]/((*recojet12_subjet_E)[0]+(*recojet12_subjet_E)[1]));
	   }else{
	     h_hist_vec[112]->Fill((*recojet12_subjet_E)[1]/((*recojet12_subjet_E)[0]+(*recojet12_subjet_E)[1]));
	     h_hist_vec[113]->Fill((*recojet12_subjet_E)[0]/((*recojet12_subjet_E)[0]+(*recojet12_subjet_E)[1]));
	   }
	   //jet1 is H jet
	   //check subjet indices 2 and 3
	   if((*recojet12_subjet_E)[2]>(*recojet12_subjet_E)[3]){
	     h_hist_vec[114]->Fill((*recojet12_subjet_E)[2]/((*recojet12_subjet_E)[2]+(*recojet12_subjet_E)[3]));
	     h_hist_vec[115]->Fill((*recojet12_subjet_E)[3]/((*recojet12_subjet_E)[2]+(*recojet12_subjet_E)[3]));
	   }else{
	     h_hist_vec[114]->Fill((*recojet12_subjet_E)[3]/((*recojet12_subjet_E)[2]+(*recojet12_subjet_E)[3]));
	     h_hist_vec[115]->Fill((*recojet12_subjet_E)[2]/((*recojet12_subjet_E)[2]+(*recojet12_subjet_E)[3]));
	   }
	 }
       }
     }//sqrtS cut

    }
    h_hist_vec_2D[0]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenTruePhP4).M());
    h_hist_vec_2D[1]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4).M());
    h_hist_vec_2D[2]->Fill(tempTotEventP4.M(),tempGenDR10JetSum.M());
    h_hist_vec_2D[3]->Fill(tempTotEventP4.M(),tempGenDR12JetSum.M());
    h_hist_vec_2D[4]->Fill(tempTotEventP4.M(),(tempTotRecoP4-tempRecoIsoPhP4).M());
    h_hist_vec_2D[5]->Fill(tempTotEventP4.M(),tempRecoDR10JetSum.M());
    h_hist_vec_2D[6]->Fill(tempTotEventP4.M(),tempRecoDR12JetSum.M());
    h_hist_vec_2D[7]->Fill((tempTotGenP4-tempGenIsoPhP4).M(),(tempTotRecoP4-tempRecoIsoPhP4).M());
    h_hist_vec_2D[8]->Fill(tempGenDR10JetSum.M(),tempRecoDR10JetSum.M());
    h_hist_vec_2D[9]->Fill(tempGenDR12JetSum.M(),tempRecoDR12JetSum.M());
    h_hist_vec_2D[10]->Fill(tempTotEventP4.M(),tempTotGenP4.M());
    h_hist_vec_2D[11]->Fill(tempTotEventP4.M(),tempTotRecoP4.M());
    h_hist_vec_2D[12]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4+tempTotInvGenP4).M());
    h_hist_vec_2D[13]->Fill(tempTotEventP4.M(),(tempTotGenP4+tempTotInvGenP4).M());
 
    if(tempTotEventP4.M()<sqrtS_low){
      h_hist_vec[4]->Fill(tempZP4.Pt());
      h_hist_vec[25]->Fill(tempZP4.Theta()*TMath::RadToDeg());
      h_hist_vec[28]->Fill(fabs(tempZP4.Theta()-tempHP4.Theta())*TMath::RadToDeg());
      h_hist_vec[31]->Fill(DeltaPhi(tempZP4.Phi(),tempHP4.Phi())*TMath::RadToDeg());
      h_hist_vec[34]->Fill(tempZP4.Angle(tempHP4.Vect())*TMath::RadToDeg());
    }else if(tempTotEventP4.M()<sqrtS_high){
      h_hist_vec[5]->Fill(tempZP4.Pt());
      h_hist_vec[26]->Fill(tempZP4.Theta()*TMath::RadToDeg());
      h_hist_vec[29]->Fill(fabs(tempZP4.Theta()-tempHP4.Theta())*TMath::RadToDeg());
      h_hist_vec[32]->Fill(DeltaPhi(tempZP4.Phi(),tempHP4.Phi())*TMath::RadToDeg());
      h_hist_vec[35]->Fill(tempZP4.Angle(tempHP4.Vect())*TMath::RadToDeg());
    }else{
      h_hist_vec[6]->Fill(tempZP4.Pt());
      h_hist_vec[27]->Fill(tempZP4.Theta()*TMath::RadToDeg());
      h_hist_vec[30]->Fill(fabs(tempZP4.Theta()-tempHP4.Theta())*TMath::RadToDeg());
      h_hist_vec[33]->Fill(DeltaPhi(tempZP4.Phi(),tempHP4.Phi())*TMath::RadToDeg());
      h_hist_vec[36]->Fill(tempZP4.Angle(tempHP4.Vect())*TMath::RadToDeg());
    }

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

  double lim_energy_low=-200;
  double lim_energy_high=3500.;

  double lim_mass_low=0;
  double lim_mass_high=500;

  TH1F* h_sqrtS_e1_e2_effective = new TH1F("h_sqrtS_e1_e2_effective","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_H_pt_sqrt_s_0_750 = new TH1F("h_H_pt_sqrt_s_0_750","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_H_pt_sqrt_s_750_2750 = new TH1F("h_H_pt_sqrt_s_750_2750","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_H_pt_sqrt_s_2750_3000 = new TH1F("h_H_pt_sqrt_s_2750_3000","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_Z_pt_sqrt_s_0_750 = new TH1F("h_Z_pt_sqrt_s_0_750","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_Z_pt_sqrt_s_750_2750 = new TH1F("h_Z_pt_sqrt_s_750_2750","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_Z_pt_sqrt_s_2750_3000 = new TH1F("h_Z_pt_sqrt_s_2750_3000","", n_bins_high, lim_energy_low,lim_energy_high);

  double lim_dalpha_low=0.;
  double lim_dalpha_high=180.;

  TH1F* h_H_dalpha_bb_sqrt_s_0_750 = new TH1F("h_H_dalpha_bb_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_H_dalpha_bb_sqrt_s_750_2750 = new TH1F("h_H_dalpha_bb_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_H_dalpha_bb_sqrt_s_2750_3000 = new TH1F("h_H_dalpha_bb_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_Z_dalpha_qq_sqrt_s_0_750 = new TH1F("h_Z_dalpha_qq_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_Z_dalpha_qq_sqrt_s_750_2750 = new TH1F("h_Z_dalpha_qq_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_Z_dalpha_qq_sqrt_s_2750_3000 = new TH1F("h_Z_dalpha_qq_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  double lim_theta_low=-180.;
  double lim_theta_high=180.;

  TH1F* h_H_theta_b_fw_sqrt_s_0_750 = new TH1F("h_H_theta_b_fw_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_H_theta_b_fw_sqrt_s_750_2750 = new TH1F("h_H_theta_b_fw_sqrt_s_750_2750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_H_theta_b_fw_sqrt_s_2750_3000 = new TH1F("h_H_theta_b_fw_sqrt_s_0_2750_3000","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_q_fw_sqrt_s_0_750 = new TH1F("h_Z_theta_q_fw_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_q_fw_sqrt_s_750_2750 = new TH1F("h_Z_theta_q_fw_sqrt_s_750_2750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_q_fw_sqrt_s_2750_3000 = new TH1F("h_Z_theta_q_fw_sqrt_s_0_2750_3000","", n_bins_high, lim_theta_low,lim_theta_high);

  TH1F* h_dalpha_min_Z_q_H_b_sqrt_s_0_750 = new TH1F("h_dalpha_min_Z_q_H_b_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dalpha_min_Z_q_H_b_sqrt_s_750_2750 = new TH1F("h_dalpha_min_Z_q_H_b_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dalpha_min_Z_q_H_b_sqrt_s_2750_3000 = new TH1F("h_dalpha_min_Z_q_H_b_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_H_theta_sqrt_s_0_750 = new TH1F("h_H_theta_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_H_theta_sqrt_s_750_2750 = new TH1F("h_H_theta_sqrt_s_750_2750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_H_theta_sqrt_s_2750_3000 = new TH1F("h_H_theta_sqrt_s_0_2750_3000","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_sqrt_s_0_750 = new TH1F("h_Z_theta_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_sqrt_s_750_2750 = new TH1F("h_Z_theta_sqrt_s_750_2750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_Z_theta_sqrt_s_2750_3000 = new TH1F("h_Z_theta_sqrt_s_0_2750_3000","", n_bins_high, lim_theta_low,lim_theta_high);

  TH1F* h_dtheta_H_Z_sqrt_s_0_750 = new TH1F("h_dtheta_H_Z_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dtheta_H_Z_sqrt_s_750_2750 = new TH1F("h_dtheta_H_Z_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dtheta_H_Z_sqrt_s_2750_3000 = new TH1F("h_dtheta_H_Z_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_dphi_H_Z_sqrt_s_0_750 = new TH1F("h_dphi_H_Z_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dphi_H_Z_sqrt_s_750_2750 = new TH1F("h_dphi_H_Z_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dphi_H_Z_sqrt_s_2750_3000 = new TH1F("h_dphi_H_Z_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_dalpha_H_Z_sqrt_s_0_750 = new TH1F("h_dalpha_H_Z_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dalpha_H_Z_sqrt_s_750_2750 = new TH1F("h_dalpha_H_Z_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_dalpha_H_Z_sqrt_s_2750_3000 = new TH1F("h_dalpha_H_Z_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_sqrtS_gen_isoPh = new TH1F("h_sqrtS_gen_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_sqrtS_reco_isoPh = new TH1F("h_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_sqrtS_HZ_effective = new TH1F("h_sqrtS_HZ_effective","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_sqrtS_gen_isoPh_inv = new TH1F("h_sqrtS_gen_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_sqrtS_reco = new TH1F("h_sqrtS_reco","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_sqrtS_gen = new TH1F("h_sqrtS_gen","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_sqrtS_gen_inv = new TH1F("h_sqrtS_gen_inv","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_sqrtS_gen_truePh_inv = new TH1F("h_sqrtS_gen_truePh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_sqrtS_reco_isoPh_inv = new TH1F("h_sqrtS_reco_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_mass1_gj10_sqrtS_2500 = new TH1F("h_mass1_gj10_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_mass2_gj10_sqrtS_2500 = new TH1F("h_mass2_gj10_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_mass1_gj12_sqrtS_2500 = new TH1F("h_mass1_gj12_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_mass2_gj12_sqrtS_2500 = new TH1F("h_mass2_gj12_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);

  TH1F* h_mass1_rj10_sqrtS_2500 = new TH1F("h_mass1_rj10_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_mass2_rj10_sqrtS_2500 = new TH1F("h_mass2_rj10_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_mass1_rj12_sqrtS_2500 = new TH1F("h_mass1_rj12_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_mass2_rj12_sqrtS_2500 = new TH1F("h_mass2_rj12_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);

  double lim_phi_low=-5.;
  double lim_phi_high=185.;

  TH1F* h_dPhi_mass1_gj10_H_sqrtS_2500 = new TH1F("h_dPhi_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass1_gj10_H_sqrtS_2500 = new TH1F("h_dTheta_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_mass2_gj10_H_sqrtS_2500 = new TH1F("h_dPhi_mass2_gj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass2_gj10_H_sqrtS_2500 = new TH1F("h_dTheta_mass2_gj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_dPhi_mass1_gj12_H_sqrtS_2500 = new TH1F("h_dPhi_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass1_gj12_H_sqrtS_2500 = new TH1F("h_dTheta_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_mass2_gj12_H_sqrtS_2500 = new TH1F("h_dPhi_mass2_gj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass2_gj12_H_sqrtS_2500 = new TH1F("h_dTheta_mass2_gj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_dPhi_mass1_rj10_H_sqrtS_2500 = new TH1F("h_dPhi_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass1_rj10_H_sqrtS_2500 = new TH1F("h_dTheta_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_mass2_rj10_H_sqrtS_2500 = new TH1F("h_dPhi_mass2_rj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass2_rj10_H_sqrtS_2500 = new TH1F("h_dTheta_mass2_rj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_dPhi_mass1_rj12_H_sqrtS_2500 = new TH1F("h_dPhi_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass1_rj12_H_sqrtS_2500 = new TH1F("h_dTheta_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_mass2_rj12_H_sqrtS_2500 = new TH1F("h_dPhi_mass2_rj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dTheta_mass2_rj12_H_sqrtS_2500 = new TH1F("h_dTheta_mass2_rj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_dAlpha_mass1_gj10_H_sqrtS_2500 = new TH1F("h_dAlpha_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dAlpha_mass2_gj10_Z_sqrtS_2500 = new TH1F("h_dAlpha_mass2_gj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_dAlpha_mass1_gj12_H_sqrtS_2500 = new TH1F("h_dAlpha_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dAlpha_mass2_gj12_Z_sqrtS_2500 = new TH1F("h_dAlpha_mass2_gj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_dAlpha_mass1_rj10_H_sqrtS_2500 = new TH1F("h_dAlpha_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dAlpha_mass2_rj10_Z_sqrtS_2500 = new TH1F("h_dAlpha_mass2_rj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_dAlpha_mass1_rj12_H_sqrtS_2500 = new TH1F("h_dAlpha_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dAlpha_mass2_rj12_Z_sqrtS_2500 = new TH1F("h_dAlpha_mass2_rj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_dPhi_MET_mass1_gj10_H_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass1_gj12_H_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass1_rj10_H_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass1_rj12_H_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  double lim_MET_low=0.;
  double lim_MET_high=250.;
  TH1F* h_genMET_sqrtS_2500 = new TH1F("h_genMET_sqrtS_2500","", n_bins_high, lim_MET_low,lim_MET_high);
  TH1F* h_recoMET_sqrtS_2500 = new TH1F("h_recoMET_sqrtS_2500","", n_bins_high, lim_MET_low,lim_MET_high);

  TH1F* h_dPhi_MET50_mass1_gj10_H_sqrtS_2500 = new TH1F("h_dPhi_MET50_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET50_mass1_gj12_H_sqrtS_2500 = new TH1F("h_dPhi_MET50_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET50_mass1_rj10_H_sqrtS_2500 = new TH1F("h_dPhi_MET50_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET50_mass1_rj12_H_sqrtS_2500 = new TH1F("h_dPhi_MET50_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_dPhi_MET100_mass1_gj10_H_sqrtS_2500 = new TH1F("h_dPhi_MET100_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET100_mass1_gj12_H_sqrtS_2500 = new TH1F("h_dPhi_MET100_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET100_mass1_rj10_H_sqrtS_2500 = new TH1F("h_dPhi_MET100_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET100_mass1_rj12_H_sqrtS_2500 = new TH1F("h_dPhi_MET100_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_dPhi_MET150_mass1_gj10_H_sqrtS_2500 = new TH1F("h_dPhi_MET150_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET150_mass1_gj12_H_sqrtS_2500 = new TH1F("h_dPhi_MET150_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET150_mass1_rj10_H_sqrtS_2500 = new TH1F("h_dPhi_MET150_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET150_mass1_rj12_H_sqrtS_2500 = new TH1F("h_dPhi_MET150_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  double lim_sj1ratio_low=0.5;
  double lim_sj1ratio_high=1.;

  double lim_sj2ratio_low=0;
  double lim_sj2ratio_high=0.5;

  TH1F* h_E_q1_over_Z = new TH1F("h_E_q1_over_Z","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_q2_over_Z = new TH1F("h_E_q2_over_Z","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_E_b1_over_H = new TH1F("h_E_b1_over_H","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_b2_over_H = new TH1F("h_E_b2_over_H","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);

  TH1F* h_E_s1_over_m2_gj10 = new TH1F("h_E_s1_over_m2_gj10","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_s2_over_m2_gj10 = new TH1F("h_E_s2_over_m2_gj10","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_E_s1_over_m1_gj10 = new TH1F("h_E_s1_over_m1_gj10","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_s2_over_m1_gj10 = new TH1F("h_E_s2_over_m1_gj10","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);

  TH1F* h_E_s1_over_m2_gj12 = new TH1F("h_E_s1_over_m2_gj12","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_s2_over_m2_gj12 = new TH1F("h_E_s2_over_m2_gj12","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_E_s1_over_m1_gj12 = new TH1F("h_E_s1_over_m1_gj12","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_s2_over_m1_gj12 = new TH1F("h_E_s2_over_m1_gj12","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);

  TH1F* h_E_s1_over_m2_rj10 = new TH1F("h_E_s1_over_m2_rj10","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_s2_over_m2_rj10 = new TH1F("h_E_s2_over_m2_rj10","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_E_s1_over_m1_rj10 = new TH1F("h_E_s1_over_m1_rj10","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_s2_over_m1_rj10 = new TH1F("h_E_s2_over_m1_rj10","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);

  TH1F* h_E_s1_over_m2_rj12 = new TH1F("h_E_s1_over_m2_rj12","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_s2_over_m2_rj12 = new TH1F("h_E_s2_over_m2_rj12","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_E_s1_over_m1_rj12 = new TH1F("h_E_s1_over_m1_rj12","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_E_s2_over_m1_rj12 = new TH1F("h_E_s2_over_m1_rj12","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);

  std::vector<TH1F*> hist_vec_HZ_parton;
  hist_vec_HZ_parton.push_back(h_sqrtS_e1_e2_effective);
  hist_vec_HZ_parton.push_back(h_H_pt_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_pt_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_H_pt_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_Z_pt_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_Z_pt_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_Z_pt_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_H_dalpha_bb_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_dalpha_bb_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_H_dalpha_bb_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_Z_dalpha_qq_sqrt_s_0_750);//10
  hist_vec_HZ_parton.push_back(h_Z_dalpha_qq_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_Z_dalpha_qq_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_H_theta_b_fw_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_theta_b_fw_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_H_theta_b_fw_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_Z_theta_q_fw_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_Z_theta_q_fw_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_Z_theta_q_fw_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_dalpha_min_Z_q_H_b_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_dalpha_min_Z_q_H_b_sqrt_s_750_2750);//20
  hist_vec_HZ_parton.push_back(h_dalpha_min_Z_q_H_b_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_H_theta_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_theta_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_H_theta_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_Z_theta_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_Z_theta_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_Z_theta_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_dtheta_H_Z_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_dtheta_H_Z_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_dtheta_H_Z_sqrt_s_2750_3000);//30
  hist_vec_HZ_parton.push_back(h_dphi_H_Z_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_dphi_H_Z_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_dphi_H_Z_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_dalpha_H_Z_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_dalpha_H_Z_sqrt_s_750_2750);
  hist_vec_HZ_parton.push_back(h_dalpha_H_Z_sqrt_s_2750_3000);
  hist_vec_HZ_parton.push_back(h_sqrtS_gen_isoPh);
  hist_vec_HZ_parton.push_back(h_sqrtS_reco_isoPh);
  hist_vec_HZ_parton.push_back(h_sqrtS_HZ_effective);
  hist_vec_HZ_parton.push_back(h_sqrtS_gen_isoPh_inv);//40
  hist_vec_HZ_parton.push_back(h_sqrtS_reco);
  hist_vec_HZ_parton.push_back(h_sqrtS_gen);
  hist_vec_HZ_parton.push_back(h_sqrtS_gen_inv);
  hist_vec_HZ_parton.push_back(h_sqrtS_gen_truePh_inv);
  hist_vec_HZ_parton.push_back(h_sqrtS_reco_isoPh_inv);
  hist_vec_HZ_parton.push_back(h_mass1_gj10_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_mass2_gj10_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_mass1_gj12_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_mass2_gj12_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_mass1_rj10_sqrtS_2500);//50
  hist_vec_HZ_parton.push_back(h_mass2_rj10_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_mass1_rj12_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_mass2_rj12_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass1_gj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass1_gj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass2_gj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass2_gj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass1_gj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass1_gj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass2_gj12_H_sqrtS_2500);//60
  hist_vec_HZ_parton.push_back(h_dTheta_mass2_gj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass1_rj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass1_rj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass2_rj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass2_rj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass1_rj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass1_rj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_mass2_rj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dTheta_mass2_rj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass1_gj10_H_sqrtS_2500);//70
  hist_vec_HZ_parton.push_back(h_dAlpha_mass2_gj10_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass1_gj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass2_gj12_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass1_rj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass2_rj10_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass1_rj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass2_rj12_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass1_gj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass1_gj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass1_rj10_H_sqrtS_2500);//80
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass1_rj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_genMET_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_recoMET_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET50_mass1_gj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET50_mass1_gj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET50_mass1_rj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET50_mass1_rj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET100_mass1_gj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET100_mass1_gj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET100_mass1_rj10_H_sqrtS_2500);//90
  hist_vec_HZ_parton.push_back(h_dPhi_MET100_mass1_rj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET150_mass1_gj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET150_mass1_gj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET150_mass1_rj10_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET150_mass1_rj12_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_E_q1_over_Z);//96
  hist_vec_HZ_parton.push_back(h_E_q2_over_Z);
  hist_vec_HZ_parton.push_back(h_E_b1_over_H);
  hist_vec_HZ_parton.push_back(h_E_b2_over_H);  
  hist_vec_HZ_parton.push_back(h_E_s1_over_m2_gj10);//100
  hist_vec_HZ_parton.push_back(h_E_s2_over_m2_gj10);
  hist_vec_HZ_parton.push_back(h_E_s1_over_m1_gj10);
  hist_vec_HZ_parton.push_back(h_E_s2_over_m1_gj10);  
  hist_vec_HZ_parton.push_back(h_E_s1_over_m2_gj12);
  hist_vec_HZ_parton.push_back(h_E_s2_over_m2_gj12);
  hist_vec_HZ_parton.push_back(h_E_s1_over_m1_gj12);
  hist_vec_HZ_parton.push_back(h_E_s2_over_m1_gj12);
  hist_vec_HZ_parton.push_back(h_E_s1_over_m2_rj10);
  hist_vec_HZ_parton.push_back(h_E_s2_over_m2_rj10);
  hist_vec_HZ_parton.push_back(h_E_s1_over_m1_rj10);//110
  hist_vec_HZ_parton.push_back(h_E_s2_over_m1_rj10);
  hist_vec_HZ_parton.push_back(h_E_s1_over_m2_rj12);
  hist_vec_HZ_parton.push_back(h_E_s2_over_m2_rj12);
  hist_vec_HZ_parton.push_back(h_E_s1_over_m1_rj12);
  hist_vec_HZ_parton.push_back(h_E_s2_over_m1_rj12);//115

  for(unsigned int i=0;i<hist_vec_HZ_parton.size();i++){
    hist_vec_HZ_parton[i]->Sumw2();
  }

 //parton vs gen true --> here we calculat all true minus the gentrue photon
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_genTrue = new TH2F("h_sqrtS_e1e2_vs_sqrtS_genTrue","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_genTrue->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_genTrue->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus iso photons
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh = new TH2F("h_sqrtS_e1e2_vs_sqrtS_gen_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
  //parton vs genJetDR10 sum
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_genJetDR10 = new TH2F("h_sqrtS_e1e2_vs_sqrtS_genJetDR10","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_genJetDR10->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_genJetDR10->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
  //parton vs genJetDR12 sum
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_genJetDR12 = new TH2F("h_sqrtS_e1e2_vs_sqrtS_genJetDR12","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_genJetDR12->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_genJetDR12->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");

 //parton vs reco true --> here we calculat all true minus iso photons
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_reco_isoPh = new TH2F("h_sqrtS_e1e2_vs_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_reco_isoPh->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_reco_isoPh->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs recoJetDR10 sum
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR10 = new TH2F("h_sqrtS_e1e2_vs_sqrtS_recoJetDR10","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR10->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR10->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs recoJetDR12 sum
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR12 = new TH2F("h_sqrtS_e1e2_vs_sqrtS_recoJetDR12","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR12->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR12->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");

 //parton vs reco true --> here we calculat all true minus iso photons
  TH2F* h_sqrtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh = new TH2F("h_sqrtS_gen_isoPh_vs_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
  h_sqrtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs recoJetDR10 sum
  TH2F* h_sqrtS_genJetDR10_eff_vs_sqrtS_recoJetDR10 = new TH2F("h_sqrtS_genJetDR10_vs_sqrtS_recoJetDR10","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_genJetDR10_eff_vs_sqrtS_recoJetDR10->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
  h_sqrtS_genJetDR10_eff_vs_sqrtS_recoJetDR10->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs recoJetDR12 sum
  TH2F* h_sqrtS_genJetDR12_eff_vs_sqrtS_recoJetDR12 = new TH2F("h_sqrtS_genJetDR12_vs_sqrtS_recoJetDR12","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_genJetDR12_eff_vs_sqrtS_recoJetDR12->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
  h_sqrtS_genJetDR12_eff_vs_sqrtS_recoJetDR12->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus the gentrue photon
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_genAll = new TH2F("h_sqrtS_e1e2_vs_sqrtS_genAll","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_genAll->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_genAll->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus the gentrue photon
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_recoAll = new TH2F("h_sqrtS_e1e2_vs_sqrtS_recoAll","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_recoAll->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_recoAll->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus iso photons
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv = new TH2F("h_sqrtS_e1e2_vs_sqrtS_gen_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus iso photons
  TH2F* h_sqrtS_e1e2_eff_vs_sqrtS_gen_inv = new TH2F("h_sqrtS_e1e2_vs_sqrtS_gen_inv","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_sqrtS_e1e2_eff_vs_sqrtS_gen_inv->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_sqrtS_e1e2_eff_vs_sqrtS_gen_inv->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");

  std::vector<TH2F*> hist_vec_HZ_2DHist;  
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_genTrue);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_genJetDR10);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_genJetDR12);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_reco_isoPh);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR10);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR12);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_genJetDR10_eff_vs_sqrtS_recoJetDR10);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_genJetDR12_eff_vs_sqrtS_recoJetDR12);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_genAll);//10
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_recoAll);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv);
  hist_vec_HZ_2DHist.push_back(h_sqrtS_e1e2_eff_vs_sqrtS_gen_inv);

 for(unsigned int i=0;i<hist_vec_HZ_2DHist.size();i++){
    hist_vec_HZ_2DHist[i]->Sumw2();
  }

 fill_HZ_histograms(file_CLIC_HZqq,hist_vec_HZ_parton,hist_vec_HZ_2DHist);


  TH1F* h_noBG_sqrtS_e1_e2_effective = new TH1F("h_noBG_sqrtS_e1_e2_effective","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_H_pt_sqrt_s_0_750 = new TH1F("h_noBG_H_pt_sqrt_s_0_750","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_H_pt_sqrt_s_750_2750 = new TH1F("h_noBG_H_pt_sqrt_s_750_2750","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_H_pt_sqrt_s_2750_3000 = new TH1F("h_noBG_H_pt_sqrt_s_2750_3000","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_Z_pt_sqrt_s_0_750 = new TH1F("h_noBG_Z_pt_sqrt_s_0_750","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_Z_pt_sqrt_s_750_2750 = new TH1F("h_noBG_Z_pt_sqrt_s_750_2750","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_Z_pt_sqrt_s_2750_3000 = new TH1F("h_noBG_Z_pt_sqrt_s_2750_3000","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_noBG_H_dalpha_bb_sqrt_s_0_750 = new TH1F("h_noBG_H_dalpha_bb_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_H_dalpha_bb_sqrt_s_750_2750 = new TH1F("h_noBG_H_dalpha_bb_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_H_dalpha_bb_sqrt_s_2750_3000 = new TH1F("h_noBG_H_dalpha_bb_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_Z_dalpha_qq_sqrt_s_0_750 = new TH1F("h_noBG_Z_dalpha_qq_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_Z_dalpha_qq_sqrt_s_750_2750 = new TH1F("h_noBG_Z_dalpha_qq_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_Z_dalpha_qq_sqrt_s_2750_3000 = new TH1F("h_noBG_Z_dalpha_qq_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_noBG_H_theta_b_fw_sqrt_s_0_750 = new TH1F("h_noBG_H_theta_b_fw_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_noBG_H_theta_b_fw_sqrt_s_750_2750 = new TH1F("h_noBG_H_theta_b_fw_sqrt_s_750_2750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_noBG_H_theta_b_fw_sqrt_s_2750_3000 = new TH1F("h_noBG_H_theta_b_fw_sqrt_s_0_2750_3000","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_noBG_Z_theta_q_fw_sqrt_s_0_750 = new TH1F("h_noBG_Z_theta_q_fw_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_noBG_Z_theta_q_fw_sqrt_s_750_2750 = new TH1F("h_noBG_Z_theta_q_fw_sqrt_s_750_2750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_noBG_Z_theta_q_fw_sqrt_s_2750_3000 = new TH1F("h_noBG_Z_theta_q_fw_sqrt_s_0_2750_3000","", n_bins_high, lim_theta_low,lim_theta_high);

  TH1F* h_noBG_dalpha_min_Z_q_H_b_sqrt_s_0_750 = new TH1F("h_noBG_dalpha_min_Z_q_H_b_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_dalpha_min_Z_q_H_b_sqrt_s_750_2750 = new TH1F("h_noBG_dalpha_min_Z_q_H_b_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_dalpha_min_Z_q_H_b_sqrt_s_2750_3000 = new TH1F("h_noBG_dalpha_min_Z_q_H_b_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_noBG_H_theta_sqrt_s_0_750 = new TH1F("h_noBG_H_theta_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_noBG_H_theta_sqrt_s_750_2750 = new TH1F("h_noBG_H_theta_sqrt_s_750_2750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_noBG_H_theta_sqrt_s_2750_3000 = new TH1F("h_noBG_H_theta_sqrt_s_0_2750_3000","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_noBG_Z_theta_sqrt_s_0_750 = new TH1F("h_noBG_Z_theta_sqrt_s_0_750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_noBG_Z_theta_sqrt_s_750_2750 = new TH1F("h_noBG_Z_theta_sqrt_s_750_2750","", n_bins_high, lim_theta_low,lim_theta_high);
  TH1F* h_noBG_Z_theta_sqrt_s_2750_3000 = new TH1F("h_noBG_Z_theta_sqrt_s_0_2750_3000","", n_bins_high, lim_theta_low,lim_theta_high);

  TH1F* h_noBG_dtheta_H_Z_sqrt_s_0_750 = new TH1F("h_noBG_dtheta_H_Z_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_dtheta_H_Z_sqrt_s_750_2750 = new TH1F("h_noBG_dtheta_H_Z_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_dtheta_H_Z_sqrt_s_2750_3000 = new TH1F("h_noBG_dtheta_H_Z_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_noBG_dphi_H_Z_sqrt_s_0_750 = new TH1F("h_noBG_dphi_H_Z_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_dphi_H_Z_sqrt_s_750_2750 = new TH1F("h_noBG_dphi_H_Z_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_dphi_H_Z_sqrt_s_2750_3000 = new TH1F("h_noBG_dphi_H_Z_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  TH1F* h_noBG_dalpha_H_Z_sqrt_s_0_750 = new TH1F("h_noBG_dalpha_H_Z_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_dalpha_H_Z_sqrt_s_750_2750 = new TH1F("h_noBG_dalpha_H_Z_sqrt_s_750_2750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_dalpha_H_Z_sqrt_s_2750_3000 = new TH1F("h_noBG_dalpha_H_Z_sqrt_s_0_2750_3000","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_noBG_sqrtS_gen_isoPh = new TH1F("h_noBG_sqrtS_gen_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_sqrtS_reco_isoPh = new TH1F("h_noBG_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_sqrtS_HZ_effective = new TH1F("h_noBG_sqrtS_HZ_effective","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_sqrtS_gen_isoPh_inv = new TH1F("h_noBG_sqrtS_gen_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_sqrtS_reco = new TH1F("h_noBG_sqrtS_reco","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_sqrtS_gen = new TH1F("h_noBG_sqrtS_gen","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_sqrtS_gen_inv = new TH1F("h_noBG_sqrtS_gen_inv","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_sqrtS_gen_truePh_inv = new TH1F("h_noBG_sqrtS_gen_truePh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_noBG_sqrtS_reco_isoPh_inv = new TH1F("h_noBG_sqrtS_reco_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_noBG_mass1_gj10_sqrtS_2500 = new TH1F("h_noBG_mass1_gj10_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_noBG_mass2_gj10_sqrtS_2500 = new TH1F("h_noBG_mass2_gj10_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_noBG_mass1_gj12_sqrtS_2500 = new TH1F("h_noBG_mass1_gj12_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_noBG_mass2_gj12_sqrtS_2500 = new TH1F("h_noBG_mass2_gj12_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);

  TH1F* h_noBG_mass1_rj10_sqrtS_2500 = new TH1F("h_noBG_mass1_rj10_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_noBG_mass2_rj10_sqrtS_2500 = new TH1F("h_noBG_mass2_rj10_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_noBG_mass1_rj12_sqrtS_2500 = new TH1F("h_noBG_mass1_rj12_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);
  TH1F* h_noBG_mass2_rj12_sqrtS_2500 = new TH1F("h_noBG_mass2_rj12_sqrtS_2500","", n_bins_high, lim_mass_low,lim_mass_high);

  TH1F* h_noBG_dPhi_mass1_gj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dTheta_mass1_gj10_H_sqrtS_2500 = new TH1F("h_noBG_dTheta_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_mass2_gj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_mass2_gj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dTheta_mass2_gj10_H_sqrtS_2500 = new TH1F("h_noBG_dTheta_mas21_gj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_noBG_dPhi_mass1_gj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dTheta_mass1_gj12_H_sqrtS_2500 = new TH1F("h_noBG_dTheta_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_mass2_gj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_mass2_gj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dTheta_mass2_gj12_H_sqrtS_2500 = new TH1F("h_noBG_dTheta_mass2_gj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_noBG_dPhi_mass1_rj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dTheta_mass1_rj10_H_sqrtS_2500 = new TH1F("h_noBG_dTheta_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_mass2_rj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_mass2_rj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dTheta_mass2_rj10_H_sqrtS_2500 = new TH1F("h_noBG_dTheta_mass2_rj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_noBG_dPhi_mass1_rj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dTheta_mass1_rj12_H_sqrtS_2500 = new TH1F("h_noBG_dTheta_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_mass2_rj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_mass2_rj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dTheta_mass2_rj12_H_sqrtS_2500 = new TH1F("h_noBG_dTheta_mass2_rj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_noBG_dAlpha_mass1_gj10_H_sqrtS_2500 = new TH1F("h_noBG_dAlpha_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dAlpha_mass2_gj10_Z_sqrtS_2500 = new TH1F("h_noBG_dAlpha_mass2_gj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_noBG_dAlpha_mass1_gj12_H_sqrtS_2500 = new TH1F("h_noBG_dAlpha_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dAlpha_mass2_gj12_Z_sqrtS_2500 = new TH1F("h_noBG_dAlpha_mass2_gj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_noBG_dAlpha_mass1_rj10_H_sqrtS_2500 = new TH1F("h_noBG_dAlpha_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dAlpha_mass2_rj10_Z_sqrtS_2500 = new TH1F("h_noBG_dAlpha_mass2_rj10_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_noBG_dAlpha_mass1_rj12_H_sqrtS_2500 = new TH1F("h_noBG_dAlpha_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dAlpha_mass2_rj12_Z_sqrtS_2500 = new TH1F("h_noBG_dAlpha_mass2_rj12_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_noBG_dPhi_MET_mass1_gj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET_mass1_gj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET_mass1_rj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET_mass1_rj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_noBG_genMET_sqrtS_2500 = new TH1F("h_noBG_genMET_sqrtS_2500","", n_bins_high, lim_MET_low,lim_MET_high);
  TH1F* h_noBG_recoMET_sqrtS_2500 = new TH1F("h_noBG_recoMET_sqrtS_2500","", n_bins_high, lim_MET_low,lim_MET_high);

  TH1F* h_noBG_dPhi_MET50_mass1_gj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET50_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET50_mass1_gj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET50_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET50_mass1_rj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET50_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET50_mass1_rj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET50_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_noBG_dPhi_MET100_mass1_gj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET100_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET100_mass1_gj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET100_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET100_mass1_rj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET100_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET100_mass1_rj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET100_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_noBG_dPhi_MET150_mass1_gj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET150_mass1_gj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET150_mass1_gj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET150_mass1_gj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET150_mass1_rj10_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET150_mass1_rj10_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_noBG_dPhi_MET150_mass1_rj12_H_sqrtS_2500 = new TH1F("h_noBG_dPhi_MET150_mass1_rj12_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_noBG_E_q1_over_Z = new TH1F("h_noBG_E_q1_over_Z","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_noBG_E_q2_over_Z = new TH1F("h_noBG_E_q2_over_Z","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_noBG_E_b1_over_H = new TH1F("h_noBG_E_b1_over_H","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_noBG_E_b2_over_H = new TH1F("h_noBG_E_b2_over_H","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);

  TH1F* h_noBG_E_s1_over_m2_gj10 = new TH1F("h_noBG_E_s1_over_m2_gj10","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_noBG_E_s2_over_m2_gj10 = new TH1F("h_noBG_E_s2_over_m2_gj10","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_noBG_E_s1_over_m1_gj10 = new TH1F("h_noBG_E_s1_over_m1_gj10","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_noBG_E_s2_over_m1_gj10 = new TH1F("h_noBG_E_s2_over_m1_gj10","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);

  TH1F* h_noBG_E_s1_over_m2_gj12 = new TH1F("h_noBG_E_s1_over_m2_gj12","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_noBG_E_s2_over_m2_gj12 = new TH1F("h_noBG_E_s2_over_m2_gj12","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_noBG_E_s1_over_m1_gj12 = new TH1F("h_noBG_E_s1_over_m1_gj12","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_noBG_E_s2_over_m1_gj12 = new TH1F("h_noBG_E_s2_over_m1_gj12","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);

  TH1F* h_noBG_E_s1_over_m2_rj10 = new TH1F("h_noBG_E_s1_over_m2_rj10","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_noBG_E_s2_over_m2_rj10 = new TH1F("h_noBG_E_s2_over_m2_rj10","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_noBG_E_s1_over_m1_rj10 = new TH1F("h_noBG_E_s1_over_m1_rj10","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_noBG_E_s2_over_m1_rj10 = new TH1F("h_noBG_E_s2_over_m1_rj10","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);

  TH1F* h_noBG_E_s1_over_m2_rj12 = new TH1F("h_noBG_E_s1_over_m2_rj12","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_noBG_E_s2_over_m2_rj12 = new TH1F("h_noBG_E_s2_over_m2_rj12","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);
  TH1F* h_noBG_E_s1_over_m1_rj12 = new TH1F("h_noBG_E_s1_over_m1_rj12","", n_bins_high, lim_sj1ratio_low,lim_sj1ratio_high);
  TH1F* h_noBG_E_s2_over_m1_rj12 = new TH1F("h_noBG_E_s2_over_m1_rj12","", n_bins_high, lim_sj2ratio_low,lim_sj2ratio_high);


  std::vector<TH1F*> hist_noBG_vec_HZ_parton;
  hist_noBG_vec_HZ_parton.push_back(h_noBG_sqrtS_e1_e2_effective);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_pt_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_pt_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_pt_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_pt_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_pt_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_pt_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_dalpha_bb_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_dalpha_bb_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_dalpha_bb_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_dalpha_qq_sqrt_s_0_750);//10
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_dalpha_qq_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_dalpha_qq_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_theta_b_fw_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_theta_b_fw_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_theta_b_fw_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_theta_q_fw_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_theta_q_fw_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_theta_q_fw_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dalpha_min_Z_q_H_b_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dalpha_min_Z_q_H_b_sqrt_s_750_2750);//20
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dalpha_min_Z_q_H_b_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_theta_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_theta_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_H_theta_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_theta_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_theta_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_Z_theta_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dtheta_H_Z_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dtheta_H_Z_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dtheta_H_Z_sqrt_s_2750_3000);//30
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dphi_H_Z_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dphi_H_Z_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dphi_H_Z_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dalpha_H_Z_sqrt_s_0_750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dalpha_H_Z_sqrt_s_750_2750);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dalpha_H_Z_sqrt_s_2750_3000);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_sqrtS_gen_isoPh);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_sqrtS_reco_isoPh);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_sqrtS_HZ_effective);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_sqrtS_gen_isoPh_inv);//40
  hist_noBG_vec_HZ_parton.push_back(h_noBG_sqrtS_reco);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_sqrtS_gen);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_sqrtS_gen_inv);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_sqrtS_gen_truePh_inv);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_sqrtS_reco_isoPh_inv);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_mass1_gj10_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_mass2_gj10_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_mass1_gj12_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_mass2_gj12_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_mass1_rj10_sqrtS_2500);//50
  hist_noBG_vec_HZ_parton.push_back(h_noBG_mass2_rj10_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_mass1_rj12_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_mass2_rj12_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_mass1_gj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dTheta_mass1_gj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_mass2_gj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dTheta_mass2_gj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_mass1_gj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dTheta_mass1_gj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_mass2_gj12_H_sqrtS_2500);//60
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dTheta_mass2_gj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_mass1_rj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dTheta_mass1_rj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_mass2_rj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dTheta_mass2_rj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_mass1_rj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dTheta_mass1_rj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_mass2_rj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dTheta_mass2_rj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dAlpha_mass1_gj10_H_sqrtS_2500);//70
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dAlpha_mass2_gj10_Z_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dAlpha_mass1_gj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dAlpha_mass2_gj12_Z_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dAlpha_mass1_rj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dAlpha_mass2_rj10_Z_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dAlpha_mass1_rj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dAlpha_mass2_rj12_Z_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET_mass1_gj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET_mass1_gj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET_mass1_rj10_H_sqrtS_2500);//80
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET_mass1_rj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_genMET_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_recoMET_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET50_mass1_gj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET50_mass1_gj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET50_mass1_rj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET50_mass1_rj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET100_mass1_gj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET100_mass1_gj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET100_mass1_rj10_H_sqrtS_2500);//90
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET100_mass1_rj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET150_mass1_gj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET150_mass1_gj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET150_mass1_rj10_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_dPhi_MET150_mass1_rj12_H_sqrtS_2500);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_q1_over_Z);//96
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_q2_over_Z);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_b1_over_H);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_b2_over_H);  
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s1_over_m2_gj10);//100
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s2_over_m2_gj10);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s1_over_m1_gj10);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s2_over_m1_gj10);  
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s1_over_m2_gj12);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s2_over_m2_gj12);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s1_over_m1_gj12);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s2_over_m1_gj12);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s1_over_m2_rj10);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s2_over_m2_rj10);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s1_over_m1_rj10);//100
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s2_over_m1_rj10);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s1_over_m2_rj12);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s2_over_m2_rj12);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s1_over_m1_rj12);
  hist_noBG_vec_HZ_parton.push_back(h_noBG_E_s2_over_m1_rj12);//115


  for(unsigned int i=0;i<hist_noBG_vec_HZ_parton.size();i++){
    hist_noBG_vec_HZ_parton[i]->Sumw2();
  }

 //parton vs gen true --> here we calculat all true minus the gentrue photon
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genTrue = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_genTrue","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genTrue->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genTrue->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus iso photons
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_gen_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
  //parton vs genJetDR10 sum
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genJetDR10 = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_genJetDR10","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genJetDR10->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genJetDR10->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
  //parton vs genJetDR12 sum
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genJetDR12 = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_genJetDR12","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genJetDR12->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genJetDR12->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");

 //parton vs reco true --> here we calculat all true minus iso photons
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_reco_isoPh = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_reco_isoPh->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_reco_isoPh->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs recoJetDR10 sum
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR10 = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_recoJetDR10","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR10->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR10->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs recoJetDR12 sum
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR12 = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_recoJetDR12","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR12->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR12->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");

 //parton vs reco true --> here we calculat all true minus iso photons
  TH2F* h_noBG_sqrtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh = new TH2F("h_noBG_sqrtS_gen_isoPh_noBG_vs_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
  h_noBG_sqrtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs recoJetDR10 sum
  TH2F* h_noBG_sqrtS_genJetDR10_eff_vs_sqrtS_recoJetDR10 = new TH2F("h_noBG_sqrtS_genJetDR10_vs_sqrtS_recoJetDR10","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_genJetDR10_eff_vs_sqrtS_recoJetDR10->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
  h_noBG_sqrtS_genJetDR10_eff_vs_sqrtS_recoJetDR10->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs recoJetDR12 sum
  TH2F* h_noBG_sqrtS_genJetDR12_eff_vs_sqrtS_recoJetDR12 = new TH2F("h_noBG_sqrtS_genJetDR12_vs_sqrtS_recoJetDR12","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_genJetDR12_eff_vs_sqrtS_recoJetDR12->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
  h_noBG_sqrtS_genJetDR12_eff_vs_sqrtS_recoJetDR12->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus the gentrue photon
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genAll = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_genAll","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genAll->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genAll->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus the gentrue photon
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoAll = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_recoAll","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoAll->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoAll->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus iso photons
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_gen_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
 //parton vs gen true --> here we calculat all true minus iso photons
  TH2F* h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_inv = new TH2F("h_noBG_sqrtS_e1e2_vs_sqrtS_gen_inv","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_inv->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
  h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_inv->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");


  std::vector<TH2F*> hist_noBG_vec_HZ_2DHist;  
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genTrue);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genJetDR10);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genJetDR12);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_reco_isoPh);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR10);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoJetDR12);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_genJetDR10_eff_vs_sqrtS_recoJetDR10);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_genJetDR12_eff_vs_sqrtS_recoJetDR12);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_genAll);//10
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_recoAll);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv);
  hist_noBG_vec_HZ_2DHist.push_back(h_noBG_sqrtS_e1e2_eff_vs_sqrtS_gen_inv);
 for(unsigned int i=0;i<hist_noBG_vec_HZ_2DHist.size();i++){
    hist_noBG_vec_HZ_2DHist[i]->Sumw2();
  }

 fill_HZ_histograms(file_CLIC_HZqq_noBGG,hist_noBG_vec_HZ_parton,hist_noBG_vec_HZ_2DHist);




  file_histogram->Write();
  file_histogram->Close();

}
