#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TColor.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TEfficiency.h"
#include "TList.h"
#include "TGraphAsymmErrors.h"

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
   gStyle->SetTitleSize(0.07,"xyz");
   gStyle->SetLabelSize(0.06,"xyz");
   /* title offset: distance between given text and axis, here x,y,z*/
   gStyle->SetLabelOffset(0.015,"xyz");
   gStyle->SetTitleOffset(1.2,"yz"); //equivalent to: gStyle->SetTitleYOffset(1.2);
   gStyle->SetTitleOffset(1.0,"x");



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
   gStyle->SetMarkerSize(1.5);
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
   gStyle->SetFuncColor(kRed);
   gStyle->SetLabelColor(kBlack,"xyz");

   //set the margins
   gStyle->SetPadBottomMargin(0.18);
   gStyle->SetPadTopMargin(0.08);
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

TCanvas* setUpperCanvas(const char* canvas_name) {
    TCanvas* c1= new TCanvas(canvas_name,canvas_name,10,50,600,500);
    c1->cd();
    gPad->SetTopMargin(0.06);
    return c1;
}


void fillHistograms(TFile* file_training, std::vector<TH1F*> t_hist_vec, std::vector<TEfficiency*> t_TEff_vec){

  TTree* tree_train= (TTree*)file_training->Get("showerData");

  vector<float> *true_Energy=0;
  vector<float> *true_Px=0;
  vector<float> *true_Py=0;
  vector<float> *true_Pz=0;
  vector<int> *true_PDGID=0;
  vector<int> *true_Status=0;

  vector<int> *true_all_6_deg_index=0;
  vector<float> *true_all_6_deg_E=0;
  vector<float> *true_all_6_deg_angle=0;
  vector<int> *true_all_6_deg_PDGID=0;//take out neutrinos, i.e. 12,14,16


  vector<float> *reco_Energy=0;
  vector<float> *reco_Px=0;
  vector<float> *reco_Py=0;
  vector<float> *reco_Pz=0;
  vector<int> *reco_PDGID=0;

  vector<int> *reco_all_6_deg_index=0;
  vector<float> *reco_all_6_deg_E=0;
  vector<float> *reco_all_6_deg_angle=0;
  vector<int> *reco_all_6_deg_PDGID=0;

  tree_train->SetBranchAddress("true_Energy", &true_Energy);
  tree_train->SetBranchAddress("true_Px", &true_Px);
  tree_train->SetBranchAddress("true_Py", &true_Py);
  tree_train->SetBranchAddress("true_Pz", &true_Pz);
  tree_train->SetBranchAddress("true_PDGID", &true_PDGID);
  tree_train->SetBranchAddress("true_GenStatus", &true_Status);

  tree_train->SetBranchAddress("true_all_6_deg_index", &true_all_6_deg_index);
  tree_train->SetBranchAddress("true_all_6_deg_E", &true_all_6_deg_E);
  tree_train->SetBranchAddress("true_all_6_deg_angle", &true_all_6_deg_angle);
  tree_train->SetBranchAddress("true_all_6_deg_PDG", &true_all_6_deg_PDGID);
 
  tree_train->SetBranchAddress("reco_Energy", &reco_Energy);
  tree_train->SetBranchAddress("reco_Px", &reco_Px);
  tree_train->SetBranchAddress("reco_Py", &reco_Py);
  tree_train->SetBranchAddress("reco_Pz", &reco_Pz);
  tree_train->SetBranchAddress("reco_PDGID", &reco_PDGID);

  tree_train->SetBranchAddress("reco_all_6_deg_index", &reco_all_6_deg_index);
  tree_train->SetBranchAddress("reco_all_6_deg_E", &reco_all_6_deg_E);
  tree_train->SetBranchAddress("reco_all_6_deg_angle", &reco_all_6_deg_angle);
  tree_train->SetBranchAddress("reco_all_6_deg_PDG", &reco_all_6_deg_PDGID);

  //for isolation values, check 2 degrees as upper limit

  //histogram vector--> GenQuantities
  //[0]:true electron-energy
  //[1]:true electron absolute isolation -->use in both cases true visible particles within 6 degrees (i.e. no additional cut)
  //[2]:true electron relative isolation
  //[3]:true electron absolute isolation -->use in both cases true visible particles within 6 degrees (i.e. no additional cut)
  //removing photons within angle of 0.010
  //[4]:true electron relative isolation
  //removing photons within angle of 0.010
  //[5]:true electron absolute isolation -->use in both cases true visible particles within 6 degrees (i.e. no additional cut)
  //removing photons within angle of 1.5 degrees
  //[6]:true electron relative isolation
  //removing photons within angle of 1.5 degrees

  //[7]:true muon-energy
  //[8]:true muon absolute isolation -->use in both cases true visible particles within 6 degrees (i.e. no additional cut)
  //[9]:true muon relative isolation
  //[10]:true muon absolute isolation -->use in both cases true visible particles within 6 degrees (i.e. no additional cut)
  //removing photons within angle of 0.010
  //[11]:true muon relative isolation
  //removing photons within angle of 0.010
  //[12]:true muon absolute isolation -->use in both cases true visible particles within 6 degrees (i.e. no additional cut)
  //removing photons within angle of 1.5 degrees
  //[13]:true muon relative isolation
  //removing photons within angle of 1.5 degrees

  //[14]: for dilepton decay, DR between both true particles 

  //[15]: DR between true and closest (spatial) reco Electron candidate --> upper limit for search 6 degrees
  //[16]: E_reco/E_true for closest (spatial) reco Electron candidate --> upper limit for search 6 degrees --> i.e response
  //[17]: E reco of closest (spatial) reco Electron candidate --> only if candidate is found within 6 degree search radius
  //[18]: DR between true and closest (energy) reco Electron candidate --> upper limit for search 2 degrees
  //[19]: E_reco/E_true for reco Electron candidate (closest in energy, but restrict search radius to within 2 degrees)
  //[20]: E_reco for reco Electron candidate (closest in energy) but restrict search radius to within 2 degrees)

  //[21]: DR between true and closest (spatial) reco Muon candidate --> upper limit for search 6 degrees
  //[22]: E_reco/E_true for closest (spatial) reco Muon candidate --> upper limit for search 6 degrees -> i.e. response
  //[23]: E reco of closest (spatial) reco Muon candidate --> only if candidate is found within 6 degree search radius
  //[24]: DR between true and closest (energy) reco Muon candidate --> upper limit for search 2 degrees
  //[25]: E_reco/E_true for reco Muon candidate (closest in energy, but restrict search radius to within 2 degrees)
  //[26]: E_reco for reco Muon candidate (closest in energy) but restrict search radius to within 2 degrees)

  //[27]: absolute isolation closest (spatial) reco Electron candidate --> check within 6 degrees
  //[28]: relative isolation closest (spatial) reco Electron candidate --> check within 6 degrees
  //[29]: absolute isolation excluding photons within 0.010 closest (spatial) reco Electron candidate --> check within 6 degrees
  //[30]: relative isolation excluding photons within 0.010 closest (spatial) reco Electron candidate --> check within 6 degrees
  //[31]: absolute isolation excluding photons within 1.5 degrees of closest (spatial) reco Electron candidate --> check within 6 degrees
  //[32]: relative isolation excluding photons within 1.5 degrees of closest (spatial) reco Electron candidate --> check within 6 degrees

  //[33]: absolute isolation closest (energy) reco Electron candidate --> check within 6 degrees
  //[34]: relative isolation closest (energy) reco Electron candidate --> check within 6 degrees
  //[35]: absolute isolation excluding photons within 0.010 closest (energy) reco Electron candidate --> check within 6 degrees
  //[36]: relative isolation excluding photons within 0.010 closest (energy) reco Electron candidate --> check within 6 degrees
  //[37]: absolute isolation excluding photons within 1.5 degree degree closest (energy) reco Electron candidate --> check within 6 degrees
  //[38]: relative isolation excluding photons within 1.5 degree degree closest (energy) reco Electron candidate --> check within 6 degrees

  //[39]: absolute isolation closest (spatial) reco Muon candidate --> check within 6 degrees
  //[40]: relative isolation closest (spatial) reco Muon candidate --> check within 6 degrees
  //[41]: absolute isolation excluding photons closest (spatial) reco Muon candidate --> check within 6 degrees
  //[42]: relative isolation excluding photons closest (spatial) reco Muon candidate --> check within 6 degrees
  //[43]: absolute isolation excluding photons within 0.010 closest (spatial) reco Muon candidate --> check within 6 degrees
  //[44]: relative isolation excluding photons within 0.010 closest (spatial) reco Muon candidate --> check within 6 degrees

  //[45]: absolute isolation closest (energy) reco Muon candidate --> check within 6 degrees
  //[46]: relative isolation closest (energy) reco Muon candidate --> check within 6 degrees
  //[47]: absolute isolation excluding photons within 0.010 closest (energy) reco Muon candidate --> check within 6 degrees
  //[48]: relative isolation excluding photons within 0.010 closest (energy) reco Muon candidate --> check within 6 degrees
  //[49]: absolute isolation excluding photons within 1.5 degree degree closest (energy) reco Muon candidate --> check within 6 degrees
  //[50]: relative isolation excluding photons within 1.5 degree degree closest (energy) reco Muon candidate --> check within 6 degrees

  //[51]: fake reco Electron candidate energy --> 
  //[52]: absolute isolation fake reco Electron candidate --> check within 6 degrees
  //[53]: relative isolation fake reco Electron candidate --> check within 6 degrees
  //[54]: absolute isolation excluding photons within 0.010 fake reco Electron candidate --> check within 6 degrees
  //[55]: relative isolation excluding photons within 0.010 fake reco Electron candidate --> check within 6 degrees
  //[56]: absolute isolation excluding photons within 1.5 degrees of fake reco Electron candidate --> check within 6 degrees
  //[57]: relative isolation excluding photons within 1.5 degrees of closest fake reco Electron candidate --> check within 6 degrees

  //[58]: fake reco Muon candidate energy --> 
  //[59]: absolute isolation fake reco Muon candidate --> check within 6 degrees
  //[60]: relative isolation fake reco Muon candidate --> check within 6 degrees
  //[61]: absolute isolation excluding photons within 0.010 fake reco Electron candidate --> check within 6 degrees
  //[62]: relative isolation excluding photons within 0.010 fake reco Muon candidate --> check within 6 degrees
  //[63]: absolute isolation excluding photons within 1.5 degrees of fake reco Muon candidate --> check within 6 degrees
  //[64]: relative isolation excluding photons within 1.5 degrees of closest fake reco Muon candidate --> check within 6 degrees

  //[65]: gen Electron energy, |cosTheta|<0.95 candidate energy --> 
  //[66]: reco candidate Electron energy, angular match of candidate within 1 degree, |cosTheta|<0.95, relative isolation 0.01
  //[67]: reco candidate Electron energy, angular match of candidate within 1 degree, |cosTheta|<0.95, relative isolation 0.02

  //[68]: gen Muon energy, |cosTheta|<0.95 candidate energy --> 
  //[69]: reco candidate Muon energy, angular match of candidate within 1 degree, |cosTheta|<0.95, relative isolation 0.01
  //[70]: reco candidate Muon energy, angular match of candidate within 1 degree, |cosTheta|<0.95, relative isolation 0.02

  float ph_angle_small=0.010;
  float ph_angle_large=1.5;
  float iso_angle_max=10.;
  float E_min_particle=0.;
  float E_match_particle=0.05;
  float Pt_match_particle=0.05;
  float dangle_match_value_max=6.;

  for (unsigned int i_entry=0;i_entry<tree_train->GetEntries();i_entry++){
    tree_train->GetEntry(i_entry);
    if(i_entry%3000==0){
      std::cout<<"i_entry "<<i_entry<<std::endl;
    }
    float true_energy=-1;
    int pdg_true1=0;
    int pdg_true2=0;
    for(unsigned int i_true=0;i_true<true_Energy->size();i_true++){
      if((*true_Energy)[i_true]<E_min_particle){
	continue;
      }
      if(i_true==0){
	pdg_true1=(*true_PDGID)[i_true];
      }else{
	pdg_true2=(*true_PDGID)[i_true];
      }
      if((*true_Status)[i_true]!=1){
	//maybe in future taus will be tackled here
	continue;
      }
      if(abs((*true_PDGID)[i_true])==11){
	//Electron
	t_hist_vec[0]->Fill((*true_Energy)[i_true]);
	float E_sum_iso=0;
	float E_sum_iso_noPh_0_010=0;
	float E_sum_iso_noPh_1_5_deg=0;
	for(unsigned int i_true_iso=0;i_true_iso<true_all_6_deg_index->size();i_true_iso++){
	  if(abs((*true_all_6_deg_PDGID)[i_true_iso])==12 || abs((*true_all_6_deg_PDGID)[i_true_iso])==14 || abs((*true_all_6_deg_PDGID)[i_true_iso])==16){
	    continue;
	  }
	  //preselection is in fact 6, so the last cut has no effect 
	  if( ((*true_all_6_deg_index)[i_true_iso]==i_true) && ((*true_all_6_deg_angle)[i_true_iso]*TMath::RadToDeg())<iso_angle_max ){
	    E_sum_iso+=(*true_all_6_deg_E)[i_true_iso];
	    //around 0.5729 degrees
	    if((*true_all_6_deg_PDGID)[i_true_iso]==22 && (*true_all_6_deg_angle)[i_true_iso]<ph_angle_small){
	      continue;
	    }
	    E_sum_iso_noPh_0_010+=(*true_all_6_deg_E)[i_true_iso];
	    if((*true_all_6_deg_PDGID)[i_true_iso]==22 && ((*true_all_6_deg_angle)[i_true_iso]*TMath::RadToDeg())<ph_angle_large){
	      continue;
	    }
	    E_sum_iso_noPh_1_5_deg+=(*true_all_6_deg_E)[i_true_iso];
	  }
	}//Electron iso loop finished
	t_hist_vec[1]->Fill(E_sum_iso);
	t_hist_vec[2]->Fill(E_sum_iso/(*true_Energy)[i_true]);
	t_hist_vec[3]->Fill(E_sum_iso_noPh_0_010);
	t_hist_vec[4]->Fill(E_sum_iso_noPh_0_010/(*true_Energy)[i_true]);
	t_hist_vec[5]->Fill(E_sum_iso_noPh_1_5_deg);
	t_hist_vec[6]->Fill(E_sum_iso_noPh_1_5_deg/(*true_Energy)[i_true]);
      }else if(abs((*true_PDGID)[i_true])==13){
	//Muons
	t_hist_vec[7]->Fill((*true_Energy)[i_true]);
	float E_sum_iso=0;
	float E_sum_iso_noPh_0_010=0;
	float E_sum_iso_noPh_1_5_deg=0;
	for(unsigned int i_true_iso=0;i_true_iso<true_all_6_deg_index->size();i_true_iso++){
	  if(abs((*true_all_6_deg_PDGID)[i_true_iso])==12 || abs((*true_all_6_deg_PDGID)[i_true_iso])==14 || abs((*true_all_6_deg_PDGID)[i_true_iso])==16){
	    continue;
	  }
	  //preselection is in fact 6, so the last cut has no effect 
	  if( ((*true_all_6_deg_index)[i_true_iso]==i_true) && ((*true_all_6_deg_angle)[i_true_iso]*TMath::RadToDeg())<iso_angle_max ){
	    E_sum_iso+=(*true_all_6_deg_E)[i_true_iso];
	    //around 0.5729 degrees
	    if((*true_all_6_deg_PDGID)[i_true_iso]==22 && (*true_all_6_deg_angle)[i_true_iso]<ph_angle_small){
	      continue;
	    }
	    E_sum_iso_noPh_0_010+=(*true_all_6_deg_E)[i_true_iso];
	    if((*true_all_6_deg_PDGID)[i_true_iso]==22 && ((*true_all_6_deg_angle)[i_true_iso]*TMath::RadToDeg())<ph_angle_large){
	      continue;
	    }
	    E_sum_iso_noPh_1_5_deg+=(*true_all_6_deg_E)[i_true_iso];
	  }
	}//Electron iso loop finished
	t_hist_vec[8]->Fill(E_sum_iso);
	t_hist_vec[9]->Fill(E_sum_iso/(*true_Energy)[i_true]);
	t_hist_vec[10]->Fill(E_sum_iso_noPh_0_010);
	t_hist_vec[11]->Fill(E_sum_iso_noPh_0_010/(*true_Energy)[i_true]);
	t_hist_vec[12]->Fill(E_sum_iso_noPh_1_5_deg);
	t_hist_vec[13]->Fill(E_sum_iso_noPh_1_5_deg/(*true_Energy)[i_true]);
      }
    }
    TLorentzVector true1(0,0,0,0);
    TLorentzVector true2(0,0,0,0);
    if(true_Energy->size()>0){
      true1.SetPxPyPzE((*true_Px)[0],(*true_Py)[0],(*true_Pz)[0],(*true_Energy)[0]);
      if(true_Energy->size()>1){
	true2.SetPxPyPzE((*true_Px)[1],(*true_Py)[1],(*true_Pz)[1],(*true_Energy)[1]);
	t_hist_vec[14]->Fill(true1.Angle(true2.Vect())*TMath::RadToDeg());
      }
    }
    int ind_reco1_Ang=-1;
    int ind_reco1_E_deg_2=-1;
    int ind_reco2_Ang=-1;
    int ind_reco2_E_deg_2=-1;
    float DAngle_min1=dangle_match_value_max;
    float DAngle_min2=dangle_match_value_max;

    float rel_ratio_E1=2.;
    float rel_ratio_E2=2.;

    //don't do threshold cut here --> allow resolution and bremstrahlung effects to not lead to artificial inefficiencies
    for(unsigned int i_reco=0;i_reco<reco_Energy->size();i_reco++){
      if(pdg_true1!=0){
	TLorentzVector tempReco(0,0,0,0);
	tempReco.SetPxPyPzE((*reco_Px)[i_reco],(*reco_Py)[i_reco],(*reco_Pz)[i_reco],(*reco_Energy)[i_reco]);
	if(abs((*reco_PDGID)[i_reco])==abs(pdg_true1) || abs((*reco_PDGID)[i_reco])==abs(pdg_true2)){
	  if(abs((*reco_PDGID)[i_reco])==abs(pdg_true1)){
	    //check deviations to true lepton 1
	    float angle1=true1.Angle(tempReco.Vect());
	    if((angle1*TMath::RadToDeg())<DAngle_min1){
	      DAngle_min1=angle1;
	      ind_reco1_Ang=i_reco;
	    }
	    if((angle1*TMath::RadToDeg())<dangle_match_value_max){
	      float delta_E1_rel=fabs(true1.Energy()-(*reco_Energy)[i_reco])/true1.Energy();
	      if(delta_E1_rel<rel_ratio_E1){
		rel_ratio_E1=delta_E1_rel;
		ind_reco1_E_deg_2=i_reco;
	      }
	    }
	  }
	  if(abs((*reco_PDGID)[i_reco])==abs(pdg_true2)){
	    //check deviations to true lepton 2
	    float angle2=true2.Angle(tempReco.Vect());
	    if((angle2*TMath::RadToDeg())<DAngle_min2){
	      DAngle_min2=angle2;
	      ind_reco2_Ang=i_reco;
	    }
	    if((angle2*TMath::RadToDeg())<dangle_match_value_max){
	      float delta_E2_rel=fabs(true2.Energy()-(*reco_Energy)[i_reco])/true2.Energy();
	      if(delta_E2_rel<rel_ratio_E2){
		rel_ratio_E2=delta_E2_rel;
		ind_reco2_E_deg_2=i_reco;
	      }
	    }
	  }
	}
      }
    }
    TLorentzVector reco1_angle(0,0,0,0);
    float E_sum_iso_rec1_angle=0;
    float E_sum_iso_rec1_angle_noph_0_010=0;
    float E_sum_iso_rec1_angle_noph_1_5_deg=0;
    TLorentzVector reco1_E_ang_deg_2(0,0,0,0);
    float E_sum_iso_rec1_E_deg2=0;
    float E_sum_iso_rec1_E_deg2_noph_0_010=0;
    float E_sum_iso_rec1_E_deg2_noph_1_5_deg=0;
    TLorentzVector reco2_angle(0,0,0,0);
    float E_sum_iso_rec2_angle=0;
    float E_sum_iso_rec2_angle_noph_0_010=0;
    float E_sum_iso_rec2_angle_noph_1_5_deg=0;
    TLorentzVector reco2_E_ang_deg_2(0,0,0,0);
    float E_sum_iso_rec2_E_deg2=0;
    float E_sum_iso_rec2_E_deg2_noph_0_010=0;
    float E_sum_iso_rec2_E_deg2_noph_1_5_deg=0;

    if(ind_reco1_Ang!=-1){
      reco1_angle.SetPxPyPzE((*reco_Px)[ind_reco1_Ang],(*reco_Py)[ind_reco1_Ang],(*reco_Pz)[ind_reco1_Ang],(*reco_Energy)[ind_reco1_Ang]);
      for(unsigned int i_rec_iso=0;i_rec_iso<reco_all_6_deg_index->size();i_rec_iso++){
	//preselection is in fact 6, so the last cut has no effect 
	if( ((*reco_all_6_deg_index)[i_rec_iso]==ind_reco1_Ang) && ((*reco_all_6_deg_angle)[i_rec_iso]*TMath::RadToDeg())<iso_angle_max ){
	  E_sum_iso_rec1_angle+=(*reco_all_6_deg_E)[i_rec_iso];
	  //around 0.5729 degrees
	  if((*reco_all_6_deg_PDGID)[i_rec_iso]==22 && (*reco_all_6_deg_angle)[i_rec_iso]<ph_angle_small){
	    continue;
	  }
	  E_sum_iso_rec1_angle_noph_0_010+=(*reco_all_6_deg_E)[i_rec_iso];
	  //additional photons cut out --> Electron check
	  if((*reco_all_6_deg_PDGID)[i_rec_iso]==22 && ((*reco_all_6_deg_angle)[i_rec_iso]*TMath::RadToDeg())<ph_angle_large){
	    continue;
	  }
	  E_sum_iso_rec1_angle_noph_1_5_deg+=(*reco_all_6_deg_E)[i_rec_iso];
	}
      }
    }
    if(ind_reco1_E_deg_2!=-1){
      reco1_E_ang_deg_2.SetPxPyPzE((*reco_Px)[ind_reco1_E_deg_2],(*reco_Py)[ind_reco1_E_deg_2],(*reco_Pz)[ind_reco1_E_deg_2],(*reco_Energy)[ind_reco1_E_deg_2]);
      for(unsigned int i_rec_iso=0;i_rec_iso<reco_all_6_deg_index->size();i_rec_iso++){
	//preselection is in fact 6, so the last cut has no effect 
	if( ((*reco_all_6_deg_index)[i_rec_iso]==ind_reco1_E_deg_2) && ((*reco_all_6_deg_angle)[i_rec_iso]*TMath::RadToDeg())<iso_angle_max ){
	  E_sum_iso_rec1_E_deg2+=(*reco_all_6_deg_E)[i_rec_iso];
	  //around 0.5729 degrees
	  if((*reco_all_6_deg_PDGID)[i_rec_iso]==22 && (*reco_all_6_deg_angle)[i_rec_iso]<ph_angle_small){
	    continue;
	  }
	  E_sum_iso_rec1_E_deg2_noph_0_010+=(*reco_all_6_deg_E)[i_rec_iso];
	  //additional photons cut out --> Electron check
	  if((*reco_all_6_deg_PDGID)[i_rec_iso]==22 && ((*reco_all_6_deg_angle)[i_rec_iso]*TMath::RadToDeg())<ph_angle_large){
	    continue;
	  }
	  E_sum_iso_rec1_E_deg2_noph_1_5_deg=(*reco_all_6_deg_E)[i_rec_iso];
	}
      }
    }



    if(ind_reco2_Ang!=-1){
      reco2_angle.SetPxPyPzE((*reco_Px)[ind_reco2_Ang],(*reco_Py)[ind_reco2_Ang],(*reco_Pz)[ind_reco2_Ang],(*reco_Energy)[ind_reco2_Ang]);
      for(unsigned int i_rec_iso=0;i_rec_iso<reco_all_6_deg_index->size();i_rec_iso++){
	//preselection is in fact 6, so the last cut has no effect 
	if( ((*reco_all_6_deg_index)[i_rec_iso]==ind_reco2_Ang) && ((*reco_all_6_deg_angle)[i_rec_iso]*TMath::RadToDeg())<iso_angle_max ){
	  E_sum_iso_rec2_angle+=(*reco_all_6_deg_E)[i_rec_iso];
	  //around 0.5729 degrees
	  if((*reco_all_6_deg_PDGID)[i_rec_iso]==22 && (*reco_all_6_deg_angle)[i_rec_iso]<ph_angle_small){
	    continue;
	  }
	  E_sum_iso_rec2_angle_noph_0_010+=(*reco_all_6_deg_E)[i_rec_iso];
	  //additional photons cut out --> Electron check
	  if((*reco_all_6_deg_PDGID)[i_rec_iso]==22 && ((*reco_all_6_deg_angle)[i_rec_iso]*TMath::RadToDeg())<ph_angle_large){
	    continue;
	  }
	  E_sum_iso_rec2_angle_noph_1_5_deg+=(*reco_all_6_deg_E)[i_rec_iso];
	}
      }
    }
    if(ind_reco2_E_deg_2!=-1){
      reco2_E_ang_deg_2.SetPxPyPzE((*reco_Px)[ind_reco2_E_deg_2],(*reco_Py)[ind_reco2_E_deg_2],(*reco_Pz)[ind_reco2_E_deg_2],(*reco_Energy)[ind_reco2_E_deg_2]);
      for(unsigned int i_rec_iso=0;i_rec_iso<reco_all_6_deg_index->size();i_rec_iso++){
	//preselection is in fact 6, so the last cut has no effect 
	if( ((*reco_all_6_deg_index)[i_rec_iso]==ind_reco2_E_deg_2) && ((*reco_all_6_deg_angle)[i_rec_iso]*TMath::RadToDeg())<iso_angle_max ){
	  E_sum_iso_rec2_E_deg2+=(*reco_all_6_deg_E)[i_rec_iso];
	  //around 0.5729 degrees
	  if((*reco_all_6_deg_PDGID)[i_rec_iso]==22 && (*reco_all_6_deg_angle)[i_rec_iso]<ph_angle_small){
	    continue;
	  }
	  E_sum_iso_rec2_E_deg2_noph_0_010+=(*reco_all_6_deg_E)[i_rec_iso];
	  //additional photons cut out --> Electron check
	  if((*reco_all_6_deg_PDGID)[i_rec_iso]==22 && ((*reco_all_6_deg_angle)[i_rec_iso]*TMath::RadToDeg())<ph_angle_large){
	    continue;
	  }
	  E_sum_iso_rec2_E_deg2_noph_1_5_deg=(*reco_all_6_deg_E)[i_rec_iso];
	}
      }
    }
    //interesting leptons have been defined at this point, now check for isolation values 
    if(abs(pdg_true1)==11){
      //leading lepton is Electron --> check performed while indices are searched for
      bool pass_teff_el_deg1_and_relIso=false;
      bool pass_teff_el_deg1=false;
      bool pass_teff_el_deg1_EMatch=false;
      bool pass_teff_el_deg1_PtMatch=false;
      //std::cout<<"true electron"<<std::endl;
      if(ind_reco1_Ang!=-1){
	if( (true1.Angle(reco1_angle.Vect())*TMath::RadToDeg())<1.0){
	//if( (true1.DeltaPhi(reco1_angle))<0.02 && fabs(true1.Theta()-reco1_angle.Theta())<0.01 && fabs((true1.Pt()-reco1_angle.Pt())/true1.Pt())<0.5){
	  pass_teff_el_deg1=true;
	  if(fabs(true1.E()-reco1_angle.E())<(E_match_particle*true1.E())){
	    pass_teff_el_deg1_EMatch=true;
	  }
	  if(fabs(true1.Pt()-reco1_angle.Pt())<(Pt_match_particle*true1.Pt())){
	    pass_teff_el_deg1_PtMatch=true;
	  }
	  if((E_sum_iso_rec1_angle/reco1_angle.E())<0.20){
	    pass_teff_el_deg1_and_relIso=true;
	  }
	}
	t_hist_vec[15]->Fill(true1.Angle(reco1_angle.Vect())*TMath::RadToDeg());
	t_hist_vec[16]->Fill(reco1_angle.E()/true1.E());
	t_hist_vec[17]->Fill(reco1_angle.E());
	t_hist_vec[27]->Fill(E_sum_iso_rec1_angle);
	t_hist_vec[28]->Fill(E_sum_iso_rec1_angle/reco1_angle.E());
	t_hist_vec[29]->Fill(E_sum_iso_rec1_angle_noph_0_010);
	t_hist_vec[30]->Fill(E_sum_iso_rec1_angle_noph_0_010/reco1_angle.E());
	t_hist_vec[31]->Fill(E_sum_iso_rec1_angle_noph_1_5_deg);
	t_hist_vec[32]->Fill(E_sum_iso_rec1_angle_noph_1_5_deg/reco1_angle.E());
      }
      if(ind_reco1_E_deg_2!=-1){
	t_hist_vec[18]->Fill(true1.Angle(reco1_E_ang_deg_2.Vect())*TMath::RadToDeg());
	t_hist_vec[19]->Fill(true1.E()/reco1_E_ang_deg_2.E());
	t_hist_vec[20]->Fill(reco1_E_ang_deg_2.E());
	t_hist_vec[33]->Fill(E_sum_iso_rec1_E_deg2);
	t_hist_vec[34]->Fill(E_sum_iso_rec1_E_deg2/reco1_E_ang_deg_2.E());
	t_hist_vec[35]->Fill(E_sum_iso_rec1_E_deg2_noph_0_010);
	t_hist_vec[36]->Fill(E_sum_iso_rec1_E_deg2_noph_0_010/reco1_E_ang_deg_2.E());
	t_hist_vec[37]->Fill(E_sum_iso_rec1_E_deg2_noph_1_5_deg);
	t_hist_vec[38]->Fill(E_sum_iso_rec1_E_deg2_noph_1_5_deg/reco1_E_ang_deg_2.E());
      }
      if(fabs(true1.CosTheta())<0.95){
	t_TEff_vec[0]->Fill(pass_teff_el_deg1,true1.E());
	t_TEff_vec[1]->Fill(pass_teff_el_deg1_and_relIso,true1.E());
      }
      if(true1.E()>=5){
	if(true1.E()<10){
	  t_TEff_vec[44]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[56]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true1.E()<30){
	  t_TEff_vec[4]->Fill(pass_teff_el_deg1,true1.CosTheta());
	  t_TEff_vec[5]->Fill(pass_teff_el_deg1_and_relIso,true1.CosTheta());
	  t_TEff_vec[6]->Fill(pass_teff_el_deg1,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[7]->Fill(pass_teff_el_deg1_and_relIso,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[45]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[57]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true1.E()<75){
	  t_TEff_vec[12]->Fill(pass_teff_el_deg1,true1.CosTheta());
	  t_TEff_vec[13]->Fill(pass_teff_el_deg1_and_relIso,true1.CosTheta());
	  t_TEff_vec[14]->Fill(pass_teff_el_deg1,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[15]->Fill(pass_teff_el_deg1_and_relIso,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[46]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[58]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true1.E()<150){
	  t_TEff_vec[20]->Fill(pass_teff_el_deg1,true1.CosTheta());
	  t_TEff_vec[21]->Fill(pass_teff_el_deg1_and_relIso,true1.CosTheta());
	  t_TEff_vec[22]->Fill(pass_teff_el_deg1,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[23]->Fill(pass_teff_el_deg1_and_relIso,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[47]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[59]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true1.E()<500){
	  t_TEff_vec[28]->Fill(pass_teff_el_deg1,true1.CosTheta());
	  t_TEff_vec[29]->Fill(pass_teff_el_deg1_and_relIso,true1.CosTheta());
	  t_TEff_vec[30]->Fill(pass_teff_el_deg1,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[31]->Fill(pass_teff_el_deg1_and_relIso,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[48]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[60]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else{
	  t_TEff_vec[36]->Fill(pass_teff_el_deg1,true1.CosTheta());
	  t_TEff_vec[37]->Fill(pass_teff_el_deg1_and_relIso,true1.CosTheta());
	  t_TEff_vec[38]->Fill(pass_teff_el_deg1,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[39]->Fill(pass_teff_el_deg1_and_relIso,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[49]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[61]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}
      }
    }else if(abs(pdg_true1)==13){
      bool pass_teff_mu_deg1_and_relIso=false;
      bool pass_teff_mu_deg1=false;
      bool pass_teff_mu_deg1_EMatch=false;
      bool pass_teff_mu_deg1_PtMatch=false;
      //std::cout<<"true muon"<<std::endl;
      //leading lepton is Muon --> check performed while indices are searched for
      if(ind_reco1_Ang!=-1){
	//if( (true1.DeltaPhi(reco1_angle))<0.02 && fabs(true1.Theta()-reco1_angle.Theta())<0.01 && fabs((true1.Pt()-reco1_angle.Pt())/true1.Pt())<0.05){
	if( (true1.Angle(reco1_angle.Vect())*TMath::RadToDeg())<1.0){
	  pass_teff_mu_deg1=true;
	  if(fabs(true1.E()-reco1_angle.E())<(E_match_particle*true1.E())){
	    pass_teff_mu_deg1_EMatch=true;
	  }
	  if(fabs(true1.Pt()-reco1_angle.Pt())<(Pt_match_particle*true1.Pt())){
	    pass_teff_mu_deg1_PtMatch=true;
	  }
	  if((E_sum_iso_rec1_angle/reco1_angle.E())<0.10){
	    pass_teff_mu_deg1_and_relIso=true;
	  }
	}
	t_hist_vec[21]->Fill(true1.Angle(reco1_angle.Vect())*TMath::RadToDeg());
	t_hist_vec[22]->Fill(reco1_angle.E()/true1.E());
	t_hist_vec[23]->Fill(reco1_angle.E());
	t_hist_vec[39]->Fill(E_sum_iso_rec1_angle);
	t_hist_vec[40]->Fill(E_sum_iso_rec1_angle/reco1_angle.E());
	t_hist_vec[41]->Fill(E_sum_iso_rec1_angle_noph_0_010);
	t_hist_vec[42]->Fill(E_sum_iso_rec1_angle_noph_0_010/reco1_angle.E());
	t_hist_vec[43]->Fill(E_sum_iso_rec1_angle_noph_1_5_deg);
	t_hist_vec[44]->Fill(E_sum_iso_rec1_angle_noph_1_5_deg/reco1_angle.E());
      }
      if(ind_reco1_E_deg_2!=-1){
	t_hist_vec[24]->Fill(true1.Angle(reco1_E_ang_deg_2.Vect())*TMath::RadToDeg());
	t_hist_vec[25]->Fill(true1.E()/reco1_E_ang_deg_2.E());
	t_hist_vec[26]->Fill(reco1_E_ang_deg_2.E());
	t_hist_vec[45]->Fill(E_sum_iso_rec1_E_deg2);
	t_hist_vec[46]->Fill(E_sum_iso_rec1_E_deg2/reco1_E_ang_deg_2.E());
	t_hist_vec[47]->Fill(E_sum_iso_rec1_E_deg2_noph_0_010);
	t_hist_vec[48]->Fill(E_sum_iso_rec1_E_deg2_noph_0_010/reco1_E_ang_deg_2.E());
	t_hist_vec[49]->Fill(E_sum_iso_rec1_E_deg2_noph_1_5_deg);
	t_hist_vec[50]->Fill(E_sum_iso_rec1_E_deg2_noph_1_5_deg/reco1_E_ang_deg_2.E());
      }
      if(fabs(true1.CosTheta())<0.95){
	t_TEff_vec[2]->Fill(pass_teff_mu_deg1,true1.E());
	t_TEff_vec[3]->Fill(pass_teff_mu_deg1_and_relIso,true1.E());
      }
      if(true1.E()>=5){
	if(true1.E()<10){
	  t_TEff_vec[50]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[62]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true1.E()<30){
	  t_TEff_vec[8]->Fill(pass_teff_mu_deg1,true1.CosTheta());
	  t_TEff_vec[9]->Fill(pass_teff_mu_deg1_and_relIso,true1.CosTheta());
	  t_TEff_vec[10]->Fill(pass_teff_mu_deg1,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[11]->Fill(pass_teff_mu_deg1_and_relIso,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[51]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[63]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true1.E()<75){
	  t_TEff_vec[16]->Fill(pass_teff_mu_deg1,true1.CosTheta());
	  t_TEff_vec[17]->Fill(pass_teff_mu_deg1_and_relIso,true1.CosTheta());
	  t_TEff_vec[18]->Fill(pass_teff_mu_deg1,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[19]->Fill(pass_teff_mu_deg1_and_relIso,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[52]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[64]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true1.E()<150){
	  t_TEff_vec[24]->Fill(pass_teff_mu_deg1,true1.CosTheta());
	  t_TEff_vec[25]->Fill(pass_teff_mu_deg1_and_relIso,true1.CosTheta());
	  t_TEff_vec[26]->Fill(pass_teff_mu_deg1,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[27]->Fill(pass_teff_mu_deg1_and_relIso,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[53]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[65]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true1.E()<500){
	  t_TEff_vec[32]->Fill(pass_teff_mu_deg1,true1.CosTheta());
	  t_TEff_vec[33]->Fill(pass_teff_mu_deg1_and_relIso,true1.CosTheta());
	  t_TEff_vec[34]->Fill(pass_teff_mu_deg1,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[35]->Fill(pass_teff_mu_deg1_and_relIso,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[54]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[66]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else{
	  t_TEff_vec[40]->Fill(pass_teff_mu_deg1,true1.CosTheta());
	  t_TEff_vec[41]->Fill(pass_teff_mu_deg1_and_relIso,true1.CosTheta());
	  t_TEff_vec[42]->Fill(pass_teff_mu_deg1,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[43]->Fill(pass_teff_mu_deg1_and_relIso,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[55]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[67]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}
      }
    }

   //interesting leptons have been defined at this point, now check for isolation values 
    if(abs(pdg_true2)==11){
      bool pass_teff_el_deg1_and_relIso=false;
      bool pass_teff_el_deg1=false;
      bool pass_teff_el_deg1_EMatch=false;
      bool pass_teff_el_deg1_PtMatch=false;
      //second leading lepton is Electron --> check performed while indices are searched for
      if(ind_reco2_Ang!=-1){
	if( (true2.Angle(reco2_angle.Vect())*TMath::RadToDeg())<1.0){
	  //if( (true2.DeltaPhi(reco2_angle))<0.02 && fabs(true2.Theta()-reco2_angle.Theta())<0.01 && fabs((true2.Pt()-reco2_angle.Pt())/true2.Pt())<0.05){
	  pass_teff_el_deg1=true;
	  if(fabs(true2.E()-reco2_angle.E())<(E_match_particle*true2.E())){
	    pass_teff_el_deg1_EMatch=true;
	  }
	  if(fabs(true2.Pt()-reco2_angle.Pt())<(Pt_match_particle*true2.Pt())){
	    pass_teff_el_deg1_PtMatch=true;
	  }
	  if((E_sum_iso_rec2_angle/reco2_angle.E())<0.20){
	    pass_teff_el_deg1_and_relIso=true;
	  }
	}
	t_hist_vec[15]->Fill(true2.Angle(reco2_angle.Vect())*TMath::RadToDeg());
	t_hist_vec[16]->Fill(reco2_angle.E()/true2.E());
	t_hist_vec[17]->Fill(reco2_angle.E());
	t_hist_vec[27]->Fill(E_sum_iso_rec2_angle);
	t_hist_vec[28]->Fill(E_sum_iso_rec2_angle/reco2_angle.E());
	t_hist_vec[29]->Fill(E_sum_iso_rec2_angle_noph_0_010);
	t_hist_vec[30]->Fill(E_sum_iso_rec2_angle_noph_0_010/reco2_angle.E());
	t_hist_vec[31]->Fill(E_sum_iso_rec2_angle_noph_1_5_deg);
	t_hist_vec[32]->Fill(E_sum_iso_rec2_angle_noph_1_5_deg/reco2_angle.E());
      }
      if(ind_reco2_E_deg_2!=-1){
	t_hist_vec[18]->Fill(true2.Angle(reco2_E_ang_deg_2.Vect())*TMath::RadToDeg());
	t_hist_vec[19]->Fill(true2.E()/reco2_E_ang_deg_2.E());
	t_hist_vec[20]->Fill(reco2_E_ang_deg_2.E());
	t_hist_vec[33]->Fill(E_sum_iso_rec2_E_deg2);
	t_hist_vec[34]->Fill(E_sum_iso_rec2_E_deg2/reco2_E_ang_deg_2.E());
	t_hist_vec[35]->Fill(E_sum_iso_rec2_E_deg2_noph_0_010);
	t_hist_vec[36]->Fill(E_sum_iso_rec2_E_deg2_noph_0_010/reco2_E_ang_deg_2.E());
	t_hist_vec[37]->Fill(E_sum_iso_rec2_E_deg2_noph_1_5_deg);
	t_hist_vec[38]->Fill(E_sum_iso_rec2_E_deg2_noph_1_5_deg/reco2_E_ang_deg_2.E());
      }
      if(fabs(true2.CosTheta())<0.95){
	t_TEff_vec[0]->Fill(pass_teff_el_deg1,true2.E());
	t_TEff_vec[1]->Fill(pass_teff_el_deg1_and_relIso,true2.E());
      }
      if(true2.E()>=5){
	if(true2.E()<10){
	  t_TEff_vec[44]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[56]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true2.E()<30){
	  t_TEff_vec[4]->Fill(pass_teff_el_deg1,true2.CosTheta());
	  t_TEff_vec[5]->Fill(pass_teff_el_deg1_and_relIso,true2.CosTheta());
	  t_TEff_vec[6]->Fill(pass_teff_el_deg1,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[7]->Fill(pass_teff_el_deg1_and_relIso,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[45]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[57]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true2.E()<75){
	  t_TEff_vec[12]->Fill(pass_teff_el_deg1,true2.CosTheta());
	  t_TEff_vec[13]->Fill(pass_teff_el_deg1_and_relIso,true2.CosTheta());
	  t_TEff_vec[14]->Fill(pass_teff_el_deg1,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[15]->Fill(pass_teff_el_deg1_and_relIso,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[46]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[58]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true2.E()<150){
	  t_TEff_vec[20]->Fill(pass_teff_el_deg1,true2.CosTheta());
	  t_TEff_vec[21]->Fill(pass_teff_el_deg1_and_relIso,true2.CosTheta());
	  t_TEff_vec[22]->Fill(pass_teff_el_deg1,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[23]->Fill(pass_teff_el_deg1_and_relIso,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[47]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[59]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true2.E()<500){
	  t_TEff_vec[28]->Fill(pass_teff_el_deg1,true2.CosTheta());
	  t_TEff_vec[29]->Fill(pass_teff_el_deg1_and_relIso,true2.CosTheta());
	  t_TEff_vec[30]->Fill(pass_teff_el_deg1,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[31]->Fill(pass_teff_el_deg1_and_relIso,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[48]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[60]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else{
	  t_TEff_vec[36]->Fill(pass_teff_el_deg1,true2.CosTheta());
	  t_TEff_vec[37]->Fill(pass_teff_el_deg1_and_relIso,true2.CosTheta());
	  t_TEff_vec[38]->Fill(pass_teff_el_deg1,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[39]->Fill(pass_teff_el_deg1_and_relIso,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[49]->Fill(pass_teff_el_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[61]->Fill(pass_teff_el_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}
      }
    }else if(abs(pdg_true2)==13){
      bool pass_teff_mu_deg1_and_relIso=false;
      bool pass_teff_mu_deg1=false;
      bool pass_teff_mu_deg1_EMatch=false;
      bool pass_teff_mu_deg1_PtMatch=false;
      //leading lepton is Muon --> check performed while indices are searched for
      if(ind_reco2_Ang!=-1){
	if( (true2.Angle(reco2_angle.Vect())*TMath::RadToDeg())<1.0){
	  //if( (true2.DeltaPhi(reco2_angle))<0.02 && fabs(true2.Theta()-reco2_angle.Theta())<0.01 && fabs((true2.Pt()-reco2_angle.Pt())/true2.Pt())<0.05){
	  pass_teff_mu_deg1=true;
	  if(fabs(true2.E()-reco2_angle.E())<(E_match_particle*true2.E())){
	    pass_teff_mu_deg1_EMatch=true;
	  }
	  if(fabs(true2.Pt()-reco2_angle.Pt())<(Pt_match_particle*true2.Pt())){
	    pass_teff_mu_deg1_PtMatch=true;
	  }
	  if((E_sum_iso_rec2_angle/reco2_angle.E())<0.10){
	    pass_teff_mu_deg1_and_relIso=true;
	  }
	}
	t_hist_vec[21]->Fill(true2.Angle(reco2_angle.Vect())*TMath::RadToDeg());
	t_hist_vec[22]->Fill(reco2_angle.E()/true2.E());
	t_hist_vec[23]->Fill(reco2_angle.E());
	t_hist_vec[39]->Fill(E_sum_iso_rec2_angle);
	t_hist_vec[40]->Fill(E_sum_iso_rec2_angle/reco2_angle.E());
	t_hist_vec[41]->Fill(E_sum_iso_rec2_angle_noph_0_010);
	t_hist_vec[42]->Fill(E_sum_iso_rec2_angle_noph_0_010/reco2_angle.E());
	t_hist_vec[43]->Fill(E_sum_iso_rec2_angle_noph_1_5_deg);
	t_hist_vec[44]->Fill(E_sum_iso_rec2_angle_noph_1_5_deg/reco2_angle.E());
      }
      if(ind_reco2_E_deg_2!=-1){
	t_hist_vec[24]->Fill(true2.Angle(reco2_E_ang_deg_2.Vect())*TMath::RadToDeg());
	t_hist_vec[25]->Fill(true2.E()/reco2_E_ang_deg_2.E());
	t_hist_vec[26]->Fill(reco2_E_ang_deg_2.E());
	t_hist_vec[45]->Fill(E_sum_iso_rec2_E_deg2);
	t_hist_vec[46]->Fill(E_sum_iso_rec2_E_deg2/reco2_E_ang_deg_2.E());
	t_hist_vec[47]->Fill(E_sum_iso_rec2_E_deg2_noph_0_010);
	t_hist_vec[48]->Fill(E_sum_iso_rec2_E_deg2_noph_0_010/reco2_E_ang_deg_2.E());
	t_hist_vec[49]->Fill(E_sum_iso_rec2_E_deg2_noph_1_5_deg);
	t_hist_vec[50]->Fill(E_sum_iso_rec2_E_deg2_noph_1_5_deg/reco2_E_ang_deg_2.E());
      }
      if(fabs(true2.CosTheta())<0.95){
	t_TEff_vec[2]->Fill(pass_teff_mu_deg1,true2.E());
	t_TEff_vec[3]->Fill(pass_teff_mu_deg1_and_relIso,true2.E());
      }    
     if(true1.E()>=5){
	if(true1.E()<10){
	  t_TEff_vec[50]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[62]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true2.E()<30){
	  t_TEff_vec[8]->Fill(pass_teff_mu_deg1,true2.CosTheta());
	  t_TEff_vec[9]->Fill(pass_teff_mu_deg1_and_relIso,true2.CosTheta());
	  t_TEff_vec[10]->Fill(pass_teff_mu_deg1,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[11]->Fill(pass_teff_mu_deg1_and_relIso,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[51]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[63]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true2.E()<75){
	  t_TEff_vec[16]->Fill(pass_teff_mu_deg1,true2.CosTheta());
	  t_TEff_vec[17]->Fill(pass_teff_mu_deg1_and_relIso,true2.CosTheta());
	  t_TEff_vec[18]->Fill(pass_teff_mu_deg1,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[19]->Fill(pass_teff_mu_deg1_and_relIso,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[52]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[64]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true2.E()<150){
	  t_TEff_vec[24]->Fill(pass_teff_mu_deg1,true2.CosTheta());
	  t_TEff_vec[25]->Fill(pass_teff_mu_deg1_and_relIso,true2.CosTheta());
	  t_TEff_vec[26]->Fill(pass_teff_mu_deg1,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[27]->Fill(pass_teff_mu_deg1_and_relIso,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[53]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[65]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else if(true2.E()<500){
	  t_TEff_vec[32]->Fill(pass_teff_mu_deg1,true2.CosTheta());
	  t_TEff_vec[33]->Fill(pass_teff_mu_deg1_and_relIso,true2.CosTheta());
	  t_TEff_vec[34]->Fill(pass_teff_mu_deg1,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[35]->Fill(pass_teff_mu_deg1_and_relIso,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[54]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[66]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}else{
	  t_TEff_vec[40]->Fill(pass_teff_mu_deg1,true2.CosTheta());
	  t_TEff_vec[41]->Fill(pass_teff_mu_deg1_and_relIso,true2.CosTheta());
	  t_TEff_vec[42]->Fill(pass_teff_mu_deg1,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[43]->Fill(pass_teff_mu_deg1_and_relIso,true2.Theta()*TMath::RadToDeg());
	  t_TEff_vec[55]->Fill(pass_teff_mu_deg1_EMatch,true1.Theta()*TMath::RadToDeg());
	  t_TEff_vec[67]->Fill(pass_teff_mu_deg1_PtMatch,true1.Theta()*TMath::RadToDeg());
	}
      }
    }

    //now check for fake Electrons and Muons -->fake means check for angles larger than 6 away from
    //W lepton generator direction and then draw corresponding histograms
    //check for certain energy of fake leptons only
    float veto_angle_min=6.0;//in degrees
    //if reco lepton within this angle from any of the true leptons, then remove them, we still should have enough statistics
    //just to be sure, don't check for misID right now
    for(unsigned int i_reco=0;i_reco<reco_Energy->size();i_reco++){
      bool is_fake_lepton = true;
      TLorentzVector tempReco(0,0,0,0);
      tempReco.SetPxPyPzE((*reco_Px)[i_reco],(*reco_Py)[i_reco],(*reco_Pz)[i_reco],(*reco_Energy)[i_reco]);
      if((*reco_Energy)[i_reco]<E_min_particle){
	continue;
      }
      if(true_Energy->size()>0){
	if((true1.Angle(tempReco.Vect())*TMath::RadToDeg())<veto_angle_min){
	  is_fake_lepton=false;
	}
	if(true_Energy->size()>1){
	  if((true2.Angle(tempReco.Vect())*TMath::RadToDeg())<veto_angle_min){
	    is_fake_lepton=false;
	  }
	}
      }
      if(is_fake_lepton){
	if(i_reco==ind_reco2_Ang || i_reco==ind_reco2_E_deg_2 || i_reco==ind_reco1_Ang || i_reco==ind_reco1_E_deg_2){
	  std::cout<<"very strange, Electron or Muon should be one of the isolated ones, WTF "<<ind_reco1_Ang<<"/"<<ind_reco1_E_deg_2<<"/"<<ind_reco2_Ang<<"/"<<ind_reco2_E_deg_2<<std::endl;
	}
	float E_sum_iso_fake=0;
	float E_sum_iso_fake_noph_0_010=0;
	float E_sum_iso_fake_noph_1_5_deg=0;
	for(unsigned int i_rec_iso=0;i_rec_iso<reco_all_6_deg_index->size();i_rec_iso++){
	  //preselection is in fact 6, so the last cut has no effect 
	  if( ((*reco_all_6_deg_index)[i_rec_iso]==i_reco) && ((*reco_all_6_deg_angle)[i_rec_iso]*TMath::RadToDeg())<iso_angle_max ){
	    E_sum_iso_fake+=(*reco_all_6_deg_E)[i_rec_iso];
	    //around 0.5729 degrees
	    if((*reco_all_6_deg_PDGID)[i_rec_iso]==22 && (*reco_all_6_deg_angle)[i_rec_iso]<ph_angle_small){
	      continue;
	    }
	    E_sum_iso_fake_noph_0_010+=(*reco_all_6_deg_E)[i_rec_iso];
	    //additional photons cut out --> Electron check
	    if((*reco_all_6_deg_PDGID)[i_rec_iso]==22 && ((*reco_all_6_deg_angle)[i_rec_iso]*TMath::RadToDeg())<ph_angle_large){
	      continue;
	    }
	    E_sum_iso_fake_noph_1_5_deg+=(*reco_all_6_deg_E)[i_rec_iso];
	  }
	}
	if(abs((*reco_PDGID)[i_reco])==11){
	  //fake Electrons
	  t_hist_vec[51]->Fill(tempReco.E());
	  t_hist_vec[52]->Fill(E_sum_iso_fake);
	  t_hist_vec[53]->Fill(E_sum_iso_fake/tempReco.E());
	  t_hist_vec[54]->Fill(E_sum_iso_fake_noph_0_010);
	  t_hist_vec[55]->Fill(E_sum_iso_fake_noph_0_010/tempReco.E());
	  t_hist_vec[56]->Fill(E_sum_iso_fake_noph_1_5_deg);
	  t_hist_vec[57]->Fill(E_sum_iso_fake_noph_1_5_deg/tempReco.E());
	}else if (abs((*reco_PDGID)[i_reco])==13){
	  //fake Muons
	  t_hist_vec[58]->Fill(tempReco.E());
	  t_hist_vec[59]->Fill(E_sum_iso_fake);
	  t_hist_vec[60]->Fill(E_sum_iso_fake/tempReco.E());
	  t_hist_vec[61]->Fill(E_sum_iso_fake_noph_0_010);
	  t_hist_vec[62]->Fill(E_sum_iso_fake_noph_0_010/tempReco.E());
	  t_hist_vec[63]->Fill(E_sum_iso_fake_noph_1_5_deg);
	  t_hist_vec[64]->Fill(E_sum_iso_fake_noph_1_5_deg/tempReco.E());
	}
      }
    }
  }
}

void plotTTBarEfficiency_Summary(){

 
  CLICdpStyle();

    //2,5,10,20,30,40,50,60,75,90,100,150,200,250,400,500,1000,1500
  //TFile* file_noOverlay=TFile::Open("/eos/user/w/weberma2/data/ttbarMacroFiles/180208_gcc62/ttbarStudy_ILC180208_gcc62_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_9549_noOverlay_PandoraPFOs_DressedLeptons.root");
  //TFile* file_Overlay=TFile::Open("/eos/user/w/weberma2/data/ttbarMacroFiles/180208_gcc62/ttbarStudy_ILC180208_gcc62_CT_FitFW_wRefit_CLIC_o3_v14_SWC_CLIC_9553_Overlay_TightSelectedPandoraPFOs_DressedLeptons.root");


  TFile* file_noOverlay=TFile::Open("/eos/user/w/weberma2/data/ttbarMacroFiles/181101_gcc62/photonStudy_tt3000_11795_noOverlay_CLIC_o3_v14_AllMCInfo_ElMu_fromW_PandoraPFOs_real.root");
  TFile* file_Overlay=TFile::Open("/eos/user/w/weberma2/data/ttbarMacroFiles/181101_gcc62/photonStudy_tt3000_11796_3TeVOverlay_CLIC_o3_v14_AllMCInfo_ElMu_fromW_TightSelectedPandoraPFOs_real.root");

  const char* final_histo_name="/eos/user/w/weberma2/data/validation181101/ttbar_11795_noOverlay_vs_11796_Overlay_E_cut_isoAngle_6_TEff_CosTrue_0_95_ang_match_1_deg_plus_rel_iso_el_0_2_mu_0_1_andDraw.root";


  int n_bins_high=200;
  double lim_energy_low=0;
  double lim_energy_high=1500.;



  double lim_rel_energy_low=0.5;
  double lim_rel_energy_high=1.5;

  double lim_E_sum_iso_low=0.;
  double lim_E_sum_iso_high=250.;

  double lim_E_sum_iso_rel_low=0.;
  double lim_E_sum_iso_rel_high=2.50;
  double lim_E_sum_iso_rel_high_fake=4.00;

  TFile* file_histogram=new TFile(final_histo_name,"recreate");

  TH1F* h_true_el_E_noOverlay = new TH1F("h_true_el_E_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_true_el_E_iso_sum_noOverlay = new TH1F("h_true_el_E_iso_sum_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_el_E_iso_rel_noOverlay = new TH1F("h_true_el_E_iso_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_true_el_E_iso_sum_noPh_0_10_noOverlay = new TH1F("h_true_el_E_iso_sum_noPh_0_10_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_el_E_iso_noPh_0_10_rel_noOverlay = new TH1F("h_true_el_E_iso_noPh_0_10_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_true_el_E_iso_sum_noPh_1_50_deg_noOverlay = new TH1F("h_true_el_E_iso_sum_noPh_1_50_deg_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_el_E_iso_noPh_1_50_deg_rel_noOverlay = new TH1F("h_true_el_E_iso_noPh_1_50_deg_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_true_mu_E_noOverlay = new TH1F("h_true_mu_E_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_true_mu_E_iso_sum_noOverlay = new TH1F("h_true_mu_E_iso_sum_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_mu_E_iso_rel_noOverlay = new TH1F("h_true_mu_E_iso_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_true_mu_E_iso_sum_noPh_0_10_noOverlay = new TH1F("h_true_mu_E_iso_sum_noPh_0_10_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_mu_E_iso_noPh_0_10_rel_noOverlay = new TH1F("h_true_mu_E_iso_noPh_0_10_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_true_mu_E_iso_sum_noPh_1_50_deg_noOverlay = new TH1F("h_true_mu_E_iso_sum_noPh_1_50_deg_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_mu_E_iso_noPh_1_50_deg_rel_noOverlay = new TH1F("h_true_mu_E_iso_noPh_1_50_deg_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_true_l1l2_dangle_noOverlay = new TH1F("h_true_l1l2_dangle_noOverlay","", n_bins_high,0.0,180.);

  TH1F* h_reco_el_true_reco_dangle_ang_match_noOverlay = new TH1F("h_reco_el_true_reco_dangle_ang_match_noOverlay","", n_bins_high,0.0,6.0);
  TH1F* h_reco_el_reco_E_over_true_E_ang_match_noOverlay = new TH1F("h_reco_el_reco_E_over_true_E_ang_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_reco_el_E_reco_ang_match_noOverlay = new TH1F("h_reco_el_E_reco_ang_match_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_reco_el_true_reco_dangle_E_match_noOverlay = new TH1F("h_reco_el_true_reco_dangle_E_match_noOverlay","", n_bins_high,0.0,6.0);
  TH1F* h_reco_el_reco_E_over_true_E_match_noOverlay = new TH1F("h_reco_el_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_reco_el_E_reco_E_match_noOverlay = new TH1F("h_reco_el_E_reco_E_match_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_reco_mu_true_reco_dangle_ang_match_noOverlay = new TH1F("h_reco_mu_true_reco_dangle_ang_match_noOverlay","", n_bins_high,0.0,6.0);
  TH1F* h_reco_mu_reco_E_over_true_E_ang_match_noOverlay = new TH1F("h_reco_mu_reco_E_over_true_E_ang_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_reco_mu_E_reco_ang_match_noOverlay = new TH1F("h_reco_mu_E_reco_ang_match_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_reco_mu_true_reco_dangle_E_match_noOverlay = new TH1F("h_reco_mu_true_reco_dangle_E_match_noOverlay","", n_bins_high,0.0,6.0);
  TH1F* h_reco_mu_reco_E_over_true_E_match_noOverlay = new TH1F("h_reco_mu_reco_E_over_true_E_match_noOverlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_reco_mu_E_reco_E_match_noOverlay = new TH1F("h_reco_mu_E_reco_E_match_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_reco_el_ang_match_E_iso_sum_noOverlay = new TH1F("h_reco_el_ang_match_E_iso_sum_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_ang_match_E_iso_rel_noOverlay = new TH1F("h_reco_el_ang_match_E_iso_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_el_ang_match_E_iso_sum_noPh_0_10_noOverlay = new TH1F("h_reco_el_ang_match_E_iso_sum_noPh_0_10_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_ang_match_E_iso_noPh_0_10_rel_noOverlay = new TH1F("h_reco_el_ang_match_E_iso_noPh_0_10_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_el_ang_match_E_iso_sum_noPh_1_50_deg_noOverlay = new TH1F("h_reco_el_ang_match_E_iso_sum_noPh_1_50_deg_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_ang_match_E_iso_noPh_1_50_deg_rel_noOverlay = new TH1F("h_reco_el_ang_match_E_iso_noPh_1_50_deg_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_reco_el_E_match_E_iso_sum_noOverlay = new TH1F("h_reco_el_E_match_E_iso_sum_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_E_match_E_iso_rel_noOverlay = new TH1F("h_reco_el_E_match_E_iso_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_el_E_match_E_iso_sum_noPh_0_10_noOverlay = new TH1F("h_reco_el_E_match_E_iso_sum_noPh_0_10_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_E_match_E_iso_noPh_0_10_rel_noOverlay = new TH1F("h_reco_el_E_match_E_iso_noPh_0_10_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_el_E_match_E_iso_sum_noPh_1_50_deg_noOverlay = new TH1F("h_reco_el_E_match_E_iso_sum_noPh_1_50_deg_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_E_match_E_iso_noPh_1_50_deg_rel_noOverlay = new TH1F("h_reco_el_E_match_E_iso_noPh_1_50_deg_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_reco_mu_ang_match_E_iso_sum_noOverlay = new TH1F("h_reco_mu_ang_match_E_iso_sum_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_ang_match_E_iso_rel_noOverlay = new TH1F("h_reco_mu_ang_match_E_iso_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_mu_ang_match_E_iso_sum_noPh_0_10_noOverlay = new TH1F("h_reco_mu_ang_match_E_iso_sum_noPh_0_10_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_ang_match_E_iso_noPh_0_10_rel_noOverlay = new TH1F("h_reco_mu_ang_match_E_iso_noPh_0_10_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_mu_ang_match_E_iso_sum_noPh_1_50_deg_noOverlay = new TH1F("h_reco_mu_ang_match_E_iso_sum_noPh_1_50_deg_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_ang_match_E_iso_noPh_1_50_deg_rel_noOverlay = new TH1F("h_reco_mu_ang_match_E_iso_noPh_1_50_deg_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_reco_mu_E_match_E_iso_sum_noOverlay = new TH1F("h_reco_mu_E_match_E_iso_sum_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_E_match_E_iso_rel_noOverlay = new TH1F("h_reco_mu_E_match_E_iso_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_mu_E_match_E_iso_sum_noPh_0_10_noOverlay = new TH1F("h_reco_mu_E_match_E_iso_sum_noPh_0_10_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_E_match_E_iso_noPh_0_10_rel_noOverlay = new TH1F("h_reco_mu_E_match_E_iso_noPh_0_10_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_mu_E_match_E_iso_sum_noPh_1_50_deg_noOverlay = new TH1F("h_reco_mu_E_match_E_iso_sum_noPh_1_50_deg_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_E_match_E_iso_noPh_1_50_deg_rel_noOverlay = new TH1F("h_reco_mu_E_match_E_iso_noPh_1_50_deg_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);



  TH1F* h_reco_fake_el_E_noOverlay = new TH1F("h_reco_fake_el_E_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_fake_el_E_iso_sum_noOverlay = new TH1F("h_reco_fake_el_E_iso_sum_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_el_E_iso_rel_noOverlay = new TH1F("h_reco_fake_el_E_iso_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);
  TH1F* h_reco_fake_el_E_iso_sum_noPh_0_10_noOverlay = new TH1F("h_reco_fake_el_E_iso_sum_noPh_0_10_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_el_E_iso_noPh_0_10_rel_noOverlay = new TH1F("h_reco_fake_el_E_iso_noPh_0_10_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);
  TH1F* h_reco_fake_el_E_iso_sum_noPh_1_50_deg_noOverlay = new TH1F("h_reco_fake_el_E_iso_sum_noPh_1_50_deg_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_el_E_iso_noPh_1_50_deg_rel_noOverlay = new TH1F("h_reco_fake_el_E_iso_noPh_1_50_deg_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);

  TH1F* h_reco_fake_mu_E_noOverlay = new TH1F("h_reco_fake_mu_E_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_fake_mu_E_iso_sum_noOverlay = new TH1F("h_reco_fake_mu_E_iso_sum_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_mu_E_iso_rel_noOverlay = new TH1F("h_reco_fake_mu_E_iso_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);
  TH1F* h_reco_fake_mu_E_iso_sum_noPh_0_10_noOverlay = new TH1F("h_reco_fake_mu_E_iso_sum_noPh_0_10_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_mu_E_iso_noPh_0_10_rel_noOverlay = new TH1F("h_reco_fake_mu_E_iso_noPh_0_10_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);
  TH1F* h_reco_fake_mu_E_iso_sum_noPh_1_50_deg_noOverlay = new TH1F("h_reco_fake_mu_E_iso_sum_noPh_1_50_deg_noOverlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_mu_E_iso_noPh_1_50_deg_rel_noOverlay = new TH1F("h_reco_fake_mu_E_iso_noPh_1_50_deg_rel_noOverlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);

  TH1F* h_true_el_E_cosTheta_0_95_noOverlay = new TH1F("h_true_el_E_cosTheta_0_95_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_el_E_reco_ang_match_rel_iso_0_01_cosTheta_0_95_noOverlay = new TH1F("h_reco_el_E_reco_ang_match_1_deg_cosTheta_0_95_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_el_E_reco_ang_match_rel_iso_0_02_cosTheta_0_95_noOverlay = new TH1F("h_reco_el_E_reco_ang_match_1_deg_rel_iso_0_02_cosTheta_0_95_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_true_mu_E_cosTheta_0_95_noOverlay = new TH1F("h_true_mu_E_cosTheta_0_95_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_mu_E_reco_ang_match_rmu_iso_0_01_cosTheta_0_95_noOverlay = new TH1F("h_reco_mu_E_reco_ang_match_1_deg_cosTheta_0_95_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_mu_E_reco_ang_match_rmu_iso_0_02_cosTheta_0_95_noOverlay = new TH1F("h_reco_mu_E_reco_ang_match_1_deg_rmu_iso_0_10_cosTheta_0_95_noOverlay","", n_bins_high, lim_energy_low,lim_energy_high);

  //[51]: fake reco Electron candidate energy --> 
  //[52]: absolute isolation fake reco Electron candidate --> check within 6 degrees
  //[53]: relative isolation fake reco Electron candidate --> check within 6 degrees
  //[54]: absolute isolation excluding photons within 0.010 fake reco Electron candidate --> check within 6 degrees
  //[55]: relative isolation excluding photons within 0.010 fake reco Electron candidate --> check within 6 degrees
  //[56]: absolute isolation excluding photons within 1.5 degrees of fake reco Electron candidate --> check within 6 degrees
  //[57]: relative isolation excluding photons within 1.5 degrees of closest fake reco Electron candidate --> check within 6 degrees

  //[58]: fake reco Muon candidate energy --> 
  //[59]: absolute isolation fake reco Muon candidate --> check within 6 degrees
  //[60]: relative isolation fake reco Muon candidate --> check within 6 degrees
  //[61]: absolute isolation excluding photons within 0.010 fake reco Electron candidate --> check within 6 degrees
  //[62]: relative isolation excluding photons within 0.010 fake reco Muon candidate --> check within 6 degrees
  //[63]: absolute isolation excluding photons within 1.5 degrees of fake reco Muon candidate --> check within 6 degrees
  //[64]: relative isolation excluding photons within 1.5 degrees of closest fake reco Muon candidate --> check within 6 degrees




  std::vector<TH1F*> hist_vector_noOverlay;

  hist_vector_noOverlay.push_back( h_true_el_E_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_el_E_iso_sum_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_el_E_iso_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_el_E_iso_sum_noPh_0_10_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_el_E_iso_noPh_0_10_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_el_E_iso_sum_noPh_1_50_deg_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_el_E_iso_noPh_1_50_deg_rel_noOverlay ); 

  hist_vector_noOverlay.push_back( h_true_mu_E_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_mu_E_iso_sum_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_mu_E_iso_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_mu_E_iso_sum_noPh_0_10_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_mu_E_iso_noPh_0_10_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_mu_E_iso_sum_noPh_1_50_deg_noOverlay ); 
  hist_vector_noOverlay.push_back( h_true_mu_E_iso_noPh_1_50_deg_rel_noOverlay ); 

  hist_vector_noOverlay.push_back( h_true_l1l2_dangle_noOverlay ); 

  hist_vector_noOverlay.push_back( h_reco_el_true_reco_dangle_ang_match_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_reco_E_over_true_E_ang_match_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_E_reco_ang_match_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_true_reco_dangle_E_match_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_reco_E_over_true_E_match_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_E_reco_E_match_noOverlay ); 

  hist_vector_noOverlay.push_back( h_reco_mu_true_reco_dangle_ang_match_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_reco_E_over_true_E_ang_match_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_E_reco_ang_match_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_true_reco_dangle_E_match_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_reco_E_over_true_E_match_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_E_reco_E_match_noOverlay ); 

  hist_vector_noOverlay.push_back( h_reco_el_ang_match_E_iso_sum_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_ang_match_E_iso_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_ang_match_E_iso_sum_noPh_0_10_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_ang_match_E_iso_noPh_0_10_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_ang_match_E_iso_sum_noPh_1_50_deg_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_ang_match_E_iso_noPh_1_50_deg_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_E_match_E_iso_sum_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_E_match_E_iso_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_E_match_E_iso_sum_noPh_0_10_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_E_match_E_iso_noPh_0_10_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_E_match_E_iso_sum_noPh_1_50_deg_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_el_E_match_E_iso_noPh_1_50_deg_rel_noOverlay ); 

  hist_vector_noOverlay.push_back( h_reco_mu_ang_match_E_iso_sum_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_ang_match_E_iso_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_ang_match_E_iso_sum_noPh_0_10_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_ang_match_E_iso_noPh_0_10_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_ang_match_E_iso_sum_noPh_1_50_deg_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_ang_match_E_iso_noPh_1_50_deg_rel_noOverlay ); 

  hist_vector_noOverlay.push_back( h_reco_mu_E_match_E_iso_sum_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_E_match_E_iso_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_E_match_E_iso_sum_noPh_0_10_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_mu_E_match_E_iso_noPh_0_10_rel_noOverlay ); 
  hist_vector_noOverlay.push_back( h_reco_mu_E_match_E_iso_sum_noPh_1_50_deg_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_mu_E_match_E_iso_noPh_1_50_deg_rel_noOverlay );

  hist_vector_noOverlay.push_back( h_reco_fake_el_E_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_el_E_iso_sum_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_el_E_iso_rel_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_el_E_iso_sum_noPh_0_10_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_el_E_iso_noPh_0_10_rel_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_el_E_iso_sum_noPh_1_50_deg_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_el_E_iso_noPh_1_50_deg_rel_noOverlay );
  
  hist_vector_noOverlay.push_back( h_reco_fake_mu_E_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_mu_E_iso_sum_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_mu_E_iso_rel_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_mu_E_iso_sum_noPh_0_10_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_mu_E_iso_noPh_0_10_rel_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_mu_E_iso_sum_noPh_1_50_deg_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_fake_mu_E_iso_noPh_1_50_deg_rel_noOverlay );

  hist_vector_noOverlay.push_back( h_true_el_E_cosTheta_0_95_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_el_E_reco_ang_match_rel_iso_0_01_cosTheta_0_95_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_el_E_reco_ang_match_rel_iso_0_02_cosTheta_0_95_noOverlay );

  hist_vector_noOverlay.push_back( h_true_mu_E_cosTheta_0_95_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_mu_E_reco_ang_match_rmu_iso_0_01_cosTheta_0_95_noOverlay );
  hist_vector_noOverlay.push_back(  h_reco_mu_E_reco_ang_match_rmu_iso_0_02_cosTheta_0_95_noOverlay ); 




  for(unsigned int i=0;i<hist_vector_noOverlay.size();i++){
    hist_vector_noOverlay[i]->Sumw2();
  }

  const unsigned int nbins_teff = 28; 
  double xbins_teff[nbins_teff+1]={0., 2.5, 5. ,7.5, 10.,20., 30., 40., 50., 60., 75., 100., 125., 150., 175., 200., 250., 300., 350., 400.,450.,500.,550.,650.,750.,850.,1000.,1150.,1500.}; 

 TEfficiency* tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay = new TEfficiency("tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay","", nbins_teff,xbins_teff);
 tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->SetTitle(";true electron energy [GeV];electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueE_AngMatch_1_deg_rel_iso_0_20_CosTheta_True_0_95_noOverlay = new TEfficiency("tEff_ElectronVsTrueE_AngMatch_1_deg_rel_iso_0_20_CosTheta_True_0_95_noOverlay","", nbins_teff,xbins_teff);
 tEff_ElectronVsTrueE_AngMatch_1_deg_rel_iso_0_20_CosTheta_True_0_95_noOverlay->SetTitle(";true electron energy [GeV];electron ID efficiency");
 TEfficiency* tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay = new TEfficiency("tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay","", nbins_teff,xbins_teff);
 tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->SetTitle(";true muon energy [GeV];muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueE_AngMatch_1_deg_rel_iso_0_10_CosTheta_True_0_95_noOverlay = new TEfficiency("tEff_MuonVsTrueE_AngMatch_1_deg_rel_iso_0_10_CosTheta_True_0_95_noOverlay","", nbins_teff,xbins_teff);
 tEff_MuonVsTrueE_AngMatch_1_deg_rel_iso_0_10_CosTheta_True_0_95_noOverlay->SetTitle(";true muon energy [GeV];muon ID efficiency");

 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay","", 20,-1.0,1.0);
 tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetTitle(";cos(#theta(e^{true}));electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay","", 20,-1.,1.);
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay","", 20,0,180);
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay","", 20,-1.0,1.0);
 tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetTitle(";cos(#theta(#mu^{true}));muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay","", 20,-1.,1.);
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay","", 20,0,180);

 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay","", 20,-1.0,1.0);
 tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetTitle(";cos(#theta(e^{true}));electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay","", 20,-1.,1.);
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay","", 20,0,180);
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay","", 20,-1.0,1.0);
 tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetTitle(";cos(#theta(#mu^{true}));muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay","", 20,-1.,1.);
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay","", 20,0,180);

 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay","", 20,-1.0,1.0);
 tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetTitle(";cos(#theta(e^{true}));electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay","", 20,-1.,1.);
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay","", 20,0,180);
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay","", 20,-1.0,1.0);
 tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetTitle(";cos(#theta(#mu^{true}));muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay","", 20,-1.,1.);
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay","", 20,0,180);

 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay","", 20,-1.0,1.0);
 tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetTitle(";cos(#theta(e^{true}));electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay","", 20,-1.,1.);
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay","", 20,0,180);
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay","", 20,-1.0,1.0);
 tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetTitle(";cos(#theta(#mu^{true}));muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay","", 20,-1.,1.);
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay","", 20,0,180);

 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay","", 20,-1.0,1.0);
 tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetTitle(";cos(#theta(e^{true}));electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay","", 20,-1.,1.);
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay","", 20,0,180);
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay","", 20,-1.0,1.0);
 tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetTitle(";cos(#theta(#mu^{true}));muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay","", 20,-1.,1.);
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay","", 20,0,180);

 //5/10/30/75/150/500/Inf
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_5_10_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_10_30_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_30_75_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_75_150_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_150_500_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_500_Inf_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");

 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_5_10_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_10_30_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_30_75_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_75_150_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_150_500_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_500_Inf_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");

 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_5_10_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_10_30_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_30_75_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_75_150_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_150_500_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_noOverlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_500_Inf_noOverlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_noOverlay->SetTitle(";#theta(e^{true});electron ID efficiency");

 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_5_10_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_10_30_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_30_75_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_75_150_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_150_500_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_noOverlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_500_Inf_noOverlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_noOverlay->SetTitle(";#theta(#mu^{true});muon ID efficiency");


  std::vector<TEfficiency*> TEff_vector_noOverlay;
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueE_AngMatch_1_deg_rel_iso_0_20_CosTheta_True_0_95_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueE_AngMatch_1_deg_rel_iso_0_10_CosTheta_True_0_95_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_noOverlay); 
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_noOverlay); 

  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_noOverlay);

  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_noOverlay);

  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_noOverlay);

  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_noOverlay);
  TEff_vector_noOverlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_noOverlay);


  fillHistograms(file_noOverlay,hist_vector_noOverlay,TEff_vector_noOverlay);

  TH1F* h_true_el_E_Overlay = new TH1F("h_true_el_E_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_true_el_E_iso_sum_Overlay = new TH1F("h_true_el_E_iso_sum_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_el_E_iso_rel_Overlay = new TH1F("h_true_el_E_iso_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_true_el_E_iso_sum_noPh_0_10_Overlay = new TH1F("h_true_el_E_iso_sum_noPh_0_10_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_el_E_iso_noPh_0_10_rel_Overlay = new TH1F("h_true_el_E_iso_noPh_0_10_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_true_el_E_iso_sum_noPh_1_50_deg_Overlay = new TH1F("h_true_el_E_iso_sum_noPh_1_50_deg_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_el_E_iso_noPh_1_50_deg_rel_Overlay = new TH1F("h_true_el_E_iso_noPh_1_50_deg_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_true_mu_E_Overlay = new TH1F("h_true_mu_E_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_true_mu_E_iso_sum_Overlay = new TH1F("h_true_mu_E_iso_sum_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_mu_E_iso_rel_Overlay = new TH1F("h_true_mu_E_iso_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_true_mu_E_iso_sum_noPh_0_10_Overlay = new TH1F("h_true_mu_E_iso_sum_noPh_0_10_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_mu_E_iso_noPh_0_10_rel_Overlay = new TH1F("h_true_mu_E_iso_noPh_0_10_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_true_mu_E_iso_sum_noPh_1_50_deg_Overlay = new TH1F("h_true_mu_E_iso_sum_noPh_1_50_deg_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_true_mu_E_iso_noPh_1_50_deg_rel_Overlay = new TH1F("h_true_mu_E_iso_noPh_1_50_deg_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_true_l1l2_dangle_Overlay = new TH1F("h_true_l1l2_dangle_Overlay","", n_bins_high,0.0,180);

  TH1F* h_reco_el_true_reco_dangle_ang_match_Overlay = new TH1F("h_reco_el_true_reco_dangle_ang_match_Overlay","", n_bins_high,0.0,6.0);
  TH1F* h_reco_el_reco_E_over_true_E_ang_match_Overlay = new TH1F("h_reco_el_reco_E_over_true_E_ang_match_Overlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_reco_el_E_reco_ang_match_Overlay = new TH1F("h_reco_el_E_reco_ang_match_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_reco_el_true_reco_dangle_E_match_Overlay = new TH1F("h_reco_el_true_reco_dangle_E_match_Overlay","", n_bins_high,0.0,6.0);
  TH1F* h_reco_el_reco_E_over_true_E_match_Overlay = new TH1F("h_reco_el_reco_E_over_true_E_match_Overlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_reco_el_E_reco_E_match_Overlay = new TH1F("h_reco_el_E_reco_E_match_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_reco_mu_true_reco_dangle_ang_match_Overlay = new TH1F("h_reco_mu_true_reco_dangle_ang_match_Overlay","", n_bins_high,0.0,6.0);
  TH1F* h_reco_mu_reco_E_over_true_E_ang_match_Overlay = new TH1F("h_reco_mu_reco_E_over_true_E_ang_match_Overlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_reco_mu_E_reco_ang_match_Overlay = new TH1F("h_reco_mu_E_reco_ang_match_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_reco_mu_true_reco_dangle_E_match_Overlay = new TH1F("h_reco_mu_true_reco_dangle_E_match_Overlay","", n_bins_high,0.0,6.0);
  TH1F* h_reco_mu_reco_E_over_true_E_match_Overlay = new TH1F("h_reco_mu_reco_E_over_true_E_match_Overlay","", n_bins_high,lim_rel_energy_low,lim_rel_energy_high);
  TH1F* h_reco_mu_E_reco_E_match_Overlay = new TH1F("h_reco_mu_E_reco_E_match_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_reco_el_ang_match_E_iso_sum_Overlay = new TH1F("h_reco_el_ang_match_E_iso_sum_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_ang_match_E_iso_rel_Overlay = new TH1F("h_reco_el_ang_match_E_iso_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_el_ang_match_E_iso_sum_noPh_0_10_Overlay = new TH1F("h_reco_el_ang_match_E_iso_sum_noPh_0_10_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_ang_match_E_iso_noPh_0_10_rel_Overlay = new TH1F("h_reco_el_ang_match_E_iso_noPh_0_10_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_el_ang_match_E_iso_sum_noPh_1_50_deg_Overlay = new TH1F("h_reco_el_ang_match_E_iso_sum_noPh_1_50_deg_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_ang_match_E_iso_noPh_1_50_deg_rel_Overlay = new TH1F("h_reco_el_ang_match_E_iso_noPh_1_50_deg_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_reco_el_E_match_E_iso_sum_Overlay = new TH1F("h_reco_el_E_match_E_iso_sum_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_E_match_E_iso_rel_Overlay = new TH1F("h_reco_el_E_match_E_iso_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_el_E_match_E_iso_sum_noPh_0_10_Overlay = new TH1F("h_reco_el_E_match_E_iso_sum_noPh_0_10_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_E_match_E_iso_noPh_0_10_rel_Overlay = new TH1F("h_reco_el_E_match_E_iso_noPh_0_10_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_el_E_match_E_iso_sum_noPh_1_50_deg_Overlay = new TH1F("h_reco_el_E_match_E_iso_sum_noPh_1_50_deg_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_el_E_match_E_iso_noPh_1_50_deg_rel_Overlay = new TH1F("h_reco_el_E_match_E_iso_noPh_1_50_deg_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_reco_mu_ang_match_E_iso_sum_Overlay = new TH1F("h_reco_mu_ang_match_E_iso_sum_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_ang_match_E_iso_rel_Overlay = new TH1F("h_reco_mu_ang_match_E_iso_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_mu_ang_match_E_iso_sum_noPh_0_10_Overlay = new TH1F("h_reco_mu_ang_match_E_iso_sum_noPh_0_10_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_ang_match_E_iso_noPh_0_10_rel_Overlay = new TH1F("h_reco_mu_ang_match_E_iso_noPh_0_10_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_mu_ang_match_E_iso_sum_noPh_1_50_deg_Overlay = new TH1F("h_reco_mu_ang_match_E_iso_sum_noPh_1_50_deg_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_ang_match_E_iso_noPh_1_50_deg_rel_Overlay = new TH1F("h_reco_mu_ang_match_E_iso_noPh_1_50_deg_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_reco_mu_E_match_E_iso_sum_Overlay = new TH1F("h_reco_mu_E_match_E_iso_sum_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_E_match_E_iso_rel_Overlay = new TH1F("h_reco_mu_E_match_E_iso_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_mu_E_match_E_iso_sum_noPh_0_10_Overlay = new TH1F("h_reco_mu_E_match_E_iso_sum_noPh_0_10_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_E_match_E_iso_noPh_0_10_rel_Overlay = new TH1F("h_reco_mu_E_match_E_iso_noPh_0_10_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);
  TH1F* h_reco_mu_E_match_E_iso_sum_noPh_1_50_deg_Overlay = new TH1F("h_reco_mu_E_match_E_iso_sum_noPh_1_50_deg_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_mu_E_match_E_iso_noPh_1_50_deg_rel_Overlay = new TH1F("h_reco_mu_E_match_E_iso_noPh_1_50_deg_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high);

  TH1F* h_reco_fake_el_E_Overlay = new TH1F("h_reco_fake_el_E_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_fake_el_E_iso_sum_Overlay = new TH1F("h_reco_fake_el_E_iso_sum_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_el_E_iso_rel_Overlay = new TH1F("h_reco_fake_el_E_iso_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);
  TH1F* h_reco_fake_el_E_iso_sum_noPh_0_10_Overlay = new TH1F("h_reco_fake_el_E_iso_sum_noPh_0_10_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_el_E_iso_noPh_0_10_rel_Overlay = new TH1F("h_reco_fake_el_E_iso_noPh_0_10_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);
  TH1F* h_reco_fake_el_E_iso_sum_noPh_1_50_deg_Overlay = new TH1F("h_reco_fake_el_E_iso_sum_noPh_1_50_deg_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_el_E_iso_noPh_1_50_deg_rel_Overlay = new TH1F("h_reco_fake_el_E_iso_noPh_1_50_deg_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);

  TH1F* h_reco_fake_mu_E_Overlay = new TH1F("h_reco_fake_m_E_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_fake_mu_E_iso_sum_Overlay = new TH1F("h_reco_fake_mu_E_iso_sum_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_mu_E_iso_rel_Overlay = new TH1F("h_reco_fake_mu_E_iso_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);
  TH1F* h_reco_fake_mu_E_iso_sum_noPh_0_10_Overlay = new TH1F("h_reco_fake_mu_E_iso_sum_noPh_0_10_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_mu_E_iso_noPh_0_10_rel_Overlay = new TH1F("h_reco_fake_mu_E_iso_noPh_0_10_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);
  TH1F* h_reco_fake_mu_E_iso_sum_noPh_1_50_deg_Overlay = new TH1F("h_reco_fake_mu_E_iso_sum_noPh_1_50_deg_Overlay","", n_bins_high, lim_E_sum_iso_low,lim_E_sum_iso_high);
  TH1F* h_reco_fake_mu_E_iso_noPh_1_50_deg_rel_Overlay = new TH1F("h_reco_fake_mu_E_iso_noPh_1_50_deg_rel_Overlay","", n_bins_high, lim_E_sum_iso_rel_low,lim_E_sum_iso_rel_high_fake);

  TH1F* h_true_el_E_cosTheta_0_95_Overlay = new TH1F("h_true_el_E_cosTheta_0_95_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_el_E_reco_ang_match_rel_iso_0_01_cosTheta_0_95_Overlay = new TH1F("h_reco_el_E_reco_ang_match_rel_iso_0_01_cosTheta_0_95_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_el_E_reco_ang_match_rel_iso_0_02_cosTheta_0_95_Overlay = new TH1F("h_reco_el_E_reco_ang_match_rel_iso_0_02_cosTheta_0_95_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);

  TH1F* h_true_mu_E_cosTheta_0_95_Overlay = new TH1F("h_true_mu_E_cosTheta_0_95_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_mu_E_reco_ang_match_rmu_iso_0_01_cosTheta_0_95_Overlay = new TH1F("h_reco_mu_E_reco_ang_match_rel_iso_0_01_cosTheta_0_95_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_reco_mu_E_reco_ang_match_rmu_iso_0_02_cosTheta_0_95_Overlay = new TH1F("h_reco_mu_E_reco_ang_match_rel_iso_0_02_cosTheta_0_95_Overlay","", n_bins_high, lim_energy_low,lim_energy_high);


  std::vector<TH1F*> hist_vector_Overlay;

  hist_vector_Overlay.push_back( h_true_el_E_Overlay ); 
  hist_vector_Overlay.push_back( h_true_el_E_iso_sum_Overlay ); 
  hist_vector_Overlay.push_back( h_true_el_E_iso_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_true_el_E_iso_sum_noPh_0_10_Overlay ); 
  hist_vector_Overlay.push_back( h_true_el_E_iso_noPh_0_10_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_true_el_E_iso_sum_noPh_1_50_deg_Overlay ); 
  hist_vector_Overlay.push_back( h_true_el_E_iso_noPh_1_50_deg_rel_Overlay ); 

  hist_vector_Overlay.push_back( h_true_mu_E_Overlay ); 
  hist_vector_Overlay.push_back( h_true_mu_E_iso_sum_Overlay ); 
  hist_vector_Overlay.push_back( h_true_mu_E_iso_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_true_mu_E_iso_sum_noPh_0_10_Overlay ); 
  hist_vector_Overlay.push_back( h_true_mu_E_iso_noPh_0_10_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_true_mu_E_iso_sum_noPh_1_50_deg_Overlay ); 
  hist_vector_Overlay.push_back( h_true_mu_E_iso_noPh_1_50_deg_rel_Overlay ); 

  hist_vector_Overlay.push_back( h_true_l1l2_dangle_Overlay ); 

  hist_vector_Overlay.push_back( h_reco_el_true_reco_dangle_ang_match_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_reco_E_over_true_E_ang_match_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_E_reco_ang_match_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_true_reco_dangle_E_match_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_reco_E_over_true_E_match_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_E_reco_E_match_Overlay ); 

  hist_vector_Overlay.push_back( h_reco_mu_true_reco_dangle_ang_match_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_reco_E_over_true_E_ang_match_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_E_reco_ang_match_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_true_reco_dangle_E_match_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_reco_E_over_true_E_match_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_E_reco_E_match_Overlay ); 

  hist_vector_Overlay.push_back( h_reco_el_ang_match_E_iso_sum_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_ang_match_E_iso_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_ang_match_E_iso_sum_noPh_0_10_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_ang_match_E_iso_noPh_0_10_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_ang_match_E_iso_sum_noPh_1_50_deg_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_ang_match_E_iso_noPh_1_50_deg_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_E_match_E_iso_sum_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_E_match_E_iso_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_E_match_E_iso_sum_noPh_0_10_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_E_match_E_iso_noPh_0_10_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_E_match_E_iso_sum_noPh_1_50_deg_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_el_E_match_E_iso_noPh_1_50_deg_rel_Overlay ); 

  hist_vector_Overlay.push_back( h_reco_mu_ang_match_E_iso_sum_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_ang_match_E_iso_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_ang_match_E_iso_sum_noPh_0_10_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_ang_match_E_iso_noPh_0_10_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_ang_match_E_iso_sum_noPh_1_50_deg_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_ang_match_E_iso_noPh_1_50_deg_rel_Overlay ); 

  hist_vector_Overlay.push_back( h_reco_mu_E_match_E_iso_sum_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_E_match_E_iso_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_E_match_E_iso_sum_noPh_0_10_Overlay );
  hist_vector_Overlay.push_back( h_reco_mu_E_match_E_iso_noPh_0_10_rel_Overlay ); 
  hist_vector_Overlay.push_back( h_reco_mu_E_match_E_iso_sum_noPh_1_50_deg_Overlay );
  hist_vector_Overlay.push_back( h_reco_mu_E_match_E_iso_noPh_1_50_deg_rel_Overlay );

  hist_vector_Overlay.push_back( h_reco_fake_el_E_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_el_E_iso_sum_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_el_E_iso_rel_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_el_E_iso_sum_noPh_0_10_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_el_E_iso_noPh_0_10_rel_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_el_E_iso_sum_noPh_1_50_deg_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_el_E_iso_noPh_1_50_deg_rel_Overlay );
  
  hist_vector_Overlay.push_back( h_reco_fake_mu_E_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_mu_E_iso_sum_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_mu_E_iso_rel_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_mu_E_iso_sum_noPh_0_10_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_mu_E_iso_noPh_0_10_rel_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_mu_E_iso_sum_noPh_1_50_deg_Overlay );
  hist_vector_Overlay.push_back( h_reco_fake_mu_E_iso_noPh_1_50_deg_rel_Overlay );

  hist_vector_noOverlay.push_back( h_true_el_E_cosTheta_0_95_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_el_E_reco_ang_match_rel_iso_0_01_cosTheta_0_95_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_el_E_reco_ang_match_rel_iso_0_02_cosTheta_0_95_noOverlay );

  hist_vector_noOverlay.push_back( h_true_mu_E_cosTheta_0_95_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_mu_E_reco_ang_match_rmu_iso_0_01_cosTheta_0_95_noOverlay );
  hist_vector_noOverlay.push_back( h_reco_mu_E_reco_ang_match_rmu_iso_0_02_cosTheta_0_95_noOverlay ); 

  for(unsigned int i=0;i<hist_vector_Overlay.size();i++){
    hist_vector_Overlay[i]->Sumw2();
  }


 TEfficiency* tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay = new TEfficiency("tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay","", nbins_teff,xbins_teff);
 tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay->SetTitle(";true electron energy [GeV];electron ID efficiency");
 TEfficiency* tEff_ElectronVsTrueE_AngMatch_1_deg_rel_iso_0_20_CosTheta_True_0_95_Overlay = new TEfficiency("tEff_ElectronVsTrueE_AngMatch_1_deg_rel_iso_0_20_CosTheta_True_0_95_Overlay","", nbins_teff,xbins_teff);
 tEff_ElectronVsTrueE_AngMatch_1_deg_rel_iso_0_20_CosTheta_True_0_95_Overlay->SetTitle(";true electron energy [GeV];electron identification efficiency");
 TEfficiency* tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay = new TEfficiency("tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay","", nbins_teff,xbins_teff);
 tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay->SetTitle(";true muon energy [GeV];muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueE_AngMatch_1_deg_rel_iso_0_10_CosTheta_True_0_95_Overlay = new TEfficiency("tEff_MuonVsTrueE_AngMatch_1_deg_rel_iso_0_10_CosTheta_True_0_95_Overlay","", nbins_teff,xbins_teff);
 tEff_MuonVsTrueE_AngMatch_1_deg_rel_iso_0_10_CosTheta_True_0_95_Overlay->SetTitle(";true muon energy [GeV];muon identification efficiency");


 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay","", 20,-1.0,1.0);
 tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay->SetTitle(";cos(#theta(e^{true}));electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay","", 20,-1.,1.);
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay","", 20,0,180);
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay","", 20,-1.0,1.0);
 tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay->SetTitle(";cos(#theta(#mu^{true}));muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay","", 20,-1.,1.);
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay","", 20,0,180);

 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay","", 20,-1.0,1.0);
 tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay->SetTitle(";cos(#theta(e^{true}));electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay","", 20,-1.,1.);
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay","", 20,0,180);
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay","", 20,-1.0,1.0);
 tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay->SetTitle(";cos(#theta(#mu^{true}));muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay","", 20,-1.,1.);
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay","", 20,0,180);

 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay","", 20,-1.0,1.0);
 tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay->SetTitle(";cos(#theta(e^{true}));electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay","", 20,-1.,1.);
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay","", 20,0,180);
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay","", 20,-1.0,1.0);
 tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay->SetTitle(";cos(#theta(#mu^{true}));muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay","", 20,-1.,1.);
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay","", 20,0,180);

 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay","", 20,-1.0,1.0);
 tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay->SetTitle(";cos(#theta(e^{true}));electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay","", 20,-1.,1.);
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay","", 20,0,180);
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay","", 20,-1.0,1.0);
 tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay->SetTitle(";cos(#theta(#mu^{true}));muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay","", 20,-1.,1.);
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay","", 20,0,180);

 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay","", 20,-1.0,1.0);
 tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay->SetTitle(";cos(#theta(e^{true}));electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay = new TEfficiency("tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay","", 20,-1.,1.);
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay","", 20,0,180);
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay","", 20,-1.0,1.0);
 tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay->SetTitle(";cos(#theta(#mu^{true}));muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay = new TEfficiency("tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay","", 20,-1.,1.);
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay","", 20,0,180);

 //5/10/30/75/150/500/Inf
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_5_10_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_10_30_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_30_75_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_75_150_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_150_500_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_EMatch_E_True_500_Inf_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");

 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_5_10_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_10_30_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_30_75_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_75_150_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_150_500_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_EMatch_E_True_500_Inf_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");

 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_5_10_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_10_30_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_30_75_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_75_150_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_150_500_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");
 TEfficiency* tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_Overlay = new TEfficiency("tEff_ElectronVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_500_Inf_Overlay","", 20,0,180);
 tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_Overlay->SetTitle(";#theta(e^{true});electron identification efficiency");

 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_5_10_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_10_30_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_30_75_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_75_150_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_150_500_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");
 TEfficiency* tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_Overlay = new TEfficiency("tEff_MuonVsTrueTheta_AngMatch_1_deg_PtMatch_E_True_500_Inf_Overlay","", 20,0,180);
 tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_Overlay->SetTitle(";#theta(#mu^{true});muon identification efficiency");

  std::vector<TEfficiency*> TEff_vector_Overlay;
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueE_AngMatch_1_deg_rel_iso_0_20_CosTheta_True_0_95_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueE_AngMatch_1_deg_rel_iso_0_10_CosTheta_True_0_95_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_rel_iso_0_20_E_True_500_Inf_Overlay);

  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_Overlay);

  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_5_10_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_EMatch_500_Inf_Overlay);

  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_Overlay);

  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_5_10_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_10_30_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_30_75_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_75_150_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_150_500_Overlay);
  TEff_vector_Overlay.push_back(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_PtMatch_500_Inf_Overlay);

  fillHistograms(file_Overlay,hist_vector_Overlay,TEff_vector_Overlay);
  /*

  TCanvas *TTEffElectronVSETrueCanvas = new TCanvas("TTEffElectronVSETrueCanvas", "TTEffElectronVSETrueCanvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSETrueCanvas->cd();
  //TTEffElectronVSETrueCanvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSETrueCanvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSETrueCanvas->SetFillColor(0);
  TTEffElectronVSETrueCanvas->SetBorderMode(0);
  TTEffElectronVSETrueCanvas->SetBorderSize(2);
  TTEffElectronVSETrueCanvas->SetRightMargin(0.0172);
  TTEffElectronVSETrueCanvas->SetTopMargin(0.05);
  TTEffElectronVSETrueCanvas->SetBottomMargin(0.150);
  TTEffElectronVSETrueCanvas->SetFrameBorderMode(0);
  TTEffElectronVSETrueCanvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->Draw();
  gPad->Update();
  tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.68,1.02);
  tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsETrue_efficiency=new TLegend(0.31,0.20,0.60,0.40);
  leg_Electron_TTVsETrue_efficiency->SetBorderSize(0);
  leg_Electron_TTVsETrue_efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsETrue_efficiency->SetTextFont(42);
  leg_Electron_TTVsETrue_efficiency->SetTextSize(0.06);
  leg_Electron_TTVsETrue_efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsETrue_efficiency->SetLineColor(1);
  leg_Electron_TTVsETrue_efficiency->SetLineStyle(1);
  leg_Electron_TTVsETrue_efficiency->SetLineWidth(1);
  leg_Electron_TTVsETrue_efficiency->SetFillColor(0);
  leg_Electron_TTVsETrue_efficiency->SetFillStyle(1001);
  leg_Electron_TTVsETrue_efficiency->SetHeader("t#bar{t}, 3 TeV, |cos#theta(e^{true})|<0.95");
  leg_Electron_TTVsETrue_efficiency->AddEntry(tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay,"no background","l");
  leg_Electron_TTVsETrue_efficiency->AddEntry( tEff_ElectronVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsETrue_efficiency->Draw();
  
  TText *clicdp_label = new TText(-1.213543,1.044,"CLICdp work in progress");
  clicdp_label->SetTextSize(0.035);
  clicdp_label->SetTextFont(61);
  //clicdp_label->SetTextFont(43);
  //clicdp_label->SetBorderSize(0);
  //clicdp_label->AddText("CLICdp Preliminary");
  //clicdp_label->Draw();
  
  
  TCanvas *TTEffMuonVSETrueCanvas = new TCanvas("TTEffMuonVSETrueCanvas", "TTEffMuonVSETrueCanvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSETrueCanvas->cd();
  //TTEffMuonVSETrueCanvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSETrueCanvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSETrueCanvas->SetFillColor(0);
  TTEffMuonVSETrueCanvas->SetBorderMode(0);
  TTEffMuonVSETrueCanvas->SetBorderSize(2);
  TTEffMuonVSETrueCanvas->SetRightMargin(0.0172);
  TTEffMuonVSETrueCanvas->SetTopMargin(0.05);
  TTEffMuonVSETrueCanvas->SetBottomMargin(0.150);
  TTEffMuonVSETrueCanvas->SetFrameBorderMode(0);
  TTEffMuonVSETrueCanvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->Draw();
  gPad->Update();
  tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.90,1.02);
  tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsETrue_efficiency=new TLegend(0.31,0.20,0.60,0.40);
  leg_Muon_TTVsETrue_efficiency->SetBorderSize(0);
  leg_Muon_TTVsETrue_efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsETrue_efficiency->SetTextFont(42);
  leg_Muon_TTVsETrue_efficiency->SetTextSize(0.06);
  leg_Muon_TTVsETrue_efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsETrue_efficiency->SetLineColor(1);
  leg_Muon_TTVsETrue_efficiency->SetLineStyle(1);
  leg_Muon_TTVsETrue_efficiency->SetLineWidth(1);
  leg_Muon_TTVsETrue_efficiency->SetFillColor(0);
  leg_Muon_TTVsETrue_efficiency->SetFillStyle(1001);
  leg_Muon_TTVsETrue_efficiency->SetHeader("t#bar{t}, 3 TeV, |cos#theta(#mu^{true})|<0.95");
  leg_Muon_TTVsETrue_efficiency->AddEntry(tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_noOverlay,"no background","l");
  leg_Muon_TTVsETrue_efficiency->AddEntry( tEff_MuonVsTrueE_AngMatch_1_deg_CosTheta_True_0_95_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsETrue_efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffElectronVSCosThetaTrue_10_30_Canvas = new TCanvas("TTEffElectronVSCosThetaTrue_10_30_Canvas", "TTEffElectronVSCosThetaTrue_10_30_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSCosThetaTrue_10_30_Canvas->cd();
  //TTEffElectronVSCosThetaTrue_10_30_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSCosThetaTrue_10_30_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSCosThetaTrue_10_30_Canvas->SetFillColor(0);
  TTEffElectronVSCosThetaTrue_10_30_Canvas->SetBorderMode(0);
  TTEffElectronVSCosThetaTrue_10_30_Canvas->SetBorderSize(2);
  TTEffElectronVSCosThetaTrue_10_30_Canvas->SetRightMargin(0.0172);
  TTEffElectronVSCosThetaTrue_10_30_Canvas->SetTopMargin(0.05);
  TTEffElectronVSCosThetaTrue_10_30_Canvas->SetBottomMargin(0.150);
  TTEffElectronVSCosThetaTrue_10_30_Canvas->SetFrameBorderMode(0);
  TTEffElectronVSCosThetaTrue_10_30_Canvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->Draw();
  gPad->Update();
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.58,1.02);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsCosThetaTrue_10_30__efficiency=new TLegend(0.20,0.18,0.49,0.38);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetBorderSize(0);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetTextFont(42);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetTextSize(0.055);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetLineColor(1);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetLineStyle(1);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetLineWidth(1);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetFillColor(0);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetFillStyle(1001);
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->SetHeader("t#bar{t}, 3 TeV, 10 GeV<E(e^{true})<30 GeV");
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->AddEntry(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay,"no background","l");
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->AddEntry(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsCosThetaTrue_10_30__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffMuonVSCosThetaTrue_10_30_Canvas = new TCanvas("TTEffMuonVSCosThetaTrue_10_30_Canvas", "TTEffMuonVSCosThetaTrue_10_30_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSCosThetaTrue_10_30_Canvas->cd();
  //TTEffMuonVSCosThetaTrue_10_30_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSCosThetaTrue_10_30_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSCosThetaTrue_10_30_Canvas->SetFillColor(0);
  TTEffMuonVSCosThetaTrue_10_30_Canvas->SetBorderMode(0);
  TTEffMuonVSCosThetaTrue_10_30_Canvas->SetBorderSize(2);
  TTEffMuonVSCosThetaTrue_10_30_Canvas->SetRightMargin(0.0172);
  TTEffMuonVSCosThetaTrue_10_30_Canvas->SetTopMargin(0.05);
  TTEffMuonVSCosThetaTrue_10_30_Canvas->SetBottomMargin(0.150);
  TTEffMuonVSCosThetaTrue_10_30_Canvas->SetFrameBorderMode(0);
  TTEffMuonVSCosThetaTrue_10_30_Canvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->Draw();
  gPad->Update();
  //tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.70,1.02);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsCosThetaTrue_10_30__efficiency=new TLegend(0.20,0.45,0.47,0.65);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetBorderSize(0);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetTextFont(42);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetTextSize(0.055);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetLineColor(1);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetLineStyle(1);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetLineWidth(1);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetFillColor(0);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetFillStyle(1001);
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->SetHeader("t#bar{t}, 3 TeV,10 GeV <E(#mu^{true})<30 GeV");
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->AddEntry(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_noOverlay,"no background","l");
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->AddEntry(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_10_30_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsCosThetaTrue_10_30__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffElectronVSCosThetaTrue_30_75_Canvas = new TCanvas("TTEffElectronVSCosThetaTrue_30_75_Canvas", "TTEffElectronVSCosThetaTrue_30_75_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSCosThetaTrue_30_75_Canvas->cd();
  //TTEffElectronVSCosThetaTrue_30_75_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSCosThetaTrue_30_75_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSCosThetaTrue_30_75_Canvas->SetFillColor(0);
  TTEffElectronVSCosThetaTrue_30_75_Canvas->SetBorderMode(0);
  TTEffElectronVSCosThetaTrue_30_75_Canvas->SetBorderSize(2);
  TTEffElectronVSCosThetaTrue_30_75_Canvas->SetRightMargin(0.0172);
  TTEffElectronVSCosThetaTrue_30_75_Canvas->SetTopMargin(0.05);
  TTEffElectronVSCosThetaTrue_30_75_Canvas->SetBottomMargin(0.150);
  TTEffElectronVSCosThetaTrue_30_75_Canvas->SetFrameBorderMode(0);
  TTEffElectronVSCosThetaTrue_30_75_Canvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->Draw();
  gPad->Update();
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.58,1.02);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsCosThetaTrue_30_75__efficiency=new TLegend(0.20,0.18,0.49,0.38);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetBorderSize(0);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetTextFont(42);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetTextSize(0.055);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetLineColor(1);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetLineStyle(1);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetLineWidth(1);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetFillColor(0);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetFillStyle(1001);
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->SetHeader("t#bar{t}, 3 TeV, 30 GeV<E(e^{true})<75 GeV");
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->AddEntry(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay,"no background","l");
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->AddEntry(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsCosThetaTrue_30_75__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffMuonVSCosThetaTrue_30_75_Canvas = new TCanvas("TTEffMuonVSCosThetaTrue_30_75_Canvas", "TTEffMuonVSCosThetaTrue_30_75_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSCosThetaTrue_30_75_Canvas->cd();
  //TTEffMuonVSCosThetaTrue_30_75_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSCosThetaTrue_30_75_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSCosThetaTrue_30_75_Canvas->SetFillColor(0);
  TTEffMuonVSCosThetaTrue_30_75_Canvas->SetBorderMode(0);
  TTEffMuonVSCosThetaTrue_30_75_Canvas->SetBorderSize(2);
  TTEffMuonVSCosThetaTrue_30_75_Canvas->SetRightMargin(0.0172);
  TTEffMuonVSCosThetaTrue_30_75_Canvas->SetTopMargin(0.05);
  TTEffMuonVSCosThetaTrue_30_75_Canvas->SetBottomMargin(0.150);
  TTEffMuonVSCosThetaTrue_30_75_Canvas->SetFrameBorderMode(0);
  TTEffMuonVSCosThetaTrue_30_75_Canvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->Draw();
  gPad->Update();
  //tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.70,1.02);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsCosThetaTrue_30_75__efficiency=new TLegend(0.20,0.45,0.47,0.65);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetBorderSize(0);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetTextFont(42);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetTextSize(0.055);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetLineColor(1);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetLineStyle(1);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetLineWidth(1);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetFillColor(0);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetFillStyle(1001);
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->SetHeader("t#bar{t}, 3 TeV,30 GeV <E(#mu^{true})<75 GeV");
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->AddEntry(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_noOverlay,"no background","l");
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->AddEntry(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_30_75_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsCosThetaTrue_30_75__efficiency->Draw();
  //clicdp_label->Draw();


  TCanvas *TTEffElectronVSCosThetaTrue_75_150_Canvas = new TCanvas("TTEffElectronVSCosThetaTrue_75_150_Canvas", "TTEffElectronVSCosThetaTrue_75_150_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSCosThetaTrue_75_150_Canvas->cd();
  //TTEffElectronVSCosThetaTrue_75_150_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSCosThetaTrue_75_150_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSCosThetaTrue_75_150_Canvas->SetFillColor(0);
  TTEffElectronVSCosThetaTrue_75_150_Canvas->SetBorderMode(0);
  TTEffElectronVSCosThetaTrue_75_150_Canvas->SetBorderSize(2);
  TTEffElectronVSCosThetaTrue_75_150_Canvas->SetRightMargin(0.0172);
  TTEffElectronVSCosThetaTrue_75_150_Canvas->SetTopMargin(0.05);
  TTEffElectronVSCosThetaTrue_75_150_Canvas->SetBottomMargin(0.150);
  TTEffElectronVSCosThetaTrue_75_150_Canvas->SetFrameBorderMode(0);
  TTEffElectronVSCosThetaTrue_75_150_Canvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->Draw();
  gPad->Update();
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.58,1.02);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsCosThetaTrue_75_150__efficiency=new TLegend(0.20,0.18,0.49,0.38);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetBorderSize(0);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetTextFont(42);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetTextSize(0.055);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetLineColor(1);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetLineStyle(1);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetLineWidth(1);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetFillColor(0);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetFillStyle(1001);
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->SetHeader("t#bar{t}, 3 TeV, 75 GeV<E(e^{true})<150 GeV");
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->AddEntry(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay,"no background","l");
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->AddEntry(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsCosThetaTrue_75_150__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffMuonVSCosThetaTrue_75_150_Canvas = new TCanvas("TTEffMuonVSCosThetaTrue_75_150_Canvas", "TTEffMuonVSCosThetaTrue_75_150_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSCosThetaTrue_75_150_Canvas->cd();
  //TTEffMuonVSCosThetaTrue_75_150_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSCosThetaTrue_75_150_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSCosThetaTrue_75_150_Canvas->SetFillColor(0);
  TTEffMuonVSCosThetaTrue_75_150_Canvas->SetBorderMode(0);
  TTEffMuonVSCosThetaTrue_75_150_Canvas->SetBorderSize(2);
  TTEffMuonVSCosThetaTrue_75_150_Canvas->SetRightMargin(0.0172);
  TTEffMuonVSCosThetaTrue_75_150_Canvas->SetTopMargin(0.05);
  TTEffMuonVSCosThetaTrue_75_150_Canvas->SetBottomMargin(0.150);
  TTEffMuonVSCosThetaTrue_75_150_Canvas->SetFrameBorderMode(0);
  TTEffMuonVSCosThetaTrue_75_150_Canvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->Draw();
  gPad->Update();
  //tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.70,1.02);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsCosThetaTrue_75_150__efficiency=new TLegend(0.20,0.45,0.47,0.65);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetBorderSize(0);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetTextFont(42);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetTextSize(0.055);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetLineColor(1);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetLineStyle(1);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetLineWidth(1);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetFillColor(0);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetFillStyle(1001);
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->SetHeader("t#bar{t}, 3 TeV,75 GeV <E(#mu^{true})<150 GeV");
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->AddEntry(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_noOverlay,"no background","l");
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->AddEntry(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_75_150_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsCosThetaTrue_75_150__efficiency->Draw();
  //clicdp_label->Draw();

 TCanvas *TTEffElectronVSCosThetaTrue_150_500_Canvas = new TCanvas("TTEffElectronVSCosThetaTrue_150_500_Canvas", "TTEffElectronVSCosThetaTrue_150_500_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSCosThetaTrue_150_500_Canvas->cd();
  //TTEffElectronVSCosThetaTrue_150_500_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSCosThetaTrue_150_500_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSCosThetaTrue_150_500_Canvas->SetFillColor(0);
  TTEffElectronVSCosThetaTrue_150_500_Canvas->SetBorderMode(0);
  TTEffElectronVSCosThetaTrue_150_500_Canvas->SetBorderSize(2);
  TTEffElectronVSCosThetaTrue_150_500_Canvas->SetRightMargin(0.0172);
  TTEffElectronVSCosThetaTrue_150_500_Canvas->SetTopMargin(0.05);
  TTEffElectronVSCosThetaTrue_150_500_Canvas->SetBottomMargin(0.150);
  TTEffElectronVSCosThetaTrue_150_500_Canvas->SetFrameBorderMode(0);
  TTEffElectronVSCosThetaTrue_150_500_Canvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->Draw();
  gPad->Update();
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.58,1.02);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsCosThetaTrue_150_500__efficiency=new TLegend(0.20,0.18,0.49,0.38);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetBorderSize(0);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetTextFont(42);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetTextSize(0.055);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetLineColor(1);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetLineStyle(1);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetLineWidth(1);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetFillColor(0);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetFillStyle(1001);
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->SetHeader("t#bar{t}, 3 TeV, 150 GeV<E(e^{true})<500 GeV");
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->AddEntry(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay,"no background","l");
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->AddEntry(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsCosThetaTrue_150_500__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffMuonVSCosThetaTrue_150_500_Canvas = new TCanvas("TTEffMuonVSCosThetaTrue_150_500_Canvas", "TTEffMuonVSCosThetaTrue_150_500_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSCosThetaTrue_150_500_Canvas->cd();
  //TTEffMuonVSCosThetaTrue_150_500_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSCosThetaTrue_150_500_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSCosThetaTrue_150_500_Canvas->SetFillColor(0);
  TTEffMuonVSCosThetaTrue_150_500_Canvas->SetBorderMode(0);
  TTEffMuonVSCosThetaTrue_150_500_Canvas->SetBorderSize(2);
  TTEffMuonVSCosThetaTrue_150_500_Canvas->SetRightMargin(0.0172);
  TTEffMuonVSCosThetaTrue_150_500_Canvas->SetTopMargin(0.05);
  TTEffMuonVSCosThetaTrue_150_500_Canvas->SetBottomMargin(0.150);
  TTEffMuonVSCosThetaTrue_150_500_Canvas->SetFrameBorderMode(0);
  TTEffMuonVSCosThetaTrue_150_500_Canvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->Draw();
  gPad->Update();
  //tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.70,1.02);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsCosThetaTrue_150_500__efficiency=new TLegend(0.20,0.45,0.47,0.65);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetBorderSize(0);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetTextFont(42);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetTextSize(0.055);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetLineColor(1);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetLineStyle(1);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetLineWidth(1);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetFillColor(0);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetFillStyle(1001);
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->SetHeader("t#bar{t}, 3 TeV,150 GeV <E(#mu^{true})<500 GeV");
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->AddEntry(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_noOverlay,"no background","l");
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->AddEntry(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_150_500_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsCosThetaTrue_150_500__efficiency->Draw();
  //clicdp_label->Draw();

 TCanvas *TTEffElectronVSCosThetaTrue_500_Inf_Canvas = new TCanvas("TTEffElectronVSCosThetaTrue_500_Inf_Canvas", "TTEffElectronVSCosThetaTrue_500_Inf_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSCosThetaTrue_500_Inf_Canvas->cd();
  //TTEffElectronVSCosThetaTrue_500_Inf_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSCosThetaTrue_500_Inf_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSCosThetaTrue_500_Inf_Canvas->SetFillColor(0);
  TTEffElectronVSCosThetaTrue_500_Inf_Canvas->SetBorderMode(0);
  TTEffElectronVSCosThetaTrue_500_Inf_Canvas->SetBorderSize(2);
  TTEffElectronVSCosThetaTrue_500_Inf_Canvas->SetRightMargin(0.0172);
  TTEffElectronVSCosThetaTrue_500_Inf_Canvas->SetTopMargin(0.05);
  TTEffElectronVSCosThetaTrue_500_Inf_Canvas->SetBottomMargin(0.150);
  TTEffElectronVSCosThetaTrue_500_Inf_Canvas->SetFrameBorderMode(0);
  TTEffElectronVSCosThetaTrue_500_Inf_Canvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->Draw();
  gPad->Update();
  //tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.58,1.02);
  tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency=new TLegend(0.20,0.18,0.49,0.38);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetBorderSize(0);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetTextFont(42);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetTextSize(0.055);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetLineColor(1);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetLineStyle(1);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetLineWidth(1);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetFillColor(0);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetFillStyle(1001);
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->SetHeader("t#bar{t}, 3 TeV, E(e^{true})>500 GeV");
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->AddEntry(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay,"no background","l");
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->AddEntry(tEff_ElectronVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsCosThetaTrue_500_Inf__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffMuonVSCosThetaTrue_500_Inf_Canvas = new TCanvas("TTEffMuonVSCosThetaTrue_500_Inf_Canvas", "TTEffMuonVSCosThetaTrue_500_Inf_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSCosThetaTrue_500_Inf_Canvas->cd();
  //TTEffMuonVSCosThetaTrue_500_Inf_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSCosThetaTrue_500_Inf_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSCosThetaTrue_500_Inf_Canvas->SetFillColor(0);
  TTEffMuonVSCosThetaTrue_500_Inf_Canvas->SetBorderMode(0);
  TTEffMuonVSCosThetaTrue_500_Inf_Canvas->SetBorderSize(2);
  TTEffMuonVSCosThetaTrue_500_Inf_Canvas->SetRightMargin(0.0172);
  TTEffMuonVSCosThetaTrue_500_Inf_Canvas->SetTopMargin(0.05);
  TTEffMuonVSCosThetaTrue_500_Inf_Canvas->SetBottomMargin(0.150);
  TTEffMuonVSCosThetaTrue_500_Inf_Canvas->SetFrameBorderMode(0);
  TTEffMuonVSCosThetaTrue_500_Inf_Canvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->Draw();
  gPad->Update();
  //tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.70,1.02);
  tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency=new TLegend(0.20,0.45,0.47,0.65);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetBorderSize(0);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetTextFont(42);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetTextSize(0.055);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetLineColor(1);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetLineStyle(1);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetLineWidth(1);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetFillColor(0);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetFillStyle(1001);
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->SetHeader("t#bar{t}, 3 TeV,E(#mu^{true})>500 GeV");
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->AddEntry(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay,"no background","l");
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->AddEntry(tEff_MuonVsTrueCosTheta_AngMatch_1_deg_E_True_500_Inf_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsCosThetaTrue_500_Inf__efficiency->Draw();
  //clicdp_label->Draw();


  TCanvas *TTEffElectronVSThetaTrue_10_30_Canvas = new TCanvas("TTEffElectronVSThetaTrue_10_30_Canvas", "TTEffElectronVSThetaTrue_10_30_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSThetaTrue_10_30_Canvas->cd();
  //TTEffElectronVSThetaTrue_10_30_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSThetaTrue_10_30_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSThetaTrue_10_30_Canvas->SetFillColor(0);
  TTEffElectronVSThetaTrue_10_30_Canvas->SetBorderMode(0);
  TTEffElectronVSThetaTrue_10_30_Canvas->SetBorderSize(2);
  TTEffElectronVSThetaTrue_10_30_Canvas->SetRightMargin(0.0172);
  TTEffElectronVSThetaTrue_10_30_Canvas->SetTopMargin(0.05);
  TTEffElectronVSThetaTrue_10_30_Canvas->SetBottomMargin(0.150);
  TTEffElectronVSThetaTrue_10_30_Canvas->SetFrameBorderMode(0);
  TTEffElectronVSThetaTrue_10_30_Canvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->Draw();
  gPad->Update();
  //tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.45,1.02);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsThetaTrue_10_30__efficiency=new TLegend(0.20,0.18,0.49,0.38);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetBorderSize(0);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetTextFont(42);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetTextSize(0.055);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetLineColor(1);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetLineStyle(1);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetLineWidth(1);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetFillColor(0);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetFillStyle(1001);
  leg_Electron_TTVsThetaTrue_10_30__efficiency->SetHeader("t#bar{t}, 3 TeV, 10 GeV<E(e^{true})<30 GeV");
  leg_Electron_TTVsThetaTrue_10_30__efficiency->AddEntry(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay,"no background","l");
  leg_Electron_TTVsThetaTrue_10_30__efficiency->AddEntry(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsThetaTrue_10_30__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffMuonVSThetaTrue_10_30_Canvas = new TCanvas("TTEffMuonVSThetaTrue_10_30_Canvas", "TTEffMuonVSThetaTrue_10_30_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSThetaTrue_10_30_Canvas->cd();
  //TTEffMuonVSThetaTrue_10_30_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSThetaTrue_10_30_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSThetaTrue_10_30_Canvas->SetFillColor(0);
  TTEffMuonVSThetaTrue_10_30_Canvas->SetBorderMode(0);
  TTEffMuonVSThetaTrue_10_30_Canvas->SetBorderSize(2);
  TTEffMuonVSThetaTrue_10_30_Canvas->SetRightMargin(0.0172);
  TTEffMuonVSThetaTrue_10_30_Canvas->SetTopMargin(0.05);
  TTEffMuonVSThetaTrue_10_30_Canvas->SetBottomMargin(0.150);
  TTEffMuonVSThetaTrue_10_30_Canvas->SetFrameBorderMode(0);
  TTEffMuonVSThetaTrue_10_30_Canvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->Draw();
  gPad->Update();
  //tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.70,1.02);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsThetaTrue_10_30__efficiency=new TLegend(0.20,0.45,0.47,0.65);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetBorderSize(0);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetTextFont(42);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetTextSize(0.055);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetLineColor(1);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetLineStyle(1);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetLineWidth(1);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetFillColor(0);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetFillStyle(1001);
  leg_Muon_TTVsThetaTrue_10_30__efficiency->SetHeader("t#bar{t}, 3 TeV,10 GeV <E(#mu^{true})<30 GeV");
  leg_Muon_TTVsThetaTrue_10_30__efficiency->AddEntry(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_noOverlay,"no background","l");
  leg_Muon_TTVsThetaTrue_10_30__efficiency->AddEntry(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_10_30_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsThetaTrue_10_30__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffElectronVSThetaTrue_30_75_Canvas = new TCanvas("TTEffElectronVSThetaTrue_30_75_Canvas", "TTEffElectronVSThetaTrue_30_75_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSThetaTrue_30_75_Canvas->cd();
  //TTEffElectronVSThetaTrue_30_75_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSThetaTrue_30_75_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSThetaTrue_30_75_Canvas->SetFillColor(0);
  TTEffElectronVSThetaTrue_30_75_Canvas->SetBorderMode(0);
  TTEffElectronVSThetaTrue_30_75_Canvas->SetBorderSize(2);
  TTEffElectronVSThetaTrue_30_75_Canvas->SetRightMargin(0.0172);
  TTEffElectronVSThetaTrue_30_75_Canvas->SetTopMargin(0.05);
  TTEffElectronVSThetaTrue_30_75_Canvas->SetBottomMargin(0.150);
  TTEffElectronVSThetaTrue_30_75_Canvas->SetFrameBorderMode(0);
  TTEffElectronVSThetaTrue_30_75_Canvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->Draw();
  gPad->Update();
  //tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.45,1.02);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsThetaTrue_30_75__efficiency=new TLegend(0.20,0.18,0.49,0.38);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetBorderSize(0);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetTextFont(42);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetTextSize(0.055);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetLineColor(1);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetLineStyle(1);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetLineWidth(1);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetFillColor(0);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetFillStyle(1001);
  leg_Electron_TTVsThetaTrue_30_75__efficiency->SetHeader("t#bar{t}, 3 TeV, 30 GeV<E(e^{true})<75 GeV");
  leg_Electron_TTVsThetaTrue_30_75__efficiency->AddEntry(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay,"no background","l");
  leg_Electron_TTVsThetaTrue_30_75__efficiency->AddEntry(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsThetaTrue_30_75__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffMuonVSThetaTrue_30_75_Canvas = new TCanvas("TTEffMuonVSThetaTrue_30_75_Canvas", "TTEffMuonVSThetaTrue_30_75_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSThetaTrue_30_75_Canvas->cd();
  //TTEffMuonVSThetaTrue_30_75_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSThetaTrue_30_75_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSThetaTrue_30_75_Canvas->SetFillColor(0);
  TTEffMuonVSThetaTrue_30_75_Canvas->SetBorderMode(0);
  TTEffMuonVSThetaTrue_30_75_Canvas->SetBorderSize(2);
  TTEffMuonVSThetaTrue_30_75_Canvas->SetRightMargin(0.0172);
  TTEffMuonVSThetaTrue_30_75_Canvas->SetTopMargin(0.05);
  TTEffMuonVSThetaTrue_30_75_Canvas->SetBottomMargin(0.150);
  TTEffMuonVSThetaTrue_30_75_Canvas->SetFrameBorderMode(0);
  TTEffMuonVSThetaTrue_30_75_Canvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->Draw();
  gPad->Update();
  //tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.70,1.02);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsThetaTrue_30_75__efficiency=new TLegend(0.20,0.45,0.47,0.65);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetBorderSize(0);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetTextFont(42);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetTextSize(0.055);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetLineColor(1);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetLineStyle(1);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetLineWidth(1);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetFillColor(0);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetFillStyle(1001);
  leg_Muon_TTVsThetaTrue_30_75__efficiency->SetHeader("t#bar{t}, 3 TeV,30 GeV <E(#mu^{true})<75 GeV");
  leg_Muon_TTVsThetaTrue_30_75__efficiency->AddEntry(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_noOverlay,"no background","l");
  leg_Muon_TTVsThetaTrue_30_75__efficiency->AddEntry(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_30_75_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsThetaTrue_30_75__efficiency->Draw();
  //clicdp_label->Draw();


  TCanvas *TTEffElectronVSThetaTrue_75_150_Canvas = new TCanvas("TTEffElectronVSThetaTrue_75_150_Canvas", "TTEffElectronVSThetaTrue_75_150_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSThetaTrue_75_150_Canvas->cd();
  //TTEffElectronVSThetaTrue_75_150_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSThetaTrue_75_150_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSThetaTrue_75_150_Canvas->SetFillColor(0);
  TTEffElectronVSThetaTrue_75_150_Canvas->SetBorderMode(0);
  TTEffElectronVSThetaTrue_75_150_Canvas->SetBorderSize(2);
  TTEffElectronVSThetaTrue_75_150_Canvas->SetRightMargin(0.0172);
  TTEffElectronVSThetaTrue_75_150_Canvas->SetTopMargin(0.05);
  TTEffElectronVSThetaTrue_75_150_Canvas->SetBottomMargin(0.150);
  TTEffElectronVSThetaTrue_75_150_Canvas->SetFrameBorderMode(0);
  TTEffElectronVSThetaTrue_75_150_Canvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->Draw();
  gPad->Update();
  //tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.58,1.02);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsThetaTrue_75_150__efficiency=new TLegend(0.20,0.18,0.49,0.38);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetBorderSize(0);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetTextFont(42);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetTextSize(0.055);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetLineColor(1);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetLineStyle(1);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetLineWidth(1);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetFillColor(0);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetFillStyle(1001);
  leg_Electron_TTVsThetaTrue_75_150__efficiency->SetHeader("t#bar{t}, 3 TeV, 75 GeV<E(e^{true})<150 GeV");
  leg_Electron_TTVsThetaTrue_75_150__efficiency->AddEntry(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay,"no background","l");
  leg_Electron_TTVsThetaTrue_75_150__efficiency->AddEntry(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsThetaTrue_75_150__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffMuonVSThetaTrue_75_150_Canvas = new TCanvas("TTEffMuonVSThetaTrue_75_150_Canvas", "TTEffMuonVSThetaTrue_75_150_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSThetaTrue_75_150_Canvas->cd();
  //TTEffMuonVSThetaTrue_75_150_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSThetaTrue_75_150_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSThetaTrue_75_150_Canvas->SetFillColor(0);
  TTEffMuonVSThetaTrue_75_150_Canvas->SetBorderMode(0);
  TTEffMuonVSThetaTrue_75_150_Canvas->SetBorderSize(2);
  TTEffMuonVSThetaTrue_75_150_Canvas->SetRightMargin(0.0172);
  TTEffMuonVSThetaTrue_75_150_Canvas->SetTopMargin(0.05);
  TTEffMuonVSThetaTrue_75_150_Canvas->SetBottomMargin(0.150);
  TTEffMuonVSThetaTrue_75_150_Canvas->SetFrameBorderMode(0);
  TTEffMuonVSThetaTrue_75_150_Canvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->Draw();
  gPad->Update();
  //tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.70,1.02);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsThetaTrue_75_150__efficiency=new TLegend(0.20,0.45,0.47,0.65);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetBorderSize(0);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetTextFont(42);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetTextSize(0.055);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetLineColor(1);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetLineStyle(1);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetLineWidth(1);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetFillColor(0);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetFillStyle(1001);
  leg_Muon_TTVsThetaTrue_75_150__efficiency->SetHeader("t#bar{t}, 3 TeV,75 GeV <E(#mu^{true})<150 GeV");
  leg_Muon_TTVsThetaTrue_75_150__efficiency->AddEntry(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_noOverlay,"no background","l");
  leg_Muon_TTVsThetaTrue_75_150__efficiency->AddEntry(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_75_150_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsThetaTrue_75_150__efficiency->Draw();
  //clicdp_label->Draw();

 TCanvas *TTEffElectronVSThetaTrue_150_500_Canvas = new TCanvas("TTEffElectronVSThetaTrue_150_500_Canvas", "TTEffElectronVSThetaTrue_150_500_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSThetaTrue_150_500_Canvas->cd();
  //TTEffElectronVSThetaTrue_150_500_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSThetaTrue_150_500_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSThetaTrue_150_500_Canvas->SetFillColor(0);
  TTEffElectronVSThetaTrue_150_500_Canvas->SetBorderMode(0);
  TTEffElectronVSThetaTrue_150_500_Canvas->SetBorderSize(2);
  TTEffElectronVSThetaTrue_150_500_Canvas->SetRightMargin(0.0172);
  TTEffElectronVSThetaTrue_150_500_Canvas->SetTopMargin(0.05);
  TTEffElectronVSThetaTrue_150_500_Canvas->SetBottomMargin(0.150);
  TTEffElectronVSThetaTrue_150_500_Canvas->SetFrameBorderMode(0);
  TTEffElectronVSThetaTrue_150_500_Canvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->Draw();
  gPad->Update();
  //tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.58,1.02);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsThetaTrue_150_500__efficiency=new TLegend(0.20,0.18,0.49,0.38);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetBorderSize(0);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetTextFont(42);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetTextSize(0.055);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetLineColor(1);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetLineStyle(1);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetLineWidth(1);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetFillColor(0);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetFillStyle(1001);
  leg_Electron_TTVsThetaTrue_150_500__efficiency->SetHeader("t#bar{t}, 3 TeV, 150 GeV<E(e^{true})<500 GeV");
  leg_Electron_TTVsThetaTrue_150_500__efficiency->AddEntry(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay,"no background","l");
  leg_Electron_TTVsThetaTrue_150_500__efficiency->AddEntry(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsThetaTrue_150_500__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffMuonVSThetaTrue_150_500_Canvas = new TCanvas("TTEffMuonVSThetaTrue_150_500_Canvas", "TTEffMuonVSThetaTrue_150_500_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSThetaTrue_150_500_Canvas->cd();
  //TTEffMuonVSThetaTrue_150_500_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSThetaTrue_150_500_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSThetaTrue_150_500_Canvas->SetFillColor(0);
  TTEffMuonVSThetaTrue_150_500_Canvas->SetBorderMode(0);
  TTEffMuonVSThetaTrue_150_500_Canvas->SetBorderSize(2);
  TTEffMuonVSThetaTrue_150_500_Canvas->SetRightMargin(0.0172);
  TTEffMuonVSThetaTrue_150_500_Canvas->SetTopMargin(0.05);
  TTEffMuonVSThetaTrue_150_500_Canvas->SetBottomMargin(0.150);
  TTEffMuonVSThetaTrue_150_500_Canvas->SetFrameBorderMode(0);
  TTEffMuonVSThetaTrue_150_500_Canvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->Draw();
  gPad->Update();
  //tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.70,1.02);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsThetaTrue_150_500__efficiency=new TLegend(0.20,0.45,0.47,0.65);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetBorderSize(0);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetTextFont(42);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetTextSize(0.055);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetLineColor(1);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetLineStyle(1);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetLineWidth(1);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetFillColor(0);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetFillStyle(1001);
  leg_Muon_TTVsThetaTrue_150_500__efficiency->SetHeader("t#bar{t}, 3 TeV,150 GeV <E(#mu^{true})<500 GeV");
  leg_Muon_TTVsThetaTrue_150_500__efficiency->AddEntry(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_noOverlay,"no background","l");
  leg_Muon_TTVsThetaTrue_150_500__efficiency->AddEntry(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_150_500_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsThetaTrue_150_500__efficiency->Draw();
  //clicdp_label->Draw();

 TCanvas *TTEffElectronVSThetaTrue_500_Inf_Canvas = new TCanvas("TTEffElectronVSThetaTrue_500_Inf_Canvas", "TTEffElectronVSThetaTrue_500_Inf_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffElectronVSThetaTrue_500_Inf_Canvas->cd();
  //TTEffElectronVSThetaTrue_500_Inf_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffElectronVSThetaTrue_500_Inf_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffElectronVSThetaTrue_500_Inf_Canvas->SetFillColor(0);
  TTEffElectronVSThetaTrue_500_Inf_Canvas->SetBorderMode(0);
  TTEffElectronVSThetaTrue_500_Inf_Canvas->SetBorderSize(2);
  TTEffElectronVSThetaTrue_500_Inf_Canvas->SetRightMargin(0.0172);
  TTEffElectronVSThetaTrue_500_Inf_Canvas->SetTopMargin(0.05);
  TTEffElectronVSThetaTrue_500_Inf_Canvas->SetBottomMargin(0.150);
  TTEffElectronVSThetaTrue_500_Inf_Canvas->SetFrameBorderMode(0);
  TTEffElectronVSThetaTrue_500_Inf_Canvas->SetFrameBorderMode(0);
  
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetFillColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetMarkerColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetLineColor(kRed);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetFillStyle(3004);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->Draw();
  gPad->Update();
  //tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.58,1.02);
  tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay->Draw("same");
  
  TLegend* leg_Electron_TTVsThetaTrue_500_Inf__efficiency=new TLegend(0.20,0.18,0.49,0.38);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetBorderSize(0);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetTextFont(42);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetTextSize(0.055);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetEntrySeparation(0.15);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetLineColor(1);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetLineStyle(1);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetLineWidth(1);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetFillColor(0);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetFillStyle(1001);
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->SetHeader("t#bar{t}, 3 TeV, E(e^{true})>500 GeV");
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->AddEntry(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay,"no background","l");
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->AddEntry(tEff_ElectronVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Electron_TTVsThetaTrue_500_Inf__efficiency->Draw();
  //clicdp_label->Draw();

  TCanvas *TTEffMuonVSThetaTrue_500_Inf_Canvas = new TCanvas("TTEffMuonVSThetaTrue_500_Inf_Canvas", "TTEffMuonVSThetaTrue_500_Inf_Canvas",0,0,800,700);
  gStyle->SetOptStat(0);
  TTEffMuonVSThetaTrue_500_Inf_Canvas->cd();
  //TTEffMuonVSThetaTrue_500_Inf_Canvas->Range(-186.894,-0.873515,1682.046,6.114605);
  //TTEffMuonVSThetaTrue_500_Inf_Canvas->Range(-1.701969,0.7303125,1.250787,1.061562);
  TTEffMuonVSThetaTrue_500_Inf_Canvas->SetFillColor(0);
  TTEffMuonVSThetaTrue_500_Inf_Canvas->SetBorderMode(0);
  TTEffMuonVSThetaTrue_500_Inf_Canvas->SetBorderSize(2);
  TTEffMuonVSThetaTrue_500_Inf_Canvas->SetRightMargin(0.0172);
  TTEffMuonVSThetaTrue_500_Inf_Canvas->SetTopMargin(0.05);
  TTEffMuonVSThetaTrue_500_Inf_Canvas->SetBottomMargin(0.150);
  TTEffMuonVSThetaTrue_500_Inf_Canvas->SetFrameBorderMode(0);
  TTEffMuonVSThetaTrue_500_Inf_Canvas->SetFrameBorderMode(0);
  
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetFillColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetMarkerColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetLineColor(kRed);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->SetFillStyle(3004);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->Draw();
  gPad->Update();
  //tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.70,1.02);
  tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay->Draw("same");
  
  TLegend* leg_Muon_TTVsThetaTrue_500_Inf__efficiency=new TLegend(0.20,0.45,0.47,0.65);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetBorderSize(0);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetFillStyle(0);
  //leg->SetBorderSize(1);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetTextFont(42);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetTextSize(0.055);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetEntrySeparation(0.15);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetLineColor(1);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetLineStyle(1);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetLineWidth(1);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetFillColor(0);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetFillStyle(1001);
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->SetHeader("t#bar{t}, 3 TeV,E(#mu^{true})>500 GeV");
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->AddEntry(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_noOverlay,"no background","l");
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->AddEntry(tEff_MuonVsTrueTheta_AngMatch_1_deg_E_True_500_Inf_Overlay,"with #gamma#gamma#rightarrow hadrons","l");
  leg_Muon_TTVsThetaTrue_500_Inf__efficiency->Draw();
  //clicdp_label->Draw();
  */

   file_histogram->Write();
   file_histogram->Close();
                                 
}
