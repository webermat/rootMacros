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
#include "TVector3.h"
#include <TVirtualFitter.h>
#include "TProfile.h"
#include "TColor.h"
#include <vector>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/IParamFunction.h"
#include "Math/WrappedFunction.h"
#include "Fit/ParameterSettings.h"

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

class  MyMETMinimizerFunction {
 public:
  double Evaluate(const double* params);
  TLorentzVector jet1; 
  TLorentzVector jet2;
  TLorentzVector isoPh; 
  bool flagValidError;
  std::vector<double> parameterValues;
  double minimumValue;

  void SetJet1Jet2IsoPh(TLorentzVector _jet1,TLorentzVector _jet2,TLorentzVector _isoPh);
  void PerformMinimization();
  std::vector<double> GetFinalParameters() const;
  bool GetValidErrorFlag() const;
  double GetMinimum() const;

  private:

};

std::vector<double> MyMETMinimizerFunction::GetFinalParameters() const
{
  return parameterValues;
}

bool MyMETMinimizerFunction::GetValidErrorFlag() const
{
  return flagValidError;
}
double MyMETMinimizerFunction::GetMinimum()  const
{
  return minimumValue;
}

void MyMETMinimizerFunction::PerformMinimization()
{
    ROOT::Math::Minimizer* pMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    pMinimizer->SetMaxFunctionCalls(100000);
    pMinimizer->SetMaxIterations(100000);
    pMinimizer->SetTolerance(0.001);
    pMinimizer->SetPrintLevel(1);

    ROOT::Math::Functor minimizerFunctor(this, &MyMETMinimizerFunction::Evaluate, 2);
    pMinimizer->SetFunction(minimizerFunctor);
    double step = 0.01;
   
    pMinimizer->SetVariable(0, "a", 0.001, step);
    //std::cout << "Initial parameter 'a' value is " << 0.001<< std::endl;
    pMinimizer->SetVariable(1, "b", 0.001, step);
    //std::cout << "Initial parameter 'b' value is " << 0.001<< std::endl;


    pMinimizer->SetVariableLimits(0, 0, 0.20);
    pMinimizer->SetVariableLimits(1, 0, 0.20);

    pMinimizer->Minimize();
    pMinimizer->Hesse();
    pMinimizer->SetValidError(true);
    const double *finalParameters(pMinimizer->X());
    //std::cout << "Reached minimum with met value of : " << pMinimizer->MinValue() << std::endl;

    ROOT::Fit::ParameterSettings parameterSettings;
    pMinimizer->GetVariableSettings(0, parameterSettings);
    //std::cout << "  => Parameter a, final value : " << finalParameters[0] << std::endl;
    pMinimizer->GetVariableSettings(1, parameterSettings);
    //std::cout << "  => Parameter b, final value : " << finalParameters[1] << std::endl;
    if(!parameterValues.empty()){
      parameterValues.clear();
    }
    parameterValues.push_back(finalParameters[0]);
    parameterValues.push_back(finalParameters[1]);
    flagValidError=pMinimizer->IsValidError();
    minimumValue=pMinimizer->MinValue();
    //std::cout<<"minimum value? "<<minimumValue<<"valid error flag?" <<flagValidError<<" params "<< parameterValues[0]<<"/"<<parameterValues[1]<<std::endl;
    //double errLow=-1;                                                                                                                                                                                              
    //double errUp=-1;                                                                                                                                                                                               
    //pMinimizer->GetMinosError(0,errLow,errUp);                                                                                                                                                                     
    //std::cout<<"minos errors 0 low up "<<errLow<<"/"<<errUp<<std::endl; 
    //pMinimizer->GetMinosError(1,errLow,errUp);                                                                                                                                                                     
    //std::cout<<"minos errors 1 low up "<<errLow<<"/"<<errUp<<std::endl; 

    delete pMinimizer;
}


void MyMETMinimizerFunction::SetJet1Jet2IsoPh(TLorentzVector _jet1,TLorentzVector _jet2,TLorentzVector _isoPh){
  jet1 = _jet1;
  jet2 = _jet2;
  isoPh = _isoPh;
}

double MyMETMinimizerFunction::Evaluate(const double* params){
  //std::cout<<"jet1 px/py/ jet 2 px/py, iso ph px/py "<<jet1.Px()<<"/"<<jet1.Py()<<"/"<<jet2.Px()<<"/"<<jet2.Py()<<"/"<<isoPh.Px()<<"/"<<isoPh.Py()<<std::endl;
     double px=(1+params[0])*jet1.Px()+(1+params[1])*jet2.Px()-isoPh.Px();
     double py=(1+params[0])*jet1.Py()+(1+params[1])*jet2.Py()-isoPh.Py();
     
     return sqrt(px*px+py*py);
}


void fill_HZ_histograms(TFile* file, std::vector<TH1F*> h_hist_parton, std::vector<TH1F*> h_hist_vec_gen, std::vector<TH1F*> h_hist_vec_reco, std::vector<TH2F*> h_hist_vec_2D, std::vector<TH1F*> hist_vec_reco_1D_reco_vs_gen_selection, std::vector<TH2F*> hist_vec_reco_2D_reco_vs_gen_selection,bool usePartonInfo ,double x_sec, bool fill_partonInfo, bool fill_genInfo){

  std::cout<<"size of histogram vectors "<<h_hist_parton.size()<<"/"<<h_hist_vec_gen.size()<<"/"<<h_hist_vec_reco.size()<<"/"<<h_hist_vec_2D.size()<<std::endl;

  MyMETMinimizerFunction metMinimizerFunctionClass;

  TTree* tree = (TTree*)file->Get("showerData");

  bool use_EMissNeutrinoProjection=true;//here use total 4 vector, or 4 vector sum of jet 1 and 2 (see flag below)- isolated photon four vector plus correction with EMiss on both jet axes
  //mass cuts are then also done after projecting the EMiss
  bool use_sqrtJets=true;//in this case use j1+j2, isolated photon is anyway NOT in jets and with upper flag still decide if EMiss projection on jets is performed
  bool use_MHMiss_over_PFOMiss = false;//if recoil jet missing energy or PFO missing energy is used in the missing energy projection
  bool performMassCuts=true;

  int counter_j1j2_gen_2000=0;
  int counter_j1j2_gen_low_2000=0;

  int counter_j1j2_met_over_gen_2000=0;
  int counter_j1j2_met_over_gen_low_2000=0;

  int counter_j1j2_reco_2000=0;
  int counter_j1j2_reco_low_2000=0;

  int counter_j1j2_met_over_reco_2000=0;
  int counter_j1j2_met_over_reco_low_2000=0;

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
  
  float cut_m_j1_c=126;
  float cut_delta_m_j1_e = 20;

  float cut_m_j2_c=95;
  float cut_delta_m_j2_e = 25;

  float cut_m_j1_c_veto=86;
  float cut_delta_m_j1_e_veto = 40;

  float cut_m_j2_c_veto=84;
  float cut_delta_m_j2_e_veto_plus = 40;
  float cut_delta_m_j2_e_veto_minus = 40;
  

  float cut_theta_window=40;//in degrees --> window around 90 degrees

  float cut_m_j1_min=0.;
  float cut_m_j1_max = 17500;
  /*
  float cut_m_j2_min=0.;
  float cut_m_j2_max = 14000;
  */
  //not applied right now, go to ellipsoid cuts
  float cut_delta_mass_low=20;
  float cut_delta_mass_high=60;
  /*
  float cut_m_j1_min=120.;
  float cut_m_j1_max = 160;
  */
  float cut_m_j2_min=55.;
  float cut_m_j2_max =135;
  
  
  //number of allowed isolated leptons (+1), isolation criteria relative 10% in cone of 10 degrees, lower energy limit on leptons and photons set to 10 GeV
  unsigned int m_cut_nLeptons = 1;

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
  vector<int> *recojet_NNH=0;
  vector<int> *recojet_Mult=0;

  vector<float> *genjet_Px=0;
  vector<float> *genjet_Py=0;
  vector<float> *genjet_Pz=0;
  vector<float> *genjet_E=0;
  vector<float> *genjet_CHFraction=0;
  vector<float> *genjet_PhFraction=0;
  vector<int> *genjet_NCH=0;
  vector<int> *genjet_NPh=0;
  vector<int> *genjet_NNH=0;
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

  vector<float> *genjet_subjet_Px=0;
  vector<float> *genjet_subjet_Py=0;
  vector<float> *genjet_subjet_Pz=0;
  vector<float> *genjet_subjet_E=0;
  vector<int> *genjet_subjet_NCH=0;
  vector<int> *genjet_subjet_jetindex=0;

  vector<float> *recojet_subjet_Px=0;
  vector<float> *recojet_subjet_Py=0;
  vector<float> *recojet_subjet_Pz=0;
  vector<float> *recojet_subjet_E=0;
  vector<int> *recojet_subjet_NCH=0;
  vector<int> *recojet_subjet_jetindex=0;

  vector<float> *recojet_subjet_jetChargeE_kappa_0_25=0;
  vector<float> *recojet_subjet_jetChargeE_kappa_0_50=0;
  vector<float> *recojet_subjet_jetChargeE_kappa_0_20=0;
  vector<float> *recojet_subjet_jetChargeE_kappa_0_30=0;
  vector<float> *recojet_subjet_jetChargePt_kappa_0_25=0;
  vector<float> *recojet_subjet_jetChargePt_kappa_0_50=0;
  vector<float> *recojet_subjet_jetChargePt_kappa_0_20=0;
  vector<float> *recojet_subjet_jetChargePt_kappa_0_30=0;

  vector<float> *genjet_subjet_jetChargeE_kappa_0_25=0;
  vector<float> *genjet_subjet_jetChargeE_kappa_0_50=0;
  vector<float> *genjet_subjet_jetChargeE_kappa_0_20=0;
  vector<float> *genjet_subjet_jetChargeE_kappa_0_30=0;
  vector<float> *genjet_subjet_jetChargePt_kappa_0_25=0;
  vector<float> *genjet_subjet_jetChargePt_kappa_0_50=0;
  vector<float> *genjet_subjet_jetChargePt_kappa_0_20=0;
  vector<float> *genjet_subjet_jetChargePt_kappa_0_30=0;


  vector<float> *trueME_Px=0;
  vector<float> *trueME_Py=0;
  vector<float> *trueME_Pz=0;
  vector<float> *trueME_E=0;
  vector<int> *trueME_PDGID=0;
 
  //if(usePartonInfo){
  tree->SetBranchAddress("trueME_E", &trueME_E);
  tree->SetBranchAddress("trueME_Px", &trueME_Px);
  tree->SetBranchAddress("trueME_Py", &trueME_Py);
  tree->SetBranchAddress("trueME_Pz", &trueME_Pz);
  tree->SetBranchAddress("trueME_PDGID", &trueME_PDGID);
  //}
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
  tree->SetBranchAddress("genjet_NNH", &genjet_NNH);

  tree->SetBranchAddress("genjet_subjet_E", &genjet_subjet_E);
  tree->SetBranchAddress("genjet_subjet_Px", &genjet_subjet_Px);
  tree->SetBranchAddress("genjet_subjet_Py", &genjet_subjet_Py);
  tree->SetBranchAddress("genjet_subjet_Pz", &genjet_subjet_Pz);
  tree->SetBranchAddress("genjet_subjet_jetindex", &genjet_subjet_jetindex);
  tree->SetBranchAddress("genjet_subjet_NCH_trackPtMin", &genjet_subjet_NCH);

  tree->SetBranchAddress("genjet_subjet_jetChargeE_kappa_0_25", &genjet_subjet_jetChargeE_kappa_0_25);
  tree->SetBranchAddress("genjet_subjet_jetChargeE_kappa_0_50", &genjet_subjet_jetChargeE_kappa_0_50);
  tree->SetBranchAddress("genjet_subjet_jetChargeE_kappa_0_20", &genjet_subjet_jetChargeE_kappa_0_20);
  tree->SetBranchAddress("genjet_subjet_jetChargeE_kappa_0_30", &genjet_subjet_jetChargeE_kappa_0_30);
  tree->SetBranchAddress("genjet_subjet_jetChargePt_kappa_0_25", &genjet_subjet_jetChargePt_kappa_0_25);
  tree->SetBranchAddress("genjet_subjet_jetChargePt_kappa_0_50", &genjet_subjet_jetChargePt_kappa_0_50);
  tree->SetBranchAddress("genjet_subjet_jetChargePt_kappa_0_20", &genjet_subjet_jetChargePt_kappa_0_20);
  tree->SetBranchAddress("genjet_subjet_jetChargePt_kappa_0_30", &genjet_subjet_jetChargePt_kappa_0_30);


  tree->SetBranchAddress("recojet_E", &recojet_E);
  tree->SetBranchAddress("recojet_Px", &recojet_Px);
  tree->SetBranchAddress("recojet_Py", &recojet_Py);
  tree->SetBranchAddress("recojet_Pz", &recojet_Pz);
  tree->SetBranchAddress("recojet_CHFraction_trackPtMin", &recojet_CHFraction);
  tree->SetBranchAddress("recojet_PhFraction", &recojet_PhFraction);
  tree->SetBranchAddress("recojet_NCH_trackPtMin", &recojet_NCH);
  tree->SetBranchAddress("recojet_NPh", &recojet_NPh);
  tree->SetBranchAddress("recojet_Mult", &recojet_Mult);
  tree->SetBranchAddress("recojet_NNH", &recojet_NNH);

  tree->SetBranchAddress("recojet_subjet_E", &recojet_subjet_E);
  tree->SetBranchAddress("recojet_subjet_Px", &recojet_subjet_Px);
  tree->SetBranchAddress("recojet_subjet_Py", &recojet_subjet_Py);
  tree->SetBranchAddress("recojet_subjet_Pz", &recojet_subjet_Pz);
  tree->SetBranchAddress("recojet_subjet_NCH", &recojet_subjet_NCH);
  tree->SetBranchAddress("recojet_subjet_jetindex", &recojet_subjet_jetindex);

  tree->SetBranchAddress("recojet_subjet_jetChargeE_kappa_0_25", &recojet_subjet_jetChargeE_kappa_0_25);
  tree->SetBranchAddress("recojet_subjet_jetChargeE_kappa_0_50", &recojet_subjet_jetChargeE_kappa_0_50);
  tree->SetBranchAddress("recojet_subjet_jetChargeE_kappa_0_20", &recojet_subjet_jetChargeE_kappa_0_20);
  tree->SetBranchAddress("recojet_subjet_jetChargeE_kappa_0_30", &recojet_subjet_jetChargeE_kappa_0_30);
  tree->SetBranchAddress("recojet_subjet_jetChargePt_kappa_0_25", &recojet_subjet_jetChargePt_kappa_0_25);
  tree->SetBranchAddress("recojet_subjet_jetChargePt_kappa_0_50", &recojet_subjet_jetChargePt_kappa_0_50);
  tree->SetBranchAddress("recojet_subjet_jetChargePt_kappa_0_20", &recojet_subjet_jetChargePt_kappa_0_20);
  tree->SetBranchAddress("recojet_subjet_jetChargePt_kappa_0_30", &recojet_subjet_jetChargePt_kappa_0_30);
 


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



  double alpha_min =200;

  double lumi = 4000.;//cross sections are in femtobarn, at 3 TeV we expect 4 inverse attobarn 

  double weight=x_sec*lumi/(double)tree->GetEntries();
  std::cout<<"weight for sample "<<weight<<std::endl;


  const char *minName = "Minuit2";
  const char *algoName = "Migrad";

  ROOT::Math::Minimizer* myMinuit = ROOT::Math::Factory::CreateMinimizer(minName,algoName);

  // set tolerance , etc...
  myMinuit->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
  myMinuit->SetMaxIterations(1000000);  // for GSL
  myMinuit->SetTolerance(0.001);
  myMinuit->SetPrintLevel(1);
  
  int npar_METFunct = 2;


  for(unsigned int i_entry=0;i_entry<tree->GetEntries();i_entry++){
    //for(unsigned int i_entry=0;i_entry<10;i_entry++){
    //fill jet energy resolution histograms
    tree->GetEntry(i_entry);

    if(i_entry%5000==0){
      std::cout<<"entry "<<i_entry<<std::endl;
    }
    bool gen_pass_mass_cuts=false;
    bool reco_pass_mass_cuts=false;

    TLorentzVector tempTotEventP4(0,0,0,0);
    TLorentzVector tempTotEventP4HZ(0,0,0,0);
    TLorentzVector tempZP4(0,0,0,0);
    TLorentzVector tempHP4(0,0,0,0);

    TLorentzVector tempZ_q1(0,0,0,0);
    TLorentzVector tempZ_q2(0,0,0,0);

    TLorentzVector tempH_q1(0,0,0,0);
    TLorentzVector tempH_q2(0,0,0,0);

    TLorentzVector tempH_b(0,0,0,0);
    TLorentzVector tempH_bbar(0,0,0,0);
    TLorentzVector tempZ_q_pos(0,0,0,0);
    TLorentzVector tempZ_q_neg(0,0,0,0);


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


    int ind_sj1_rj1=-1;
    float E_sj1_rj1=-1;//jets forced into two subjets, thus only one check necessary
    int ind_sj2_rj1=-1;
    int ind_sj1_rj2=-1;
    float E_sj1_rj2=-1;
    int ind_sj2_rj2=-1;

    double angle_min_sj_rj_Z_q_pos=TMath::Pi();
    int ind_Z_q_pos_sj_rj_angle=-1;
    double angle_min_sj_rj_Z_q_neg=TMath::Pi();
    int ind_Z_q_neg_sj_rj_angle=-1;

    double angle_min_sj_gj_Z_q_pos=TMath::Pi();
    double angle_min_sj_gj_Z_q_neg=TMath::Pi();
    int ind_Z_q_pos_sj_gj_angle=-1;
    int ind_Z_q_neg_sj_gj_angle=-1;

    double angle_min_sj_gj_H_bbar=TMath::Pi();
    int ind_H_bbar_sj_gj_angle=-1;
    double angle_min_sj_gj_H_b=TMath::Pi();
    int ind_H_b_sj_gj_angle=-1;
    
    double angle_min_sj_rj_H_bbar=TMath::Pi();
    int ind_H_bbar_sj_rj_angle=-1;
    double angle_min_sj_rj_H_b=TMath::Pi();
    int ind_H_b_sj_rj_angle=-1;

    TLorentzVector temp_sj1_rj1(0,0,0,0);
    TLorentzVector temp_sj2_rj1(0,0,0,0);
    TLorentzVector temp_sj1_rj2(0,0,0,0);
    TLorentzVector temp_sj2_rj2(0,0,0,0);
    

    bool H_decays_bbar=true;
    if(usePartonInfo){
      for(unsigned int i=0;i<trueME_E->size();i++){
	TLorentzVector temp(0,0,0,0);
	temp.SetPxPyPzE((*trueME_Px)[i],(*trueME_Py)[i],(*trueME_Pz)[i],(*trueME_E)[i]);
	if(i<2){
	  //determined by outgoing e+e- pair
	  tempTotEventP4+=temp;
	}
      }
      if(fabs((*trueME_PDGID)[8])!=5 || fabs((*trueME_PDGID)[9])!=5){
	H_decays_bbar=false;
	//std::cout<<"entry "<<i_entry <<" is NO bbar "<<(*trueME_PDGID)[7]<<"/"<<(*trueME_PDGID)[8]<<std::endl;
	if((*trueME_PDGID)[8]>0){
	  tempH_b.SetPxPyPzE((*trueME_Px)[8],(*trueME_Py)[8],(*trueME_Pz)[8],(*trueME_E)[8]);
	  tempH_bbar.SetPxPyPzE((*trueME_Px)[9],(*trueME_Py)[9],(*trueME_Pz)[9],(*trueME_E)[9]);
	}else{
	  tempH_bbar.SetPxPyPzE((*trueME_Px)[8],(*trueME_Py)[8],(*trueME_Pz)[7],(*trueME_E)[8]);
	  tempH_b.SetPxPyPzE((*trueME_Px)[9],(*trueME_Py)[9],(*trueME_Pz)[9],(*trueME_E)[9]);
	}
      }else{
	if((*trueME_PDGID)[7]==5){
	  tempH_b.SetPxPyPzE((*trueME_Px)[8],(*trueME_Py)[8],(*trueME_Pz)[8],(*trueME_E)[8]);
	  tempH_bbar.SetPxPyPzE((*trueME_Px)[9],(*trueME_Py)[9],(*trueME_Pz)[9],(*trueME_E)[9]);
	}else{
	  tempH_bbar.SetPxPyPzE((*trueME_Px)[8],(*trueME_Py)[8],(*trueME_Pz)[8],(*trueME_E)[8]);
	  tempH_b.SetPxPyPzE((*trueME_Px)[9],(*trueME_Py)[9],(*trueME_Pz)[9],(*trueME_E)[9]);
	}
      }
    }
    if(!H_decays_bbar){
      continue;
    }
    
    if(usePartonInfo){//d=1,u=2,s=3,c=4,b=5,t=6
      int ind_Z_q_pos=4;
      int ind_Z_q_neg=5;
      //positive charge, aka up type quarks or down bar type quarks
      if((*trueME_PDGID)[6]==-1 || (*trueME_PDGID)[6]==2 || (*trueME_PDGID)[6]==-3 || (*trueME_PDGID)[6]==4 || (*trueME_PDGID)[6]==-5 || (*trueME_PDGID)[6]==6){
	tempZ_q_pos.SetPxPyPzE((*trueME_Px)[6],(*trueME_Py)[6],(*trueME_Pz)[6],(*trueME_E)[6]);
	tempZ_q_neg.SetPxPyPzE((*trueME_Px)[7],(*trueME_Py)[7],(*trueME_Pz)[7],(*trueME_E)[7]);
      }else{//quark index 4 is negatively charged
	tempZ_q_neg.SetPxPyPzE((*trueME_Px)[6],(*trueME_Py)[6],(*trueME_Pz)[6],(*trueME_E)[6]);
	tempZ_q_pos.SetPxPyPzE((*trueME_Px)[7],(*trueME_Py)[7],(*trueME_Pz)[7],(*trueME_E)[7]);
      }
      //6 and 7 are Z daugthers
      if((*trueME_E)[6]>(*trueME_E)[7]){
	tempZ_q1.SetPxPyPzE((*trueME_Px)[6],(*trueME_Py)[6],(*trueME_Pz)[6],(*trueME_E)[6]);
	tempZ_q2.SetPxPyPzE((*trueME_Px)[7],(*trueME_Py)[7],(*trueME_Pz)[7],(*trueME_E)[7]);
      }else{
	tempZ_q1.SetPxPyPzE((*trueME_Px)[7],(*trueME_Py)[7],(*trueME_Pz)[7],(*trueME_E)[7]);
	tempZ_q2.SetPxPyPzE((*trueME_Px)[6],(*trueME_Py)[6],(*trueME_Pz)[6],(*trueME_E)[6]);
      }
      //8 and 9 are H daugthers
      if((*trueME_E)[8]>(*trueME_E)[9]){
	tempH_q1.SetPxPyPzE((*trueME_Px)[8],(*trueME_Py)[8],(*trueME_Pz)[8],(*trueME_E)[8]);
	tempH_q2.SetPxPyPzE((*trueME_Px)[9],(*trueME_Py)[9],(*trueME_Pz)[9],(*trueME_E)[9]);
      }else{
	tempH_q1.SetPxPyPzE((*trueME_Px)[9],(*trueME_Py)[9],(*trueME_Pz)[9],(*trueME_E)[9]);
	tempH_q2.SetPxPyPzE((*trueME_Px)[8],(*trueME_Py)[8],(*trueME_Pz)[8],(*trueME_E)[8]);
      }
      tempHP4=tempH_q1+tempH_q2;
      tempZP4=tempZ_q1+tempZ_q2;
      //std::cout<<"before H p/bbar p/b p/sum bbbar p/z p/zbbbar p "<<tempHP4.P()<<"/"<<tempH_b.P()<<"/"<<tempH_bbar.P()<<"/"<<(tempH_b+tempH_bbar).P()<<"/"<<tempZP4.P()<<"/"<<(tempH_b+tempH_bbar+tempZP4).P()<<std::endl;
      TVector3 boostH_COM = -tempHP4.BoostVector();
      TVector3 boostZ_COM = -tempZP4.BoostVector();
      TLorentzVector tempHP4_boostH_COM=tempHP4;
      tempHP4_boostH_COM.Boost(boostH_COM);              
      //std::cout<<"H px/py/pz/ boost pz/py/pz "<<tempHP4.Px()<<"/"<<tempHP4.Py()<<"/"<<tempHP4.Pz()<<" "<<boostH_COM.Px()<<"/"<<boostH_COM.Py()<<"/"<<boostH_COM.Pz()<<std::endl;
      TLorentzVector tempH_b_boostH_COM=tempH_b;
      tempH_b_boostH_COM.Boost(boostH_COM);
      TLorentzVector tempH_bbar_boostH_COM=tempH_bbar;
      tempH_bbar_boostH_COM.Boost(boostH_COM);
      TLorentzVector tempZP4_boostZ_COM=tempZP4;
      tempZP4_boostZ_COM.Boost(boostZ_COM);
      TLorentzVector tempH_bbar_boostZ_COM=tempH_bbar;
      tempH_bbar_boostZ_COM.Boost(boostZ_COM);
      TLorentzVector tempH_b_boostZ_COM=tempH_b;
      tempH_b_boostZ_COM.Boost(boostZ_COM);
      TLorentzVector tempZP4_boostH_COM=tempZP4;
      tempZP4_boostH_COM.Boost(boostH_COM);
      TLorentzVector tempZ_q1_boostZ_COM=tempZ_q1;
      tempZ_q1_boostZ_COM.Boost(boostZ_COM);
      TLorentzVector tempZ_q2_boostZ_COM=tempZ_q2;
      tempZ_q2_boostZ_COM.Boost(boostZ_COM);
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
	 h_hist_parton[58]->Fill(tempH_b_boostH_COM.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[61]->Fill(tempH_bbar_boostH_COM.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[64]->Fill(cos(tempH_b_boostH_COM.Angle(tempHP4.Vect())),weight);
	 h_hist_parton[67]->Fill(cos(tempH_bbar_boostH_COM.Angle(tempHP4.Vect())),weight);
	 h_hist_parton[70]->Fill(cos(min(tempH_b_boostH_COM.Angle(tempHP4.Vect()),tempH_bbar_boostH_COM.Angle(tempHP4.Vect()))),weight);
	 h_hist_parton[73]->Fill(cos(max(tempH_b_boostH_COM.Angle(tempHP4.Vect()),tempH_bbar_boostH_COM.Angle(tempHP4.Vect()))),weight);
	 h_hist_parton[76]->Fill(cos(tempH_bbar_boostZ_COM.Angle(tempZP4.Vect())),weight);
	 h_hist_parton[79]->Fill(cos(tempH_b_boostZ_COM.Angle(tempZP4.Vect())),weight);
	 h_hist_parton[82]->Fill(cos(min(tempH_bbar_boostZ_COM.Angle(tempZP4.Vect()),tempH_b_boostZ_COM.Angle(tempZP4.Vect()))),weight);
	 h_hist_parton[85]->Fill(cos(max(tempH_bbar_boostZ_COM.Angle(tempZP4.Vect()),tempH_b_boostZ_COM.Angle(tempZP4.Vect()))),weight);
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
	 h_hist_parton[59]->Fill(tempH_b_boostH_COM.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[62]->Fill(tempH_bbar_boostH_COM.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);	 
	 h_hist_parton[65]->Fill(cos(tempH_b_boostH_COM.Angle(tempHP4.Vect())),weight);
	 h_hist_parton[68]->Fill(cos(tempH_bbar_boostH_COM.Angle(tempHP4.Vect())),weight);
	 h_hist_parton[71]->Fill(cos(min(tempH_b_boostH_COM.Angle(tempHP4.Vect()),tempH_bbar_boostH_COM.Angle(tempHP4.Vect()))),weight);
	 h_hist_parton[74]->Fill(cos(max(tempH_b_boostH_COM.Angle(tempHP4.Vect()),tempH_bbar_boostH_COM.Angle(tempHP4.Vect()))),weight);
	 h_hist_parton[77]->Fill(cos(tempH_bbar_boostZ_COM.Angle(tempZP4.Vect())),weight);
	 h_hist_parton[80]->Fill(cos(tempH_b_boostZ_COM.Angle(tempZP4.Vect())),weight);
	 h_hist_parton[83]->Fill(cos(min(tempH_bbar_boostZ_COM.Angle(tempZP4.Vect()),tempH_b_boostZ_COM.Angle(tempZP4.Vect()))),weight);
	 h_hist_parton[86]->Fill(cos(max(tempH_bbar_boostZ_COM.Angle(tempZP4.Vect()),tempH_b_boostZ_COM.Angle(tempZP4.Vect()))),weight);
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
	 h_hist_parton[60]->Fill(tempH_b_boostH_COM.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[63]->Fill(tempH_bbar_boostH_COM.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	 h_hist_parton[66]->Fill(cos(tempH_b_boostH_COM.Angle(tempHP4.Vect())),weight);
	 h_hist_parton[69]->Fill(cos(tempH_bbar_boostH_COM.Angle(tempHP4.Vect())),weight);
	 h_hist_parton[72]->Fill(cos(min(tempH_b_boostH_COM.Angle(tempHP4.Vect()),tempH_bbar_boostH_COM.Angle(tempHP4.Vect()))),weight);
	 h_hist_parton[75]->Fill(cos(max(tempH_b_boostH_COM.Angle(tempHP4.Vect()),tempH_bbar_boostH_COM.Angle(tempHP4.Vect()))),weight);
	 h_hist_parton[78]->Fill(cos(tempH_bbar_boostZ_COM.Angle(tempZP4.Vect())),weight);
	 h_hist_parton[81]->Fill(cos(tempH_b_boostZ_COM.Angle(tempZP4.Vect())),weight);
	 h_hist_parton[84]->Fill(cos(min(tempH_bbar_boostZ_COM.Angle(tempZP4.Vect()),tempH_b_boostZ_COM.Angle(tempZP4.Vect()))),weight);
	 h_hist_parton[87]->Fill(cos(max(tempH_bbar_boostZ_COM.Angle(tempZP4.Vect()),tempH_b_boostZ_COM.Angle(tempZP4.Vect()))),weight);
	 h_hist_parton[109]->Fill(tempHP4.E(),weight);
	 h_hist_parton[110]->Fill(tempZP4.E(),weight);
	 h_hist_parton[111]->Fill(tempHP4.E()-tempZP4.E(),weight);
	 h_hist_parton[112]->Fill((tempHP4.E()-tempZP4.E())/(tempHP4.E()+tempZP4.E()),weight);
       }
       h_hist_parton[108]->Fill((tempHP4+tempZP4).E(),weight);
     }//if parton info is to be filled      
    }

    TLorentzVector tempGenEMissCorrP4(0,0,0,0);
    TLorentzVector tempRecoEMissCorrP4(0,0,0,0);


    TLorentzVector tempRecoMHMissCorrP4(0,0,0,0);
    TLorentzVector tempGenMHMissCorrP4(0,0,0,0);

    TLorentzVector tempTotGenP4(0,0,0,
0);
    tempTotGenP4.SetPxPyPzE(true_Px,true_Py,true_Pz,true_E);

    TLorentzVector tempTotInvGenP4(0,0,0,0);
    //tempTotInvGenP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);
    //want to fake the reco value --> beam recoil changes it very slightly from true_inv vector
    tempTotInvGenP4.SetPxPyPzE(-true_Px,-true_Py,-true_Pz,sqrt(pow(true_Px,2)+pow(true_Py,2)+pow(true_Pz,2)));

    TLorentzVector tempRecoEMissP4(0,0,0,0);
    tempRecoEMissP4.SetPxPyPzE(-totPFO_Px,-totPFO_Py,-totPFO_Pz,sqrt(pow(totPFO_Px,2)+pow(totPFO_Py,2)+pow(totPFO_Pz,2)));

    TLorentzVector tempGenMETP4(0,0,0,0);
    tempGenMETP4.SetPxPyPzE(tempTotInvGenP4.Px(),tempTotInvGenP4.Py(),0,tempTotInvGenP4.Pt());
    TLorentzVector tempRecoMETP4(0,0,0,0);
    tempRecoMETP4.SetPxPyPzE(tempRecoEMissP4.Px(),tempRecoEMissP4.Py(),0,tempRecoEMissP4.Pt());

    TLorentzVector tempTotRecoP4(0,0,0,0);
    tempTotRecoP4.SetPxPyPzE(totPFO_Px,totPFO_Py,totPFO_Pz,totPFO_E);

    h_hist_vec_reco[83]->Fill(tempRecoEMissP4.Pt(),weight);

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
    double sqrtS_eff_gen=(tempTotGenP4-tempGenIsoPhP4).M();
    
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
    double sqrtS_eff_reco=(tempTotRecoP4-tempRecoIsoPhP4).M();
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
    TLorentzVector gj_m1_orig(0,0,0,0);
    TLorentzVector gj_m2_orig(0,0,0,0);
    TLorentzVector gj1_EMiss(0,0,0,0);
    TLorentzVector gj2_EMiss(0,0,0,0);
    TLorentzVector gj1_MHMiss(0,0,0,0);
    TLorentzVector gj2_MHMiss(0,0,0,0);
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

      //double METMinimizer(TLorentzVector jet1, TLorentzVector jet2, TLorentzVector isoPh, const double* par){
      //metMinimizerFunctionClass.SetJet1Jet2IsoPh(gj_m1,gj_m2, tempGenIsoPhP4);
      //double* par_gen;

      //double gj1_corr=0;
      //double gj2_corr=0;


      /*
      double met = (gj_m1+gj_m2-tempRecoIsoPhP4).Pt();
      if((gj_m1.E()+gj_m2.E())>2000 && (met/(gj_m1.E()+gj_m2.E()))<0.30){
	//std::cout<<"met gen value "<<met<<std::endl;
	metMinimizerFunctionClass.PerformMinimization();
	if(metMinimizerFunctionClass.GetFinalParameters()[0]<(0.20-1.e-3) && metMinimizerFunctionClass.GetFinalParameters()[1]<(0.20-1.e-3)){
	  gj1_corr=metMinimizerFunctionClass.GetFinalParameters()[0];
	  gj2_corr=metMinimizerFunctionClass.GetFinalParameters()[1];
	}
	//if(usePartonInfo && fill_partonInfo && fill_genInfo){
	if(genjet_E->size()==2 && n_IsoLep_gen<m_cut_nLeptons){
	  if(tempTotEventP4.M()>sqrtS_high){
	    h_hist_vec_2D[3]->Fill(metMinimizerFunctionClass.GetFinalParameters()[0],metMinimizerFunctionClass.GetFinalParameters()[1],weight);
	  }
	}
	  //}
      //h_hist_vec_2D[3]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenTruePhP4).M(),weight);
      //h_hist_vec_2D[4]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4).M(),weight);
	//if(metMinimizerFunctionClass.GetValidErrorFlag()){
	//std::cout<<"gen values for a and b "<<metMinimizerFunctionClass.GetFinalParameters()[0]<<"/"<<metMinimizerFunctionClass.GetFinalParameters()[1]<<" with min value "<<metMinimizerFunctionClass.GetMinimum()<<std::endl;
	//}else{//fill negative values
	//std::cout<<"gen values for a and b no valid error though "<<metMinimizerFunctionClass.GetFinalParameters()[0]<<"/"<<metMinimizerFunctionClass.GetFinalParameters()[1]<<" with min value "<<metMinimizerFunctionClass.GetMinimum()<<std::endl;
	//}
      }
 


      //ROOT::Minuit2::Minuit2Minimizer
 


      // create function wrapper for minmizer
      //ROOT::Math::Functor f(&MyMETMinimizerFunction::Evaluate,npar_METFunct);
      //ROOT::Math::Functor f(&MetFunctionTemplate,npar_METFunct);

      //myMinuit->SetFunction(f);
      */

      unsigned int ind_gj1_orig=ind_gj1;
      unsigned int ind_gj2_orig=ind_gj2;
      if(gj_m2.M()>=gj_m1.M()){
	std::cout<<"gj mass order should have not been the case nowadays "<<gj_m2.M()<<"/"<<gj_m1.M()<<std::endl;
      }
      //recoil missing energy based on both jets
      TLorentzVector tempGenMHMissP4=-gj_m1-gj_m2;
      TLorentzVector tempGenMHTP4(0,0,0,0);
      tempGenMHTP4.SetPxPyPzE(tempGenMHMissP4.Px(),tempGenMHMissP4.Py(),0.,tempGenMHMissP4.Pt());

      gj_m1_orig=gj_m1;
      gj_m2_orig=gj_m2;
      if((gj_m1_orig.E()+gj_m2_orig.E())>2000){
	counter_j1j2_gen_2000+=1;
      }else{
	counter_j1j2_gen_low_2000+=1;
      }
      TLorentzVector gj1_METProjVecProp(0,0,0,0);
      TLorentzVector gj1_EMissProjVecProp(0,0,0,0);
      TLorentzVector gj1_MHMissProjVecProp(0,0,0,0);
      TLorentzVector gj2_METProjVecProp(0,0,0,0);
      TLorentzVector gj2_EMissProjVecProp(0,0,0,0);
      TLorentzVector gj2_MHMissProjVecProp(0,0,0,0);
      TVector3 genmet_proj_sum(0,0,0);
      if(tempGenMETP4.Pt()>0){
	//TLorentzVector tempGenMETP4_gj1_Dir(0,0,0,0);
	double gj1_METProj=tempGenMETP4.Vect().Dot(gj_m1_orig.Vect().Unit());
	TLorentzVector gj1_METProjVec(0,0,0,0);
	gj1_METProjVec.SetPxPyPzE(gj1_METProj*gj_m1_orig.Px()/gj_m1_orig.P(),gj1_METProj*gj_m1_orig.Py()/gj_m1_orig.P(),0,gj1_METProj);
	//tempTotInvGenP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);
	TVector3 gj1_METProjProp=(tempGenMETP4.Vect().Dot(gj_m1_orig.Vect().Unit()))*gj_m1_orig.Vect().Unit();
	if(tempGenMETP4.Vect().Unit().Dot(gj_m1_orig.Vect().Unit())>0){//gj1 and MET in same hemisphere
	  genmet_proj_sum+=gj1_METProjProp;
	  gj1_METProjVecProp.SetPxPyPzE(gj1_METProjProp.Px(),gj1_METProjProp.Py(),0,gj1_METProjProp.Pt());
	  gj1_EMissProjVecProp.SetPxPyPzE(gj1_METProjProp.Px(),gj1_METProjProp.Py(),gj1_METProjProp.Pt()*gj_m1_orig.Pz()/gj_m1_orig.Pt(),gj1_METProjProp.Pt()*gj_m1_orig.P()/gj_m1_orig.Pt());
	}

	double gj2_METProj=tempGenMETP4.Vect().Dot(gj_m2_orig.Vect().Unit());
	TLorentzVector gj2_METProjVec(0,0,0,0);
	gj2_METProjVec.SetPxPyPzE(gj2_METProj*gj_m2_orig.Px()/gj_m2_orig.P(),gj2_METProj*gj_m2_orig.Py()/gj_m2_orig.P(),0,gj2_METProj);
	//tempTotInvGenP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);
	TVector3 gj2_METProjProp=(tempGenMETP4.Vect().Dot(gj_m2_orig.Vect().Unit()))*gj_m2_orig.Vect().Unit();
	if(tempGenMETP4.Vect().Unit().Dot(gj_m2_orig.Vect().Unit())>0){//gj2 and MET in same hemisphere
	  genmet_proj_sum+=gj2_METProjProp;
	  gj2_METProjVecProp.SetPxPyPzE(gj2_METProjProp.Px(),gj2_METProjProp.Py(),0,gj2_METProjProp.Pt());
	  gj2_EMissProjVecProp.SetPxPyPzE(gj2_METProjProp.Px(),gj2_METProjProp.Py(),gj2_METProjProp.Pt()*gj_m2_orig.Pz()/gj_m2_orig.Pt(),gj2_METProjProp.Pt()*gj_m2_orig.P()/gj_m2_orig.Pt());
	}
	if(genmet_proj_sum.Pt()>tempGenMETP4.Pt()){
	  //std::cout<<"danger, we have to do something with met projection gen "<<genmet_proj_sum.Pt()<<"/"<<tempGenMETP4.Pt()<<"/"<<sqrt(gj1_METProjProp.Pz()*gj1_METProjProp.Pz()+gj1_METProjProp.Pt()*gj1_METProjProp.Pt())<<"/"<<sqrt(gj2_METProjProp.Pz()*gj2_METProjProp.Pz()+gj2_METProjProp.Pt()*gj2_METProjProp.Pt())<<"/"<<gj1_METProjProp.Pt()<<"/"<<gj2_METProjProp.Pt()<<"/"<<acos(tempGenMETP4.Vect().Unit().Dot(gj_m1_orig.Vect().Unit()))*TMath::RadToDeg()<<"/"<<acos(tempGenMETP4.Vect().Unit().Dot(gj_m2_orig.Vect().Unit()))*TMath::RadToDeg()<<std::endl;
	}
      }


      if(genmet_proj_sum.Pt()>tempGenMETP4.Pt()){
	if((gj_m1_orig.E()+gj_m2_orig.E())>2000){
	  counter_j1j2_met_over_gen_2000+=1;
	}else{
	  counter_j1j2_met_over_gen_low_2000+=1;
	}
      }

      TLorentzVector gj1_MET=gj_m1_orig+gj1_METProjVecProp;
      TLorentzVector gj2_MET=gj_m2_orig+gj2_METProjVecProp;
      if(gj2_MET.M()>gj1_MET.M()){
	TLorentzVector temp=gj1_MET;
	gj1_MET=gj2_MET;
	gj2_MET=temp;
      }
      gj1_EMiss=gj_m1_orig+gj1_EMissProjVecProp;
      gj2_EMiss=gj_m2_orig+gj2_EMissProjVecProp;
      h_hist_vec_2D[3]->Fill(gj1_EMissProjVecProp.Pt()/gj_m1_orig.Pt(),gj2_EMissProjVecProp.Pt()/gj_m2_orig.Pt(),weight);
      unsigned int ind_gj1_EMiss=ind_gj1;
      unsigned int ind_gj2_EMiss=ind_gj2;
      if(gj2_EMiss.M()>gj1_EMiss.M()){
	//original jet ordering switches
	TLorentzVector temp=gj1_EMiss;
	gj1_EMiss=gj2_EMiss;
	gj2_EMiss=temp;
	ind_gj2_EMiss=ind_gj1;
	ind_gj1_EMiss=ind_gj2;
	//std::cout<<"original jet ordering switches before/after lead "<<ind_gj1_orig<<"/"<<ind_gj1_EMiss<<std::endl;
      }
      tempGenEMissCorrP4=gj1_EMiss+gj2_EMiss-gj_m1_orig-gj_m2_orig;

      TLorentzVector genMET2D(0,0,0,0);
      genMET2D.SetPxPyPzE(tempGenMETP4.Px(),tempGenMETP4.Py(),0,tempGenMETP4.Pt());
      TLorentzVector genMET2D_EMiss(0,0,0,0);
      genMET2D_EMiss.SetPxPyPzE((tempGenMETP4-tempGenEMissCorrP4).Px(), (tempGenMETP4-tempGenEMissCorrP4).Py(),0,(tempGenMETP4-tempGenEMissCorrP4).Pt());
      
      if(genMET2D_EMiss.Vect().Unit().Dot(genMET2D.Vect().Unit())<0){
	//std::cout<<"seems met and met corr change direction "<<acos(genMET2D_EMiss.Vect().Unit().Dot(genMET2D.Vect().Unit()))*TMath::RadToDeg()<<" orig/new "<< genMET2D.Pt()<<"/"<<genMET2D_EMiss.Pt()<<std::endl;
      }
      
      if(tempGenMHTP4.Pt()>0){
	TVector3 genmht_proj_sum(0,0,0);
	double gj1_MHTProj=tempGenMHTP4.Vect().Dot(gj_m1_orig.Vect().Unit());
	TLorentzVector gj1_MHTProjVec(0,0,0,0);
	gj1_MHTProjVec.SetPxPyPzE(gj1_MHTProj*gj_m1_orig.Px()/gj_m1_orig.P(),gj1_MHTProj*gj_m1_orig.Py()/gj_m1_orig.P(),0,gj1_MHTProj);
	//tempTotInvGenP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);
	TVector3 gj1_MHTProjProp=(tempGenMHTP4.Vect().Dot(gj_m1_orig.Vect().Unit()))*gj_m1_orig.Vect().Unit();
	if(tempGenMHTP4.Vect().Unit().Dot(gj_m1_orig.Vect().Unit())>0){//gj1 and MHT in same hemisphere
	  genmht_proj_sum+=gj1_MHTProjProp;
	  //scale pz with values obtained by MHT projection for px and py along the jet z momentum axis, assume massless neutrino vector addition, aka scale MHT with p/pt ratio
	  gj1_MHMissProjVecProp.SetPxPyPzE(gj1_MHTProjProp.Px(),gj1_MHTProjProp.Py(),gj1_MHTProjProp.Pt()*gj_m1_orig.Pz()/gj_m1_orig.Pt(),gj1_MHTProjProp.Pt()*gj_m1_orig.P()/gj_m1_orig.Pt());
	}
	double gj2_MHTProj=tempGenMHTP4.Vect().Dot(gj_m2_orig.Vect().Unit());
	TLorentzVector gj2_MHTProjVec(0,0,0,0);
	gj2_MHTProjVec.SetPxPyPzE(gj2_MHTProj*gj_m2_orig.Px()/gj_m2_orig.P(),gj2_MHTProj*gj_m2_orig.Py()/gj_m2_orig.P(),0,gj2_MHTProj);
	//tempTotInvGenP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);
	TVector3 gj2_MHTProjProp=(tempGenMHTP4.Vect().Dot(gj_m2_orig.Vect().Unit()))*gj_m2_orig.Vect().Unit();
	if(tempGenMHTP4.Vect().Unit().Dot(gj_m2_orig.Vect().Unit())>0){//gj2 and MHT in same hemisphere
	  genmht_proj_sum+=gj2_MHTProjProp;
	  gj2_MHMissProjVecProp.SetPxPyPzE(gj2_MHTProjProp.Px(),gj2_MHTProjProp.Py(),gj2_MHTProjProp.Pt()*gj_m2_orig.Pz()/gj_m2_orig.Pt(),gj2_MHTProjProp.Pt()*gj_m2_orig.P()/gj_m2_orig.Pt());
	}
	if(genmht_proj_sum.Pt()>tempGenMHTP4.Pt()){
	  //std::cout<<"danger, we have to do something with mht projection gen "<<genmht_proj_sum.Pt()<<"/"<<tempGenMHTP4.Pt()<<std::endl;
	}
      }
      gj1_MHMiss=gj_m1_orig+gj1_MHMissProjVecProp;
      gj2_MHMiss=gj_m2_orig+gj2_MHMissProjVecProp;
      //std::cout<<"minimizer/calcMH gj1 calc "<<gj1_corr<<"/"<<gj1_MHMissProjVecProp.Px()/gj_m1_orig.Px()<<"/"<<gj1_MHMissProjVecProp.Py()/gj_m1_orig.Py()<<std::endl;
      //std::cout<<"minimizer/calcMH gj2 calc "<<gj2_corr<<"/"<<gj2_MHMissProjVecProp.Px()/gj_m2_orig.Px()<<"/"<<gj2_MHMissProjVecProp.Py()/gj_m2_orig.Py()<<std::endl;
      unsigned int ind_gj1_MHMiss=ind_gj1;
      unsigned int ind_gj2_MHMiss=ind_gj2;
      if(gj2_MHMiss.M()>gj1_MHMiss.M()){
	TLorentzVector temp=gj1_MHMiss;
	gj1_MHMiss=gj2_MHMiss;
	gj2_MHMiss=temp;
	ind_gj1_MHMiss=ind_gj2;
	ind_gj2_MHMiss=ind_gj1;
      }
      tempGenMHMissCorrP4=gj1_MHMiss+gj2_MHMiss-gj_m1_orig-gj_m2_orig;
      
      if(use_sqrtJets){//for jets don't need to subtract photons, already done in preselection of jet filling
	sqrtS_eff_gen=(gj_m1_orig+gj_m2_orig).M();
      }
      

      if(use_EMissNeutrinoProjection){
	gj_m1=gj1_EMiss;
	gj_m2=gj2_EMiss;
	ind_gj1=ind_gj1_EMiss;
	ind_gj2=ind_gj2_EMiss;
	if(use_sqrtJets){//isolated photons and leptons NOT included in jet clustering
	  sqrtS_eff_gen=(gj_m1+gj_m2+tempGenEMissCorrP4).M();
	}else{//ALL stable particles included in total sum
	  sqrtS_eff_gen=(tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).M();
	}
      }
      if(use_MHMiss_over_PFOMiss){
	gj_m1=gj1_MHMiss;
	gj_m2=gj2_MHMiss;
	ind_gj1=ind_gj1_MHMiss;
	ind_gj2=ind_gj2_MHMiss;
	if(use_sqrtJets){
	  sqrtS_eff_gen=(gj_m1+gj_m2+tempGenMHMissCorrP4).M();
	}else{
	  sqrtS_eff_gen=(tempTotGenP4+tempGenMHMissCorrP4).M();
	}
      }
      //if(gj_m1.M()==gj_m1_orig.M()){
      //std::cout<<"mass before after the same "<<gj_m1_orig.M()<<std::endl;
      //}else{
      //std::cout<<"mass change there before/after "<<gj_m1_orig.M()<<"/"<<gj_m1.M()<<std::endl;
      //}
      
      if(usePartonInfo && fill_partonInfo && fill_genInfo){
	h_hist_vec_2D[17]->Fill(tempTotEventP4.M(),(gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).M(),weight);
	h_hist_vec_2D[18]->Fill(tempTotEventP4.M(),(gj1_EMiss+gj2_EMiss-tempRecoIsoPhP4).M(),weight);	
	if(fill_genInfo && genjet_E->size()==2 && n_IsoLep_gen<m_cut_nLeptons){
	  if(tempTotEventP4.M()<sqrtS_low){
	    h_hist_parton[128]->Fill((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).Pt()/(tempTotGenP4-tempGenIsoPhP4).E(),weight);
	    h_hist_parton[134]->Fill((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).Pt()/(gj_m1_orig.E()+gj_m2_orig.E()),weight);
	    h_hist_parton[140]->Fill((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).Pt(),weight);
	    h_hist_parton[146]->Fill((tempTotGenP4-tempGenIsoPhP4).E(),weight);
	    h_hist_parton[152]->Fill(gj_m1_orig.E()+gj_m2_orig.E(),weight);	     
	  }else if(tempTotEventP4.M()<sqrtS_high){
	    h_hist_parton[129]->Fill((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).Pt()/(tempTotGenP4-tempGenIsoPhP4).E(),weight);
	    h_hist_parton[135]->Fill((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).Pt()/(gj_m1_orig.E()+gj_m2_orig.E()),weight);
	    h_hist_parton[141]->Fill((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).Pt(),weight);
	    h_hist_parton[147]->Fill((tempTotGenP4-tempGenIsoPhP4).E(),weight);
	    h_hist_parton[153]->Fill(gj_m1_orig.E()+gj_m2_orig.E(),weight);	    
	  }else{
	    h_hist_parton[130]->Fill((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).Pt()/(tempTotGenP4-tempGenIsoPhP4).E(),weight);
	    h_hist_parton[136]->Fill((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).Pt()/(gj_m1_orig.E()+gj_m2_orig.E()),weight);
	    h_hist_parton[142]->Fill((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).Pt(),weight);
	    h_hist_parton[148]->Fill((tempTotGenP4-tempGenIsoPhP4).E(),weight);
	    h_hist_parton[154]->Fill(gj_m1_orig.E()+gj_m2_orig.E(),weight);	    
	  }
	}
	if(sqrtS_eff_gen>sqrtS_high){
	  h_hist_parton[88]->Fill(DeltaPhi(gj_m1.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	  h_hist_parton[90]->Fill(DeltaPhi(gj_m2.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	  if((gj_m1.Angle(tempHP4.Vect())*TMath::RadToDeg())<25.){//jet 1 and H spatially matched
	    h_hist_parton[92]->Fill(DeltaPhi(gj_m1.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	  }else if((gj_m1.Angle(tempHP4.Vect())*TMath::RadToDeg())>155.){//jet 1 and H spatially unmatched
	    h_hist_parton[94]->Fill(DeltaPhi(gj_m1.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	  }
	  if((gj_m2.Angle(tempZP4.Vect())*TMath::RadToDeg())<25.){//jet 2 and Z spatially matched
	    h_hist_parton[96]->Fill(DeltaPhi(gj_m2.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	  }else if((gj_m2.Angle(tempZP4.Vect())*TMath::RadToDeg())>155.){//jet 2 and Z spatially unmatched
	    h_hist_parton[98]->Fill(DeltaPhi(gj_m2.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	  }
	  if((gj_m1.Angle(tempHP4.Vect())*TMath::RadToDeg())<25.){//jet 1 and H spatially matched
	    h_hist_parton[100]->Fill(DeltaPhi(gj_m1.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	  }
	  if((gj_m2.Angle(tempHP4.Vect())*TMath::RadToDeg())<25.){//jet 2 and H spatially matched
	    h_hist_parton[100]->Fill(DeltaPhi(gj_m2.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	  }
	  if((gj_m1.Angle(tempZP4.Vect())*TMath::RadToDeg())<25.){//jet 1 and Z spatially matched
	    h_hist_parton[102]->Fill(DeltaPhi(gj_m1.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	  }
	  if((gj_m2.Angle(tempZP4.Vect())*TMath::RadToDeg())<25.){//jet 2 and Z spatially matched
	    h_hist_parton[102]->Fill(DeltaPhi(gj_m2.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	  }
	}
	if(tempTotEventP4.M()>sqrtS_high){
	  h_hist_parton[38]->Fill(gj_m1.M(),weight);
	  h_hist_parton[39]->Fill(gj_m1.M(),weight);
	  h_hist_parton[42]->Fill(DeltaPhi(gj_m1.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	  h_hist_parton[43]->Fill(fabs(gj_m1.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	  h_hist_parton[44]->Fill(DeltaPhi(gj_m1.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	  h_hist_parton[45]->Fill(fabs(gj_m1.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	  h_hist_parton[50]->Fill(gj_m1.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	  h_hist_parton[51]->Fill(gj_m2.Angle(tempZP4.Vect())*TMath::RadToDeg(),weight);
	  h_hist_parton[104]->Fill(gj1_EMiss.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	  h_hist_parton[105]->Fill(gj2_EMiss.Angle(tempZP4.Vect())*TMath::RadToDeg(),weight);
	}
      }

      
      //if too few components in jet, then jet index NOT found -->i.e. remains at -1
      for(unsigned int i=0;i<genjet_subjet_E->size();i++){
	if((*genjet_subjet_jetindex)[i]==ind_gj1){
	  if((*genjet_subjet_E)[i]>E_sj1_gj1){
	    ind_sj2_gj1=ind_sj1_gj1;
	    ind_sj1_gj1=i;
	    E_sj1_gj1=(*genjet_subjet_E)[i];
	  }else{
	    ind_sj2_gj1=i;
	  }
	}
	if((*genjet_subjet_jetindex)[i]==ind_gj2){
	  if((*genjet_subjet_E)[i]>E_sj1_gj2){
	    ind_sj2_gj2=ind_sj1_gj2;
	    ind_sj1_gj2=i;
	    E_sj1_gj2=(*genjet_subjet_E)[i];
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

      
      TVector3 boost_gj1_COM = -gj_m1_orig.BoostVector();
      TVector3 boost_gj2_COM = -gj_m2_orig.BoostVector();
      if(use_EMissNeutrinoProjection && ind_gj1!=ind_gj1_orig){
	boost_gj1_COM = -gj_m2_orig.BoostVector();
	boost_gj2_COM = -gj_m1_orig.BoostVector();
      }
      if(use_MHMiss_over_PFOMiss && ind_gj1!=ind_gj1_orig){
	boost_gj1_COM = -gj_m2_orig.BoostVector();
	boost_gj2_COM = -gj_m1_orig.BoostVector();
      }
      TLorentzVector temp_sj1_gj1_boost_gj1_COM=temp_sj1_gj1;
      temp_sj1_gj1_boost_gj1_COM.Boost(boost_gj1_COM);
      TLorentzVector temp_sj2_gj1_boost_gj1_COM=temp_sj2_gj1;
      temp_sj2_gj1_boost_gj1_COM.Boost(boost_gj1_COM);
      TLorentzVector temp_sj1_gj2_boost_gj2_COM=temp_sj1_gj2;
      temp_sj1_gj2_boost_gj2_COM.Boost(boost_gj2_COM);
      TLorentzVector temp_sj2_gj2_boost_gj2_COM=temp_sj2_gj2;
      temp_sj2_gj2_boost_gj2_COM.Boost(boost_gj2_COM);
      
      TLorentzVector temp_gj1_boost_gj1_COM=gj_m1;
      temp_gj1_boost_gj1_COM.Boost(boost_gj1_COM);
      TLorentzVector temp_gj2_boost_gj2_COM=gj_m2;
      temp_gj2_boost_gj2_COM.Boost(boost_gj2_COM);
      TLorentzVector temp_gj1_boost_gj2_COM=gj_m1;
      temp_gj1_boost_gj2_COM.Boost(boost_gj2_COM);
      TLorentzVector temp_gj2_boost_gj1_COM=gj_m2;
      temp_gj2_boost_gj1_COM.Boost(boost_gj1_COM);
      

      bool veto_qqqq_mass_ellipse_gen = false;
      if( (gj_m2.M()<cut_m_j2_c_veto && ((pow((gj_m2.M()-cut_m_j2_c_veto)/cut_delta_m_j2_e_veto_minus,2)+pow((gj_m1.M()-cut_m_j1_c_veto)/cut_delta_m_j1_e_veto,2))<1.)) || (gj_m2.M()>cut_m_j2_c_veto && ((pow((gj_m2.M()-cut_m_j2_c_veto)/cut_delta_m_j2_e_veto_plus,2)+pow((gj_m1.M()-cut_m_j1_c_veto)/cut_delta_m_j1_e_veto,2))<1.))){
	veto_qqqq_mass_ellipse_gen =true;
	//std::cout<<"veto is hit "<<std::endl;
      }
      if( (pow((gj_m2.M()-cut_m_j2_c)/cut_delta_m_j2_e,2)+pow((gj_m1.M()-cut_m_j1_c)/cut_delta_m_j1_e,2))<1.   && !veto_qqqq_mass_ellipse_gen && ((fabs(gj_m2.Theta()-0.5*TMath::Pi())*TMath::RadToDeg())<cut_theta_window && (fabs(gj_m1.Theta()-0.5*TMath::Pi())*TMath::RadToDeg())<cut_theta_window)){
	gen_pass_mass_cuts=true;
      }
      if(!performMassCuts){
	gen_pass_mass_cuts=true;
      }

      
      if(fill_genInfo && genjet_E->size()==2 && n_IsoLep_gen<m_cut_nLeptons && gen_pass_mass_cuts ){//effectively no cut here, but exactly two genjets
	TLorentzVector temp_sj;	
	if(usePartonInfo && tempTotEventP4.M()>sqrtS_low && tempTotEventP4.M()>sqrtS_high ){
	  if(ind_sj1_gj1!=-1){	    
	    temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj1],(*genjet_subjet_Py)[ind_sj1_gj1],(*genjet_subjet_Pz)[ind_sj1_gj1],(*genjet_subjet_E)[ind_sj1_gj1]);
	    if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
	      angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	      ind_Z_q_pos_sj_gj_angle=ind_sj1_gj1;
	    }
	    if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
	      angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	      ind_Z_q_neg_sj_gj_angle=ind_sj1_gj1;
	    }
	  }
	  if(ind_sj2_gj1!=-1){
	    temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj1],(*genjet_subjet_Py)[ind_sj2_gj1],(*genjet_subjet_Pz)[ind_sj2_gj1],(*genjet_subjet_E)[ind_sj2_gj1]);
	    if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
	      angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	      ind_Z_q_pos_sj_gj_angle=ind_sj2_gj1;
	    }
	    if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
	      angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	      ind_Z_q_neg_sj_gj_angle=ind_sj2_gj1;
	    }
	  }
	  if(ind_sj1_gj2!=-1){
	    temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj2],(*genjet_subjet_Py)[ind_sj1_gj2],(*genjet_subjet_Pz)[ind_sj1_gj2],(*genjet_subjet_E)[ind_sj1_gj2]);
	    if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
	      angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	      ind_Z_q_pos_sj_gj_angle=ind_sj1_gj2;
	    }
	    if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
	      angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	      ind_Z_q_neg_sj_gj_angle=ind_sj1_gj2;
	    }
	  }
	  if(ind_sj2_gj2!=-1){
	    temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj2],(*genjet_subjet_Py)[ind_sj2_gj2],(*genjet_subjet_Pz)[ind_sj2_gj2],(*genjet_subjet_E)[ind_sj2_gj2]);
	    if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
	      angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	      ind_Z_q_pos_sj_gj_angle=ind_sj2_gj2;
	    }
	    if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
	      angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	      ind_Z_q_neg_sj_gj_angle=ind_sj2_gj2;
	    }
	  }
	  if(ind_Z_q_pos_sj_gj_angle==ind_Z_q_neg_sj_gj_angle){
	    //decide which angle is the more propable
	    if(angle_min_sj_gj_Z_q_pos<angle_min_sj_gj_Z_q_neg){
	      if(ind_sj1_gj1!=-1 && ind_sj1_gj1!=ind_Z_q_pos_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj1],(*genjet_subjet_Py)[ind_sj1_gj1],(*genjet_subjet_Pz)[ind_sj1_gj1],(*genjet_subjet_E)[ind_sj1_gj1]);
		if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
		  angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
		  ind_Z_q_neg_sj_gj_angle=ind_sj1_gj1;
		}
	      }
	      if(ind_sj2_gj1!=-1 && ind_sj2_gj1!=ind_Z_q_pos_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj1],(*genjet_subjet_Py)[ind_sj2_gj1],(*genjet_subjet_Pz)[ind_sj2_gj1],(*genjet_subjet_E)[ind_sj2_gj1]);
		if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
		  angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
		  ind_Z_q_neg_sj_gj_angle=ind_sj2_gj1;
		}
	      }
	      if(ind_sj1_gj2!=-1 && ind_sj1_gj2!=ind_Z_q_pos_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj2],(*genjet_subjet_Py)[ind_sj1_gj2],(*genjet_subjet_Pz)[ind_sj1_gj2],(*genjet_subjet_E)[ind_sj1_gj2]);
		if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
		  angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
		  ind_Z_q_neg_sj_gj_angle=ind_sj1_gj2;
		}
	      }
	      if(ind_sj2_gj2!=-1 && ind_sj2_gj2!=ind_Z_q_pos_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj2],(*genjet_subjet_Py)[ind_sj2_gj2],(*genjet_subjet_Pz)[ind_sj2_gj2],(*genjet_subjet_E)[ind_sj2_gj2]);
		if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
		  angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
		  ind_Z_q_neg_sj_gj_angle=ind_sj2_gj2;
		}
	      }
	    }else{//closer to b
	      if(ind_sj1_gj1!=-1 && ind_sj1_gj1!=ind_Z_q_neg_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj1],(*genjet_subjet_Py)[ind_sj1_gj1],(*genjet_subjet_Pz)[ind_sj1_gj1],(*genjet_subjet_E)[ind_sj1_gj1]);
		if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
		  angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
		  ind_Z_q_pos_sj_gj_angle=ind_sj1_gj1;
		}
	      }
	      if(ind_sj2_gj1!=-1 && ind_sj2_gj1!=ind_Z_q_neg_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj1],(*genjet_subjet_Py)[ind_sj2_gj1],(*genjet_subjet_Pz)[ind_sj2_gj1],(*genjet_subjet_E)[ind_sj2_gj1]);
		if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
		  angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
		  ind_Z_q_pos_sj_gj_angle=ind_sj2_gj1;
		}
	      }
	      if(ind_sj1_gj2!=-1 && ind_sj1_gj2!=ind_Z_q_neg_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj2],(*genjet_subjet_Py)[ind_sj1_gj2],(*genjet_subjet_Pz)[ind_sj1_gj2],(*genjet_subjet_E)[ind_sj1_gj2]);
		if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
		  angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
		  ind_Z_q_pos_sj_gj_angle=ind_sj1_gj2;
		}
	      }
	      if(ind_sj2_gj2!=-1 && ind_sj2_gj2!=ind_Z_q_neg_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj2],(*genjet_subjet_Py)[ind_sj2_gj2],(*genjet_subjet_Pz)[ind_sj2_gj2],(*genjet_subjet_E)[ind_sj2_gj2]);
		if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
		  angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
		  ind_Z_q_pos_sj_gj_angle=ind_sj2_gj2;
		}
	      }
	    }
	  }
	}
	if(usePartonInfo && tempTotEventP4.M()>sqrtS_low && tempTotEventP4.M()>sqrtS_high){
	  if(ind_sj1_gj1!=-1){	    
	    temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj1],(*genjet_subjet_Py)[ind_sj1_gj1],(*genjet_subjet_Pz)[ind_sj1_gj1],(*genjet_subjet_E)[ind_sj1_gj1]);
	    if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_gj_H_bbar){
	      angle_min_sj_gj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
	      ind_H_bbar_sj_gj_angle=ind_sj1_gj1;
	    }
	    if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_gj_H_b){
	      angle_min_sj_gj_H_b=temp_sj.Angle(tempH_b.Vect());
	      ind_H_b_sj_gj_angle=ind_sj1_gj1;
	    }
	  }
	  if(ind_sj2_gj1!=-1){
	    temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj1],(*genjet_subjet_Py)[ind_sj2_gj1],(*genjet_subjet_Pz)[ind_sj2_gj1],(*genjet_subjet_E)[ind_sj2_gj1]);
	    if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_gj_H_bbar){
	      angle_min_sj_gj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
	      ind_H_bbar_sj_gj_angle=ind_sj2_gj1;
	    }
	    if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_gj_H_b){
	      angle_min_sj_gj_H_b=temp_sj.Angle(tempH_b.Vect());
	      ind_H_b_sj_gj_angle=ind_sj2_gj1;
	    }
	  }
	  if(ind_sj1_gj2!=-1){
	    temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj2],(*genjet_subjet_Py)[ind_sj1_gj2],(*genjet_subjet_Pz)[ind_sj1_gj2],(*genjet_subjet_E)[ind_sj1_gj2]);
	    if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_gj_H_bbar){
	      angle_min_sj_gj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
	      ind_H_bbar_sj_gj_angle=ind_sj1_gj2;
	    }
	    if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_gj_H_b){
	      angle_min_sj_gj_H_b=temp_sj.Angle(tempH_b.Vect());
	      ind_H_b_sj_gj_angle=ind_sj1_gj2;
	    }
	  }
	  if(ind_sj2_gj2!=-1){
	    temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj2],(*genjet_subjet_Py)[ind_sj2_gj2],(*genjet_subjet_Pz)[ind_sj2_gj2],(*genjet_subjet_E)[ind_sj2_gj2]);
	    if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_gj_H_bbar){
	      angle_min_sj_gj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
	      ind_H_bbar_sj_gj_angle=ind_sj2_gj2;
	    }
	    if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_gj_H_b){
	      angle_min_sj_gj_H_b=temp_sj.Angle(tempH_b.Vect());
	      ind_H_b_sj_gj_angle=ind_sj2_gj2;
	    }
	  }
	  if(ind_H_bbar_sj_gj_angle==ind_H_b_sj_gj_angle){
	    //decide which angle is the more propable
	    if(angle_min_sj_gj_H_bbar<angle_min_sj_gj_H_b){
	      if(ind_sj1_gj1!=-1 && ind_sj1_gj1!=ind_H_bbar_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj1],(*genjet_subjet_Py)[ind_sj1_gj1],(*genjet_subjet_Pz)[ind_sj1_gj1],(*genjet_subjet_E)[ind_sj1_gj1]);
		if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_gj_H_b){
		  angle_min_sj_gj_H_b=temp_sj.Angle(tempH_b.Vect());
		  ind_H_b_sj_gj_angle=ind_sj1_gj1;
		}
	      }
	      if(ind_sj2_gj1!=-1 && ind_sj2_gj1!=ind_H_bbar_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj1],(*genjet_subjet_Py)[ind_sj2_gj1],(*genjet_subjet_Pz)[ind_sj2_gj1],(*genjet_subjet_E)[ind_sj2_gj1]);
		if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_gj_H_b){
		  angle_min_sj_gj_H_b=temp_sj.Angle(tempH_b.Vect());
		  ind_H_b_sj_gj_angle=ind_sj2_gj1;
		}
	      }
	      if(ind_sj1_gj2!=-1 && ind_sj1_gj2!=ind_H_bbar_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj2],(*genjet_subjet_Py)[ind_sj1_gj2],(*genjet_subjet_Pz)[ind_sj1_gj2],(*genjet_subjet_E)[ind_sj1_gj2]);
		if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_gj_H_b){
		  angle_min_sj_gj_H_b=temp_sj.Angle(tempH_b.Vect());
		  ind_H_b_sj_gj_angle=ind_sj1_gj2;
		}
	      }
	      if(ind_sj2_gj2!=-1 && ind_sj2_gj2!=ind_H_bbar_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj2],(*genjet_subjet_Py)[ind_sj2_gj2],(*genjet_subjet_Pz)[ind_sj2_gj2],(*genjet_subjet_E)[ind_sj2_gj2]);
		if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_gj_H_b){
		  angle_min_sj_gj_H_b=temp_sj.Angle(tempH_b.Vect());
		  ind_H_b_sj_gj_angle=ind_sj2_gj2;
		}
	      }
	    }else{//closer to b
	      if(ind_sj1_gj1!=-1 && ind_sj1_gj1!=ind_H_b_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj1],(*genjet_subjet_Py)[ind_sj1_gj1],(*genjet_subjet_Pz)[ind_sj1_gj1],(*genjet_subjet_E)[ind_sj1_gj1]);
		if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_gj_H_bbar){
		  angle_min_sj_gj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
		  ind_H_bbar_sj_gj_angle=ind_sj1_gj1;
		}
	      }
	      if(ind_sj2_gj1!=-1 && ind_sj2_gj1!=ind_H_b_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj1],(*genjet_subjet_Py)[ind_sj2_gj1],(*genjet_subjet_Pz)[ind_sj2_gj1],(*genjet_subjet_E)[ind_sj2_gj1]);
		if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_gj_H_bbar){
		  angle_min_sj_gj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
		  ind_H_bbar_sj_gj_angle=ind_sj2_gj1;
		}
	      }
	      if(ind_sj1_gj2!=-1 && ind_sj1_gj2!=ind_H_b_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj2],(*genjet_subjet_Py)[ind_sj1_gj2],(*genjet_subjet_Pz)[ind_sj1_gj2],(*genjet_subjet_E)[ind_sj1_gj2]);
		if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_gj_H_bbar){
		  angle_min_sj_gj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
		  ind_H_bbar_sj_gj_angle=ind_sj1_gj2;
		}
	      }
	      if(ind_sj2_gj2!=-1 && ind_sj2_gj2!=ind_H_b_sj_gj_angle){	    
		temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj2],(*genjet_subjet_Py)[ind_sj2_gj2],(*genjet_subjet_Pz)[ind_sj2_gj2],(*genjet_subjet_E)[ind_sj2_gj2]);
		if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_gj_H_bbar){
		  angle_min_sj_gj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
		  ind_H_bbar_sj_gj_angle=ind_sj2_gj2;
		}
	      }
	    }
	  }
	  if(ind_Z_q_pos_sj_gj_angle!=-1/* && (angle_min_sj_gj_Z_q_pos*TMath::RadToDeg())<15.0*/){
	    h_hist_vec_gen[368]->Fill((*genjet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_gj_angle],weight);
	    h_hist_vec_gen[369]->Fill((*genjet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_pos_sj_gj_angle],weight);
	    h_hist_vec_gen[374]->Fill((*genjet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_pos_sj_gj_angle],weight);
	    h_hist_vec_gen[375]->Fill((*genjet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_pos_sj_gj_angle],weight);
	    h_hist_vec_gen[376]->Fill((*genjet_subjet_jetChargePt_kappa_0_25)[ind_Z_q_pos_sj_gj_angle],weight);
	    h_hist_vec_gen[377]->Fill((*genjet_subjet_jetChargePt_kappa_0_50)[ind_Z_q_pos_sj_gj_angle],weight);
	    h_hist_vec_gen[382]->Fill((*genjet_subjet_jetChargePt_kappa_0_20)[ind_Z_q_pos_sj_gj_angle],weight);
	    h_hist_vec_gen[383]->Fill((*genjet_subjet_jetChargePt_kappa_0_30)[ind_Z_q_pos_sj_gj_angle],weight);
	  }
	  if(ind_Z_q_neg_sj_gj_angle!=-1/* && (angle_min_sj_gj_Z_q_neg*TMath::RadToDeg())<15.0*/){
	    h_hist_vec_gen[384]->Fill((*genjet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_gj_angle],weight);
	    h_hist_vec_gen[385]->Fill((*genjet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_neg_sj_gj_angle],weight);
	    h_hist_vec_gen[390]->Fill((*genjet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_neg_sj_gj_angle],weight);
	    h_hist_vec_gen[391]->Fill((*genjet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_neg_sj_gj_angle],weight);
	    h_hist_vec_gen[392]->Fill((*genjet_subjet_jetChargePt_kappa_0_25)[ind_Z_q_neg_sj_gj_angle],weight);
	    h_hist_vec_gen[393]->Fill((*genjet_subjet_jetChargePt_kappa_0_50)[ind_Z_q_neg_sj_gj_angle],weight);
	    h_hist_vec_gen[398]->Fill((*genjet_subjet_jetChargePt_kappa_0_20)[ind_Z_q_neg_sj_gj_angle],weight);
	    h_hist_vec_gen[399]->Fill((*genjet_subjet_jetChargePt_kappa_0_30)[ind_Z_q_neg_sj_gj_angle],weight);
	  }
	}//finish cut on MC truth info

	if(sqrtS_eff_gen<sqrtS_low){
	  h_hist_vec_gen[0]->Fill(gj_m1.Angle(gj_m2.Vect())*TMath::RadToDeg(),weight);
	  h_hist_vec_gen[3]->Fill(DeltaPhi(gj_m1.Phi(),gj_m2.Phi())*TMath::RadToDeg(),weight);
	  h_hist_vec_gen[6]->Fill(fabs(gj_m1.Theta()-gj_m2.Theta())*TMath::RadToDeg(),weight);
	  h_hist_vec_gen[9]->Fill(gj_m1.Theta()*TMath::RadToDeg(),weight);
	  h_hist_vec_gen[12]->Fill(gj_m2.Theta()*TMath::RadToDeg(),weight);
	  h_hist_vec_gen[15]->Fill(gj_m1.M(),weight);
	  h_hist_vec_gen[18]->Fill(gj_m2.M(),weight);

	 if(ind_sj1_gj1!=-1){
	   h_hist_vec_gen[237]->Fill((*genjet_subjet_E)[ind_sj1_gj1],weight);
	   h_hist_vec_gen[249]->Fill((*genjet_subjet_E)[ind_sj1_gj1]/(*genjet_E)[ind_gj1],weight);
	   if(ind_sj2_gj1!=-1){
	     h_hist_vec_gen[255]->Fill(temp_sj1_gj1.Angle(temp_sj2_gj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj1!=-1){
	   h_hist_vec_gen[240]->Fill((*genjet_subjet_E)[ind_sj2_gj1],weight);
	 }
	 if(ind_sj1_gj2!=-1){
	   h_hist_vec_gen[243]->Fill((*genjet_subjet_E)[ind_sj1_gj2],weight);
	   h_hist_vec_gen[252]->Fill((*genjet_subjet_E)[ind_sj1_gj2]/(*genjet_E)[ind_gj2],weight);
	   if(ind_sj2_gj2!=-1){
	     h_hist_vec_gen[258]->Fill(temp_sj1_gj2.Angle(temp_sj2_gj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj2!=-1){
	   h_hist_vec_gen[246]->Fill((*genjet_subjet_E)[ind_sj2_gj2],weight);
	 }
	 h_hist_vec_gen[266]->Fill(gj1_MET.M(),weight);
	 h_hist_vec_gen[269]->Fill(gj2_MET.M(),weight);
	 h_hist_vec_gen[272]->Fill(gj1_EMiss.M(),weight);
	 h_hist_vec_gen[275]->Fill(gj2_EMiss.M(),weight);
	 h_hist_vec_gen[293]->Fill(gj1_MHMiss.M(),weight);
	 h_hist_vec_gen[296]->Fill(gj2_MHMiss.M(),weight);
	 //while the boost does change with energy scaling based on MET, the original direction of the vector does not
	 h_hist_vec_gen[314]->Fill(cos(min(temp_sj1_gj1_boost_gj1_COM.Angle(gj_m1.Vect()),temp_sj2_gj1_boost_gj1_COM.Angle(gj_m1.Vect()))),weight);
	 h_hist_vec_gen[317]->Fill(cos(max(temp_sj1_gj1_boost_gj1_COM.Angle(gj_m1.Vect()),temp_sj2_gj1_boost_gj1_COM.Angle(gj_m1.Vect()))),weight);
	 h_hist_vec_gen[320]->Fill(cos(min(temp_sj1_gj2_boost_gj2_COM.Angle(gj_m2.Vect()),temp_sj2_gj2_boost_gj2_COM.Angle(gj_m2.Vect()))),weight);
	 h_hist_vec_gen[323]->Fill(cos(max(temp_sj1_gj2_boost_gj2_COM.Angle(gj_m2.Vect()),temp_sj2_gj2_boost_gj2_COM.Angle(gj_m2.Vect()))),weight);
	 h_hist_vec_gen[330]->Fill(gj_m1_orig.M()-gj_m2_orig.M(),weight);
	 h_hist_vec_gen[333]->Fill(gj1_EMiss.M()-gj2_EMiss.M(),weight);
      }else if(sqrtS_eff_gen<sqrtS_high_reco){
	h_hist_vec_gen[1]->Fill(gj_m1.Angle(gj_m2.Vect())*TMath::RadToDeg(),weight);
	h_hist_vec_gen[4]->Fill(DeltaPhi(gj_m1.Phi(),gj_m2.Phi())*TMath::RadToDeg(),weight);
	h_hist_vec_gen[7]->Fill(fabs(gj_m1.Theta()-gj_m2.Theta())*TMath::RadToDeg(),weight);
	h_hist_vec_gen[10]->Fill(gj_m1.Theta()*TMath::RadToDeg(),weight);
	h_hist_vec_gen[13]->Fill(gj_m2.Theta()*TMath::RadToDeg(),weight);
	h_hist_vec_gen[16]->Fill(gj_m1.M(),weight);
	h_hist_vec_gen[19]->Fill(gj_m2.M(),weight);

	 if(ind_sj1_gj1!=-1){
	   h_hist_vec_gen[238]->Fill((*genjet_subjet_E)[ind_sj1_gj1],weight);
	   h_hist_vec_gen[250]->Fill((*genjet_subjet_E)[ind_sj1_gj1]/(*genjet_E)[ind_gj1],weight);
	   if(ind_sj2_gj1!=-1){
	     h_hist_vec_gen[256]->Fill(temp_sj1_gj1.Angle(temp_sj2_gj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj1!=-1){
	   h_hist_vec_gen[241]->Fill((*genjet_subjet_E)[ind_sj2_gj1],weight);
	 }
	 if(ind_sj1_gj2!=-1){
	   h_hist_vec_gen[244]->Fill((*genjet_subjet_E)[ind_sj1_gj2],weight);
	   h_hist_vec_gen[253]->Fill((*genjet_subjet_E)[ind_sj1_gj2]/(*genjet_E)[ind_gj2],weight);
	   if(ind_sj2_gj2!=-1){
	     h_hist_vec_gen[259]->Fill(temp_sj1_gj2.Angle(temp_sj2_gj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj2!=-1){
	   h_hist_vec_gen[247]->Fill((*genjet_subjet_E)[ind_sj2_gj2],weight);
	 }
	 h_hist_vec_gen[267]->Fill(gj1_MET.M(),weight);
	 h_hist_vec_gen[270]->Fill(gj2_MET.M(),weight);
	 h_hist_vec_gen[273]->Fill(gj1_EMiss.M(),weight);
	 h_hist_vec_gen[276]->Fill(gj2_EMiss.M(),weight);
	 h_hist_vec_gen[294]->Fill(gj1_MHMiss.M(),weight);
	 h_hist_vec_gen[297]->Fill(gj2_MHMiss.M(),weight);
	 //while the boost does change with energy scaling based on MET, the original direction of the vector does not
	 h_hist_vec_gen[315]->Fill(cos(min(temp_sj1_gj1_boost_gj1_COM.Angle(gj_m1.Vect()),temp_sj2_gj1_boost_gj1_COM.Angle(gj_m1.Vect()))),weight);
	 h_hist_vec_gen[318]->Fill(cos(max(temp_sj1_gj1_boost_gj1_COM.Angle(gj_m1.Vect()),temp_sj2_gj1_boost_gj1_COM.Angle(gj_m1.Vect()))),weight);
	 h_hist_vec_gen[321]->Fill(cos(min(temp_sj1_gj2_boost_gj2_COM.Angle(gj_m2.Vect()),temp_sj2_gj2_boost_gj2_COM.Angle(gj_m2.Vect()))),weight);
	 h_hist_vec_gen[324]->Fill(cos(max(temp_sj1_gj2_boost_gj2_COM.Angle(gj_m2.Vect()),temp_sj2_gj2_boost_gj2_COM.Angle(gj_m2.Vect()))),weight);
	 h_hist_vec_gen[331]->Fill(gj_m1_orig.M()-gj_m2_orig.M(),weight);
	 h_hist_vec_gen[334]->Fill(gj1_EMiss.M()-gj2_EMiss.M(),weight);
	 
       }else{
	 h_hist_vec_gen[2]->Fill(gj_m1.Angle(gj_m2.Vect())*TMath::RadToDeg(),weight);    
	 h_hist_vec_gen[5]->Fill(DeltaPhi(gj_m1.Phi(),gj_m2.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[8]->Fill(fabs(gj_m1.Theta()-gj_m2.Theta())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[11]->Fill(gj_m1.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[14]->Fill(gj_m2.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[17]->Fill(gj_m1.M(),weight);
	 h_hist_vec_gen[20]->Fill(gj_m2.M(),weight);
	 if(ind_sj1_gj1!=-1){
	   if(usePartonInfo){
	     temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj1],(*genjet_subjet_Py)[ind_sj1_gj1],(*genjet_subjet_Pz)[ind_sj1_gj1],(*genjet_subjet_E)[ind_sj1_gj1]);
	     if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
	       angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	       ind_Z_q_pos_sj_gj_angle=ind_sj1_gj1;
	     }
	     if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
	       angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	       ind_Z_q_neg_sj_gj_angle=ind_sj1_gj1;
	     }
	   }
	   h_hist_vec_gen[239]->Fill((*genjet_subjet_E)[ind_sj1_gj1],weight);
	   h_hist_vec_gen[251]->Fill((*genjet_subjet_E)[ind_sj1_gj1]/(*genjet_E)[ind_gj1],weight);
	   if(ind_sj2_gj1!=-1){
	     h_hist_vec_gen[257]->Fill(temp_sj1_gj1.Angle(temp_sj2_gj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj1!=-1){
	   if(usePartonInfo){
	     temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj1],(*genjet_subjet_Py)[ind_sj2_gj1],(*genjet_subjet_Pz)[ind_sj2_gj1],(*genjet_subjet_E)[ind_sj2_gj1]);
	     if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
	       angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	       ind_Z_q_pos_sj_gj_angle=ind_sj2_gj1;
	     }
	     if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
	       angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	       ind_Z_q_neg_sj_gj_angle=ind_sj2_gj1;
	     }
	   }
	   h_hist_vec_gen[242]->Fill((*genjet_subjet_E)[ind_sj2_gj1],weight);
	 }
	 if(ind_sj1_gj2!=-1){
	   if(usePartonInfo){
	     temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj1_gj2],(*genjet_subjet_Py)[ind_sj1_gj2],(*genjet_subjet_Pz)[ind_sj1_gj2],(*genjet_subjet_E)[ind_sj1_gj2]);
	     if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
	       angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	       ind_Z_q_pos_sj_gj_angle=ind_sj1_gj2;
	     }
	     if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
	       angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	       ind_Z_q_neg_sj_gj_angle=ind_sj1_gj2;
	     }
	   }
	   h_hist_vec_gen[245]->Fill((*genjet_subjet_E)[ind_sj1_gj2],weight);
	   h_hist_vec_gen[254]->Fill((*genjet_subjet_E)[ind_sj1_gj2]/(*genjet_E)[ind_gj2],weight);
	   if(ind_sj2_gj2!=-1){
	     h_hist_vec_gen[260]->Fill(temp_sj1_gj2.Angle(temp_sj2_gj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_gj2!=-1){
	   if(usePartonInfo){
	     temp_sj.SetPxPyPzE((*genjet_subjet_Px)[ind_sj2_gj2],(*genjet_subjet_Py)[ind_sj2_gj2],(*genjet_subjet_Pz)[ind_sj2_gj2],(*genjet_subjet_E)[ind_sj2_gj2]);
	     if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
	       angle_min_sj_gj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	       ind_Z_q_pos_sj_gj_angle=ind_sj2_gj2;
	     }
	     if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_gj_Z_q_neg){
	       angle_min_sj_gj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	       ind_Z_q_neg_sj_gj_angle=ind_sj2_gj2;
	    }
	   }
	   h_hist_vec_gen[248]->Fill((*genjet_subjet_E)[ind_sj2_gj2],weight);
	 }
	 if(usePartonInfo && i_entry%500==0){
	   std::cout<<"min angle/index "<<angle_min_sj_gj_Z_q_pos*TMath::RadToDeg()<<"/"<<ind_Z_q_pos_sj_gj_angle<<std::endl;
	 }
	 if(usePartonInfo && fill_partonInfo && (ind_Z_q_pos_sj_gj_angle!=-1 && (angle_min_sj_gj_Z_q_pos*TMath::RadToDeg())<15.0) && ( ind_Z_q_neg_sj_gj_angle!=-1 && (angle_min_sj_gj_Z_q_neg*TMath::RadToDeg())<15.0)){
	   h_hist_vec_2D[42]->Fill((*genjet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_pos_sj_gj_angle],(*genjet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_neg_sj_gj_angle],weight); 
	   h_hist_vec_2D[43]->Fill((*genjet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_gj_angle],(*genjet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_gj_angle],weight); 
	   h_hist_vec_2D[44]->Fill((*genjet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_pos_sj_gj_angle],(*genjet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_neg_sj_gj_angle],weight); 
	   if(((*genjet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_pos_sj_gj_angle])!=0 && ((*genjet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_neg_sj_gj_angle])!=0){
	     h_hist_vec_gen[401]->Fill((*genjet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_pos_sj_gj_angle]-(*genjet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_neg_sj_gj_angle],weight);
	   }
	   if(((*genjet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_gj_angle])!=0 && ((*genjet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_gj_angle])!=0){
	     h_hist_vec_gen[402]->Fill((*genjet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_gj_angle]-(*genjet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_gj_angle],weight);
	   }
	   if(((*genjet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_pos_sj_gj_angle])!=0 && ((*genjet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_neg_sj_gj_angle])!=0){
	     h_hist_vec_gen[403]->Fill((*genjet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_pos_sj_gj_angle]-(*genjet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_neg_sj_gj_angle],weight);
	   }
	 }
	 if(usePartonInfo && fill_partonInfo && (ind_H_bbar_sj_gj_angle!=-1 && (angle_min_sj_gj_H_bbar*TMath::RadToDeg())<15.0) && ( ind_H_b_sj_gj_angle!=-1 && (angle_min_sj_gj_H_b*TMath::RadToDeg())<15.0)){
	   h_hist_vec_2D[50]->Fill((*genjet_subjet_jetChargeE_kappa_0_20)[ind_H_bbar_sj_gj_angle],(*genjet_subjet_jetChargeE_kappa_0_20)[ind_H_b_sj_gj_angle],weight); 
	   h_hist_vec_2D[51]->Fill((*genjet_subjet_jetChargeE_kappa_0_25)[ind_H_bbar_sj_gj_angle],(*genjet_subjet_jetChargeE_kappa_0_25)[ind_H_b_sj_gj_angle],weight); 
	   h_hist_vec_2D[52]->Fill((*genjet_subjet_jetChargeE_kappa_0_30)[ind_H_bbar_sj_gj_angle],(*genjet_subjet_jetChargeE_kappa_0_30)[ind_H_b_sj_gj_angle],weight); 
	   if(((*genjet_subjet_jetChargeE_kappa_0_20)[ind_H_bbar_sj_gj_angle])!=0 && ((*genjet_subjet_jetChargeE_kappa_0_20)[ind_H_b_sj_gj_angle])!=0){
	     h_hist_vec_gen[405]->Fill((*genjet_subjet_jetChargeE_kappa_0_20)[ind_H_bbar_sj_gj_angle]-(*genjet_subjet_jetChargeE_kappa_0_20)[ind_H_b_sj_gj_angle],weight);
	   }
	   if(((*genjet_subjet_jetChargeE_kappa_0_25)[ind_H_bbar_sj_gj_angle])!=0 && ((*genjet_subjet_jetChargeE_kappa_0_25)[ind_H_b_sj_gj_angle])!=0){
	     h_hist_vec_gen[406]->Fill((*genjet_subjet_jetChargeE_kappa_0_25)[ind_H_bbar_sj_gj_angle]-(*genjet_subjet_jetChargeE_kappa_0_25)[ind_H_b_sj_gj_angle],weight);
	   }
	   if(((*genjet_subjet_jetChargeE_kappa_0_30)[ind_H_bbar_sj_gj_angle])!=0 && ((*genjet_subjet_jetChargeE_kappa_0_30)[ind_H_b_sj_gj_angle])!=0){
	     h_hist_vec_gen[407]->Fill((*genjet_subjet_jetChargeE_kappa_0_30)[ind_H_bbar_sj_gj_angle]-(*genjet_subjet_jetChargeE_kappa_0_30)[ind_H_b_sj_gj_angle],weight);
	   }
	 }
	 if(usePartonInfo && ind_Z_q_pos_sj_gj_angle!=-1 && (angle_min_sj_gj_Z_q_pos*TMath::RadToDeg())<15.0){
	   h_hist_vec_gen[336]->Fill((*genjet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_gj_angle],weight);
	   h_hist_vec_gen[337]->Fill((*genjet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_pos_sj_gj_angle],weight);
	   h_hist_vec_gen[342]->Fill((*genjet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_pos_sj_gj_angle],weight);
	   h_hist_vec_gen[343]->Fill((*genjet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_pos_sj_gj_angle],weight);
	   h_hist_vec_gen[344]->Fill((*genjet_subjet_jetChargePt_kappa_0_25)[ind_Z_q_pos_sj_gj_angle],weight);
	   h_hist_vec_gen[345]->Fill((*genjet_subjet_jetChargePt_kappa_0_50)[ind_Z_q_pos_sj_gj_angle],weight);
	   h_hist_vec_gen[350]->Fill((*genjet_subjet_jetChargePt_kappa_0_20)[ind_Z_q_pos_sj_gj_angle],weight);
	   h_hist_vec_gen[351]->Fill((*genjet_subjet_jetChargePt_kappa_0_30)[ind_Z_q_pos_sj_gj_angle],weight);
	 }
	 if(usePartonInfo && ind_Z_q_neg_sj_gj_angle!=-1 && (angle_min_sj_gj_Z_q_neg*TMath::RadToDeg())<15.0){
	   h_hist_vec_gen[352]->Fill((*genjet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_gj_angle],weight);
	   h_hist_vec_gen[353]->Fill((*genjet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_neg_sj_gj_angle],weight);
	   h_hist_vec_gen[358]->Fill((*genjet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_neg_sj_gj_angle],weight);
	   h_hist_vec_gen[359]->Fill((*genjet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_neg_sj_gj_angle],weight);
	   h_hist_vec_gen[360]->Fill((*genjet_subjet_jetChargePt_kappa_0_25)[ind_Z_q_neg_sj_gj_angle],weight);
	   h_hist_vec_gen[361]->Fill((*genjet_subjet_jetChargePt_kappa_0_50)[ind_Z_q_neg_sj_gj_angle],weight);
	   h_hist_vec_gen[366]->Fill((*genjet_subjet_jetChargePt_kappa_0_20)[ind_Z_q_neg_sj_gj_angle],weight);
	   h_hist_vec_gen[367]->Fill((*genjet_subjet_jetChargePt_kappa_0_30)[ind_Z_q_neg_sj_gj_angle],weight);
	 }
	 h_hist_vec_gen[268]->Fill(gj1_MET.M(),weight);
	 h_hist_vec_gen[271]->Fill(gj2_MET.M(),weight);
	 h_hist_vec_gen[274]->Fill(gj1_EMiss.M(),weight);
	 h_hist_vec_gen[277]->Fill(gj2_EMiss.M(),weight);
	 h_hist_vec_gen[295]->Fill(gj1_MHMiss.M(),weight);
	 h_hist_vec_gen[298]->Fill(gj2_MHMiss.M(),weight);
	 h_hist_vec_gen[281]->Fill(gj_m1_orig.E()+gj_m2_orig.E(),weight);
	 h_hist_vec_gen[282]->Fill(gj1_EMiss.E()+gj2_EMiss.E(),weight);
	 h_hist_vec_gen[283]->Fill((tempTotGenP4-tempGenIsoPhP4).E(),weight);
	 h_hist_vec_gen[284]->Fill((tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).E(),weight);
	 h_hist_vec_gen[285]->Fill(gj_m1_orig.E(),weight);
	 h_hist_vec_gen[286]->Fill(gj_m2_orig.E(),weight);
	 h_hist_vec_gen[287]->Fill(gj1_EMiss.E(),weight);
	 h_hist_vec_gen[288]->Fill(gj2_EMiss.E(),weight);
	 h_hist_vec_gen[289]->Fill(gj_m1_orig.E()-gj_m2_orig.E(),weight);
	 h_hist_vec_gen[290]->Fill(gj1_EMiss.E()-gj2_EMiss.E(),weight);
	 h_hist_vec_gen[291]->Fill((gj_m1_orig.E()-gj_m2_orig.E())/(gj_m1_orig.E()+gj_m2_orig.E()),weight);
	 h_hist_vec_gen[292]->Fill((gj1_EMiss.E()-gj2_EMiss.E())/(gj1_EMiss.E()+gj2_EMiss.E()),weight);

	 h_hist_vec_gen[299]->Fill(gj1_MHMiss.E()+gj2_MHMiss.E(),weight);
	 h_hist_vec_gen[300]->Fill((tempTotGenP4-tempGenIsoPhP4+tempGenMHMissCorrP4).E(),weight);
	 h_hist_vec_gen[301]->Fill(gj1_MHMiss.E(),weight);
	 h_hist_vec_gen[302]->Fill(gj2_MHMiss.E(),weight);
	 h_hist_vec_gen[303]->Fill((gj1_MHMiss.E()-gj2_MHMiss.E())/(gj1_MHMiss.E()+gj2_MHMiss.E()),weight);
	 h_hist_vec_gen[304]->Fill((tempGenEMissCorrP4.Pt()-tempTotInvGenP4.Pt())/tempTotInvGenP4.Pt(),weight);
	 h_hist_vec_gen[305]->Fill((tempGenMHMissCorrP4.Pt()-tempTotInvGenP4.Pt())/tempTotInvGenP4.Pt(),weight);
	 h_hist_vec_gen[306]->Fill(((-gj_m1_orig-gj_m2_orig).Pt()-tempTotInvGenP4.Pt())/tempTotInvGenP4.Pt(),weight);
	 h_hist_vec_gen[307]->Fill((tempGenEMissCorrP4.E()-tempTotInvGenP4.E())/tempTotInvGenP4.E(),weight);
	 h_hist_vec_gen[308]->Fill((tempGenMHMissCorrP4.E()-tempTotInvGenP4.E())/tempTotInvGenP4.E(),weight);
	 h_hist_vec_gen[309]->Fill(DeltaPhi(tempGenEMissCorrP4.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[310]->Fill(DeltaPhi(tempGenMHMissCorrP4.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[311]->Fill(DeltaPhi((-gj_m1_orig-gj_m2_orig).Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[312]->Fill(tempGenEMissCorrP4.Angle(tempTotInvGenP4.Vect())*TMath::RadToDeg(),weight);
	 h_hist_vec_gen[313]->Fill(tempGenMHMissCorrP4.Angle(tempTotInvGenP4.Vect())*TMath::RadToDeg(),weight);
	 //while the boost does change with energy scaling based on MET, the original direction of the vector does not
	 h_hist_vec_gen[316]->Fill(cos(min(temp_sj1_gj1_boost_gj1_COM.Angle(gj_m1.Vect()),temp_sj2_gj1_boost_gj1_COM.Angle(gj_m1.Vect()))),weight);
	 h_hist_vec_gen[319]->Fill(cos(max(temp_sj1_gj1_boost_gj1_COM.Angle(gj_m1.Vect()),temp_sj2_gj1_boost_gj1_COM.Angle(gj_m1.Vect()))),weight);
	 h_hist_vec_gen[322]->Fill(cos(min(temp_sj1_gj2_boost_gj2_COM.Angle(gj_m2.Vect()),temp_sj2_gj2_boost_gj2_COM.Angle(gj_m2.Vect()))),weight);
	 h_hist_vec_gen[325]->Fill(cos(max(temp_sj1_gj2_boost_gj2_COM.Angle(gj_m2.Vect()),temp_sj2_gj2_boost_gj2_COM.Angle(gj_m2.Vect()))),weight);

	 h_hist_vec_gen[326]->Fill(tempTotInvGenP4.Pt()/(tempTotGenP4-tempGenIsoPhP4).E(),weight);
	 h_hist_vec_gen[327]->Fill((-gj_m1_orig-gj_m2_orig).Pt()/(tempTotGenP4-tempGenIsoPhP4).E(),weight);
	 h_hist_vec_gen[328]->Fill(tempTotInvGenP4.Pt()/(gj_m1_orig+gj_m2_orig).E(),weight);
	 h_hist_vec_gen[329]->Fill((-gj_m1_orig-gj_m2_orig).Pt()/(gj_m1_orig+gj_m2_orig).E(),weight);

	 h_hist_vec_gen[332]->Fill(gj_m1_orig.M()-gj_m2_orig.M(),weight);
	 h_hist_vec_gen[335]->Fill(gj1_EMiss.M()-gj2_EMiss.M(),weight);
      }//selection on highest sqrtS closed
     }//bracket closed on two genjets plus isolated lepton cut+mass cut+ fill gen info
    //hilst SqrtS for all events including two genjets
     if(fill_genInfo){
       h_hist_vec_gen[261]->Fill((tempTotGenP4-tempGenIsoPhP4).M(),weight);
       h_hist_vec_gen[262]->Fill((tempTotGenP4-tempGenIsoPhP4+tempTotInvGenP4).M(),weight);
       h_hist_vec_gen[263]->Fill(tempTotGenP4.M(),weight);
       h_hist_vec_gen[264]->Fill((tempTotGenP4+tempTotInvGenP4).M(),weight);
       h_hist_vec_gen[265]->Fill((tempTotGenP4-tempGenIsoPhP4+tempTotInvGenP4).M(),weight);
       h_hist_vec_gen[278]->Fill((tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).M(),weight);
       h_hist_vec_gen[279]->Fill((gj_m1_orig+gj_m2_orig).M(),weight);
       h_hist_vec_gen[280]->Fill((gj1_EMiss+gj2_EMiss).M(),weight);
       if(fill_partonInfo && tempTotEventP4.M()>sqrtS_high_reco){
	 h_hist_parton[115]->Fill(((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
	 h_hist_parton[116]->Fill(((gj1_EMiss+gj2_EMiss-tempGenIsoPhP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
	 h_hist_parton[119]->Fill(((tempTotGenP4-tempGenIsoPhP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
	 h_hist_parton[120]->Fill(((tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
	 h_hist_parton[122]->Fill(((gj1_MHMiss+gj2_MHMiss-tempGenIsoPhP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
	 h_hist_parton[124]->Fill(((tempTotGenP4-tempGenIsoPhP4+tempGenMHMissCorrP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
       }
     }
    }
    //two genjet loop AGAIN 
    
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
      TLorentzVector tempRecoMHMissP4=-rj_m1-rj_m2;
      TLorentzVector tempRecoMHTP4(0,0,0,0);
      tempRecoMHTP4.SetPxPyPzE(tempRecoMHMissP4.Px(),tempRecoMHMissP4.Py(),0.,tempRecoMHMissP4.Pt());

      //needed in case we switch default mass jets to EMiss rescaled once, then subjet becomes a "raw" subjet energy, which has to be compared
      //to original jet energy properties etc
      TLorentzVector rj_m1_orig=rj_m1;
      TLorentzVector rj_m2_orig=rj_m2;
      unsigned int ind_rj1_orig=ind_rj1;
      unsigned int ind_rj2_orig=ind_rj2;
      if((rj_m1_orig.E()+rj_m2_orig.E())>2000){
	counter_j1j2_reco_2000+=1;
      }else{
	counter_j1j2_reco_low_2000+=1;
      }


      TLorentzVector rj1_METProjVecProp(0,0,0,0);
      TLorentzVector rj1_EMissProjVecProp(0,0,0,0);
      TLorentzVector rj1_MHMissProjVecProp(0,0,0,0);
      TLorentzVector rj2_METProjVecProp(0,0,0,0);
      TLorentzVector rj2_EMissProjVecProp(0,0,0,0);
      TLorentzVector rj2_MHMissProjVecProp(0,0,0,0);
      TVector3 recomet_proj_sum(0,0,0);
      if(tempRecoMETP4.Pt()>0){
	//TLorentzVector tempRecoMETP4_rj1_Dir(0,0,0,0);
	double rj1_METProj=tempRecoMETP4.Vect().Dot(rj_m1_orig.Vect().Unit());
	TLorentzVector rj1_METProjVec(0,0,0,0);
	rj1_METProjVec.SetPxPyPzE(rj1_METProj*rj_m1_orig.Px()/rj_m1_orig.P(),rj1_METProj*rj_m1_orig.Py()/rj_m1_orig.P(),0,rj1_METProj);
	//std::cout<<"proj r1 MET "<<rj1_METProj<<" x/y/z/E/pt "<<rj1_METProjVec.Px()<<"/"<<rj1_METProjVec.Py()<<"/"<<rj1_METProjVec.Pz()<<"/"<<rj1_METProjVec.E()<<"/"<<rj1_METProjVec.Pt()<<std::endl;
	//tempTotInvRecoP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);
	TVector3 rj1_METProjProp=(tempRecoMETP4.Vect().Dot(rj_m1_orig.Vect().Unit()))*rj_m1_orig.Vect().Unit();
	if(tempRecoMETP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit())>0){//rj1 and MET in same hemisphere
	  rj1_METProjVecProp.SetPxPyPzE(rj1_METProjProp.Px(),rj1_METProjProp.Py(),0,rj1_METProjProp.Pt());
	  rj1_EMissProjVecProp.SetPxPyPzE(rj1_METProjProp.Px(),rj1_METProjProp.Py(),rj1_METProjProp.Pt()*rj_m1_orig.Pz()/rj_m1_orig.Pt(),rj1_METProjProp.Pt()*rj_m1_orig.P()/rj_m1_orig.Pt());
	  recomet_proj_sum+=rj1_METProjProp;
	  //std::cout<<"right hemisphere rj1 MET "<<(tempRecoMETP4.Vect().Dot(rj_m1_orig.Vect().Unit()))<<"/"<<i_entry<<" "<<acos(tempRecoMETP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit()))*TMath::RadToDeg()<<std::endl;
	}else{
	  if(tempRecoMETP4.Vect().Dot(rj_m1_orig.Vect().Unit())>0){
	    std::cout<<"wrong hemisphere rj1 MET"<<(tempRecoMETP4.Vect().Dot(rj_m1_orig.Vect().Unit()))<<"/"<<i_entry<<" "<<acos(tempRecoMETP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit()))*TMath::RadToDeg()<<std::endl;
	  }
	}
	double rj2_METProj=tempRecoMETP4.Vect().Dot(rj_m2_orig.Vect().Unit());
	TLorentzVector rj2_METProjVec(0,0,0,0);
	rj2_METProjVec.SetPxPyPzE(rj2_METProj*rj_m2_orig.Px()/rj_m2_orig.P(),rj2_METProj*rj_m2_orig.Py()/rj_m2_orig.P(),0,rj2_METProj);
	TVector3 rj2_METProjProp=(tempRecoMETP4.Vect().Dot(rj_m2_orig.Vect().Unit()))*rj_m2_orig.Vect().Unit();
	if(tempRecoMETP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit())>0){//rj2 and MET in same hemisphere
	  rj2_METProjVecProp.SetPxPyPzE(rj2_METProjProp.Px(),rj2_METProjProp.Py(),0,rj2_METProjProp.Pt());
	  rj2_EMissProjVecProp.SetPxPyPzE(rj2_METProjProp.Px(),rj2_METProjProp.Py(),rj2_METProjProp.Pt()*rj_m2_orig.Pz()/rj_m2_orig.Pt(),rj2_METProjProp.Pt()*rj_m2_orig.P()/rj_m2_orig.Pt());
	  recomet_proj_sum+=rj2_METProjProp;
	}else{
	  if(tempRecoMETP4.Vect().Dot(rj_m2_orig.Vect().Unit())>0){
	    std::cout<<"wrong hemisphere rj1 MET"<<(tempRecoMETP4.Vect().Dot(rj_m2_orig.Vect().Unit()))<<"/"<<i_entry<<" "<<acos(tempRecoMETP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit()))*TMath::RadToDeg()<<std::endl;
	  }
	  //std::cout<<"wrong hemisphere rj2 MET "<<(tempRecoMETP4.Vect().Dot(rj_m2_orig.Vect().Unit()))<<"/"<<i_entry<<" "<<acos(tempRecoMETP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit()))*TMath::RadToDeg()<<std::endl;
	}
	if(recomet_proj_sum.Pt()>tempRecoMETP4.Pt()){
	  //std::cout<<"danger, we have to do something with met projection reco "<<recomet_proj_sum.Pt()<<"/"<<tempRecoMETP4.Pt()<<std::endl;
	}
      }
      if(recomet_proj_sum.Pt()>tempRecoMETP4.Pt()){
	if((rj_m1_orig.E()+rj_m2_orig.E())>2000){
	  counter_j1j2_met_over_reco_2000+=1;
	}else{
	  counter_j1j2_met_over_reco_low_2000+=1;
	}
      }


      TLorentzVector rj1_MET=rj_m1_orig+rj1_METProjVecProp;
      TLorentzVector rj2_MET=rj_m2_orig+rj2_METProjVecProp;
      if(rj2_MET.M()>rj1_MET.M()){
	TLorentzVector temp=rj1_MET;
	rj1_MET=rj2_MET;
	rj2_MET=temp;
      }
      if((rj1_EMissProjVecProp.M()!=0 || rj2_EMissProjVecProp.M()!=0) && (rj1_EMissProjVecProp.M()>1e-3 || rj2_EMissProjVecProp.M()>1e-3)){
	std::cout<<"mass was supposed to be 0 "<<rj1_EMissProjVecProp.M()<<"/"<<rj2_EMissProjVecProp.M()<<std::endl;
      }
      TLorentzVector rj1_EMiss=rj_m1_orig+rj1_EMissProjVecProp;
      TLorentzVector rj2_EMiss=rj_m2_orig+rj2_EMissProjVecProp;
      unsigned int ind_rj1_EMiss=ind_rj1;
      unsigned int ind_rj2_EMiss=ind_rj2;

      h_hist_vec_2D[4]->Fill(rj1_EMissProjVecProp.Pt()/rj_m1_orig.Pt(),rj2_EMissProjVecProp.Pt()/rj_m2_orig.Pt(),weight);


      if(rj2_EMiss.M()>rj1_EMiss.M()){
	TLorentzVector temp=rj1_EMiss;
	rj1_EMiss=rj2_EMiss;
	rj2_EMiss=temp;
	ind_rj1_EMiss=ind_rj2;
	ind_rj2_EMiss=ind_rj1;
      }
      tempRecoEMissCorrP4=rj1_EMiss+rj2_EMiss-rj_m1_orig-rj_m2_orig;
      
      if(tempRecoMHTP4.Pt()>0){
	double rj1_MHTProj=tempRecoMHTP4.Vect().Dot(rj_m1_orig.Vect().Unit());
	TLorentzVector rj1_MHTProjVec(0,0,0,0);
	rj1_MHTProjVec.SetPxPyPzE(rj1_MHTProj*rj_m1_orig.Px()/rj_m1_orig.P(),rj1_MHTProj*rj_m1_orig.Py()/rj_m1_orig.P(),0,rj1_MHTProj);
	//tempTotInvRecoP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);
	TVector3 recomht_proj_sum(0,0,0);
	TVector3 rj1_MHTProjProp=(tempRecoMHTP4.Vect().Dot(rj_m1_orig.Vect().Unit()))*rj_m1_orig.Vect().Unit();
	if(tempRecoMHTP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit())>0){//rj1 and MHT in same hemisphere
	  //scale pz with values obtained by MHT projection for px and py along the jet z momentum axis, assume massless neutrino vector addition, aka scale MHT with p/pt ratio
	  rj1_MHMissProjVecProp.SetPxPyPzE(rj1_MHTProjProp.Px(),rj1_MHTProjProp.Py(),rj1_MHTProjProp.Pt()*rj_m1_orig.Pz()/rj_m1_orig.Pt(),rj1_MHTProjProp.Pt()*rj_m1_orig.P()/rj_m1_orig.Pt());
	  recomht_proj_sum+=rj1_MHTProjProp;
	}
	double rj2_MHTProj=tempRecoMHTP4.Vect().Dot(rj_m2_orig.Vect().Unit());
	TLorentzVector rj2_MHTProjVec(0,0,0,0);
	rj2_MHTProjVec.SetPxPyPzE(rj2_MHTProj*rj_m2_orig.Px()/rj_m2_orig.P(),rj2_MHTProj*rj_m2_orig.Py()/rj_m2_orig.P(),0,rj2_MHTProj);
	//tempTotInvRecoP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);
	TVector3 rj2_MHTProjProp=(tempRecoMHTP4.Vect().Dot(rj_m2_orig.Vect().Unit()))*rj_m2_orig.Vect().Unit();
	if(tempRecoMHTP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit())>0){//rj2 and MHT in same hemisphere
	  rj2_MHMissProjVecProp.SetPxPyPzE(rj2_MHTProjProp.Px(),rj2_MHTProjProp.Py(),rj2_MHTProjProp.Pt()*rj_m2_orig.Pz()/rj_m2_orig.Pt(),rj2_MHTProjProp.Pt()*rj_m2_orig.P()/rj_m2_orig.Pt());
	  recomht_proj_sum+=rj2_MHTProjProp;
	}
	if(recomht_proj_sum.Pt()>tempRecoMHTP4.Pt()){
	  std::cout<<"danger, we have to do something with mht projection reco "<<recomht_proj_sum.Pt()<<"/"<<tempRecoMHTP4.Pt()<<std::endl;
	}
      }
      //std::cout<<"minimizer/calcMH rj1 calc "<<rj1_corr<<"/"<<rj1_MHMissProjVecProp.Px()/rj_m1_orig.Px()<<"/"<<rj1_MHMissProjVecProp.Py()/rj_m1_orig.Py()<<std::endl;
      //std::cout<<"minimizer/calcMH rj2 calc "<<rj2_corr<<"/"<<rj2_MHMissProjVecProp.Px()/rj_m2_orig.Px()<<"/"<<rj2_MHMissProjVecProp.Py()/rj_m2_orig.Py()<<std::endl;
      TLorentzVector rj1_MHMiss=rj_m1_orig+rj1_MHMissProjVecProp;
      TLorentzVector rj2_MHMiss=rj_m2_orig+rj2_MHMissProjVecProp;
      unsigned int ind_rj1_MHMiss=ind_rj1;
      unsigned int ind_rj2_MHMiss=ind_rj2;
      if(rj2_MHMiss.M()>rj1_MHMiss.M()){
	TLorentzVector temp=rj1_MHMiss;
	rj1_MHMiss=rj2_MHMiss;
	rj2_MHMiss=temp;
	ind_rj1_MHMiss=ind_rj2;
	ind_rj2_MHMiss=ind_rj1;
      }
      tempRecoMHMissCorrP4=rj1_MHMiss+rj2_MHMiss-rj_m1_orig-rj_m2_orig;
      
      //double px1=(1+rj1_corr)*rj_m1_orig.Px()+(1+rj2_corr)*rj_m2_orig.Px()-tempRecoIsoPhP4.Px();
      //double py2=(1+rj1_corr)*rj_m1_orig.Py()+(1+rj2_corr)*rj_m2_orig.Py()-tempRecoIsoPhP4.Py();
      
      //double met_reco =sqrt(px1*px1+py2*py2);
      
      //std::cout<<"before/minimizer MHMiss "<<(rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).Pt()<<"/"<<met_reco<<" EMiss/MHMiss "<<(rj1_EMiss+rj2_EMiss-tempRecoIsoPhP4).Pt()<<"/"<<(rj_m1_orig+rj1_MHMissProjVecProp+rj_m2_orig+rj2_MHMissProjVecProp-tempRecoIsoPhP4).Pt() <<std::endl;
      
      if(tempRecoMHMissCorrP4.Pt()==0){
	//std::cout<<"reco MC is 0, but RecoEMissCorr is not "<<tempRecoEMissCorrP4.Pt()<<"/"<<tempRecoEMissCorrP4.Pt()/(rj_m1_orig+rj_m2_orig).Pt()<<" WTF "<<tempRecoMHTP4.Pt()<<"/"<<tempRecoMHTP4.Px()<<"/"<<tempRecoMHTP4.Py()<<"/"<<(rj_m1_orig+rj_m2_orig).Px()<<"/"<< (rj_m1_orig+rj_m2_orig).Py()<<std::endl;
	double rj1_MHTProj=tempRecoMHTP4.Vect().Dot(rj_m1_orig.Vect().Unit());
	TLorentzVector rj1_MHTProjVec(0,0,0,0);
	rj1_MHTProjVec.SetPxPyPzE(rj1_MHTProj*rj_m1_orig.Px()/rj_m1_orig.P(),rj1_MHTProj*rj_m1_orig.Py()/rj_m1_orig.P(),0,rj1_MHTProj);
	//std::cout<<"proj r1 "<<rj1_MHTProj<<" x/y/z/E/pt "<<rj1_MHTProjVec.Px()<<"/"<<rj1_MHTProjVec.Py()<<"/"<<rj1_MHTProjVec.Pz()<<"/"<<rj1_MHTProjVec.E()<<"/"<<rj1_MHTProjVec.Pt()<<std::endl;
	//tempTotInvRecoP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);
	TVector3 rj1_MHTProjProp=(tempRecoMHTP4.Vect().Dot(rj_m1_orig.Vect().Unit()))*rj_m1_orig.Vect().Unit();
	if(tempRecoMHTP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit())>0){//rj1 and MHT in same hemisphere
	  //scale pz with values obtained by MHT projection for px and py along the jet z momentum axis, assume massless neutrino vector addition, aka scale MHT with p/pt ratio
	  rj1_MHMissProjVecProp.SetPxPyPzE(rj1_MHTProjProp.Px(),rj1_MHTProjProp.Py(),rj1_MHTProjProp.Pt()*rj_m1_orig.Pz()/rj_m1_orig.Pt(),rj1_MHTProjProp.Pt()*rj_m1_orig.P()/rj_m1_orig.Pt());
	}else{
	  if((tempRecoMHTP4.Vect().Dot(rj_m1_orig.Vect().Unit()))>0){
	    std::cout<<"what is happening rj1 MHT "<<i_entry<<"/"<<tempRecoMHTP4.Vect().Dot(rj_m1_orig.Vect().Unit())<<std::endl;
	  }
	  //std::cout<<"wrong hemisphere rj1 "<<(tempRecoMHTP4.Vect().Dot(rj_m1_orig.Vect().Unit()))<<"/"<<acos(tempRecoMHTP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit()))*TMath::RadToDeg()<<" val "<<rj1_MHTProj<<"/"<<i_entry<<std::endl;
	}
	//if(tempRecoMHTP4.Vect().Dot(rj_m1_orig.Vect().Unit())<0){
	//std::cout<<"what is happening "<<i_entry<<"/"<<tempRecoMHTP4.Vect().Dot(rj_m1_orig.Vect().Unit())<<"/"<<acos(tempRecoMHTP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit()))*TMath::RadToDeg()<<std::endl;
	//}
	double rj2_MHTProj=tempRecoMHTP4.Vect().Dot(rj_m2_orig.Vect().Unit());
	TLorentzVector rj2_MHTProjVec(0,0,0,0);
	rj2_MHTProjVec.SetPxPyPzE(rj2_MHTProj*rj_m2_orig.Px()/rj_m2_orig.P(),rj2_MHTProj*rj_m2_orig.Py()/rj_m2_orig.P(),0,rj2_MHTProj);
	//std::cout<<"proj r2 "<<rj2_MHTProj<<" x/y/z/E/pt "<<rj2_MHTProjVec.Px()<<"/"<<rj2_MHTProjVec.Py()<<"/"<<rj2_MHTProjVec.Pz()<<"/"<<rj2_MHTProjVec.E()<<"/"<<rj2_MHTProjVec.Pt()<<std::endl;
	//tempTotInvRecoP4.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E);
	TVector3 rj2_MHTProjProp=(tempRecoMHTP4.Vect().Dot(rj_m2_orig.Vect().Unit()))*rj_m2_orig.Vect().Unit();
	//if(tempRecoMHTP4.Vect().Dot(rj_m2_orig.Vect().Unit())<0){
	//std::cout<<"what is happening "<<i_entry<<"/"<<tempRecoMHTP4.Vect().Dot(rj_m2_orig.Vect().Unit())<<"/"<<acos(tempRecoMHTP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit()))*TMath::RadToDeg()<<std::endl;
	//}
	if(tempRecoMHTP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit())>0){//rj2 and MHT in same hemisphere
	  rj2_MHMissProjVecProp.SetPxPyPzE(rj2_MHTProjProp.Px(),rj2_MHTProjProp.Py(),rj2_MHTProjProp.Pt()*rj_m2_orig.Pz()/rj_m2_orig.Pt(),rj2_MHTProjProp.Pt()*rj_m2_orig.P()/rj_m2_orig.Pt());
	  //std::cout<<"right hemisphere rj2 "<<acos(tempRecoMHTP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit()))*TMath::RadToDeg()<<" val "<<rj2_MHTProj<<"/"<<i_entry<<std::endl;
	}else{
	  if((tempRecoMHTP4.Vect().Dot(rj_m2_orig.Vect().Unit()))>0){
	    std::cout<<"what is happening rj2 MHT "<<i_entry<<"/"<<tempRecoMHTP4.Vect().Dot(rj_m2_orig.Vect().Unit())<<std::endl;
	  }
	  //std::cout<<"wrong hemisphere rj2 "<<(tempRecoMHTP4.Vect().Dot(rj_m2_orig.Vect().Unit()))<<"/"<<acos(tempRecoMHTP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit()))*TMath::RadToDeg()<<" val "<<rj2_MHTProj<<"/"<<i_entry<<std::endl;
	}
      }

      //h_hist_vec_2D[3]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenTruePhP4).M(),weight);
 


      if(use_sqrtJets){
	 sqrtS_eff_reco=(rj_m1_orig+rj_m2_orig).M();
       }
      //if EMiss neutrino projection is used, default for jet masses goes to EMiss corrected jets
      //sqrtS calculation is updated with neutrino projected four momenta too
      if(use_EMissNeutrinoProjection){
	rj_m1=rj1_EMiss;
	rj_m2=rj2_EMiss;
	ind_rj1=ind_rj1_EMiss;
	ind_rj2=ind_rj2_EMiss;
	if(use_sqrtJets){
	  sqrtS_eff_reco=(rj_m1+rj_m2+tempRecoEMissCorrP4).M();
	}else{
	  sqrtS_eff_reco=(tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M();
	}
      }
      if(use_MHMiss_over_PFOMiss){
	rj_m1=rj1_MHMiss;
	rj_m2=rj2_MHMiss;
	ind_rj1=ind_rj1_MHMiss;
	ind_rj2=ind_rj2_MHMiss;
	if(use_sqrtJets){
	  sqrtS_eff_reco=(rj_m1+rj_m2+tempRecoMHMissCorrP4).M();
	}else{
	  sqrtS_eff_reco=(tempTotRecoP4-tempRecoIsoPhP4+tempRecoMHMissCorrP4).M();
	}
      }
      if(usePartonInfo && fill_partonInfo){
	if(recojet_E->size()==2 && n_IsoLep_reco<m_cut_nLeptons){
	  if(tempTotEventP4.M()<sqrtS_low){
	    h_hist_parton[125]->Fill((rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).Pt()/(tempTotRecoP4-tempRecoIsoPhP4).E(),weight);
	    h_hist_parton[131]->Fill((rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).Pt()/(rj_m1_orig.E()+rj_m2_orig.E()),weight);
	    h_hist_parton[137]->Fill((rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).Pt(),weight);
	    h_hist_parton[143]->Fill((tempTotRecoP4-tempRecoIsoPhP4).E(),weight);
	    h_hist_parton[149]->Fill(rj_m1_orig.E()+rj_m2_orig.E(),weight);	     
	  }else if(tempTotEventP4.M()<sqrtS_high){
	    h_hist_parton[126]->Fill((rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).Pt()/(tempTotRecoP4-tempRecoIsoPhP4).E(),weight);
	    h_hist_parton[132]->Fill((rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).Pt()/(rj_m1_orig.E()+rj_m2_orig.E()),weight);
	    h_hist_parton[138]->Fill((rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).Pt(),weight);
	    h_hist_parton[144]->Fill((tempTotRecoP4-tempRecoIsoPhP4).E(),weight);
	    h_hist_parton[150]->Fill(rj_m1_orig.E()+rj_m2_orig.E(),weight);	    
	  }else{
	    h_hist_parton[127]->Fill((rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).Pt()/(tempTotRecoP4-tempRecoIsoPhP4).E(),weight);
	    h_hist_parton[133]->Fill((rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).Pt()/(rj_m1_orig.E()+rj_m2_orig.E()),weight);
	    h_hist_parton[139]->Fill((rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).Pt(),weight);
	    h_hist_parton[145]->Fill((tempTotRecoP4-tempRecoIsoPhP4).E(),weight);
	    h_hist_parton[151]->Fill(rj_m1_orig.E()+rj_m2_orig.E(),weight);	    
	  }
	}
	if(sqrtS_eff_reco>sqrtS_high){
	  h_hist_parton[89]->Fill(DeltaPhi(rj_m1.Phi(),tempRecoMETP4.Phi())*TMath::RadToDeg(),weight);
	  h_hist_parton[91]->Fill(DeltaPhi(rj_m2.Phi(),tempRecoMETP4.Phi())*TMath::RadToDeg(),weight);
	  if((rj_m1.Angle(tempHP4.Vect())*TMath::RadToDeg())<25.){//jet 1 and H spatially matched
	    h_hist_parton[93]->Fill(DeltaPhi(rj_m1.Phi(),tempRecoMETP4.Phi())*TMath::RadToDeg(),weight);
	  }else if((rj_m1.Angle(tempHP4.Vect())*TMath::RadToDeg())>155.){//jet 1 and H spatially unmatched
	    h_hist_parton[95]->Fill(DeltaPhi(rj_m1.Phi(),tempRecoMETP4.Phi())*TMath::RadToDeg(),weight);
	  }
	  if((rj_m2.Angle(tempZP4.Vect())*TMath::RadToDeg())<25.){//jet 2 and Z spatially matched
	    h_hist_parton[97]->Fill(DeltaPhi(rj_m2.Phi(),tempRecoMETP4.Phi())*TMath::RadToDeg(),weight);
	  }else if((rj_m2.Angle(tempZP4.Vect())*TMath::RadToDeg())>155.){//jet 2 and Z spatially unmatched
	    h_hist_parton[99]->Fill(DeltaPhi(rj_m2.Phi(),tempRecoMETP4.Phi())*TMath::RadToDeg(),weight);
	  }
	  if((rj_m1.Angle(tempHP4.Vect())*TMath::RadToDeg())<25.){//jet 1 and H spatially matched
	    h_hist_parton[101]->Fill(DeltaPhi(rj_m1.Phi(),tempRecoMETP4.Phi())*TMath::RadToDeg(),weight);
	  }
	  if((rj_m2.Angle(tempHP4.Vect())*TMath::RadToDeg())<25.){//jet 2 and H spatially matched
	    h_hist_parton[101]->Fill(DeltaPhi(rj_m2.Phi(),tempRecoMETP4.Phi())*TMath::RadToDeg(),weight);
	  }
	  if((rj_m1.Angle(tempZP4.Vect())*TMath::RadToDeg())<25.){//jet 1 and Z spatially matched
	    h_hist_parton[103]->Fill(DeltaPhi(rj_m1.Phi(),tempRecoMETP4.Phi())*TMath::RadToDeg(),weight);
	  }
	  if((rj_m2.Angle(tempZP4.Vect())*TMath::RadToDeg())<25.){//jet 2 and Z spatially matched
	    h_hist_parton[103]->Fill(DeltaPhi(rj_m2.Phi(),tempRecoMETP4.Phi())*TMath::RadToDeg(),weight);
	}
	}
	if(tempTotEventP4.M()>sqrtS_high){
	  h_hist_parton[40]->Fill(rj_m1.M(),weight);
	  h_hist_parton[41]->Fill(rj_m2.M(),weight);
	  h_hist_parton[46]->Fill(DeltaPhi(rj_m1.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	  h_hist_parton[47]->Fill(fabs(rj_m1.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	  h_hist_parton[48]->Fill(DeltaPhi(rj_m1.Phi(),tempHP4.Phi())*TMath::RadToDeg(),weight);
	  h_hist_parton[49]->Fill(fabs(rj_m1.Theta()-tempHP4.Theta())*TMath::RadToDeg(),weight);
	  h_hist_parton[52]->Fill(rj_m1.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	  h_hist_parton[53]->Fill(rj_m2.Angle(tempZP4.Vect())*TMath::RadToDeg(),weight);
	  h_hist_parton[106]->Fill(rj1_EMiss.Angle(tempHP4.Vect())*TMath::RadToDeg(),weight);
	  h_hist_parton[107]->Fill(rj2_EMiss.Angle(tempZP4.Vect())*TMath::RadToDeg(),weight);
	}
      }
      
      //if too few components in jet, then jet index NOT found -->i.e. remains at -1
      for(unsigned int i=0;i<recojet_subjet_E->size();i++){
	if((*recojet_subjet_jetindex)[i]==ind_rj1){
	  if((*recojet_subjet_E)[i]>E_sj1_rj1){
	    ind_sj2_rj1=ind_sj1_rj1;
	    ind_sj1_rj1=i;
	    E_sj1_rj1=(*recojet_subjet_E)[i];
	  }else{
	    ind_sj2_rj1=i;
	  }
	}
	if((*recojet_subjet_jetindex)[i]==ind_rj2){
	  if((*recojet_subjet_E)[i]>E_sj1_rj2){
	    ind_sj2_rj2=ind_sj1_rj2;
	    ind_sj1_rj2=i;
	    E_sj1_rj2=(*recojet_subjet_E)[i];
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
      
      //check if something is wrong here --> somehow boosts etc look weird
      TVector3 boost_rj1_COM = -rj_m1_orig.BoostVector();
      TVector3 boost_rj2_COM = -rj_m2_orig.BoostVector();
      if(use_EMissNeutrinoProjection && ind_rj1!=ind_rj1_orig){
	boost_rj1_COM = -rj_m2_orig.BoostVector();
	boost_rj2_COM = -rj_m1_orig.BoostVector();
      }
      if(use_MHMiss_over_PFOMiss && ind_rj1!=ind_rj1_orig){
	boost_rj1_COM = -rj_m2_orig.BoostVector();
	boost_rj2_COM = -rj_m1_orig.BoostVector();
      }
      TLorentzVector temp_sj1_rj1_boost_rj1_COM=temp_sj1_rj1;
      temp_sj1_rj1_boost_rj1_COM.Boost(boost_rj1_COM);
      TLorentzVector temp_sj2_rj1_boost_rj1_COM=temp_sj2_rj1;
      temp_sj2_rj1_boost_rj1_COM.Boost(boost_rj1_COM);
      TLorentzVector temp_sj1_rj2_boost_rj2_COM=temp_sj1_rj2;
      temp_sj1_rj2_boost_rj2_COM.Boost(boost_rj2_COM);
      TLorentzVector temp_sj2_rj2_boost_rj2_COM=temp_sj2_rj2;
      temp_sj2_rj2_boost_rj2_COM.Boost(boost_rj2_COM);

 
      bool veto_qqqq_mass_ellipse_reco = false;
      if( (rj_m2.M()<cut_m_j2_c_veto && ((pow((rj_m2.M()-cut_m_j2_c_veto)/cut_delta_m_j2_e_veto_minus,2)+pow((rj_m1.M()-cut_m_j1_c_veto)/cut_delta_m_j1_e_veto,2))<1.)) || (rj_m2.M()>cut_m_j2_c_veto && ((pow((rj_m2.M()-cut_m_j2_c_veto)/cut_delta_m_j2_e_veto_plus,2)+pow((rj_m1.M()-cut_m_j1_c_veto)/cut_delta_m_j1_e_veto,2))<1.))){
	veto_qqqq_mass_ellipse_reco =true;
	//std::cout<<"veto reco is hit"<<std::endl;
      }

      if( (pow((rj_m2.M()-cut_m_j2_c)/cut_delta_m_j2_e,2)+pow((rj_m1.M()-cut_m_j1_c)/cut_delta_m_j1_e,2))<1. && !veto_qqqq_mass_ellipse_reco  && (fabs((rj_m2.Theta()-0.5*TMath::Pi())*TMath::RadToDeg())<cut_theta_window && (fabs(rj_m1.Theta()-0.5*TMath::Pi())*TMath::RadToDeg())<cut_theta_window)){
	  reco_pass_mass_cuts=true;
	}

      //if((rj_m2.M()> cut_m_j2_min && rj_m2.M()< cut_m_j2_max) && (rj_m1.M()> cut_m_j1_min && rj_m1.M()< cut_m_j1_max) && ((rj_m1.M()-rj_m2.M())>cut_delta_mass_low && (rj_m1.M()-rj_m2.M())<cut_delta_mass_high)){
      //reco_pass_mass_cuts=true;
      //}
 
      if(!performMassCuts){
	reco_pass_mass_cuts=true;
      } //if(ind_sj1_rj1==ind_sj1_rj2 || ind_sj1_rj1==ind_sj2_rj1 || ind_sj1_rj1==ind_sj2_rj2 || ind_sj2_rj1==ind_sj1_rj2 || ind_sj2_rj1==ind_sj2_rj2 ||ind_sj1_rj2==ind_sj2_rj2){  
      //std::cout<<"what the hell is going on with reco subjet indices "<<ind_sj1_rj1<<"/"<<ind_sj2_rj1<<"/"<<ind_sj1_rj2<<"/"<<ind_sj2_rj2<<std::endl;
      //}
 


      if(recojet_E->size()==2 && n_IsoLep_reco<m_cut_nLeptons && reco_pass_mass_cuts){//effectively no cut here, but exactly two recojets
	TLorentzVector temp_sj;
	if(usePartonInfo && tempTotEventP4.M()>sqrtS_low && tempTotEventP4.M()>sqrtS_high){
	  if(ind_sj1_rj1!=-1){	    
	    temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj1],(*recojet_subjet_Py)[ind_sj1_rj1],(*recojet_subjet_Pz)[ind_sj1_rj1],(*recojet_subjet_E)[ind_sj1_rj1]);
	    if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_gj_Z_q_pos){
	      angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	      ind_Z_q_pos_sj_rj_angle=ind_sj1_rj1;
	    }
	    if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
	      angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	      ind_Z_q_neg_sj_rj_angle=ind_sj1_rj1;
	    }
	  }
	  if(ind_sj2_rj1!=-1){
	    temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj1],(*recojet_subjet_Py)[ind_sj2_rj1],(*recojet_subjet_Pz)[ind_sj2_rj1],(*recojet_subjet_E)[ind_sj2_rj1]);
	    if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
	      angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	      ind_Z_q_pos_sj_rj_angle=ind_sj2_rj1;
	    }
	    if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
	      angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	      ind_Z_q_neg_sj_rj_angle=ind_sj2_rj1;
	    }
	  }
	  if(ind_sj1_rj2!=-1){
	    temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj2],(*recojet_subjet_Py)[ind_sj1_rj2],(*recojet_subjet_Pz)[ind_sj1_rj2],(*recojet_subjet_E)[ind_sj1_rj2]);
	    if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
	      angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	      ind_Z_q_pos_sj_rj_angle=ind_sj1_rj2;
	    }
	    if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
	      angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	      ind_Z_q_neg_sj_rj_angle=ind_sj1_rj2;
	    }
	  }
	  if(ind_sj2_rj2!=-1){
	    temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj2],(*recojet_subjet_Py)[ind_sj2_rj2],(*recojet_subjet_Pz)[ind_sj2_rj2],(*recojet_subjet_E)[ind_sj2_rj2]);
	    if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
	      angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	      ind_Z_q_pos_sj_rj_angle=ind_sj2_rj2;
	    }
	    if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
	      angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	      ind_Z_q_neg_sj_rj_angle=ind_sj2_rj2;
	    }
	  }
	  if(ind_Z_q_pos_sj_rj_angle==ind_Z_q_neg_sj_rj_angle){
	    //decide which angle is the more propable
	    if(angle_min_sj_rj_Z_q_pos<angle_min_sj_rj_Z_q_neg){
	      if(ind_sj1_rj1!=-1 && ind_sj1_rj1!=ind_Z_q_pos_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj1],(*recojet_subjet_Py)[ind_sj1_rj1],(*recojet_subjet_Pz)[ind_sj1_rj1],(*recojet_subjet_E)[ind_sj1_rj1]);
		if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
		  angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
		  ind_Z_q_neg_sj_rj_angle=ind_sj1_rj1;
		}
	      }
	      if(ind_sj2_rj1!=-1 && ind_sj2_rj1!=ind_Z_q_pos_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj1],(*recojet_subjet_Py)[ind_sj2_rj1],(*recojet_subjet_Pz)[ind_sj2_rj1],(*recojet_subjet_E)[ind_sj2_rj1]);
		if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
		  angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
		  ind_Z_q_neg_sj_rj_angle=ind_sj2_rj1;
		}
	      }
	      if(ind_sj1_rj2!=-1 && ind_sj1_rj2!=ind_Z_q_pos_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj2],(*recojet_subjet_Py)[ind_sj1_rj2],(*recojet_subjet_Pz)[ind_sj1_rj2],(*recojet_subjet_E)[ind_sj1_rj2]);
		if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
		  angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
		  ind_Z_q_neg_sj_rj_angle=ind_sj1_rj2;
		}
	      }
	      if(ind_sj2_rj2!=-1 && ind_sj2_rj2!=ind_Z_q_pos_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj2],(*recojet_subjet_Py)[ind_sj2_rj2],(*recojet_subjet_Pz)[ind_sj2_rj2],(*recojet_subjet_E)[ind_sj2_rj2]);
		if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
		  angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
		  ind_Z_q_neg_sj_rj_angle=ind_sj2_rj2;
		}
	      }
	    }else{//closer to b
	      if(ind_sj1_rj1!=-1 && ind_sj1_rj1!=ind_Z_q_neg_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj1],(*recojet_subjet_Py)[ind_sj1_rj1],(*recojet_subjet_Pz)[ind_sj1_rj1],(*recojet_subjet_E)[ind_sj1_rj1]);
		if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
		  angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
		  ind_Z_q_pos_sj_rj_angle=ind_sj1_rj1;
		}
	      }
	      if(ind_sj2_rj1!=-1 && ind_sj2_rj1!=ind_Z_q_neg_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj1],(*recojet_subjet_Py)[ind_sj2_rj1],(*recojet_subjet_Pz)[ind_sj2_rj1],(*recojet_subjet_E)[ind_sj2_rj1]);
		if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
		  angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
		  ind_Z_q_pos_sj_rj_angle=ind_sj2_rj1;
		}
	      }
	      if(ind_sj1_rj2!=-1 && ind_sj1_rj2!=ind_Z_q_neg_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj2],(*recojet_subjet_Py)[ind_sj1_rj2],(*recojet_subjet_Pz)[ind_sj1_rj2],(*recojet_subjet_E)[ind_sj1_rj2]);
		if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
		  angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
		  ind_Z_q_pos_sj_rj_angle=ind_sj1_rj2;
		}
	      }
	      if(ind_sj2_rj2!=-1 && ind_sj2_rj2!=ind_Z_q_neg_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj2],(*recojet_subjet_Py)[ind_sj2_rj2],(*recojet_subjet_Pz)[ind_sj2_rj2],(*recojet_subjet_E)[ind_sj2_rj2]);
		if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
		  angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
		  ind_Z_q_pos_sj_rj_angle=ind_sj2_rj2;
		}
	      }
	    }
	  }
	}
	if(usePartonInfo && tempTotEventP4.M()>sqrtS_low && tempTotEventP4.M()>sqrtS_high){
	  if(ind_sj1_rj1!=-1){	    
	    temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj1],(*recojet_subjet_Py)[ind_sj1_rj1],(*recojet_subjet_Pz)[ind_sj1_rj1],(*recojet_subjet_E)[ind_sj1_rj1]);
	    if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_rj_H_bbar){
	      angle_min_sj_rj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
	      ind_H_bbar_sj_rj_angle=ind_sj1_rj1;
	    }
	    if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_rj_H_b){
	      angle_min_sj_rj_H_b=temp_sj.Angle(tempH_b.Vect());
	      ind_H_b_sj_rj_angle=ind_sj1_rj1;
	    }
	  }
	  if(ind_sj2_rj1!=-1){
	    temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj1],(*recojet_subjet_Py)[ind_sj2_rj1],(*recojet_subjet_Pz)[ind_sj2_rj1],(*recojet_subjet_E)[ind_sj2_rj1]);
	    if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_rj_H_bbar){
	      angle_min_sj_rj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
	      ind_H_bbar_sj_rj_angle=ind_sj2_rj1;
	    }
	    if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_rj_H_b){
	      angle_min_sj_rj_H_b=temp_sj.Angle(tempH_b.Vect());
	      ind_H_b_sj_rj_angle=ind_sj2_rj1;
	    }
	  }
	  if(ind_sj1_rj2!=-1){
	    temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj2],(*recojet_subjet_Py)[ind_sj1_rj2],(*recojet_subjet_Pz)[ind_sj1_rj2],(*recojet_subjet_E)[ind_sj1_rj2]);
	    if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_rj_H_bbar){
	      angle_min_sj_rj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
	      ind_H_bbar_sj_rj_angle=ind_sj1_rj2;
	    }
	    if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_rj_H_b){
	      angle_min_sj_rj_H_b=temp_sj.Angle(tempH_b.Vect());
	      ind_H_b_sj_rj_angle=ind_sj1_rj2;
	    }
	  }
	  if(ind_sj2_rj2!=-1){
	    temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj2],(*recojet_subjet_Py)[ind_sj2_rj2],(*recojet_subjet_Pz)[ind_sj2_rj2],(*recojet_subjet_E)[ind_sj2_rj2]);
	    if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_rj_H_bbar){
	      angle_min_sj_rj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
	      ind_H_bbar_sj_rj_angle=ind_sj2_rj2;
	    }
	    if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_rj_H_b){
	      angle_min_sj_rj_H_b=temp_sj.Angle(tempH_b.Vect());
	      ind_H_b_sj_rj_angle=ind_sj2_rj2;
	    }
	  }
	  if(ind_H_bbar_sj_rj_angle==ind_H_b_sj_rj_angle){
	    //decide which angle is the more propable
	    if(angle_min_sj_rj_H_bbar<angle_min_sj_rj_H_b){
	      if(ind_sj1_rj1!=-1 && ind_sj1_rj1!=ind_H_bbar_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj1],(*recojet_subjet_Py)[ind_sj1_rj1],(*recojet_subjet_Pz)[ind_sj1_rj1],(*recojet_subjet_E)[ind_sj1_rj1]);
		if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_rj_H_b){
		  angle_min_sj_rj_H_b=temp_sj.Angle(tempH_b.Vect());
		  ind_H_b_sj_rj_angle=ind_sj1_rj1;
		}
	      }
	      if(ind_sj2_rj1!=-1 && ind_sj2_rj1!=ind_H_bbar_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj1],(*recojet_subjet_Py)[ind_sj2_rj1],(*recojet_subjet_Pz)[ind_sj2_rj1],(*recojet_subjet_E)[ind_sj2_rj1]);
		if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_rj_H_b){
		  angle_min_sj_rj_H_b=temp_sj.Angle(tempH_b.Vect());
		  ind_H_b_sj_rj_angle=ind_sj2_rj1;
		}
	      }
	      if(ind_sj1_rj2!=-1 && ind_sj1_rj2!=ind_H_bbar_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj2],(*recojet_subjet_Py)[ind_sj1_rj2],(*recojet_subjet_Pz)[ind_sj1_rj2],(*recojet_subjet_E)[ind_sj1_rj2]);
		if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_rj_H_b){
		  angle_min_sj_rj_H_b=temp_sj.Angle(tempH_b.Vect());
		  ind_H_b_sj_rj_angle=ind_sj1_rj2;
		}
	      }
	      if(ind_sj2_rj2!=-1 && ind_sj2_rj2!=ind_H_bbar_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj2],(*recojet_subjet_Py)[ind_sj2_rj2],(*recojet_subjet_Pz)[ind_sj2_rj2],(*recojet_subjet_E)[ind_sj2_rj2]);
		if(temp_sj.Angle(tempH_b.Vect())<angle_min_sj_rj_H_b){
		  angle_min_sj_rj_H_b=temp_sj.Angle(tempH_b.Vect());
		  ind_H_b_sj_rj_angle=ind_sj2_rj2;
		}
	      }
	    }else{//closer to b
	      if(ind_sj1_rj1!=-1 && ind_sj1_rj1!=ind_H_b_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj1],(*recojet_subjet_Py)[ind_sj1_rj1],(*recojet_subjet_Pz)[ind_sj1_rj1],(*recojet_subjet_E)[ind_sj1_rj1]);
		if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_rj_H_bbar){
		  angle_min_sj_rj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
		  ind_H_bbar_sj_rj_angle=ind_sj1_rj1;
		}
	      }
	      if(ind_sj2_rj1!=-1 && ind_sj2_rj1!=ind_H_b_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj1],(*recojet_subjet_Py)[ind_sj2_rj1],(*recojet_subjet_Pz)[ind_sj2_rj1],(*recojet_subjet_E)[ind_sj2_rj1]);
		if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_rj_H_bbar){
		  angle_min_sj_rj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
		  ind_H_bbar_sj_rj_angle=ind_sj2_rj1;
		}
	      }
	      if(ind_sj1_rj2!=-1 && ind_sj1_rj2!=ind_H_b_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj2],(*recojet_subjet_Py)[ind_sj1_rj2],(*recojet_subjet_Pz)[ind_sj1_rj2],(*recojet_subjet_E)[ind_sj1_rj2]);
		if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_rj_H_bbar){
		  angle_min_sj_rj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
		  ind_H_bbar_sj_rj_angle=ind_sj1_rj2;
		}
	      }
	      if(ind_sj2_rj2!=-1 && ind_sj2_rj2!=ind_H_b_sj_rj_angle){	    
		temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj2],(*recojet_subjet_Py)[ind_sj2_rj2],(*recojet_subjet_Pz)[ind_sj2_rj2],(*recojet_subjet_E)[ind_sj2_rj2]);
		if(temp_sj.Angle(tempH_bbar.Vect())<angle_min_sj_rj_H_bbar){
		  angle_min_sj_rj_H_bbar=temp_sj.Angle(tempH_bbar.Vect());
		  ind_H_bbar_sj_rj_angle=ind_sj2_rj2;
		}
	      }
	    }
	  }
	  if(tempZ_q_neg.Angle(tempZ_q_pos.Vect())*TMath::RadToDeg()<30.){
	    if(ind_Z_q_pos_sj_rj_angle!=-1/* && (angle_min_sj_rj_Z_q_pos*TMath::RadToDeg())<15.0*/){
	      if((*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_rj_angle]==0){
		//std::cout<<"event/ind/inangle "<<i_entry<<"/"<<ind_Z_q_pos_sj_rj_angle<<"/"<<angle_min_sj_rj_Z_q_pos*TMath::RadToDeg()<<std::endl;
	      }
	      h_hist_vec_reco[391]->Fill((*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_rj_angle],weight);
	      h_hist_vec_reco[392]->Fill((*recojet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_pos_sj_rj_angle],weight);
	      h_hist_vec_reco[397]->Fill((*recojet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_pos_sj_rj_angle],weight);
	      h_hist_vec_reco[398]->Fill((*recojet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_pos_sj_rj_angle],weight);
	      h_hist_vec_reco[399]->Fill((*recojet_subjet_jetChargePt_kappa_0_25)[ind_Z_q_pos_sj_rj_angle],weight);
	      h_hist_vec_reco[400]->Fill((*recojet_subjet_jetChargePt_kappa_0_50)[ind_Z_q_pos_sj_rj_angle],weight);
	      h_hist_vec_reco[405]->Fill((*recojet_subjet_jetChargePt_kappa_0_20)[ind_Z_q_pos_sj_rj_angle],weight);
	      h_hist_vec_reco[406]->Fill((*recojet_subjet_jetChargePt_kappa_0_30)[ind_Z_q_pos_sj_rj_angle],weight);
	    }
	    if(ind_Z_q_neg_sj_rj_angle!=-1/* && (angle_min_sj_rj_Z_q_neg*TMath::RadToDeg())<15.0*/){
	      if((*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_rj_angle]==0){
		//std::cout<<"event/ind/inangle Z_q_pos "<<i_entry<<"/"<<ind_Z_q_neg_sj_rj_angle<<"/"<<angle_min_sj_rj_Z_q_neg*TMath::RadToDeg()<<std::endl;
	      }
	      h_hist_vec_reco[407]->Fill((*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_rj_angle],weight);
	      h_hist_vec_reco[408]->Fill((*recojet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_neg_sj_rj_angle],weight);
	      h_hist_vec_reco[413]->Fill((*recojet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_neg_sj_rj_angle],weight);
	      h_hist_vec_reco[414]->Fill((*recojet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_neg_sj_rj_angle],weight);
	      h_hist_vec_reco[415]->Fill((*recojet_subjet_jetChargePt_kappa_0_25)[ind_Z_q_neg_sj_rj_angle],weight);
	      h_hist_vec_reco[416]->Fill((*recojet_subjet_jetChargePt_kappa_0_50)[ind_Z_q_neg_sj_rj_angle],weight);
	      h_hist_vec_reco[421]->Fill((*recojet_subjet_jetChargePt_kappa_0_20)[ind_Z_q_neg_sj_rj_angle],weight);
	      h_hist_vec_reco[422]->Fill((*recojet_subjet_jetChargePt_kappa_0_30)[ind_Z_q_neg_sj_rj_angle],weight);
	    }
	  }
	  //parton matching used here too
	  if(fill_genInfo && genjet_E->size()==2 && n_IsoLep_gen<m_cut_nLeptons && gen_pass_mass_cuts){
	    int ind_HJ_reco=-1;
	    int ind_ZJ_reco=-1;
	    int ind_HJ_gen=-1;
	    int ind_ZJ_gen=-1;
	    if(rj_m1.Angle(tempHP4.Vect())<rj_m2.Angle(tempHP4.Vect())){//r1 matched to higgs
	      ind_HJ_reco=ind_rj1;
	    }else{
	      ind_HJ_reco=ind_rj2;
	    }
	    if(rj_m1.Angle(tempZP4.Vect())<rj_m2.Angle(tempZP4.Vect())){//r1 matched to higgs
	      ind_ZJ_reco=ind_rj1;
	    }else{
	      ind_ZJ_reco=ind_rj2;
	    }
	    if(gj_m1.Angle(tempHP4.Vect())<gj_m2.Angle(tempHP4.Vect())){//r1 matched to higgs
	      ind_HJ_gen=ind_gj1;
	    }else{
	      ind_HJ_gen=ind_gj2;
	    }
	    if(gj_m1.Angle(tempZP4.Vect())<gj_m2.Angle(tempZP4.Vect())){//r1 matched to higgs
	      ind_ZJ_gen=ind_gj1;
	    }else{
	      ind_ZJ_gen=ind_gj2;
	    }
	    if(ind_HJ_reco!=-1 && ind_HJ_gen!=-1){
	      if(sqrtS_eff_gen<sqrtS_low){
		h_hist_vec_2D[23]->Fill((*genjet_Mult)[ind_HJ_gen]-(*genjet_NPh)[ind_HJ_gen]-(*genjet_NNH)[ind_HJ_gen],(*recojet_Mult)[ind_HJ_reco]-(*recojet_NPh)[ind_HJ_reco]-(*recojet_NNH)[ind_HJ_reco],weight);
	      }else if(sqrtS_eff_gen<sqrtS_high_reco){
		h_hist_vec_2D[24]->Fill((*genjet_Mult)[ind_HJ_gen]-(*genjet_NPh)[ind_HJ_gen]-(*genjet_NNH)[ind_HJ_gen],(*recojet_Mult)[ind_HJ_reco]-(*recojet_NPh)[ind_HJ_reco]-(*recojet_NNH)[ind_HJ_reco],weight);
	      }else{
		h_hist_vec_2D[25]->Fill((*genjet_Mult)[ind_HJ_gen]-(*genjet_NPh)[ind_HJ_gen]-(*genjet_NNH)[ind_HJ_gen],(*recojet_Mult)[ind_HJ_reco]-(*recojet_NPh)[ind_HJ_reco]-(*recojet_NNH)[ind_HJ_reco],weight);
	      }
	    }
	    if(ind_ZJ_reco!=-1 && ind_ZJ_gen!=-1){
	      if(sqrtS_eff_gen<sqrtS_low){
		h_hist_vec_2D[26]->Fill((*genjet_Mult)[ind_ZJ_gen]-(*genjet_NPh)[ind_ZJ_gen]-(*genjet_NNH)[ind_ZJ_gen],(*recojet_Mult)[ind_ZJ_reco]-(*recojet_NPh)[ind_ZJ_reco]-(*recojet_NNH)[ind_ZJ_reco],weight);
	      }else if(sqrtS_eff_gen<sqrtS_high_reco){
		h_hist_vec_2D[27]->Fill((*genjet_Mult)[ind_ZJ_gen]-(*genjet_NPh)[ind_ZJ_gen]-(*genjet_NNH)[ind_ZJ_gen],(*recojet_Mult)[ind_ZJ_reco]-(*recojet_NPh)[ind_ZJ_reco]-(*recojet_NNH)[ind_ZJ_reco],weight);
	      }else{
		h_hist_vec_2D[28]->Fill((*genjet_Mult)[ind_ZJ_gen]-(*genjet_NPh)[ind_ZJ_gen]-(*genjet_NNH)[ind_ZJ_gen],(*recojet_Mult)[ind_ZJ_reco]-(*recojet_NPh)[ind_ZJ_reco]-(*recojet_NNH)[ind_ZJ_reco],weight);
	      }
	    }
	  }
	}
	//parton S cut
	h_hist_vec_reco[279]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	h_hist_vec_reco[280]->Fill((tempTotRecoP4).M(),weight);
	h_hist_vec_reco[281]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissP4).M(),weight);
	h_hist_vec_reco[294]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
	if(sqrtS_eff_reco<sqrtS_low){
	  h_hist_vec_2D[0]->Fill(rj_m1.M(),rj_m2.M(),weight);
	}else if (sqrtS_eff_reco<sqrtS_high_reco){
	  h_hist_vec_2D[1]->Fill(rj_m1.M(),rj_m2.M(),weight);
	}else{
	  h_hist_vec_2D[2]->Fill(rj_m1.M(),rj_m2.M(),weight);    
	}
       if(sqrtS_eff_reco<sqrtS_low){
	 h_hist_vec_reco[0]->Fill(rj_m1.Angle(rj_m2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[3]->Fill(DeltaPhi(rj_m1.Phi(),rj_m2.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[6]->Fill(fabs(rj_m1.Theta()-rj_m2.Theta())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[9]->Fill(rj_m1.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[12]->Fill(rj_m2.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[15]->Fill(rj_m1.M(),weight);
	 h_hist_vec_reco[18]->Fill(rj_m2.M(),weight);
	 if(ind_sj1_rj1!=-1){
	   h_hist_vec_reco[237]->Fill((*recojet_subjet_E)[ind_sj1_rj1],weight);
	   h_hist_vec_reco[249]->Fill((*recojet_subjet_E)[ind_sj1_rj1]/(*recojet_E)[ind_rj1],weight);
	   if(ind_sj2_rj1!=-1){
	     h_hist_vec_reco[255]->Fill(temp_sj1_rj1.Angle(temp_sj2_rj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj1!=-1){
	   h_hist_vec_reco[240]->Fill((*recojet_subjet_E)[ind_sj2_rj1],weight);
	 }
	 if(ind_sj1_rj2!=-1){
	   h_hist_vec_reco[243]->Fill((*recojet_subjet_E)[ind_sj1_rj2],weight);
	   h_hist_vec_reco[252]->Fill((*recojet_subjet_E)[ind_sj1_rj2]/(*recojet_E)[ind_rj2],weight);
	   if(ind_sj2_rj2!=-1){
	     h_hist_vec_reco[258]->Fill(temp_sj1_rj2.Angle(temp_sj2_rj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj2!=-1){
	   h_hist_vec_reco[246]->Fill((*recojet_subjet_E)[ind_sj2_rj2],weight);
	 }
	 //std::cout<<"jet 1 first/second angle/sum "<<min(temp_sj1_rj1_boost_rj1_COM.Angle(rj_m1.Vect()),temp_sj2_rj1_boost_rj1_COM.Angle(rj_m1.Vect()))*TMath::RadToDeg()<<"/"<<max(temp_sj1_rj1_boost_rj1_COM.Angle(rj_m1.Vect()),temp_sj2_rj1_boost_rj1_COM.Angle(rj_m1.Vect()))*TMath::RadToDeg()<<"/"<<(min(temp_sj1_rj1_boost_rj1_COM.Angle(rj_m1.Vect()),temp_sj2_rj1_boost_rj1_COM.Angle(rj_m1.Vect()))+max(temp_sj1_rj1_boost_rj1_COM.Angle(rj_m1.Vect()),temp_sj2_rj1_boost_rj1_COM.Angle(rj_m1.Vect())))*TMath::RadToDeg()<<std::endl;
	 //std::cout<<"jet 2 first/second angle/sum "<<min(temp_sj1_rj2_boost_rj2_COM.Angle(rj_m1.Vect()),temp_sj2_rj2_boost_rj2_COM.Angle(rj_m1.Vect()))*TMath::RadToDeg()<<"/"<<max(temp_sj1_rj2_boost_rj2_COM.Angle(rj_m1.Vect()),temp_sj2_rj2_boost_rj2_COM.Angle(rj_m1.Vect()))*TMath::RadToDeg()<<"/"<<(min(temp_sj1_rj2_boost_rj2_COM.Angle(rj_m1.Vect()),temp_sj2_rj2_boost_rj2_COM.Angle(rj_m1.Vect()))+max(temp_sj1_rj2_boost_rj2_COM.Angle(rj_m1.Vect()),temp_sj2_rj2_boost_rj2_COM.Angle(rj_m1.Vect())))*TMath::RadToDeg()<<std::endl;
	 h_hist_vec_reco[282]->Fill(rj1_MET.M(),weight);
	 h_hist_vec_reco[285]->Fill(rj2_MET.M(),weight);
	 h_hist_vec_reco[288]->Fill(rj1_EMiss.M(),weight);
	 h_hist_vec_reco[291]->Fill(rj2_EMiss.M(),weight);
	 //while the boost does change with energy scaling based on MET, the original direction of the vector does not
	 h_hist_vec_reco[337]->Fill(cos(min(temp_sj1_rj1_boost_rj1_COM.Angle(rj_m1.Vect()),temp_sj2_rj1_boost_rj1_COM.Angle(rj_m1.Vect()))),weight);
	 h_hist_vec_reco[340]->Fill(cos(max(temp_sj1_rj1_boost_rj1_COM.Angle(rj_m1.Vect()),temp_sj2_rj1_boost_rj1_COM.Angle(rj_m1.Vect()))),weight);
	 h_hist_vec_reco[343]->Fill(cos(min(temp_sj1_rj2_boost_rj2_COM.Angle(rj_m2.Vect()),temp_sj2_rj2_boost_rj2_COM.Angle(rj_m2.Vect()))),weight);
	 h_hist_vec_reco[346]->Fill(cos(max(temp_sj1_rj2_boost_rj2_COM.Angle(rj_m2.Vect()),temp_sj2_rj2_boost_rj2_COM.Angle(rj_m2.Vect()))),weight);
	 h_hist_vec_reco[353]->Fill(rj_m1_orig.M()-rj_m2_orig.M(),weight);
	 h_hist_vec_reco[356]->Fill(rj1_EMiss.M()-rj2_EMiss.M(),weight);
       }else if(sqrtS_eff_reco<sqrtS_high_reco){
	 h_hist_vec_reco[1]->Fill(rj_m1.Angle(rj_m2.Vect())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[4]->Fill(DeltaPhi(rj_m1.Phi(),rj_m2.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[7]->Fill(fabs(rj_m1.Theta()-rj_m2.Theta())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[10]->Fill(rj_m1.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[13]->Fill(rj_m2.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[16]->Fill(rj_m1.M(),weight);
	 h_hist_vec_reco[19]->Fill(rj_m2.M(),weight);
	 if(ind_sj1_rj1!=-1){
	   h_hist_vec_reco[238]->Fill((*recojet_subjet_E)[ind_sj1_rj1],weight);
	   h_hist_vec_reco[250]->Fill((*recojet_subjet_E)[ind_sj1_rj1]/(*recojet_E)[ind_rj1],weight);
	   if(ind_sj2_rj1!=-1){
	     h_hist_vec_reco[256]->Fill(temp_sj1_rj1.Angle(temp_sj2_rj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj1!=-1){
	   h_hist_vec_reco[241]->Fill((*recojet_subjet_E)[ind_sj2_rj1],weight);
	 }
	 if(ind_sj1_rj2!=-1){
	   h_hist_vec_reco[244]->Fill((*recojet_subjet_E)[ind_sj1_rj2],weight);
	   h_hist_vec_reco[253]->Fill((*recojet_subjet_E)[ind_sj1_rj2]/(*recojet_E)[ind_rj2],weight);
	   if(ind_sj2_rj2!=-1){
	     h_hist_vec_reco[259]->Fill(temp_sj1_rj2.Angle(temp_sj2_rj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj2!=-1){
	   h_hist_vec_reco[247]->Fill((*recojet_subjet_E)[ind_sj2_rj2],weight);
	 }
	 h_hist_vec_reco[283]->Fill(rj1_MET.M(),weight);
	 h_hist_vec_reco[286]->Fill(rj2_MET.M(),weight);
	 h_hist_vec_reco[289]->Fill(rj1_EMiss.M(),weight);
	 h_hist_vec_reco[292]->Fill(rj2_EMiss.M(),weight);
	 //while the boost does change with energy scaling based on MET, the original direction of the vector does not
	 h_hist_vec_reco[338]->Fill(cos(min(temp_sj1_rj1_boost_rj1_COM.Angle(rj_m1.Vect()),temp_sj2_rj1_boost_rj1_COM.Angle(rj_m1.Vect()))),weight);
	 h_hist_vec_reco[341]->Fill(cos(max(temp_sj1_rj1_boost_rj1_COM.Angle(rj_m1.Vect()),temp_sj2_rj1_boost_rj1_COM.Angle(rj_m1.Vect()))),weight);
	 h_hist_vec_reco[344]->Fill(cos(min(temp_sj1_rj2_boost_rj2_COM.Angle(rj_m2.Vect()),temp_sj2_rj2_boost_rj2_COM.Angle(rj_m2.Vect()))),weight);
	 h_hist_vec_reco[347]->Fill(cos(max(temp_sj1_rj2_boost_rj2_COM.Angle(rj_m2.Vect()),temp_sj2_rj2_boost_rj2_COM.Angle(rj_m2.Vect()))),weight);
	 h_hist_vec_reco[354]->Fill(rj_m1_orig.M()-rj_m2_orig.M(),weight);
	 h_hist_vec_reco[357]->Fill(rj1_EMiss.M()-rj2_EMiss.M(),weight);
       }else{
	 h_hist_vec_reco[2]->Fill(rj_m1.Angle(rj_m2.Vect())*TMath::RadToDeg(),weight);    
	 h_hist_vec_reco[5]->Fill(DeltaPhi(rj_m1.Phi(),rj_m2.Phi())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[8]->Fill(fabs(rj_m1.Theta()-rj_m2.Theta())*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[11]->Fill(rj_m1.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[14]->Fill(rj_m2.Theta()*TMath::RadToDeg(),weight);
	 h_hist_vec_reco[17]->Fill(rj_m1.M(),weight);
	 h_hist_vec_reco[20]->Fill(rj_m2.M(),weight);
	 TLorentzVector temp_sj;
	 if(ind_sj1_rj1!=-1){
	   if(usePartonInfo){
	     temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj1],(*recojet_subjet_Py)[ind_sj1_rj1],(*recojet_subjet_Pz)[ind_sj1_rj1],(*recojet_subjet_E)[ind_sj1_rj1]);
	     if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
	       angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	       ind_Z_q_pos_sj_rj_angle=ind_sj1_rj1;
	     }
	     if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
	       angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	       ind_Z_q_neg_sj_rj_angle=ind_sj1_rj1;
	     }
	   }
	   h_hist_vec_reco[239]->Fill((*recojet_subjet_E)[ind_sj1_rj1],weight);
	   h_hist_vec_reco[251]->Fill((*recojet_subjet_E)[ind_sj1_rj1]/(*recojet_E)[ind_rj1],weight);
	   if(ind_sj2_rj1!=-1){
	     h_hist_vec_reco[257]->Fill(temp_sj1_rj1.Angle(temp_sj2_rj1.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj1!=-1){
	   if(usePartonInfo){
	     temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj1],(*recojet_subjet_Py)[ind_sj2_rj1],(*recojet_subjet_Pz)[ind_sj2_rj1],(*recojet_subjet_E)[ind_sj2_rj1]);
	     if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
	       angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	       ind_Z_q_pos_sj_rj_angle=ind_sj2_rj1;
	     }
	     if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
	       angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	       ind_Z_q_neg_sj_rj_angle=ind_sj2_rj1;
	     }
	   }
	   h_hist_vec_reco[242]->Fill((*recojet_subjet_E)[ind_sj2_rj1],weight);
	 }
	 if(ind_sj1_rj2!=-1){
	   if(usePartonInfo){
	     temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj1_rj2],(*recojet_subjet_Py)[ind_sj1_rj2],(*recojet_subjet_Pz)[ind_sj1_rj2],(*recojet_subjet_E)[ind_sj1_rj2]);
	     if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
	       angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	       ind_Z_q_pos_sj_rj_angle=ind_sj1_rj2;
	     }
	     if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
	       angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	       ind_Z_q_neg_sj_rj_angle=ind_sj1_rj2;
	     }
	   }
	   h_hist_vec_reco[245]->Fill((*recojet_subjet_E)[ind_sj1_rj2],weight);
	   h_hist_vec_reco[254]->Fill((*recojet_subjet_E)[ind_sj1_rj2]/(*recojet_E)[ind_rj2],weight);
	   if(ind_sj2_rj2!=-1){
	     h_hist_vec_reco[260]->Fill(temp_sj1_rj2.Angle(temp_sj2_rj2.Vect())*TMath::RadToDeg(),weight);
	   }
	 }
	 if(ind_sj2_rj2!=-1){
	   if(usePartonInfo){
	     temp_sj.SetPxPyPzE((*recojet_subjet_Px)[ind_sj2_rj2],(*recojet_subjet_Py)[ind_sj2_rj2],(*recojet_subjet_Pz)[ind_sj2_rj2],(*recojet_subjet_E)[ind_sj2_rj2]);
	     if(temp_sj.Angle(tempZ_q_pos.Vect())<angle_min_sj_rj_Z_q_pos){
	       angle_min_sj_rj_Z_q_pos=temp_sj.Angle(tempZ_q_pos.Vect());
	       ind_Z_q_pos_sj_rj_angle=ind_sj2_rj2;
	     }
	     if(temp_sj.Angle(tempZ_q_neg.Vect())<angle_min_sj_rj_Z_q_neg){
	       angle_min_sj_rj_Z_q_neg=temp_sj.Angle(tempZ_q_neg.Vect());
	       ind_Z_q_neg_sj_rj_angle=ind_sj2_rj2;
	     }
	   }
	   h_hist_vec_reco[248]->Fill((*recojet_subjet_E)[ind_sj2_rj2],weight);
	 }
	 if(usePartonInfo && fill_partonInfo && (ind_Z_q_pos_sj_rj_angle!=-1 && (angle_min_sj_rj_Z_q_pos*TMath::RadToDeg())<15.0) && ( ind_Z_q_neg_sj_rj_angle!=-1 && (angle_min_sj_rj_Z_q_neg*TMath::RadToDeg())<15.0)){
	   h_hist_vec_2D[46]->Fill((*recojet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_pos_sj_rj_angle],(*recojet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_neg_sj_rj_angle],weight); 
	   h_hist_vec_2D[47]->Fill((*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_rj_angle],(*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_rj_angle],weight); 
	   h_hist_vec_2D[48]->Fill((*recojet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_pos_sj_rj_angle],(*recojet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_neg_sj_rj_angle],weight); 
	   if(((*recojet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_pos_sj_rj_angle])!=0 && ((*recojet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_neg_sj_rj_angle])!=0){
	     h_hist_vec_reco[424]->Fill((*recojet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_pos_sj_rj_angle]-(*recojet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_neg_sj_rj_angle],weight);
	   }
	   if(((*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_rj_angle])!=0 && ((*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_rj_angle])!=0){
	     h_hist_vec_reco[425]->Fill((*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_rj_angle]-(*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_rj_angle],weight);
	   }
	   if(((*recojet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_pos_sj_rj_angle])!=0 && ((*recojet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_neg_sj_rj_angle])!=0){
	     h_hist_vec_reco[426]->Fill((*recojet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_pos_sj_rj_angle]-(*recojet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_neg_sj_rj_angle],weight);
	   }
	   if(((*recojet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_pos_sj_rj_angle])!=0 && ((*recojet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_neg_sj_rj_angle])!=0){
	     h_hist_vec_reco[427]->Fill((*recojet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_pos_sj_rj_angle]-(*recojet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_neg_sj_rj_angle],weight);
	   }
	 }
	 if(usePartonInfo && fill_partonInfo && (ind_H_bbar_sj_rj_angle!=-1 && (angle_min_sj_rj_H_bbar*TMath::RadToDeg())<15.0) && ( ind_H_b_sj_rj_angle!=-1 && (angle_min_sj_rj_H_b*TMath::RadToDeg())<15.0)){
	   h_hist_vec_2D[54]->Fill((*recojet_subjet_jetChargeE_kappa_0_20)[ind_H_bbar_sj_rj_angle],(*recojet_subjet_jetChargeE_kappa_0_20)[ind_H_b_sj_rj_angle],weight); 
	   h_hist_vec_2D[55]->Fill((*recojet_subjet_jetChargeE_kappa_0_25)[ind_H_bbar_sj_rj_angle],(*recojet_subjet_jetChargeE_kappa_0_25)[ind_H_b_sj_rj_angle],weight); 
	   h_hist_vec_2D[56]->Fill((*recojet_subjet_jetChargeE_kappa_0_30)[ind_H_bbar_sj_rj_angle],(*recojet_subjet_jetChargeE_kappa_0_30)[ind_H_b_sj_rj_angle],weight); 
	   if(((*recojet_subjet_jetChargeE_kappa_0_20)[ind_H_bbar_sj_rj_angle])!=0 && ((*recojet_subjet_jetChargeE_kappa_0_20)[ind_H_b_sj_rj_angle])!=0){
	     h_hist_vec_reco[433]->Fill((*recojet_subjet_jetChargeE_kappa_0_20)[ind_H_bbar_sj_rj_angle]-(*recojet_subjet_jetChargeE_kappa_0_20)[ind_H_b_sj_rj_angle],weight);
	   }
	   if(((*recojet_subjet_jetChargeE_kappa_0_25)[ind_H_bbar_sj_rj_angle])!=0 && ((*recojet_subjet_jetChargeE_kappa_0_25)[ind_H_b_sj_rj_angle])!=0){
	     h_hist_vec_reco[434]->Fill((*recojet_subjet_jetChargeE_kappa_0_25)[ind_H_bbar_sj_rj_angle]-(*recojet_subjet_jetChargeE_kappa_0_25)[ind_H_b_sj_rj_angle],weight);
	   }
	   if(((*recojet_subjet_jetChargeE_kappa_0_30)[ind_H_bbar_sj_rj_angle])!=0 && ((*recojet_subjet_jetChargeE_kappa_0_30)[ind_H_b_sj_rj_angle])!=0){
	     h_hist_vec_reco[435]->Fill((*recojet_subjet_jetChargeE_kappa_0_30)[ind_H_bbar_sj_rj_angle]-(*recojet_subjet_jetChargeE_kappa_0_30)[ind_H_b_sj_rj_angle],weight);
	   }
	   if(((*recojet_subjet_jetChargeE_kappa_0_50)[ind_H_bbar_sj_rj_angle])!=0 && ((*recojet_subjet_jetChargeE_kappa_0_50)[ind_H_b_sj_rj_angle])!=0){
	     h_hist_vec_reco[436]->Fill((*recojet_subjet_jetChargeE_kappa_0_50)[ind_H_bbar_sj_rj_angle]-(*recojet_subjet_jetChargeE_kappa_0_50)[ind_H_b_sj_rj_angle],weight);
	   }
	 }

	 if(usePartonInfo && ind_Z_q_pos_sj_rj_angle!=-1 && (angle_min_sj_rj_Z_q_pos*TMath::RadToDeg())<15.0){
	   h_hist_vec_reco[359]->Fill((*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_pos_sj_rj_angle],weight);
	   h_hist_vec_reco[360]->Fill((*recojet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_pos_sj_rj_angle],weight);
	   h_hist_vec_reco[365]->Fill((*recojet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_pos_sj_rj_angle],weight);
	   h_hist_vec_reco[366]->Fill((*recojet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_pos_sj_rj_angle],weight);
	   h_hist_vec_reco[367]->Fill((*recojet_subjet_jetChargePt_kappa_0_25)[ind_Z_q_pos_sj_rj_angle],weight);
	   h_hist_vec_reco[373]->Fill((*recojet_subjet_jetChargePt_kappa_0_20)[ind_Z_q_pos_sj_rj_angle],weight);
	   h_hist_vec_reco[374]->Fill((*recojet_subjet_jetChargePt_kappa_0_30)[ind_Z_q_pos_sj_rj_angle],weight);
	 }
	 if(usePartonInfo && ind_Z_q_neg_sj_rj_angle!=-1 && (angle_min_sj_rj_Z_q_neg*TMath::RadToDeg())<15.0){
	   h_hist_vec_reco[375]->Fill((*recojet_subjet_jetChargeE_kappa_0_25)[ind_Z_q_neg_sj_rj_angle],weight);
	   h_hist_vec_reco[376]->Fill((*recojet_subjet_jetChargeE_kappa_0_50)[ind_Z_q_neg_sj_rj_angle],weight);
	   h_hist_vec_reco[381]->Fill((*recojet_subjet_jetChargeE_kappa_0_20)[ind_Z_q_neg_sj_rj_angle],weight);
	   h_hist_vec_reco[382]->Fill((*recojet_subjet_jetChargeE_kappa_0_30)[ind_Z_q_neg_sj_rj_angle],weight);
	   h_hist_vec_reco[383]->Fill((*recojet_subjet_jetChargePt_kappa_0_25)[ind_Z_q_neg_sj_rj_angle],weight);
	   h_hist_vec_reco[384]->Fill((*recojet_subjet_jetChargePt_kappa_0_50)[ind_Z_q_neg_sj_rj_angle],weight);
	   h_hist_vec_reco[389]->Fill((*recojet_subjet_jetChargePt_kappa_0_20)[ind_Z_q_neg_sj_rj_angle],weight);
	   h_hist_vec_reco[390]->Fill((*recojet_subjet_jetChargePt_kappa_0_30)[ind_Z_q_neg_sj_rj_angle],weight);
	 }
	 h_hist_vec_reco[284]->Fill(rj1_MET.M(),weight);
	 h_hist_vec_reco[287]->Fill(rj2_MET.M(),weight);
	 h_hist_vec_reco[290]->Fill(rj1_EMiss.M(),weight);
	 h_hist_vec_reco[293]->Fill(rj2_EMiss.M(),weight);
	 h_hist_vec_reco[315]->Fill(rj_m1_orig.E()+rj_m2_orig.E(),weight);
	 h_hist_vec_reco[316]->Fill(rj1_EMiss.E()+rj2_EMiss.E(),weight);
	 h_hist_vec_reco[317]->Fill((tempTotRecoP4-tempRecoIsoPhP4).E(),weight);
	 h_hist_vec_reco[318]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).E(),weight);
	 h_hist_vec_reco[319]->Fill(rj_m1_orig.E(),weight);
	 h_hist_vec_reco[320]->Fill(rj_m2_orig.E(),weight);
	 h_hist_vec_reco[321]->Fill(rj1_EMiss.E(),weight);
	 h_hist_vec_reco[322]->Fill(rj2_EMiss.E(),weight);
	 h_hist_vec_reco[323]->Fill(rj_m1_orig.E()-rj_m2_orig.E(),weight);
	 h_hist_vec_reco[324]->Fill(rj1_EMiss.E()-rj2_EMiss.E(),weight);
	 h_hist_vec_reco[325]->Fill((rj_m1_orig.E()-rj_m2_orig.E())/(rj_m1_orig.E()+rj_m2_orig.E()),weight);
	 h_hist_vec_reco[326]->Fill((rj1_EMiss.E()-rj2_EMiss.E())/(rj1_EMiss.E()+rj2_EMiss.E()),weight);
	 if(fill_genInfo){//reco met projection compared to true MET
	   h_hist_vec_reco[327]->Fill((tempRecoEMissCorrP4.Pt()-tempTotInvGenP4.Pt())/tempTotInvGenP4.Pt(),weight);
	   h_hist_vec_reco[328]->Fill((tempRecoMHMissCorrP4.Pt()-tempTotInvGenP4.Pt())/tempTotInvGenP4.Pt(),weight);
	   h_hist_vec_reco[329]->Fill(((-rj_m1_orig-rj_m2_orig).Pt()-tempTotInvGenP4.Pt())/tempTotInvGenP4.Pt(),weight);
	   h_hist_vec_reco[330]->Fill((tempRecoEMissCorrP4.E()-tempTotInvGenP4.E())/tempTotInvGenP4.E(),weight);
	   h_hist_vec_reco[331]->Fill((tempRecoMHMissCorrP4.E()-tempTotInvGenP4.E())/tempTotInvGenP4.E(),weight);
	   h_hist_vec_reco[332]->Fill(DeltaPhi(tempRecoEMissCorrP4.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	   h_hist_vec_reco[333]->Fill(DeltaPhi(tempRecoMHMissCorrP4.Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	   h_hist_vec_reco[334]->Fill(DeltaPhi((-rj_m1_orig-rj_m2_orig).Phi(),tempTotInvGenP4.Phi())*TMath::RadToDeg(),weight);
	   h_hist_vec_reco[335]->Fill(tempRecoEMissCorrP4.Angle(tempTotInvGenP4.Vect())*TMath::RadToDeg(),weight);
	   h_hist_vec_reco[336]->Fill(tempRecoMHMissCorrP4.Angle(tempTotInvGenP4.Vect())*TMath::RadToDeg(),weight);
	 }
	 //while the boost does change with energy scaling based on MET, the original direction of the vector does not
 	 h_hist_vec_reco[339]->Fill(cos(min(temp_sj1_rj1_boost_rj1_COM.Angle(rj_m1.Vect()),temp_sj2_rj1_boost_rj1_COM.Angle(rj_m1.Vect()))),weight);
	 h_hist_vec_reco[342]->Fill(cos(max(temp_sj1_rj1_boost_rj1_COM.Angle(rj_m1.Vect()),temp_sj2_rj1_boost_rj1_COM.Angle(rj_m1.Vect()))),weight);
	 h_hist_vec_reco[345]->Fill(cos(min(temp_sj1_rj2_boost_rj2_COM.Angle(rj_m2.Vect()),temp_sj2_rj2_boost_rj2_COM.Angle(rj_m2.Vect()))),weight);
	 h_hist_vec_reco[348]->Fill(cos(max(temp_sj1_rj2_boost_rj2_COM.Angle(rj_m2.Vect()),temp_sj2_rj2_boost_rj2_COM.Angle(rj_m2.Vect()))),weight);

	 h_hist_vec_reco[349]->Fill(tempRecoMETP4.Pt()/(tempTotRecoP4-tempRecoIsoPhP4).E());
	 h_hist_vec_reco[350]->Fill((-rj_m1_orig-rj_m2_orig).Pt()/(tempTotRecoP4-tempRecoIsoPhP4).E());
	 h_hist_vec_reco[351]->Fill(tempRecoMETP4.Pt()/(rj_m1_orig+rj_m2_orig).E());
	 h_hist_vec_reco[352]->Fill((-rj_m1_orig-rj_m2_orig).Pt()/(rj_m1_orig+rj_m2_orig).E());

	 h_hist_vec_reco[355]->Fill(rj_m1_orig.M()-rj_m2_orig.M(),weight);
	 h_hist_vec_reco[358]->Fill(rj1_EMiss.M()-rj2_EMiss.M(),weight);

       }//largest sqrtS bin
       if(usePartonInfo&& fill_partonInfo){
	 h_hist_vec_2D[19]->Fill(tempTotEventP4.M(),(rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).M(),weight);
	 h_hist_vec_2D[20]->Fill(tempTotEventP4.M(),(rj1_EMiss+rj2_EMiss-tempRecoIsoPhP4).M(),weight);
       }
       if(fill_genInfo){
	 h_hist_vec_2D[21]->Fill((gj_m1_orig+gj_m2_orig-tempGenIsoPhP4).M(),(rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).M(),weight);
	 h_hist_vec_2D[22]->Fill((gj1_EMiss+gj2_EMiss-tempRecoIsoPhP4).M(),(rj1_EMiss+rj2_EMiss-tempRecoIsoPhP4).M(),weight);
       }
      }//bracket closed on two recojets plus isolated lepton cut
      h_hist_vec_reco[313]->Fill((rj_m1_orig+rj_m2_orig).M(),weight);
      h_hist_vec_reco[314]->Fill((rj1_EMiss+rj2_EMiss).M(),weight);
      if(fill_partonInfo && tempTotEventP4.M()>sqrtS_high_reco){
	h_hist_parton[113]->Fill(((rj_m1_orig+rj_m2_orig-tempRecoIsoPhP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
	h_hist_parton[114]->Fill(((rj1_EMiss+rj2_EMiss-tempRecoIsoPhP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
	h_hist_parton[117]->Fill(((tempTotRecoP4-tempRecoIsoPhP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
	h_hist_parton[118]->Fill(((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
	h_hist_parton[121]->Fill(((rj1_MHMiss+rj2_MHMiss-tempRecoIsoPhP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
	h_hist_parton[123]->Fill(((tempTotRecoP4-tempRecoIsoPhP4+tempRecoMHMissCorrP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight);
       }
    }//recojets (==2_ are closed

   if((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M()<sqrtS_low){
     if((tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).M()<sqrtS_low){
       h_hist_vec_reco[295]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
     }else if((tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).M()<sqrtS_high_reco) {
       h_hist_vec_reco[296]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
     }else{
       h_hist_vec_reco[297]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
     }
     if(fill_partonInfo){
       if(tempTotEventP4.M()<sqrtS_low){
	 h_hist_vec_reco[298]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
       }else if(tempTotEventP4.M()<sqrtS_high_reco) {
	 h_hist_vec_reco[299]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
       }else{
	 h_hist_vec_reco[300]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
       }
     }
   }else if((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M()<sqrtS_high_reco){
     if((tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).M()<sqrtS_low){
       h_hist_vec_reco[301]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
     }else if((tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).M()<sqrtS_high_reco) {
       h_hist_vec_reco[302]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
     }else{
       h_hist_vec_reco[303]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
     }
     if(fill_partonInfo){
       if(tempTotEventP4.M()<sqrtS_low){
	 h_hist_vec_reco[304]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
       }else if(tempTotEventP4.M()<sqrtS_high_reco) {
	 h_hist_vec_reco[305]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
       }else{
	 h_hist_vec_reco[306]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
       }
     }
   }else{
     if((tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).M()<sqrtS_low){
       h_hist_vec_reco[307]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
     }else if((tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).M()<sqrtS_high_reco) {
       h_hist_vec_reco[308]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
     }else{
       h_hist_vec_reco[309]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
     }
     if(fill_partonInfo){
       if(tempTotEventP4.M()<sqrtS_low){
	 h_hist_vec_reco[310]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
       }else if(tempTotEventP4.M()<sqrtS_high_reco) {
	 h_hist_vec_reco[311]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
       }else{
	 h_hist_vec_reco[312]->Fill((tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
       }
     }
   }
   //gen info should always be the same order, parton info however is very process dependent
    if(sqrtS_eff_reco<sqrtS_low){
      //if(fill_genInfo){
      if(sqrtS_eff_gen<sqrtS_low){
	h_hist_vec_reco[261]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      }else if(sqrtS_eff_gen<sqrtS_high_reco) {
	h_hist_vec_reco[262]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      }else{
	h_hist_vec_reco[263]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      }
      //}
      if(fill_partonInfo){
	if(tempTotEventP4.M()<sqrtS_low){
	  h_hist_vec_reco[264]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else if(tempTotEventP4.M()<sqrtS_high_reco) {
	  h_hist_vec_reco[265]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}else{
	  h_hist_vec_reco[266]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
	}
      }
    }else if(sqrtS_eff_reco<sqrtS_high_reco){
      //if(fill_genInfo){
      if(sqrtS_eff_gen<sqrtS_low){
	h_hist_vec_reco[267]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      }else if(sqrtS_eff_gen<sqrtS_high_reco) {
	h_hist_vec_reco[268]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      }else{
	h_hist_vec_reco[269]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      }
      //}
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
      //if(fill_genInfo){
      if(sqrtS_eff_gen<sqrtS_low){
	h_hist_vec_reco[273]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      }else if(sqrtS_eff_gen<sqrtS_high_reco) {
	h_hist_vec_reco[274]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      }else{
	h_hist_vec_reco[275]->Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      }
      //}
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
      //h_hist_vec_2D[3]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenTruePhP4).M(),weight);
      //h_hist_vec_2D[4]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4).M(),weight);
      h_hist_vec_2D[5]->Fill(tempTotEventP4.M(),tempGenJetSum.M(),weight);    
      h_hist_vec_2D[6]->Fill(tempTotEventP4.M(),(tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      h_hist_vec_2D[7]->Fill(tempTotEventP4.M(),tempRecoJetSum.M(),weight);
      h_hist_vec_2D[8]->Fill((tempTotGenP4-tempGenIsoPhP4).M(),(tempTotRecoP4-tempRecoIsoPhP4).M(),weight);
      h_hist_vec_2D[9]->Fill(tempGenJetSum.M(),tempRecoJetSum.M(),weight);
      h_hist_vec_2D[10]->Fill(tempTotEventP4.M(),tempTotGenP4.M(),weight);
      h_hist_vec_2D[11]->Fill(tempTotEventP4.M(),tempTotRecoP4.M(),weight);
      h_hist_vec_2D[12]->Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4+tempTotInvGenP4).M(),weight);
      h_hist_vec_2D[13]->Fill(tempTotEventP4.M(),(tempTotGenP4+tempTotInvGenP4).M(),weight);
      h_hist_vec_2D[14]->Fill(tempTotEventP4.M(),(tempGenEMissCorrP4+tempTotGenP4-tempGenIsoPhP4).M(),weight);
      h_hist_vec_2D[15]->Fill(tempTotEventP4.M(),(tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
      h_hist_vec_2D[16]->Fill((tempGenEMissCorrP4+tempTotGenP4-tempGenIsoPhP4).M(), (tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M(),weight);
    }


    if(hist_vec_reco_1D_reco_vs_gen_selection.size()!=0 && usePartonInfo){
      //delta phi between gen and reco subjets
      if(ind_sj1_rj1!=-1 && ind_sj1_gj1!=-1){
	hist_vec_reco_1D_reco_vs_gen_selection[0]->Fill(DeltaPhi(temp_sj1_rj1.Phi(),temp_sj1_gj1.Phi()),weight);
	hist_vec_reco_1D_reco_vs_gen_selection[4]->Fill(temp_sj1_rj1.Theta()-temp_sj1_gj1.Theta(),weight);
	hist_vec_reco_1D_reco_vs_gen_selection[8]->Fill((temp_sj1_rj1.E()-temp_sj1_gj1.E())/temp_sj1_gj1.E(),weight);
      }
      if(ind_sj2_rj1!=-1 && ind_sj2_gj1!=-1){
	hist_vec_reco_1D_reco_vs_gen_selection[1]->Fill(DeltaPhi(temp_sj2_rj1.Phi(),temp_sj2_gj1.Phi()),weight);
	hist_vec_reco_1D_reco_vs_gen_selection[5]->Fill(temp_sj2_rj1.Theta()-temp_sj2_gj1.Theta(),weight);
	hist_vec_reco_1D_reco_vs_gen_selection[9]->Fill((temp_sj1_rj1.E()-temp_sj1_gj1.E())/temp_sj1_gj1.E(),weight);
      }
      if(ind_sj1_rj2!=-1 && ind_sj1_gj2!=-1){
	hist_vec_reco_1D_reco_vs_gen_selection[2]->Fill(DeltaPhi(temp_sj1_rj2.Phi(),temp_sj1_gj2.Phi()),weight);
	hist_vec_reco_1D_reco_vs_gen_selection[6]->Fill(temp_sj1_rj2.Theta()-temp_sj1_gj2.Theta(),weight);
	hist_vec_reco_1D_reco_vs_gen_selection[11]->Fill((temp_sj1_rj1.E()-temp_sj1_gj1.E())/temp_sj1_gj1.E(),weight);
      }
      if(ind_sj2_rj2!=-1 && ind_sj2_gj2!=-1){
	hist_vec_reco_1D_reco_vs_gen_selection[3]->Fill(DeltaPhi(temp_sj2_rj2.Phi(),temp_sj2_gj2.Phi()),weight);
	hist_vec_reco_1D_reco_vs_gen_selection[7]->Fill(temp_sj2_rj2.Theta()-temp_sj2_gj2.Theta(),weight);
	hist_vec_reco_1D_reco_vs_gen_selection[12]->Fill((temp_sj1_rj1.E()-temp_sj1_gj1.E())/temp_sj1_gj1.E(),weight);
      }
      //if(ind_sj1_rj1!=-1 && ind_sj1_gj1!=1 && ind_sj2_rj1!=-1 && ind_sj2_gj1!=1)){
      //hist_vec_reco_1D_reco_vs_gen_selection[8]->Fill(DeltaPhi(temp_sj1_rj1.Phi(),temp_sj1_gj1.Phi()),weight);
      //hist_vec_reco_1D_reco_vs_gen_selection[13]->Fill(temp_sj1_rj1.Theta()-temp_sj1_gj1.Theta(),weight);
      //}
    }
    if(hist_vec_reco_2D_reco_vs_gen_selection.size()!=0){

    }



  }//loop over tree

  std::cout<<"end of tree, b tagging survived by "<<((float)h_hist_vec_2D[0]->GetEntries()+(float)h_hist_vec_2D[1]->GetEntries()+(float)h_hist_vec_2D[2]->GetEntries())/(float)tree->GetEntries()<<std::endl;

  std::cout<<"counter_gen_met_over/totreco2000/rat_gen/reco_met_over/totreco2000/rat_reco "<<counter_j1j2_met_over_gen_2000<<"/"<<counter_j1j2_gen_2000<<"/"<<(float)counter_j1j2_met_over_gen_2000/(float)counter_j1j2_gen_2000<<"/"<<counter_j1j2_met_over_reco_2000<<"/"<<counter_j1j2_reco_2000<<"/"<<(float)counter_j1j2_met_over_reco_2000/(float)counter_j1j2_reco_2000<<" out of "<<tree->GetEntries()<<std::endl;


  std::cout<<"counter_gen_met_over/totrecolow_2000/rat_gen/reco_met_over/totrecolow_2000/rat_reco "<<counter_j1j2_met_over_gen_low_2000<<"/"<<counter_j1j2_gen_low_2000<<"/"<<(float)counter_j1j2_met_over_gen_low_2000/(float)counter_j1j2_gen_low_2000<<"/"<<counter_j1j2_met_over_reco_low_2000<<"/"<<counter_j1j2_reco_low_2000<<"/"<<(float)counter_j1j2_met_over_reco_low_2000/(float)counter_j1j2_reco_low_2000<<" out of "<<tree->GetEntries()<<std::endl;
}

  


void HZAnalyzerFull(){

  CLICdpStyle();

  gROOT->ProcessLine("#include <vector>");

  bool is_polp=false;


  const char* final_histo_name="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/HZqq_and_BG_1339x_noMassCuts_lepveto_190417_SignalAndBackground_DR7_sqrtS_j1_j2_EMissProj_noMHMiss_allHZEvts_polm_MassCuts.root";
  
  TFile* file_CLIC_HZqq = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/HZStudy_hz_qq_13391_polm80_3TeV_wO_CLIC_o3_v14_DR7.root");
  TFile* file_CLIC_tt  = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/HZStudy_qqqq_13394_to_13397_polm80_3TeV_wO_CLIC_o3_v14_DR7.root");
  TFile* file_CLIC_ee_qq  = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/HZStudy_ee_qq_13399_to_13402_polm80_3TeV_wO_CLIC_o3_v14_DR7.root");
  TFile* file_CLIC_qqqq  = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/HZStudy_qqqq_13394_to_13397_polm80_3TeV_wO_CLIC_o3_v14_DR7.root");
  if(is_polp){
    file_CLIC_HZqq = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/HZStudy_hz_qq_13392_polp80_3TeV_wO_CLIC_o3_v14_DR7.root");
    file_CLIC_tt  = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/HZStudy_qqqq_13393_polp80_3TeV_wO_CLIC_o3_v14_DR7.root");
    file_CLIC_ee_qq  = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/HZStudy_ee_qq_13398_polp80_3TeV_wO_CLIC_o3_v14_DR7.root");
    file_CLIC_qqqq  = TFile::Open("/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/HZStudy_qqqq_13393_polp80_3TeV_wO_CLIC_o3_v14_DR7.root");
  }

  //polm numbers
  double xsec_qqqq=902.;//+-2.02 fb, --lumi 0.1823fb-1)
  double xsec_ee_qq=1269.;//+-48.4 --lumi 0.03509 fb-1)
  double xsec_hz_qq=3.83;//+-0.0035 fb, --lumi 27.25 fb-1
  double xsec_tt=52.569;//+-0.031, lumi 0.9511 fb-1)

  //polp numbers
  if(is_polp){
    float lumi_factor=0.25;
    xsec_qqqq=120.*lumi_factor;//+-2.02 fb, --lumi 0.1823fb-1)
    xsec_ee_qq=786.*lumi_factor;//+-48.4 --lumi 0.03509 fb-1)
    xsec_hz_qq=2.67*lumi_factor;//+-0.0035 fb, --lumi 27.25 fb-1
    xsec_tt=52.569*lumi_factor;//+-0.031, lumi 0.9511 fb-1)
  }
  bool usePartonInfo = true;
  bool fillPartonInfo = true;
  bool fillGenInfo = true;


  TFile* file_histogram=new TFile(final_histo_name,"recreate");

  int n_bins_high=200;
  int n_bins_low=50;

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

  TH1F* h_H_b_HelicityAngle_sqrt_s_0_750 = new TH1F("h_H_b_HelicityAngle_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_H_b_HelicityAngle_sqrt_s_750_2500 = new TH1F("h_H_b_HelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_H_b_HelicityAngle_sqrt_s_2500 = new TH1F("h_H_b_HelicityAngle_sqrt_s_2500","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_H_Z_q_pos_HelicityAngle_sqrt_s_0_750 = new TH1F("h_H_Z_q_pos_HelicityAngle_sqrt_s_0_750","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_H_Z_q_pos_HelicityAngle_sqrt_s_750_2500 = new TH1F("h_H_Z_q_pos_HelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
  TH1F* h_H_Z_q_pos_HelicityAngle_sqrt_s_2500 = new TH1F("h_H_Z_q_pos_HelicityAngle_sqrt_s_2500","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

  double lim_cosdalpha_low=-1.;
  double lim_cosdalpha_high=1.;

  TH1F* h_H_b_CosHelicityAngle_sqrt_s_0_750 = new TH1F("h_H_b_CosHelicityAngle_sqrt_s_0_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_H_b_CosHelicityAngle_sqrt_s_750_2500 = new TH1F("h_H_b_CosHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_H_b_CosHelicityAngle_sqrt_s_2500 = new TH1F("h_H_b_CosHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_H_Z_q_pos_CosHelicityAngle_sqrt_s_0_750 = new TH1F("h_H_Z_q_pos_CosHelicityAngle_sqrt_s_0_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_H_Z_q_pos_CosHelicityAngle_sqrt_s_750_2500 = new TH1F("h_H_Z_q_pos_CosHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_H_Z_q_pos_CosHelicityAngle_sqrt_s_2500 = new TH1F("h_H_Z_q_pos_CosHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);


  TH1F* h_Z_q_CosHelicityAngle_sqrt_s_0_750 = new TH1F("h_Z_q_CosHelicityAngle_sqrt_s_0_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_Z_q_CosHelicityAngle_sqrt_s_750_2500 = new TH1F("h_Z_q_CosHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_Z_q_CosHelicityAngle_sqrt_s_2500 = new TH1F("h_Z_q_CosHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_Z_qbar_CosHelicityAngle_sqrt_s_0_750 = new TH1F("h_Z_qbar_CosHelicityAngle_sqrt_s_0_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_Z_qbar_CosHelicityAngle_sqrt_s_750_2500 = new TH1F("h_Z_qbar_CosHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_Z_qbar_CosHelicityAngle_sqrt_s_2500 = new TH1F("h_Z_qbar_CosHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  
  TH1F* h_H_bZ_q_pos_CosMinHelicityAngle_sqrt_s_0_750 = new TH1F("h_H_bZ_q_pos_CosMinHelicityAngle_sqrt_s_0_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_H_bZ_q_pos_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_H_bZ_q_pos_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_H_bZ_q_pos_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_H_bZ_q_pos_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_H_bZ_q_pos_CosMaxHelicityAngle_sqrt_s_0_750 = new TH1F("h_H_bZ_q_pos_CosMaxHelicityAngle_sqrt_s_0_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_H_bZ_q_pos_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_H_bZ_q_pos_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_H_bZ_q_pos_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_H_bZ_q_pos_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);

  TH1F* h_Z_qqbar_CosMinHelicityAngle_sqrt_s_0_750 = new TH1F("h_Z_qqbar_CosMinHelicityAngle_sqrt_s_0_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_Z_qqbar_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_Z_qqbar_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_Z_qqbar_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_Z_qqbar_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_Z_qqbar_CosMaxHelicityAngle_sqrt_s_0_750 = new TH1F("h_Z_qqbar_CosMaxHelicityAngle_sqrt_s_0_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_Z_qqbar_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_Z_qqbar_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  TH1F* h_Z_qqbar_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_Z_qqbar_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);

  TH1F* h_dPhi_MET_mass1_gj_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_gj_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass1_rj_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_rj_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass2_gj_sqrtS_2500 = new TH1F("h_dPhi_MET_mass2_gj_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass2_rj_sqrtS_2500 = new TH1F("h_dPhi_MET_mass2_rj_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  //dalpha j1/H < 25 degrees
  TH1F* h_dPhi_MET_mass1_gj_matched_H_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_gj_matched_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass1_rj_matched_H_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_rj_matched_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
 //dalpha j1/H > 155 degrees
  TH1F* h_dPhi_MET_mass1_gj_unmatched_H_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_gj_unmatched_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass1_rj_unmatched_H_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_rj_unmatched_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  //dalpha j2/Z < 25 degrees
  TH1F* h_dPhi_MET_mass2_gj_matched_Z_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_gj_matched_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass2_rj_matched_Z_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_rj_matched_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
 //dalpha j2/Z > 155 degrees
  TH1F* h_dPhi_MET_mass2_gj_unmatched_Z_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_gj_unmatched_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_mass2_rj_unmatched_Z_sqrtS_2500 = new TH1F("h_dPhi_MET_mass1_rj_unmatched_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  //dalpha j/H < 25 degrees
  TH1F* h_dPhi_MET_gj_matched_H_sqrtS_2500 = new TH1F("h_dPhi_MET_gj_matched_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_rj_matched_H_sqrtS_2500 = new TH1F("h_dPhi_MET_rj_matched_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  //dalpha j/Z < 25 degrees
  TH1F* h_dPhi_MET_gj_matched_Z_sqrtS_2500 = new TH1F("h_dPhi_MET_gj_matched_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dPhi_MET_rj_matched_Z_sqrtS_2500 = new TH1F("h_dPhi_MET_rj_matched_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_dAlpha_mass1_gj_EMiss_H_sqrtS_2500 = new TH1F("h_dAlpha_mass1_gj_EMiss_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dAlpha_mass2_gj_EMiss_Z_sqrtS_2500 = new TH1F("h_dAlpha_mass2_gj_EMiss_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  
  TH1F* h_dAlpha_mass1_rj_EMiss_H_sqrtS_2500 = new TH1F("h_dAlpha_mass1_rj_EMiss_H_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);
  TH1F* h_dAlpha_mass2_rj_EMiss_Z_sqrtS_2500 = new TH1F("h_dAlpha_mass2_rj_EMiss_Z_sqrtS_2500","", n_bins_high, lim_phi_low,lim_phi_high);

  TH1F* h_sqrtS_E_tot_H_Z = new TH1F("h_sqrtS_E_tot_H_Z","", n_bins_high, lim_energy_low,lim_energy_high);

  double lim_energy_jet_low=800;
  double lim_energy_jet_high=1700.;
  TH1F* h_E_H_sqrtS_2500 = new TH1F("h_sqrtS_E_H_sqrtS_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
  TH1F* h_E_Z_sqrtS_2500 = new TH1F("h_sqrtS_E_Z_sqrtS_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);

  double lim_delta_energy_jet_low=-200.;
  double lim_delta_energy_jet_high=200.;
  TH1F* h_E_H_min_Z_sqrtS_2500 = new TH1F("h_E_H_min_Z_sqrtS_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);

  double lim_delta_energy_jet_rel_low=-0.35;
  double lim_delta_energy_jet_rel_high=0.35;
  TH1F* h_E_H_min_Z_over_E_tot_H_Z_sqrtS_2500 = new TH1F("h_E_H_min_Z_over_E_tot_H_Z_sqrtS_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);


  double lim_delta_sqrtS_rel_low=-0.30;
  double lim_delta_sqrtS_rel_high=0.30;
  TH1F* h_HZ_delta_sqrtS_reco_j1_j2_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_reco_j1_j2_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_reco_j1_j2_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{reco}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");
  TH1F* h_HZ_delta_sqrtS_reco_j1_j2_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_reco_j1_j2_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_reco_j1_j2_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{reco}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");

  TH1F* h_HZ_delta_sqrtS_gen_j1_j2_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_gen_j1_j2_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_gen_j1_j2_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{gen}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");
  TH1F* h_HZ_delta_sqrtS_gen_j1_j2_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_gen_j1_j2_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_gen_j1_j2_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{gen}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");

  TH1F* h_HZ_delta_sqrtS_reco_Etot_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_reco_Etot_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_reco_Etot_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{reco}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");
  TH1F* h_HZ_delta_sqrtS_reco_Etot_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_reco_Etot_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_reco_Etot_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{reco}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");

  TH1F* h_HZ_delta_sqrtS_gen_Etot_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_gen_Etot_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_gen_Etot_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{gen}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");
  TH1F* h_HZ_delta_sqrtS_gen_Etot_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_gen_Etot_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_gen_Etot_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{gen}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");

  TH1F* h_HZ_delta_sqrtS_reco_j1_j2_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_reco_j1_j2_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_reco_j1_j2_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{reco}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");
  TH1F* h_HZ_delta_sqrtS_gen_j1_j2_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_gen_j1_j2_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_gen_j1_j2_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{gen}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");
  TH1F* h_HZ_delta_sqrtS_reco_Etot_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_reco_Etot_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_reco_Etot_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{reco}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");
  TH1F* h_HZ_delta_sqrtS_gen_Etot_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500 = new TH1F("h_HZ_delta_sqrtS_gen_Etot_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500","", n_bins_high, lim_delta_sqrtS_rel_low,lim_delta_sqrtS_rel_high);
  h_HZ_delta_sqrtS_gen_Etot_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500->GetXaxis()->SetTitle("(#sqrt{s_{eff}^{gen}}-#sqrt{s_{eff}^{parton}})/#sqrt{s_{eff}^{parton}}");

  double lim_met_min_rel=0;
  double lim_met_high_rel=0.30;
  //in the following what we call met for this particular histograms refers to met calculated by the dijet sum minus the isolated photon sum
  TH1F* h_HZ_reco_met_over_E_tot_sqrt_part_s_0_750 = new TH1F("h_HZ_reco_met_over_E_tot_sqrt_part_s_0_750","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_reco_met_over_E_tot_sqrt_part_s_0_750->GetXaxis()->SetTitle("missing p_{T}/E_{tot}");
  TH1F* h_HZ_reco_met_over_E_tot_sqrt_part_s_750_2500 = new TH1F("h_HZ_reco_met_over_E_tot_sqrt_part_s_750_2500","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_reco_met_over_E_tot_sqrt_part_s_750_2500->GetXaxis()->SetTitle("missing p_{T}/E_{tot}");
  TH1F* h_HZ_reco_met_over_E_tot_sqrt_part_s_2500 = new TH1F("h_HZ_reco_met_over_E_tot_sqrt_part_s_2500","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_reco_met_over_E_tot_sqrt_part_s_2500 ->GetXaxis()->SetTitle("missing p_{T}/E_{tot}");

  TH1F* h_HZ_gen_met_over_E_tot_sqrt_part_s_0_750 = new TH1F("h_HZ_gen_met_over_E_tot_sqrt_part_s_0_750","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_gen_met_over_E_tot_sqrt_part_s_0_750->GetXaxis()->SetTitle("missing p_{T}/E_{tot}");
  TH1F* h_HZ_gen_met_over_E_tot_sqrt_part_s_750_2500 = new TH1F("h_HZ_gen_met_over_E_tot_sqrt_part_s_750_2500","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_gen_met_over_E_tot_sqrt_part_s_750_2500->GetXaxis()->SetTitle("missing p_{T}/E_{tot}");
  TH1F* h_HZ_gen_met_over_E_tot_sqrt_part_s_2500 = new TH1F("h_HZ_gen_met_over_E_tot_sqrt_part_s_2500","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_gen_met_over_E_tot_sqrt_part_s_2500 ->GetXaxis()->SetTitle("missing p_{T}/E_{tot}");

  TH1F* h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_0_750 = new TH1F("h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_0_750","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_0_750->GetXaxis()->SetTitle("missing p_{T}/(E_{j1}+E_{j2})");
  TH1F* h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_750_2500 = new TH1F("h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_750_2500","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_750_2500->GetXaxis()->SetTitle("missing p_{T}/(E_{j1}+E_{j2})");
  TH1F* h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_2500 = new TH1F("h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_2500","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_2500 ->GetXaxis()->SetTitle("missing p_{T}/(E_{j1}+E_{j2})");

  TH1F* h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_0_750 = new TH1F("h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_0_750","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_0_750->GetXaxis()->SetTitle("missing p_{T}/(E_{j1}+E_{j2})");
  TH1F* h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_750_2500 = new TH1F("h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_750_2500","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_750_2500->GetXaxis()->SetTitle("missing p_{T}/(E_{j1}+E_{j2})");
  TH1F* h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_2500 = new TH1F("h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_2500","", n_bins_high, lim_met_min_rel,lim_met_high_rel);
  h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_2500 ->GetXaxis()->SetTitle("missing p_{T}/(E_{j1}+E_{j2})");


  double lim_met_min=0;
  double lim_met_high=600;
  TH1F* h_HZ_reco_met_sqrt_part_s_0_750 = new TH1F("h_HZ_reco_met_sqrt_part_s_0_750","", n_bins_high, lim_met_min,lim_met_high);
  h_HZ_reco_met_sqrt_part_s_0_750->GetXaxis()->SetTitle("missing p_{T} [GeV]");
  TH1F* h_HZ_reco_met_sqrt_part_s_750_2500 = new TH1F("h_HZ_reco_met_sqrt_part_s_750_2500","", n_bins_high, lim_met_min,lim_met_high);
  h_HZ_reco_met_sqrt_part_s_750_2500->GetXaxis()->SetTitle("missing p_{T} [GeV]");
  TH1F* h_HZ_reco_met_sqrt_part_s_2500 = new TH1F("h_HZ_reco_met_sqrt_part_s_2500","", n_bins_high, lim_met_min,lim_met_high);
  h_HZ_reco_met_sqrt_part_s_2500 ->GetXaxis()->SetTitle("missing p_{T} [GeV]");

  TH1F* h_HZ_gen_met_sqrt_part_s_0_750 = new TH1F("h_HZ_gen_met_sqrt_part_s_0_750","", n_bins_high, lim_met_min,lim_met_high);
  h_HZ_gen_met_sqrt_part_s_0_750->GetXaxis()->SetTitle("missing p_{T} [GeV]");
  TH1F* h_HZ_gen_met_sqrt_part_s_750_2500 = new TH1F("h_HZ_gen_met_sqrt_part_s_750_2500","", n_bins_high, lim_met_min,lim_met_high);
  h_HZ_gen_met_sqrt_part_s_750_2500->GetXaxis()->SetTitle("missing p_{T} [GeV]");
  TH1F* h_HZ_gen_met_sqrt_part_s_2500 = new TH1F("h_HZ_gen_met_sqrt_part_s_2500","", n_bins_high, lim_met_min,lim_met_high);
  h_HZ_gen_met_sqrt_part_s_2500 ->GetXaxis()->SetTitle("missing p_{T} [GeV]");


  TH1F* h_HZ_reco_E_tot_sqrt_part_s_0_750 = new TH1F("h_HZ_reco_E_tot_sqrt_part_s_0_750","", n_bins_high, 0,1000);
  h_HZ_reco_E_tot_sqrt_part_s_0_750->GetXaxis()->SetTitle("E_{tot} [GeV]");
  TH1F* h_HZ_reco_E_tot_sqrt_part_s_750_2500 = new TH1F("h_HZ_reco_E_tot_sqrt_part_s_750_2500","", n_bins_high,500,3000);
  h_HZ_reco_E_tot_sqrt_part_s_750_2500->GetXaxis()->SetTitle("E_{tot} [GeV]");
  TH1F* h_HZ_reco_E_tot_sqrt_part_s_2500 = new TH1F("h_HZ_reco_E_tot_sqrt_part_s_2500","", n_bins_high, 2000,4000);
  h_HZ_reco_E_tot_sqrt_part_s_2500 ->GetXaxis()->SetTitle("E_{tot} [GeV]");

  TH1F* h_HZ_gen_E_tot_sqrt_part_s_0_750 = new TH1F("h_HZ_gen_E_tot_sqrt_part_s_0_750","", n_bins_high, 0,1000);
  h_HZ_gen_E_tot_sqrt_part_s_0_750->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
  TH1F* h_HZ_gen_E_tot_sqrt_part_s_750_2500 = new TH1F("h_HZ_gen_E_tot_sqrt_part_s_750_2500","", n_bins_high, 500,3000);
  h_HZ_gen_E_tot_sqrt_part_s_750_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
  TH1F* h_HZ_gen_E_tot_sqrt_part_s_2500 = new TH1F("h_HZ_gen_E_tot_sqrt_part_s_2500","", n_bins_high, 2000,4000);
  h_HZ_gen_E_tot_sqrt_part_s_2500 ->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");

  TH1F* h_HZ_reco_E_j1_j2_sqrt_part_s_0_750 = new TH1F("h_HZ_reco_E_j1_j2_sqrt_part_s_0_750","", n_bins_high, 0,1000);
  h_HZ_reco_E_j1_j2_sqrt_part_s_0_750->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
  TH1F* h_HZ_reco_E_j1_j2_sqrt_part_s_750_2500 = new TH1F("h_HZ_reco_E_j1_j2_sqrt_part_s_750_2500","", n_bins_high,500,3000);
  h_HZ_reco_E_j1_j2_sqrt_part_s_750_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
  TH1F* h_HZ_reco_E_j1_j2_sqrt_part_s_2500 = new TH1F("h_HZ_reco_E_j1_j2_sqrt_part_s_2500","", n_bins_high, 2000,4000);
  h_HZ_reco_E_j1_j2_sqrt_part_s_2500 ->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");

  TH1F* h_HZ_gen_E_j1_j2_sqrt_part_s_0_750 = new TH1F("h_HZ_gen_E_j1_j2_sqrt_part_s_0_750","", n_bins_high, 0,1000);
  h_HZ_gen_E_j1_j2_sqrt_part_s_0_750->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
  TH1F* h_HZ_gen_E_j1_j2_sqrt_part_s_750_2500 = new TH1F("h_HZ_gen_E_j1_j2_sqrt_part_s_750_2500","", n_bins_high, 500,3000);
  h_HZ_gen_E_j1_j2_sqrt_part_s_750_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
  TH1F* h_HZ_gen_E_j1_j2_sqrt_part_s_2500 = new TH1F("h_HZ_gen_E_j1_j2_sqrt_part_s_2500","", n_bins_high, 2000,4000);
  h_HZ_gen_E_j1_j2_sqrt_part_s_2500 ->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");


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
  hist_vec_HZ_parton.push_back(h_H_b_HelicityAngle_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_b_HelicityAngle_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_H_b_HelicityAngle_sqrt_s_2500);//60
  hist_vec_HZ_parton.push_back(h_H_Z_q_pos_HelicityAngle_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_Z_q_pos_HelicityAngle_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_H_Z_q_pos_HelicityAngle_sqrt_s_2500);//63
  hist_vec_HZ_parton.push_back(h_H_b_CosHelicityAngle_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_b_CosHelicityAngle_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_H_b_CosHelicityAngle_sqrt_s_2500);//66
  hist_vec_HZ_parton.push_back(h_H_Z_q_pos_CosHelicityAngle_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_Z_q_pos_CosHelicityAngle_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_H_Z_q_pos_CosHelicityAngle_sqrt_s_2500);//69
  hist_vec_HZ_parton.push_back(h_H_bZ_q_pos_CosMinHelicityAngle_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_bZ_q_pos_CosMinHelicityAngle_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_H_bZ_q_pos_CosMinHelicityAngle_sqrt_s_2500);
  hist_vec_HZ_parton.push_back(h_H_bZ_q_pos_CosMaxHelicityAngle_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_H_bZ_q_pos_CosMaxHelicityAngle_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_H_bZ_q_pos_CosMaxHelicityAngle_sqrt_s_2500);//75
  hist_vec_HZ_parton.push_back(h_Z_q_CosHelicityAngle_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_Z_q_CosHelicityAngle_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_Z_q_CosHelicityAngle_sqrt_s_2500);//78
  hist_vec_HZ_parton.push_back(h_Z_qbar_CosHelicityAngle_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_Z_qbar_CosHelicityAngle_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_Z_qbar_CosHelicityAngle_sqrt_s_2500);//81
  hist_vec_HZ_parton.push_back(h_Z_qqbar_CosMinHelicityAngle_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_Z_qqbar_CosMinHelicityAngle_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_Z_qqbar_CosMinHelicityAngle_sqrt_s_2500);
  hist_vec_HZ_parton.push_back(h_Z_qqbar_CosMaxHelicityAngle_sqrt_s_0_750);
  hist_vec_HZ_parton.push_back(h_Z_qqbar_CosMaxHelicityAngle_sqrt_s_750_2500);
  hist_vec_HZ_parton.push_back(h_Z_qqbar_CosMaxHelicityAngle_sqrt_s_2500);//87
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass1_gj_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass1_rj_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass2_gj_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass2_rj_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass1_gj_matched_H_sqrtS_2500);//92 //dalpha j1/H < 25 degrees
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass1_rj_matched_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass1_gj_unmatched_H_sqrtS_2500); //dalpha j1/H > 155 degrees
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass1_rj_unmatched_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass2_gj_matched_Z_sqrtS_2500);  //dalpha j2/Z < 25 degrees
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass2_rj_matched_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass2_gj_unmatched_Z_sqrtS_2500); //dalpha j2/Z > 155 degrees
  hist_vec_HZ_parton.push_back(h_dPhi_MET_mass2_rj_unmatched_Z_sqrtS_2500);//99
  hist_vec_HZ_parton.push_back(h_dPhi_MET_gj_matched_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_rj_matched_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_gj_matched_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dPhi_MET_rj_matched_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass1_gj_EMiss_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass2_gj_EMiss_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass1_rj_EMiss_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_dAlpha_mass2_rj_EMiss_Z_sqrtS_2500);//107
  hist_vec_HZ_parton.push_back(h_sqrtS_E_tot_H_Z); 
  hist_vec_HZ_parton.push_back(h_E_H_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_E_Z_sqrtS_2500); 
  hist_vec_HZ_parton.push_back(h_E_H_min_Z_sqrtS_2500);
  hist_vec_HZ_parton.push_back(h_E_H_min_Z_over_E_tot_H_Z_sqrtS_2500);//112
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_reco_j1_j2_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_reco_j1_j2_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_gen_j1_j2_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);//115
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_gen_j1_j2_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_reco_Etot_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_reco_Etot_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);//118
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_gen_Etot_isoPh_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_gen_Etot_isoPh_EMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);//120

  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_reco_j1_j2_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_gen_j1_j2_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_reco_Etot_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_delta_sqrtS_gen_Etot_isoPh_MHMiss_min_sqrtS_e1_e2_over_sqrtS_e1_e2_sqrt_s_2500); //124

  hist_vec_HZ_parton.push_back(h_HZ_reco_met_over_E_tot_sqrt_part_s_0_750);//125
  hist_vec_HZ_parton.push_back(h_HZ_reco_met_over_E_tot_sqrt_part_s_750_2500);
  hist_vec_HZ_parton.push_back(h_HZ_reco_met_over_E_tot_sqrt_part_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_gen_met_over_E_tot_sqrt_part_s_0_750);
  hist_vec_HZ_parton.push_back(h_HZ_gen_met_over_E_tot_sqrt_part_s_750_2500);
  hist_vec_HZ_parton.push_back(h_HZ_gen_met_over_E_tot_sqrt_part_s_2500);//130

  hist_vec_HZ_parton.push_back(h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_0_750);
  hist_vec_HZ_parton.push_back(h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_750_2500);
  hist_vec_HZ_parton.push_back(h_HZ_reco_met_over_E_j1_j2_sqrt_part_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_0_750);
  hist_vec_HZ_parton.push_back(h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_750_2500);
  hist_vec_HZ_parton.push_back(h_HZ_gen_met_over_E_j1_j2_sqrt_part_s_2500);//136

  hist_vec_HZ_parton.push_back(h_HZ_reco_met_sqrt_part_s_0_750);
  hist_vec_HZ_parton.push_back(h_HZ_reco_met_sqrt_part_s_750_2500);
  hist_vec_HZ_parton.push_back(h_HZ_reco_met_sqrt_part_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_gen_met_sqrt_part_s_0_750);
  hist_vec_HZ_parton.push_back(h_HZ_gen_met_sqrt_part_s_750_2500);
  hist_vec_HZ_parton.push_back(h_HZ_gen_met_sqrt_part_s_2500);//142

  hist_vec_HZ_parton.push_back(h_HZ_reco_E_tot_sqrt_part_s_0_750);
  hist_vec_HZ_parton.push_back(h_HZ_reco_E_tot_sqrt_part_s_750_2500);
  hist_vec_HZ_parton.push_back(h_HZ_reco_E_tot_sqrt_part_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_gen_E_tot_sqrt_part_s_0_750);
  hist_vec_HZ_parton.push_back(h_HZ_gen_E_tot_sqrt_part_s_750_2500);
  hist_vec_HZ_parton.push_back(h_HZ_gen_E_tot_sqrt_part_s_2500);//148

  hist_vec_HZ_parton.push_back(h_HZ_reco_E_j1_j2_sqrt_part_s_0_750);
  hist_vec_HZ_parton.push_back(h_HZ_reco_E_j1_j2_sqrt_part_s_750_2500);
  hist_vec_HZ_parton.push_back(h_HZ_reco_E_j1_j2_sqrt_part_s_2500);
  hist_vec_HZ_parton.push_back(h_HZ_gen_E_j1_j2_sqrt_part_s_0_750);
  hist_vec_HZ_parton.push_back(h_HZ_gen_E_j1_j2_sqrt_part_s_750_2500);
  hist_vec_HZ_parton.push_back(h_HZ_gen_E_j1_j2_sqrt_part_s_2500);//154


  
  for(unsigned int i=0;i<hist_vec_HZ_parton.size();i++){
    hist_vec_HZ_parton[i]->Sumw2();
    hist_vec_HZ_parton[i]->SetLineWidth(2);
    hist_vec_HZ_parton[i]->GetYaxis()->SetTitle("Events");
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


 TH1F* h_HZ_sqrtS_gen_isoPh = new TH1F("h_HZ_sqrtS_gen_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);//default
 TH1F* h_HZ_sqrtS_gen_isoPh_inv = new TH1F("h_HZ_sqrtS_gen_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
 TH1F* h_HZ_sqrtS_gen = new TH1F("h_HZ_sqrtS_gen","", n_bins_high, lim_energy_low,lim_energy_high);
 TH1F* h_HZ_sqrtS_gen_inv = new TH1F("h_HZ_sqrtS_gen_inv","", n_bins_high, lim_energy_low,lim_energy_high);
 TH1F* h_HZ_sqrtS_gen_truePh_inv = new TH1F("h_HZ_sqrtS_gen_truePh_inv","", n_bins_high, lim_energy_low,lim_energy_high);

 h_HZ_sqrtS_gen_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_gen_isoPh_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_gen->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_gen_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_gen_truePh_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_HZ_mass_j1_METProj_gen_sqrt_s_0_750 = new TH1F("h_HZ_mass_j1_METProj_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_METProj_gen_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j1_METProj_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_METProj_gen_sqrt_s_2500 = new TH1F("h_HZ_mass_j1_METProj_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 
 TH1F* h_HZ_mass_j2_METProj_gen_sqrt_s_0_750 = new TH1F("h_HZ_mass_j2_METProj_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_METProj_gen_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j2_METProj_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_METProj_gen_sqrt_s_2500 = new TH1F("h_HZ_mass_j2_METProj_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);

 TH1F* h_HZ_mass_j1_EMissProj_gen_sqrt_s_0_750 = new TH1F("h_HZ_mass_j1_EMissProj_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_EMissProj_gen_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j1_EMissProj_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_EMissProj_gen_sqrt_s_2500 = new TH1F("h_HZ_mass_j1_EMissProj_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 
 TH1F* h_HZ_mass_j2_EMissProj_gen_sqrt_s_0_750 = new TH1F("h_HZ_mass_j2_EMissProj_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_EMissProj_gen_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j2_EMissProj_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_EMissProj_gen_sqrt_s_2500 = new TH1F("h_HZ_mass_j2_EMissProj_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);

 TH1F* h_HZ_sqrtS_gen_isoPh_EMissCorr = new TH1F("h_HZ_sqrtS_gen_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_sqrtS_gen_isoPh_EMissCorr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 TH1F* h_HZ_sqrtS_gen_j1_j2_isoPh = new TH1F("h_HZ_sqrtS_gen_j1_j2_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_sqrtS_gen_j1_j2_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 TH1F* h_HZ_sqrtS_gen_j1_j2_isoPh_EMissCorr = new TH1F("h_HZ_sqrtS_gen_j1_j2_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_sqrtS_gen_j1_j2_isoPh_EMissCorr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 
 double lim_jet_energy_low_sqrt_s_2500=2000;
 TH1F* h_HZ_E_tot_j1_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_E_tot_j1_j2_gen_sqrt_s_2500","", n_bins_high, lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_HZ_E_tot_j1_j2_gen_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
 TH1F* h_HZ_E_tot_j1_j2_EMiss_gen_sqrt_s_2500 = new TH1F("h_HZ_E_tot_j1_j2_EMiss_gen_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_HZ_E_tot_j1_j2_EMiss_gen_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");

 TH1F* h_HZ_E_tot_isoPh_gen_sqrt_s_2500 = new TH1F("h_HZ_E_tot_isoPh_gen_sqrt_s_2500","", n_bins_high, lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_HZ_E_tot_isoPh_gen_sqrt_s_2500->GetXaxis()->SetTitle("E_{tot} [GeV]");
 TH1F* h_HZ_E_tot_isoPh_EMiss_gen_sqrt_s_2500 = new TH1F("h_HZ_E_tot_isoPh_EMiss_gen_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_HZ_E_tot_isoPh_EMiss_gen_sqrt_s_2500->GetXaxis()->SetTitle("(E_{tot}+E_{miss}^{jetProj}) [GeV]");
 
 TH1F* h_HZ_E_j1_gen_sqrt_s_2500 = new TH1F("h_HZ_E_j1_gen_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_HZ_E_j1_gen_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_HZ_E_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_E_j2_gen_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_HZ_E_j2_gen_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");
 TH1F* h_HZ_E_j1_gen_EMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j1_gen_EMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_HZ_E_j1_gen_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_HZ_E_j2_gen_EMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j2_gen_EMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_HZ_E_j2_gen_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");

 TH1F* h_HZ_E_j1_min_E_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_E_j1_min_E_j2_gen_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);
 h_HZ_E_j1_min_E_j2_gen_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2}) [GeV]");
 TH1F* h_HZ_E_j1_min_E_j2_gen_EMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j1_min_E_j2_gen_EMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);
 h_HZ_E_j1_min_E_j2_gen_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2}) [GeV]");

 TH1F* h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_sqrt_s_2500 = new TH1F("h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");
 TH1F* h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_EMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_EMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");


 TH1F* h_HZ_mass_j1_MHMissProj_gen_sqrt_s_0_750 = new TH1F("h_HZ_mass_j1_MHMissProj_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_MHMissProj_gen_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j1_MHMissProj_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_MHMissProj_gen_sqrt_s_2500 = new TH1F("h_HZ_mass_j1_MHMissProj_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_MHMissProj_gen_sqrt_s_0_750 = new TH1F("h_HZ_mass_j2_MHMissProj_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_MHMissProj_gen_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j2_MHMissProj_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_MHMissProj_gen_sqrt_s_2500 = new TH1F("h_HZ_mass_j2_MHMissProj_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_E_tot_j1_j2_MHMiss_gen_sqrt_s_2500 = new TH1F("h_HZ_E_tot_j1_j2_MHMiss_gen_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_HZ_E_tot_j1_j2_MHMiss_gen_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
 TH1F* h_HZ_E_tot_isoPh_MHMiss_gen_sqrt_s_2500 = new TH1F("h_HZ_E_tot_isoPh_MHMiss_gen_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_HZ_E_tot_isoPh_MHMiss_gen_sqrt_s_2500->GetXaxis()->SetTitle("(E_{tot}+E_{miss}^{jetProj}) [GeV]");
 TH1F* h_HZ_E_j1_gen_MHMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j1_gen_MHMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_HZ_E_j1_gen_MHMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_HZ_E_j2_gen_MHMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j2_gen_MHMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_HZ_E_j2_gen_MHMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");
 TH1F* h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_MHMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_MHMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_MHMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");

 float lim_delta_energy_met_rel_low=-1.0;
 float lim_delta_energy_met_rel_high=1.0; 


 TH1F* h_HZ_delta_EMissProj_gen_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_HZ_delta_EMissProj_gen_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_HZ_delta_EMissProj_gen_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MET_{jet-proj,gen}-MET_{gen})/MET_{gen}");
 TH1F* h_HZ_delta_MHMissProj_gen_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_HZ_delta_MHMissProj_gen_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_HZ_delta_MHMissProj_gen_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MHT_{jet-proj}-MET_{gen})/MET_{gen}");
 TH1F* h_HZ_delta_MHT_gen_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_HZ_delta_MHT_gen_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_HZ_delta_MHT_gen_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MHT_{jet-proj,gen}-MET_{gen})/MET_{gen}");
 
 TH1F* h_HZ_delta_EMiss_E_gen_MET_E_gen_over_MET_E_gen_sqrtS_2500 = new TH1F("h_HZ_delta_EMiss_E_gen_MET_E_gen_over_MET_E_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_HZ_delta_EMiss_E_gen_MET_E_gen_over_MET_E_gen_sqrtS_2500->GetXaxis()->SetTitle("(E_{miss,jet-proj,gen}-E_{miss,gen}/E_{miss,gen}");
 TH1F* h_HZ_delta_MHMiss_E_gen_MET_E_gen_over_MET_E_gen_sqrtS_2500 = new TH1F("h_HZ_delta_MHMiss_E_gen_MET_E_gen_over_MET_E_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_HZ_delta_MHMiss_E_gen_MET_E_gen_over_MET_E_gen_sqrtS_2500->GetXaxis()->SetTitle("(E^{MHT}_{miss,jet-proj,gen}-E_{miss,gen})/E_{miss,gen}");
 
 TH1F* h_HZ_dPhi_EMissProj_gen_MET_gen_sqrtS_2500 = new TH1F("h_HZ_dPhi_EMissProj_gen_MET_gen_sqrtS_2500","", n_bins_high_gen,lim_dalpha_low, lim_dalpha_high);
 h_HZ_dPhi_EMissProj_gen_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MET_{jet-proj,gen},MET_{gen}) [#circ]");
 TH1F* h_HZ_dPhi_MHMissProj_gen_MET_gen_sqrtS_2500 = new TH1F("h_HZ_dPhi_MHMissProj_gen_MET_gen_sqrtS_2500","", n_bins_high_gen,lim_dalpha_low, lim_dalpha_high);
 h_HZ_dPhi_MHMissProj_gen_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MET_{jet-proj,gen},MET_{gen}) [#circ]");
 TH1F* h_HZ_dPhi_MHT_gen_MET_gen_sqrtS_2500 = new TH1F("h_HZ_dPhi_MHT_gen_MET_gen_sqrtS_2500","", n_bins_high_gen,lim_dalpha_low, lim_dalpha_high);
 h_HZ_dPhi_MHT_gen_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MHT_{gen},MET_{gen}) [#circ]");
 
 TH1F* h_HZ_dAlpha_EMissProj_gen_MET_gen_sqrtS_2500 = new TH1F("h_HZ_dAlpha_EMissProj_gen_MET_gen_sqrtS_2500","", n_bins_high_gen,lim_dalpha_low, lim_dalpha_high);
 h_HZ_dAlpha_EMissProj_gen_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#alpha(MET_{jet-proj,gen},MET_{gen}) [#circ]");
 TH1F* h_HZ_dAlpha_MHMissProj_gen_MET_gen_sqrtS_2500 = new TH1F("h_HZ_dAlpha_MHMissProj_gen_MET_gen_sqrtS_2500","", n_bins_high_gen,lim_dalpha_low, lim_dalpha_high);
 h_HZ_dAlpha_MHMissProj_gen_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#alpha(MET_{jet-proj,gen},MET_{gen}) [#circ]");

  TH1F* h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_750 = new TH1F("h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_750 = new TH1F("h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");

  TH1F* h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_750 = new TH1F("h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_750 = new TH1F("h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");

  TH1F* h_HZ_MET_gen_over_E_tot_gen_sqrtS_2500 = new TH1F("h_HZ_MET_over_E_tot_gen_sqrtS_2500","", n_bins_high_gen,0.,lim_delta_energy_jet_rel_high);
  h_HZ_MET_gen_over_E_tot_gen_sqrtS_2500->GetXaxis()->SetTitle("MET_{gen}/E_{tot}");
  TH1F* h_HZ_MHT_gen_over_E_tot_gen_sqrtS_2500 = new TH1F("h_HZ_MHT_over_E_tot_gen_sqrtS_2500","", n_bins_high_gen,0.,lim_delta_energy_jet_rel_high);
  h_HZ_MHT_gen_over_E_tot_gen_sqrtS_2500->GetXaxis()->SetTitle("MHT_{gen}/E_{tot}");
  TH1F* h_HZ_MET_gen_over_E_j1_j2_gen_sqrtS_2500 = new TH1F("h_HZ_MET_over_E_j1_j2_gen_sqrtS_2500","", n_bins_high_gen,0.,lim_delta_energy_jet_rel_high);
  h_HZ_MET_gen_over_E_j1_j2_gen_sqrtS_2500->GetXaxis()->SetTitle("MET_{gen}/(E_{j1}+E_{j2})");
  TH1F* h_HZ_MHT_gen_over_E_j1_j2_gen_sqrtS_2500 = new TH1F("h_HZ_MHT_over_E_j1_j2_gen_sqrtS_2500","", n_bins_high_gen,0.,lim_delta_energy_jet_rel_high);
  h_HZ_MHT_gen_over_E_j1_j2_gen_sqrtS_2500->GetXaxis()->SetTitle("MHT_{gen}/(E_{j1}+E_{j2})");

  float lim_delta_mass_low=0;
  float lim_delta_mass_high=100;

  TH1F* h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_750 = new TH1F("h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_750","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_750->GetXaxis()->SetTitle("#Delta mass(gj1,gj2)");
  TH1F* h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_750_2500 = new TH1F("h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_700_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta mass(gj1,gj2)");
  TH1F* h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_2500 = new TH1F("h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_2500->GetXaxis()->SetTitle("#Delta mass(gj1,gj2)");

  TH1F* h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_750 = new TH1F("h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_750","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_750->GetXaxis()->SetTitle("#Delta mass(gj1,gj2)");
  TH1F* h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_750_2500 = new TH1F("h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_700_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta mass(gj1,gj2)");
  TH1F* h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_2500 = new TH1F("h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_2500->GetXaxis()->SetTitle("#Delta mass(gj1,gj2)");

  float lim_jet_charge_low=-2.0;
  float lim_jet_charge_high=2.0;

  TH1F* h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");


  TH1F* h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  //now sqrtS based on partons

  TH1F* h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");


  TH1F* h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");

  float delta_lim_jet_charge_low=-2.5;
  float delta_lim_jet_charge_high=3.5;

  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b}(H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b}(H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b}(H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500","", n_bins_high_gen, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b}(H)-matched)-subjet-charge_pt (b(H)-matched)");

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

 hist_vec_gen_HZ_1D.push_back(h_HZ_sqrtS_gen_isoPh);
 hist_vec_gen_HZ_1D.push_back(h_HZ_sqrtS_gen_isoPh_inv);
 hist_vec_gen_HZ_1D.push_back(h_HZ_sqrtS_gen);
 hist_vec_gen_HZ_1D.push_back(h_HZ_sqrtS_gen_inv);
 hist_vec_gen_HZ_1D.push_back(h_HZ_sqrtS_gen_truePh_inv);//265

 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_METProj_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_METProj_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_METProj_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_METProj_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_METProj_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_METProj_gen_sqrt_s_2500);//271
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_EMissProj_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_EMissProj_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_EMissProj_gen_sqrt_s_2500);
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_EMissProj_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_EMissProj_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_EMissProj_gen_sqrt_s_2500); //277
 hist_vec_gen_HZ_1D.push_back(h_HZ_sqrtS_gen_isoPh_EMissCorr);
 hist_vec_gen_HZ_1D.push_back(h_HZ_sqrtS_gen_j1_j2_isoPh);
 hist_vec_gen_HZ_1D.push_back(h_HZ_sqrtS_gen_j1_j2_isoPh_EMissCorr);//280
 
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_tot_j1_j2_gen_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_tot_j1_j2_EMiss_gen_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_tot_isoPh_gen_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_tot_isoPh_EMiss_gen_sqrt_s_2500); 
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j1_gen_sqrt_s_2500);//285
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j2_gen_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j1_gen_EMiss_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j2_gen_EMiss_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j1_min_E_j2_gen_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j1_min_E_j2_gen_EMiss_sqrt_s_2500);//290
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_EMiss_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_MHMissProj_gen_sqrt_s_0_750);//293
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_MHMissProj_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j1_MHMissProj_gen_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_MHMissProj_gen_sqrt_s_0_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_MHMissProj_gen_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_mass_j2_MHMissProj_gen_sqrt_s_2500);//298
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_tot_j1_j2_MHMiss_gen_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_tot_isoPh_MHMiss_gen_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j1_gen_MHMiss_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j2_gen_MHMiss_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_gen_MHMiss_sqrt_s_2500);//303

 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_EMissProj_gen_MET_gen_over_MET_gen_sqrtS_2500);//304
 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_MHMissProj_gen_MET_gen_over_MET_gen_sqrtS_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_MHT_gen_MET_gen_over_MET_gen_sqrtS_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_EMiss_E_gen_MET_E_gen_over_MET_E_gen_sqrtS_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_MHMiss_E_gen_MET_E_gen_over_MET_E_gen_sqrtS_2500);//308
 hist_vec_gen_HZ_1D.push_back(h_HZ_dPhi_EMissProj_gen_MET_gen_sqrtS_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dPhi_MHMissProj_gen_MET_gen_sqrtS_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dPhi_MHT_gen_MET_gen_sqrtS_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_EMissProj_gen_MET_gen_sqrtS_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_dAlpha_MHMissProj_gen_MET_gen_sqrtS_2500);//313
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj1_sj_CosMinHelicityAngle_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500);//318
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj1_sj_CosMaxHelicityAngle_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj2_sj_CosMinHelicityAngle_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_750);//323
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_gj2_sj_CosMaxHelicityAngle_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_MET_gen_over_E_tot_gen_sqrtS_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_MHT_gen_over_E_tot_gen_sqrtS_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_MET_gen_over_E_j1_j2_gen_sqrtS_2500);//328
 hist_vec_gen_HZ_1D.push_back(h_HZ_MHT_gen_over_E_j1_j2_gen_sqrtS_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_mass_gen_gj1_gj2_sqrt_s_2500);//332

 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_750);
 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_750_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_delta_mass_EMiss_gen_gj1_gj2_sqrt_s_2500);//335

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);//340
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);//345
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);//350
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);//352
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);//357
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);//362
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);//367


 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);//372
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);//377
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);//382
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_matched_H_sqrt_s_part_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);//384
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);//389
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_50_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_75_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);//394
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_00_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_10_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);//399

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);//400
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_Z_q_pos_gj_min_matched_Z_q_neg_gj_matched_H_sqrt_s_part_2500);

 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500);//404
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500);
 hist_vec_gen_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_gen_sj_matched_H_bbar_gj_min_matched_H_b_gj_matched_H_sqrt_s_part_2500);

 
 for(unsigned int i=0;i<hist_vec_gen_HZ_1D.size();i++){
   hist_vec_gen_HZ_1D[i]->Sumw2();
   hist_vec_gen_HZ_1D[i]->SetLineColor(kBlack);
   hist_vec_gen_HZ_1D[i]->SetLineWidth(2);
   hist_vec_gen_HZ_1D[i]->GetYaxis()->SetTitle("Events");
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
 
 TH1F* h_HZ_mass_j1_METProj_reco_sqrt_s_0_750 = new TH1F("h_HZ_mass_j1_METProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_METProj_reco_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j1_METProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_METProj_reco_sqrt_s_2500 = new TH1F("h_HZ_mass_j1_METProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_HZ_mass_j1_METProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_HZ_mass_j1_METProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_HZ_mass_j1_METProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 
 TH1F* h_HZ_mass_j2_METProj_reco_sqrt_s_0_750 = new TH1F("h_HZ_mass_j2_METProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_METProj_reco_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j2_METProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_METProj_reco_sqrt_s_2500 = new TH1F("h_HZ_mass_j2_METProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_HZ_mass_j2_METProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_HZ_mass_j2_METProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_HZ_mass_j2_METProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_HZ_mass_j1_EMissProj_reco_sqrt_s_0_750 = new TH1F("h_HZ_mass_j1_EMissProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_EMissProj_reco_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j1_EMissProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j1_EMissProj_reco_sqrt_s_2500 = new TH1F("h_HZ_mass_j1_EMissProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_HZ_mass_j1_EMissProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_HZ_mass_j1_EMissProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_HZ_mass_j1_EMissProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 
 TH1F* h_HZ_mass_j2_EMissProj_reco_sqrt_s_0_750 = new TH1F("h_HZ_mass_j2_EMissProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_EMissProj_reco_sqrt_s_750_2500 = new TH1F("h_HZ_mass_j2_EMissProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_HZ_mass_j2_EMissProj_reco_sqrt_s_2500 = new TH1F("h_HZ_mass_j2_EMissProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_HZ_mass_j2_EMissProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_HZ_mass_j2_EMissProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_HZ_mass_j2_EMissProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 TH1F* h_HZ_sqrtS_reco_isoPh_EMissCorr = new TH1F("h_HZ_sqrtS_reco_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_sqrtS_reco_isoPh_EMissCorr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500 = new TH1F("h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_HZ_sqrtS_reco_j1_j2_isoPh = new TH1F("h_HZ_sqrtS_reco_j1_j2_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_sqrtS_reco_j1_j2_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 TH1F* h_HZ_sqrtS_reco_j1_j2_isoPh_EMissCorr = new TH1F("h_HZ_sqrtS_reco_j1_j2_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_sqrtS_reco_j1_j2_isoPh_EMissCorr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 
 TH1F* h_HZ_E_tot_j1_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_E_tot_j1_j2_reco_sqrt_s_2500","", n_bins_high, lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_HZ_E_tot_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
 TH1F* h_HZ_E_tot_j1_j2_EMiss_reco_sqrt_s_2500 = new TH1F("h_HZ_E_tot_j1_j2_EMiss_reco_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_HZ_E_tot_j1_j2_EMiss_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");

 TH1F* h_HZ_E_tot_isoPh_reco_sqrt_s_2500 = new TH1F("h_HZ_E_tot_isoPh_reco_sqrt_s_2500","", n_bins_high, lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_HZ_E_tot_isoPh_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{tot} [GeV]");
 TH1F* h_HZ_E_tot_isoPh_EMiss_reco_sqrt_s_2500 = new TH1F("h_HZ_E_tot_isoPh_EMiss_reco_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_HZ_E_tot_isoPh_EMiss_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{tot}+E_{miss}^{jetProj}) [GeV]");
 
 TH1F* h_HZ_E_j1_reco_sqrt_s_2500 = new TH1F("h_HZ_E_j1_reco_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_HZ_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_HZ_E_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_HZ_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");
 TH1F* h_HZ_E_j1_reco_EMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j1_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_HZ_E_j1_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_HZ_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_HZ_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");

 TH1F* h_HZ_E_j1_min_E_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_E_j1_min_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);
 h_HZ_E_j1_min_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2}) [GeV]");
 TH1F* h_HZ_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);
 h_HZ_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2}) [GeV]");

 TH1F* h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500 = new TH1F("h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");
 TH1F* h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");


 //
 TH1F* h_HZ_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_HZ_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_HZ_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MET_{jet-proj,reco}-MET_{gen})/MET_{gen}");
 TH1F* h_HZ_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_HZ_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_HZ_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MHT_{jet-proj}-MET_{reco})/MET_{gen}");
 TH1F* h_HZ_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_HZ_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_HZ_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MHT_{jet-proj,reco}-MET_{gen})/MET_{gen}");
 
 TH1F* h_HZ_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500 = new TH1F("h_HZ_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_HZ_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500->GetXaxis()->SetTitle("(E_{miss,jet-proj,reco}-E_{miss,gen}/E_{miss,gen}");
 TH1F* h_HZ_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500 = new TH1F("h_HZ_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_HZ_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500->GetXaxis()->SetTitle("(E^{MHT}_{miss,jet-proj,reco}-E_{miss,gen})/E_{miss,gen}");
 
 TH1F* h_HZ_dPhi_EMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_HZ_dPhi_EMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_HZ_dPhi_EMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_HZ_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_HZ_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_HZ_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_HZ_dPhi_MHT_reco_MET_gen_sqrtS_2500 = new TH1F("h_HZ_dPhi_MHT_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_HZ_dPhi_MHT_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MHT_{reco},MET_{gen}) [#circ]");
 
 TH1F* h_HZ_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_HZ_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_HZ_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#alpha(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_HZ_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_HZ_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_HZ_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#alpha(MET_{jet-proj,reco},MET_{gen}) [#circ]");


  TH1F* h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_750 = new TH1F("h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_750 = new TH1F("h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");

  TH1F* h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_750 = new TH1F("h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_750 = new TH1F("h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");

  TH1F* h_HZ_MET_reco_over_E_tot_reco_sqrtS_2500 = new TH1F("h_HZ_MET_reco_over_E_tot_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_HZ_MET_reco_over_E_tot_reco_sqrtS_2500->GetXaxis()->SetTitle("MET_{reco}/E_{tot}");
  TH1F* h_HZ_MHT_reco_over_E_tot_reco_sqrtS_2500 = new TH1F("h_HZ_MHT_reco_over_E_tot_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_HZ_MHT_reco_over_E_tot_reco_sqrtS_2500->GetXaxis()->SetTitle("MHT_{reco}/E_{tot}");
  TH1F* h_HZ_MET_reco_over_E_j1_j2_reco_sqrtS_2500 = new TH1F("h_HZ_MET_reco_over_E_j1_j2_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_HZ_MET_reco_over_E_j1_j2_reco_sqrtS_2500->GetXaxis()->SetTitle("MET_{reco}/(E_{j1}+E_{j2})");
  TH1F* h_HZ_MHT_reco_over_E_j1_j2_reco_sqrtS_2500 = new TH1F("h_HZ_MHT_reco_over_E_j1_j2_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_HZ_MHT_reco_over_E_j1_j2_reco_sqrtS_2500->GetXaxis()->SetTitle("MHT_{reco}/(E_{j1}+E_{j2})");

  TH1F* h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_750 = new TH1F("h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_750","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_750->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_750_2500 = new TH1F("h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_700_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_2500 = new TH1F("h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");

  TH1F* h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750 = new TH1F("h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500 = new TH1F("h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_700_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500 = new TH1F("h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");


  TH1F* h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");


  TH1F* h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  //now but on parton sqrtS
  TH1F* h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q+(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)");


  TH1F* h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_E (q-(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, lim_jet_charge_low,lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q-(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");

  TH1F* h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (q+(Z)-matched)-subjet-charge_pt (q-(Z)-matched)");




  TH1F* h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b} (H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b} (H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b} (H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b} (H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b} (H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b} (H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b} (H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_25_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_25_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_25_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b} (H)-matched)-subjet-charge_pt (b(H)-matched)");
  TH1F* h_HZ_subjet_charge_pt_kappa_1_50_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500 = new TH1F("h_HZ_subjet_charge_pt_kappa_1_50_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500","", n_bins_high_reco, delta_lim_jet_charge_low,delta_lim_jet_charge_high);
  h_HZ_subjet_charge_pt_kappa_1_50_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500->GetXaxis()->SetTitle("subjet-charge_pt (#bar{b} (H)-matched)-subjet-charge_pt (b(H)-matched)");



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

 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j1_METProj_reco_sqrt_s_0_750);//282
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j1_METProj_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j1_METProj_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j2_METProj_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j2_METProj_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j2_METProj_reco_sqrt_s_2500);//287
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j1_EMissProj_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j1_EMissProj_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j1_EMissProj_reco_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j2_EMissProj_reco_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j2_EMissProj_reco_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_mass_j2_EMissProj_reco_sqrt_s_2500); //293

 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_isoPh_EMissCorr);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500);//297
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500);//303
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500);//309
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500);//312
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_j1_j2_isoPh);
 hist_vec_reco_HZ_1D.push_back(h_HZ_sqrtS_reco_j1_j2_isoPh_EMissCorr);//314

 hist_vec_reco_HZ_1D.push_back(h_HZ_E_tot_j1_j2_reco_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_tot_j1_j2_EMiss_reco_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_tot_isoPh_reco_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_tot_isoPh_EMiss_reco_sqrt_s_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_j1_reco_sqrt_s_2500);//319
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_j2_reco_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_j1_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_j2_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_j1_min_E_j2_reco_sqrt_s_2500);//323
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500);

 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500);//327
 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_dPhi_EMissProj_reco_MET_gen_sqrtS_2500);//332
 hist_vec_reco_HZ_1D.push_back(h_HZ_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dPhi_MHT_reco_MET_gen_sqrtS_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_750);//337
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj1_sj_CosMinHelicityAngle_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500);//342
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj2_sj_CosMinHelicityAngle_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500);//347
 hist_vec_reco_HZ_1D.push_back(h_HZ_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_MET_reco_over_E_tot_reco_sqrtS_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_MHT_reco_over_E_tot_reco_sqrtS_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_MET_reco_over_E_j1_j2_reco_sqrtS_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_MHT_reco_over_E_j1_j2_reco_sqrtS_2500);//352

 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_2500);//355

 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750);
 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500);//358


 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);//359
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);//364
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);//369
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_2500);//374

 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);//379
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);//384
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);//389
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
 //parton S cut
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);//391
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);//396
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);//401
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_matched_H_sqrt_s_part_2500);//406

 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);//411
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);//416
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_10_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);//421
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);//423
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500); //427
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500);//431 

 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_15_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500);//432
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500);
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500);//436
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_00_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_25_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500); 
 hist_vec_reco_HZ_1D.push_back(h_HZ_subjet_charge_pt_kappa_1_50_reco_sj_matched_H_bbar_rj_min_matched_H_b_rj_matched_H_sqrt_s_part_2500); 




 for(unsigned int i=0;i<hist_vec_reco_HZ_1D.size();i++){
   hist_vec_reco_HZ_1D[i]->Sumw2();
   hist_vec_reco_HZ_1D[i]->SetLineColor(kBlack);
   hist_vec_reco_HZ_1D[i]->SetLineWidth(2);
   hist_vec_reco_HZ_1D[i]->GetYaxis()->SetTitle("Events");
 }

 
 TH2F* h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_0_750 = new TH2F("h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_0_750->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500 = new TH2F("h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_2500 = new TH2F("h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_2500->GetYaxis()->SetTitle("jet2 mass [GeV]");

 float n_bins_low_METCorr=-0.01;
 float n_bins_high_METCorr=0.21;

 TH2F* h_HZ_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500 = new TH2F("h_HZ_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500","", n_bins_high,n_bins_low_METCorr,n_bins_high_METCorr , n_bins_high,n_bins_low_METCorr, n_bins_high_METCorr);
 h_HZ_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetXaxis()->SetTitle("gj1 METCorr");
 h_HZ_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetYaxis()->SetTitle("gj2 METCorr");

 TH2F* h_HZ_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500 = new TH2F("h_HZ_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500","", n_bins_high,n_bins_low_METCorr,n_bins_high_METCorr , n_bins_high,n_bins_low_METCorr, n_bins_high_METCorr);
 h_HZ_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetXaxis()->SetTitle("rj1 METCorr");
 h_HZ_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetYaxis()->SetTitle("rj2 METCorr");

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

  //parton vs gen true --> here we calculat all true minus iso photons
 TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh_EMissCorr = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_gen_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh_EMissCorr->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh_EMissCorr->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
 //parton vs reco true --> here we calculat all true minus iso photons
 TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_isoPh_EMissCorr = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_reco_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_isoPh_EMissCorr->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_isoPh_EMissCorr->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
 //parton vs reco true --> here we calculat all true minus iso photons
 TH2F* h_HZ_srqtS_gen_isoPh_EMissCorr_eff_vs_sqrtS_reco_isoPh_EMissCorr = new TH2F("h_HZ_srqtS_gen_isoPh_EMissCorr_vs_sqrtS_reco_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_srqtS_gen_isoPh_EMissCorr_eff_vs_sqrtS_reco_isoPh_EMissCorr->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
 h_HZ_srqtS_gen_isoPh_EMissCorr_eff_vs_sqrtS_reco_isoPh_EMissCorr->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");

  //parton vs gen true --> here we calculat all true minus iso photons
 TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_j1_j2_isoPh = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_gen_j1_j2_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_j1_j2_isoPh->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_j1_j2_isoPh->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
  //parton vs gen true --> here we calculat all true minus iso photons
 TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_j1_j2_isoPh_EMissCorr = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_gen_j1_j2_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_j1_j2_isoPh_EMissCorr->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_j1_j2_isoPh_EMissCorr->GetYaxis()->SetTitle("gen #sqrt{s}_{eff}");
  //parton vs reco true --> here we calculat all true minus iso photons
 TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_j1_j2_isoPh = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_reco_j1_j2_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_j1_j2_isoPh->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_j1_j2_isoPh->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs reco true --> here we calculat all true minus iso photons
 TH2F* h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_j1_j2_isoPh_EMissCorr = new TH2F("h_HZ_srqtS_e1e2_vs_sqrtS_reco_j1_j2_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_j1_j2_isoPh_EMissCorr->GetXaxis()->SetTitle("parton #sqrt{s}_{eff}");
 h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_j1_j2_isoPh_EMissCorr->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs reco true --> here we calculat all true minus iso photons
 TH2F* h_HZ_sqrtS_gen_j1_j2_isoPh_vs_sqrtS_reco_j1_j2_isoPh = new TH2F("h_HZ_sqrtS_gen_j1_j2_isoPh_vs_sqrtS_reco_j1_j2_isoPh","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_sqrtS_gen_j1_j2_isoPh_vs_sqrtS_reco_j1_j2_isoPh->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
 h_HZ_sqrtS_gen_j1_j2_isoPh_vs_sqrtS_reco_j1_j2_isoPh->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");
  //parton vs reco true --> here we calculat all true minus iso photons
 TH2F* h_HZ_sqrtS_gen_j1_j2_isoPh_EMissCorr_eff_vs_sqrtS_reco_j1_j2_isoPh_EMissCorr = new TH2F("h_HZ_sqrtS_gen_j1_j2_isoPh_EMissCorr_vs_sqrtS_reco_j1_j2_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high);
 h_HZ_sqrtS_gen_j1_j2_isoPh_EMissCorr_eff_vs_sqrtS_reco_j1_j2_isoPh_EMissCorr->GetXaxis()->SetTitle("gen #sqrt{s}_{eff}");
 h_HZ_sqrtS_gen_j1_j2_isoPh_EMissCorr_eff_vs_sqrtS_reco_j1_j2_isoPh_EMissCorr->GetYaxis()->SetTitle("reco #sqrt{s}_{eff}");

 int n_bins_mult=101;
 float n_bins_low_mul=0.5;
 float n_bins_high_mult=100.5;

 TH2F* h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_0_750 = new TH2F("h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_0_750","", n_bins_mult,n_bins_low_mul,n_bins_high_mult, n_bins_mult,n_bins_low_mul,n_bins_high_mult);
 h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_0_750->GetXaxis()->SetTitle("N_{tracks}^{gen}");
 h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_0_750->GetYaxis()->SetTitle("N_{tracks}^{reco}");

 TH2F* h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_750_2500 = new TH2F("h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_750_2500","", n_bins_mult,n_bins_low_mul,n_bins_high_mult, n_bins_mult,n_bins_low_mul,n_bins_high_mult);
 h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{tracks}^{gen}");
 h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_750_2500->GetYaxis()->SetTitle("N_{tracks}^{reco}");

 TH2F* h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_2500 = new TH2F("h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_2500","", n_bins_mult,n_bins_low_mul,n_bins_high_mult, n_bins_mult,n_bins_low_mul,n_bins_high_mult);
 h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_2500->GetXaxis()->SetTitle("N_{tracks}^{gen}");
 h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_2500->GetYaxis()->SetTitle("N_{tracks}^{reco}");

 TH2F* h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_0_750 = new TH2F("h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_0_750","", n_bins_mult,n_bins_low_mul,n_bins_high_mult, n_bins_mult,n_bins_low_mul,n_bins_high_mult);
 h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_0_750->GetXaxis()->SetTitle("N_{tracks}^{gen}");
 h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_0_750->GetYaxis()->SetTitle("N_{tracks}^{reco}");

 TH2F* h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_750_2500 = new TH2F("h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_750_2500","", n_bins_mult,n_bins_low_mul,n_bins_high_mult, n_bins_mult,n_bins_low_mul,n_bins_high_mult);
 h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{tracks}^{gen}");
 h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_750_2500->GetYaxis()->SetTitle("N_{tracks}^{reco}");

 TH2F* h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_2500 = new TH2F("h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_2500","", n_bins_mult,n_bins_low_mul,n_bins_high_mult, n_bins_mult,n_bins_low_mul,n_bins_high_mult);
 h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_2500->GetXaxis()->SetTitle("N_{tracks}^{gen}");
 h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_2500->GetYaxis()->SetTitle("N_{tracks}^{reco}");

 float n_bins_low_dalpha=0;
 float n_bins_high_dalpha=15;


 TH2F* h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_0_750 = new TH2F("h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},b)");
 h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_0_750->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},b)");

 TH2F* h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_750_2500 = new TH2F("h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},b)");
 h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_750_2500->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},b)");

 TH2F* h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_2500 = new TH2F("h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},b)");
 h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_2500->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},b)");

 TH2F* h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_0_750 = new TH2F("h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},bbar)");
 h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_0_750->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},bbar)");

 TH2F* h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_750_2500 = new TH2F("h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},bbar)");
 h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_750_2500->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},bbar)");

 TH2F* h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_2500 = new TH2F("h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},bbar)");
 h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_2500->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},bbar)");




 TH2F* h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_0_750 = new TH2F("h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},q^{-})");
 h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_0_750->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},q^{-})");

 TH2F* h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_750_2500 = new TH2F("h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},q^{-})");
 h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_750_2500->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},q^{-})");

 TH2F* h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_2500 = new TH2F("h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},q^{-})");
 h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_2500->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},q^{-})");

 TH2F* h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_0_750 = new TH2F("h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},q^{+})");
 h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_0_750->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},q^{+})");

 TH2F* h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_750_2500 = new TH2F("h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},q^{+})");
 h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_750_2500->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},q^{+})");

 TH2F* h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_2500 = new TH2F("h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha, n_bins_high_gen_mass,n_bins_low_dalpha,n_bins_high_dalpha);
 h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj_{g},q^{+})");
 h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_2500->GetYaxis()->SetTitle("#Delta#alpha(sj_{r},q^{+})");

  TH2F* h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);

  TH2F* h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);


  TH2F* h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);

  TH2F* h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);
  TH2F* h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500 = new TH2F("h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500","", n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high, n_bins_high_gen, lim_jet_charge_low,lim_jet_charge_high);


  std::vector<TH2F*> hist_vec_HZ_2DHist;  
  hist_vec_HZ_2DHist.push_back(h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_0_750);
  hist_vec_HZ_2DHist.push_back(h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_mass_j1_vs_mass_j2_reco_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500);
  //hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_genTrue);
  //hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_genJet);//5
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_isoPh);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_recoJet);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_gen_isoPh_eff_vs_sqrtS_reco_isoPh);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_genJet_eff_vs_sqrtS_recoJet);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_genAll);//10
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_recoAll);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh_inv);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_inv);//13
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_isoPh_EMissCorr);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_isoPh_EMissCorr);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_gen_isoPh_EMissCorr_eff_vs_sqrtS_reco_isoPh_EMissCorr);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_j1_j2_isoPh);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_gen_j1_j2_isoPh_EMissCorr);//18
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_j1_j2_isoPh);
  hist_vec_HZ_2DHist.push_back(h_HZ_srqtS_e1e2_eff_vs_sqrtS_reco_j1_j2_isoPh_EMissCorr);
  hist_vec_HZ_2DHist.push_back(h_HZ_sqrtS_gen_j1_j2_isoPh_vs_sqrtS_reco_j1_j2_isoPh);
  hist_vec_HZ_2DHist.push_back(h_HZ_sqrtS_gen_j1_j2_isoPh_EMissCorr_eff_vs_sqrtS_reco_j1_j2_isoPh_EMissCorr);//22
  hist_vec_HZ_2DHist.push_back(h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_0_750);//23
  hist_vec_HZ_2DHist.push_back(h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_750_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_NTrack_gen_vs_reco_J_H_matched_gen_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_0_750);
  hist_vec_HZ_2DHist.push_back(h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_750_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_NTrack_gen_vs_reco_J_Z_matched_gen_sqrt_s_2500);//28
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_0_750);
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_750_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_b_vs_sj_reco_b_J_H_matched_gen_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_0_750);
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_750_2500);//33
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_bbar_vs_sj_reco_bbar_J_H_matched_gen_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_0_750);
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_750_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_q_min_vs_sj_reco_q_min_J_Z_matched_gen_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_0_750);//38
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_750_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_dAlpha_sj_gen_q_plus_vs_sj_reco_q_plus_J_Z_matched_gen_sqrt_s_2500);//40
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);//41
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_vs_matched_Z_q_neg_gj_matched_H_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);//45
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_Z_q_pos_vs_matched_Z_q_neg_rj_matched_H_sqrt_s_2500);

  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_15_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500);//49
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_20_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_25_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_30_gen_sj_matched_H_bbar_vs_matched_H_b_gj_matched_H_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_15_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500);//53
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_20_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_25_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500);
  hist_vec_HZ_2DHist.push_back(h_HZ_subjet_charge_E_kappa_0_30_reco_sj_matched_H_bbar_vs_matched_H_b_rj_matched_H_sqrt_s_2500);



  for(unsigned int i=0;i<hist_vec_HZ_2DHist.size();i++){
    hist_vec_HZ_2DHist[i]->Sumw2();
  }

 std::vector<TH1F*> hist_vec_reco_HZ_1D_reco_vs_gen_selection;  
 std::vector<TH2F*> hist_vec_reco_HZ_2D_reco_vs_gen_selection;
 std::cout<<"before filling HZ"<<std::endl;
 fill_HZ_histograms(file_CLIC_HZqq, hist_vec_HZ_parton, hist_vec_gen_HZ_1D, hist_vec_reco_HZ_1D, hist_vec_HZ_2DHist, hist_vec_reco_HZ_1D_reco_vs_gen_selection, hist_vec_reco_HZ_2D_reco_vs_gen_selection, usePartonInfo ,xsec_hz_qq,fillPartonInfo,fillGenInfo);
 std::cout<<"after filling HZ"<<std::endl;
 int index0_cut=-1;
 for(int i=0;i<(h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1);i++){
   if(index0_cut==-1 && h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetBinLowEdge(i)>=0){
     index0_cut=i;
     break;
   }
 }
 std::cout<<"ratio kappa 0_20 "<<h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(0,index0_cut)/h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(0,h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1)<<" events left "<<h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(index0_cut,h_HZ_subjet_charge_pt_kappa_0_20_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1)<<std::endl;
 index0_cut=-1;
 for(int i=0;i<(h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1);i++){
   if(index0_cut==-1 && h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetBinLowEdge(i)>=0){
     index0_cut=i;
     break;
   }
 }
 std::cout<<"ratio kappa 0_25 "<<h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(0,index0_cut)/h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(0,h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1)<<" events left "<<h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(index0_cut,h_HZ_subjet_charge_pt_kappa_0_25_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1)<<std::endl;
 index0_cut=-1;
 for(int i=0;i<(h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1);i++){
   if(index0_cut==-1 && h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetBinLowEdge(i)>=0){
     index0_cut=i;
     break;
   }
 }
 std::cout<<"ratio kappa 0_30 "<<h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(0,index0_cut)/h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(0,h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1)<<" events left "<<h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(index0_cut,h_HZ_subjet_charge_pt_kappa_0_30_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1)<<std::endl;
 index0_cut=-1;
 for(int i=0;i<(h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1);i++){
   if(index0_cut==-1 && h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetBinLowEdge(i)>=0){
     index0_cut=i;
     break;
   }
 }
 std::cout<<"ratio kappa 0_50 "<<h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(0,index0_cut)/h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(0,h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1)<<" events left "<<h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->Integral(index0_cut,h_HZ_subjet_charge_pt_kappa_0_50_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1)<<std::endl;
 index0_cut=-1;
 for(int i=0;i<(h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetNbinsX()+1);i++){
   if(index0_cut==-1 && h_HZ_subjet_charge_pt_kappa_0_75_reco_sj_matched_Z_q_pos_rj_min_matched_Z_q_neg_rj_matched_H_sqrt_s_part_2500->GetBinLowEdge(i)>=0){
     index0_cut=i;
     break;
   }
 }
  /*

  //NOW _tt

 TH1F* h_tt_dAlpha_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_dAlpha_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_tt_dAlpha_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_dAlpha_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_tt_dAlpha_j1_j2_reco_sqrt_s_2500 = new TH1F("h_tt_dAlpha_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_tt_dAlpha_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");
 h_tt_dAlpha_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");
 h_tt_dAlpha_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");

 TH1F* h_tt_dPhi_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_dPhi_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_tt_dPhi_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_dPhi_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_tt_dPhi_j1_j2_reco_sqrt_s_2500 = new TH1F("h_tt_dPhi_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_tt_dPhi_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");
 h_tt_dPhi_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");
 h_tt_dPhi_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");

 TH1F* h_tt_dTheta_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_dTheta_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_tt_dTheta_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_dTheta_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_tt_dTheta_j1_j2_reco_sqrt_s_2500 = new TH1F("h_tt_dTheta_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_tt_dTheta_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");
 h_tt_dTheta_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");
 h_tt_dTheta_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");

 TH1F* h_tt_Theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_Theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_tt_Theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_Theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_tt_Theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_Theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_tt_Theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#theta(j1)[#circ]");
 h_tt_Theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#theta(j1)[#circ]");
 h_tt_Theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("#theta(j1)[#circ]");

 TH1F* h_tt_Theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_Theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_tt_Theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_Theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_tt_Theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_Theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_tt_Theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 h_tt_Theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 h_tt_Theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 
 TH1F* h_tt_mass_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_mass_j1_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_mass_j1_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j1_reco_sqrt_s_2500 = new TH1F("h_tt_mass_j1_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 
 h_tt_mass_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_tt_mass_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_tt_mass_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");

 TH1F* h_tt_mass_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_mass_j2_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_mass_j2_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j2_reco_sqrt_s_2500 = new TH1F("h_tt_mass_j2_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);

 h_tt_mass_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_tt_mass_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_tt_mass_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_tt_tau21_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_tau21_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau21_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_tau21_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau21_j1_reco_sqrt_s_2500 = new TH1F("h_tt_tau21_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_tt_tau21_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_tt_tau21_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_tt_tau21_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j1)");

 TH1F* h_tt_tau21_lrz_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_tau21_lrz_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau21_lrz_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_tau21_lrz_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau21_lrz_j1_reco_sqrt_s_2500 = new TH1F("h_tt_tau21_lrz_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_tt_tau21_lrz_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_tt_tau21_lrz_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_tt_tau21_lrz_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j1)");

 TH1F* h_tt_tau32_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_tau32_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau32_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_tau32_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau32_j1_reco_sqrt_s_2500 = new TH1F("h_tt_tau32_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_tt_tau32_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_tt_tau32_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_tt_tau32_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j1)");

 TH1F* h_tt_tau32_lrz_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_tau32_lrz_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau32_lrz_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_tau32_lrz_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau32_lrz_j1_reco_sqrt_s_2500 = new TH1F("h_tt_tau32_lrz_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_tt_tau32_lrz_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_tt_tau32_lrz_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_tt_tau32_lrz_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j1)");

 TH1F* h_tt_tau21_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_tau21_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau21_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_tau21_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau21_j2_reco_sqrt_s_2500 = new TH1F("h_tt_tau21_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_tt_tau21_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_tt_tau21_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_tt_tau21_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j2)");

 TH1F* h_tt_tau21_lrz_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_tau21_lrz_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau21_lrz_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_tau21_lrz_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau21_lrz_j2_reco_sqrt_s_2500 = new TH1F("h_tt_tau21_lrz_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_tt_tau21_lrz_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_tt_tau21_lrz_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_tt_tau21_lrz_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j2)");

 TH1F* h_tt_tau32_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_tau32_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau32_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_tau32_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau32_j2_reco_sqrt_s_2500 = new TH1F("h_tt_tau32_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_tt_tau32_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_tt_tau32_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_tt_tau32_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j2)");

 TH1F* h_tt_tau32_lrz_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_tau32_lrz_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau32_lrz_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_tau32_lrz_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_tt_tau32_lrz_j2_reco_sqrt_s_2500 = new TH1F("h_tt_tau32_lrz_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_tt_tau32_lrz_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_tt_tau32_lrz_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_tt_tau32_lrz_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 //N2 j1
 TH1F* h_tt_N2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_tt_N2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_tt_N2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");

 TH1F* h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");

 TH1F* h_tt_N2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_tt_N2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_tt_N2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");

 TH1F* h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");

 TH1F* h_tt_N2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_tt_N2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_tt_N2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");

 TH1F* h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 //N2 j2
 TH1F* h_tt_N2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_tt_N2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_tt_N2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");

 TH1F* h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");

 TH1F* h_tt_N2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_tt_N2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_tt_N2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");

 TH1F* h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");

 TH1F* h_tt_N2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_tt_N2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_tt_N2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");

 TH1F* h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");



 //N3 j1
 TH1F* h_tt_N3_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_tt_N3_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_tt_N3_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");

 TH1F* h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");

 TH1F* h_tt_N3_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_tt_N3_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_tt_N3_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");

 TH1F* h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");

 TH1F* h_tt_N3_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_tt_N3_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_tt_N3_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");

 TH1F* h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 //N3 j2
 TH1F* h_tt_N3_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_tt_N3_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_tt_N3_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");

 TH1F* h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");

 TH1F* h_tt_N3_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_tt_N3_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_tt_N3_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");

 TH1F* h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");

 TH1F* h_tt_N3_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_tt_N3_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_tt_N3_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");

 TH1F* h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");


 //now the C2 and C3 series
 //C2 j1
 TH1F* h_tt_C2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_tt_C2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_tt_C2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");

 TH1F* h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");

 TH1F* h_tt_C2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_tt_C2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_tt_C2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");

 TH1F* h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");

 TH1F* h_tt_C2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_tt_C2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_tt_C2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");

 TH1F* h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 //C2 j2
 TH1F* h_tt_C2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_tt_C2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_tt_C2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");

 TH1F* h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");

 TH1F* h_tt_C2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_tt_C2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_tt_C2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");

 TH1F* h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");

 TH1F* h_tt_C2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_tt_C2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_tt_C2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");

 TH1F* h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");



 //C3 j1
 TH1F* h_tt_C3_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_tt_C3_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_tt_C3_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");

 TH1F* h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");

 TH1F* h_tt_C3_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_tt_C3_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_tt_C3_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");

 TH1F* h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");

 TH1F* h_tt_C3_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_tt_C3_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_tt_C3_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");

 TH1F* h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 //C3 j2
 TH1F* h_tt_C3_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_tt_C3_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_tt_C3_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");

 TH1F* h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");

 TH1F* h_tt_C3_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_tt_C3_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_tt_C3_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");

 TH1F* h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");

 TH1F* h_tt_C3_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_tt_C3_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_tt_C3_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");

 TH1F* h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");

 //D2 series
 //D2 j1
 TH1F* h_tt_D2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_tt_D2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_tt_D2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");

 TH1F* h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");

 TH1F* h_tt_D2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_tt_D2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_tt_D2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");

 TH1F* h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");

 TH1F* h_tt_D2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_tt_D2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_tt_D2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");

 TH1F* h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 //D2 j2
 TH1F* h_tt_D2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_tt_D2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_tt_D2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");

 TH1F* h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");

 TH1F* h_tt_D2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_tt_D2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_tt_D2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");

 TH1F* h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");

 TH1F* h_tt_D2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_tt_D2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_tt_D2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");

 TH1F* h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");




 //D2 (1,2) series
 //D2 (1,2) j1
 TH1F* h_tt_D2_1_2_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_1_2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_tt_D2_1_2_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_1_2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_tt_D2_1_2_j1_reco_sqrt_s_2500 = new TH1F("h_tt_D2_1_2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_tt_D2_1_2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_tt_D2_1_2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_tt_D2_1_2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");

 TH1F* h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");

 //D2_1_2 (1,2) j2
 TH1F* h_tt_D2_1_2_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_1_2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_tt_D2_1_2_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_1_2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_tt_D2_1_2_j2_reco_sqrt_s_2500 = new TH1F("h_tt_D2_1_2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_tt_D2_1_2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_tt_D2_1_2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_tt_D2_1_2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");

 TH1F* h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");



 TH1F* h_tt_subjet1_E_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_subjet1_E_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_tt_subjet1_E_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_subjet1_E_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_tt_subjet1_E_j1_reco_sqrt_s_2500 = new TH1F("h_tt_subjet1_E_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_tt_subjet1_E_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");
 h_tt_subjet1_E_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");
 h_tt_subjet1_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");

 TH1F* h_tt_subjet2_E_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_subjet2_E_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_tt_subjet2_E_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_subjet2_E_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_tt_subjet2_E_j1_reco_sqrt_s_2500 = new TH1F("h_tt_subjet2_E_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_tt_subjet2_E_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");
 h_tt_subjet2_E_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");
 h_tt_subjet2_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");

 TH1F* h_tt_subjet1_E_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_subjet1_E_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_tt_subjet1_E_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_subjet1_E_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_tt_subjet1_E_j2_reco_sqrt_s_2500 = new TH1F("h_tt_subjet1_E_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_tt_subjet1_E_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");
 h_tt_subjet1_E_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");
 h_tt_subjet1_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");

 TH1F* h_tt_subjet2_E_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_subjet2_E_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_tt_subjet2_E_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_subjet2_E_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_tt_subjet2_E_j2_reco_sqrt_s_2500 = new TH1F("h_tt_subjet2_E_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_tt_subjet2_E_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");
 h_tt_subjet2_E_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");
 h_tt_subjet2_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");

 TH1F* h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_2500 = new TH1F("h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}/E (j1)");
 h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}/E (j1)");
 h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}/E (j1)");

 TH1F* h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_2500 = new TH1F("h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}/E (j2)");
 h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}/E (j2)");
 h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}/E (j2)");


 TH1F* h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750 = new TH1F("h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500 = new TH1F("h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500 = new TH1F("h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1 [#circ]");

 TH1F* h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750 = new TH1F("h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500 = new TH1F("h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500 = new TH1F("h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");


 TH1F* h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 TH1F* h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 TH1F* h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 TH1F* h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

  TH1F* h_tt_sqrtS_reco_isoPh = new TH1F("h_tt_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);//default
  TH1F* h_tt_sqrtS_reco = new TH1F("h_tt_sqrtS_reco","", n_bins_high, lim_energy_low,lim_energy_high);
  TH1F* h_tt_sqrtS_reco_isoPh_inv = new TH1F("h_tt_sqrtS_reco_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
 h_tt_sqrtS_reco_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_isoPh_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_tt_mass_j1_METProj_reco_sqrt_s_0_750 = new TH1F("h_tt_mass_j1_METProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j1_METProj_reco_sqrt_s_750_2500 = new TH1F("h_tt_mass_j1_METProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j1_METProj_reco_sqrt_s_2500 = new TH1F("h_tt_mass_j1_METProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_tt_mass_j1_METProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_tt_mass_j1_METProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_tt_mass_j1_METProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 
 TH1F* h_tt_mass_j2_METProj_reco_sqrt_s_0_750 = new TH1F("h_tt_mass_j2_METProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j2_METProj_reco_sqrt_s_750_2500 = new TH1F("h_tt_mass_j2_METProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j2_METProj_reco_sqrt_s_2500 = new TH1F("h_tt_mass_j2_METProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_tt_mass_j2_METProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_tt_mass_j2_METProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_tt_mass_j2_METProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_tt_mass_j1_EMissProj_reco_sqrt_s_0_750 = new TH1F("h_tt_mass_j1_EMissProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j1_EMissProj_reco_sqrt_s_750_2500 = new TH1F("h_tt_mass_j1_EMissProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j1_EMissProj_reco_sqrt_s_2500 = new TH1F("h_tt_mass_j1_EMissProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_tt_mass_j1_EMissProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_tt_mass_j1_EMissProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_tt_mass_j1_EMissProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 
 TH1F* h_tt_mass_j2_EMissProj_reco_sqrt_s_0_750 = new TH1F("h_tt_mass_j2_EMissProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j2_EMissProj_reco_sqrt_s_750_2500 = new TH1F("h_tt_mass_j2_EMissProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_tt_mass_j2_EMissProj_reco_sqrt_s_2500 = new TH1F("h_tt_mass_j2_EMissProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_tt_mass_j2_EMissProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_tt_mass_j2_EMissProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_tt_mass_j2_EMissProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_tt_sqrtS_reco_isoPh_EMissCorr = new TH1F("h_tt_sqrtS_reco_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high);
 h_tt_sqrtS_reco_isoPh_EMissCorr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500 = new TH1F("h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_tt_sqrtS_reco_j1_j2_isoPh = new TH1F("h_tt_sqrtS_reco_j1_j2_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);
 h_tt_sqrtS_reco_j1_j2_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 TH1F* h_tt_sqrtS_reco_j1_j2_isoPh_EMissCorr = new TH1F("h_tt_sqrtS_reco_j1_j2_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high);
 h_tt_sqrtS_reco_j1_j2_isoPh_EMissCorr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 
 TH1F* h_tt_E_tot_j1_j2_reco_sqrt_s_2500 = new TH1F("h_tt_E_tot_j1_j2_reco_sqrt_s_2500","", n_bins_high, lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_tt_E_tot_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
 TH1F* h_tt_E_tot_j1_j2_EMiss_reco_sqrt_s_2500 = new TH1F("h_tt_E_tot_j1_j2_EMiss_reco_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_tt_E_tot_j1_j2_EMiss_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");

 TH1F* h_tt_E_tot_isoPh_reco_sqrt_s_2500 = new TH1F("h_tt_E_tot_isoPh_reco_sqrt_s_2500","", n_bins_high, lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_tt_E_tot_isoPh_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{tot} [GeV]");
 TH1F* h_tt_E_tot_isoPh_EMiss_reco_sqrt_s_2500 = new TH1F("h_tt_E_tot_isoPh_EMiss_reco_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_tt_E_tot_isoPh_EMiss_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{tot}+E_{miss}^{jetProj}) [GeV]");
 
 TH1F* h_tt_E_j1_reco_sqrt_s_2500 = new TH1F("h_tt_E_j1_reco_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_tt_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_tt_E_j2_reco_sqrt_s_2500 = new TH1F("h_tt_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_tt_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");
 TH1F* h_tt_E_j1_reco_EMiss_sqrt_s_2500 = new TH1F("h_tt_E_j1_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_tt_E_j1_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_tt_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_tt_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_tt_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");

 TH1F* h_tt_E_j1_min_E_j2_reco_sqrt_s_2500 = new TH1F("h_tt_E_j1_min_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);
 h_tt_E_j1_min_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2}) [GeV]");
 TH1F* h_tt_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_tt_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);
 h_tt_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2}) [GeV]");

 TH1F* h_tt_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500 = new TH1F("h_tt_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_tt_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");
 TH1F* h_tt_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_tt_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_tt_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");

 TH1F* h_tt_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_tt_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_tt_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MET_{jet-proj,reco}-MET_{gen})/MET_{gen}");
 TH1F* h_tt_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_tt_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_tt_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MHT_{jet-proj}-MET_{reco})/MET_{gen}");
 TH1F* h_tt_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_tt_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_tt_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MHT_{jet-proj,reco}-MET_{gen})/MET_{gen}");
 
 TH1F* h_tt_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500 = new TH1F("h_tt_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_tt_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500->GetXaxis()->SetTitle("(E_{miss,jet-proj,reco}-E_{miss,gen}/E_{miss,gen}");
 TH1F* h_tt_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500 = new TH1F("h_tt_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_tt_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500->GetXaxis()->SetTitle("(E^{MHT}_{miss,jet-proj,reco}-E_{miss,gen})/E_{miss,gen}");
 
 TH1F* h_tt_dPhi_EMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_tt_dPhi_EMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_tt_dPhi_EMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_tt_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_tt_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_tt_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_tt_dPhi_MHT_reco_MET_gen_sqrtS_2500 = new TH1F("h_tt_dPhi_MHT_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_tt_dPhi_MHT_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MHT_{reco},MET_{gen}) [#circ]");
 
 TH1F* h_tt_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_tt_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_tt_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#alpha(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_tt_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_tt_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_tt_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#alpha(MET_{jet-proj,reco},MET_{gen}) [#circ]");


  TH1F* h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_750 = new TH1F("h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_750 = new TH1F("h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");

  TH1F* h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_750 = new TH1F("h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_750 = new TH1F("h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");

  TH1F* h_tt_MET_reco_over_E_tot_reco_sqrtS_2500 = new TH1F("h_tt_MET_reco_over_E_tot_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_tt_MET_reco_over_E_tot_reco_sqrtS_2500->GetXaxis()->SetTitle("MET_{reco}/E_{tot}");
  TH1F* h_tt_MHT_reco_over_E_tot_reco_sqrtS_2500 = new TH1F("h_tt_MHT_reco_over_E_tot_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_tt_MHT_reco_over_E_tot_reco_sqrtS_2500->GetXaxis()->SetTitle("MHT_{reco}/E_{tot}");
  TH1F* h_tt_MET_reco_over_E_j1_j2_reco_sqrtS_2500 = new TH1F("h_tt_MET_reco_over_E_j1_j2_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_tt_MET_reco_over_E_j1_j2_reco_sqrtS_2500->GetXaxis()->SetTitle("MET_{reco}/(E_{j1}+E_{j2})");
  TH1F* h_tt_MHT_reco_over_E_j1_j2_reco_sqrtS_2500 = new TH1F("h_tt_MHT_reco_over_E_j1_j2_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_tt_MHT_reco_over_E_j1_j2_reco_sqrtS_2500->GetXaxis()->SetTitle("MHT_{reco}/(E_{j1}+E_{j2})");

  TH1F* h_tt_delta_mass_reco_rj1_rj2_sqrt_s_750 = new TH1F("h_tt_delta_mass_reco_rj1_rj2_sqrt_s_750","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_tt_delta_mass_reco_rj1_rj2_sqrt_s_750->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_tt_delta_mass_reco_rj1_rj2_sqrt_s_750_2500 = new TH1F("h_tt_delta_mass_reco_rj1_rj2_sqrt_s_700_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_tt_delta_mass_reco_rj1_rj2_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_tt_delta_mass_reco_rj1_rj2_sqrt_s_2500 = new TH1F("h_tt_delta_mass_reco_rj1_rj2_sqrt_s_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_tt_delta_mass_reco_rj1_rj2_sqrt_s_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");

  TH1F* h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750 = new TH1F("h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500 = new TH1F("h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_700_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500 = new TH1F("h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");

 std::vector<TH1F*> hist_vec_reco_tt_1D;  
 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_j1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_dPhi_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_dPhi_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_dPhi_j1_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_dTheta_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_dTheta_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_dTheta_j1_j2_reco_sqrt_s_2500);//8
 
 hist_vec_reco_tt_1D.push_back(h_tt_Theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_Theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_Theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_Theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_Theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_Theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_mass_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j1_reco_sqrt_s_2500);//17
 
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_lrz_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_lrz_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_lrz_j1_reco_sqrt_s_2500);//26
 
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_lrz_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_lrz_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_lrz_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_tau21_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_j2_reco_sqrt_s_2500);//35
 
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_lrz_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_lrz_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_tau21_lrz_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_lrz_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_lrz_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_tau32_lrz_j2_reco_sqrt_s_2500);//44
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_j1_reco_sqrt_s_2500);//53
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//62

 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta1_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_j2_reco_sqrt_s_2500);//71
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta2_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//80
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_j1_reco_sqrt_s_2500);//89
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500);//98
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta1_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_j2_reco_sqrt_s_2500);//107
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta2_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500);//116
 
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_j1_reco_sqrt_s_2500);//125

 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//134

 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_j2_reco_sqrt_s_2500);//143

 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//152
 
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_j1_reco_sqrt_s_2500);//161

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500);//170

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_j2_reco_sqrt_s_2500);//179

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500);//188

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_j1_reco_sqrt_s_2500);//197

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//206

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_j2_reco_sqrt_s_2500);//215

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//224

 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_E_theta_j1_reco_sqrt_s_2500);


 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_j2_reco_sqrt_s_2500);//233

 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_D2_1_2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_subjet2_E_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet2_E_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet2_E_j1_reco_sqrt_s_2500);//242

 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_subjet2_E_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet2_E_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet2_E_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_over_jetE_j1_reco_sqrt_s_2500);//251

 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_subjet1_E_over_jetE_j2_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500);//260

 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500);//269

 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500);//278

 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_isoPh);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_isoPh_inv);//281

 hist_vec_reco_tt_1D.push_back(h_tt_mass_j1_METProj_reco_sqrt_s_0_750);//282
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j1_METProj_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j1_METProj_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j2_METProj_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j2_METProj_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j2_METProj_reco_sqrt_s_2500);//287
 
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j1_EMissProj_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j1_EMissProj_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j1_EMissProj_reco_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j2_EMissProj_reco_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j2_EMissProj_reco_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_mass_j2_EMissProj_reco_sqrt_s_2500); //293

 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_isoPh_EMissCorr);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500);//297
 
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500);//303
 
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500);
 
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500);//309
 
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500);//312
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_j1_j2_isoPh);
 hist_vec_reco_tt_1D.push_back(h_tt_sqrtS_reco_j1_j2_isoPh_EMissCorr);//314

 hist_vec_reco_tt_1D.push_back(h_tt_E_tot_j1_j2_reco_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_E_tot_j1_j2_EMiss_reco_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_E_tot_isoPh_reco_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_E_tot_isoPh_EMiss_reco_sqrt_s_2500); 
 hist_vec_reco_tt_1D.push_back(h_tt_E_j1_reco_sqrt_s_2500);//319
 hist_vec_reco_tt_1D.push_back(h_tt_E_j2_reco_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_E_j1_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_E_j2_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_E_j1_min_E_j2_reco_sqrt_s_2500);//323
 hist_vec_reco_tt_1D.push_back(h_tt_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500);

 hist_vec_reco_tt_1D.push_back(h_tt_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500);//327
 hist_vec_reco_tt_1D.push_back(h_tt_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500); 
 hist_vec_reco_tt_1D.push_back(h_tt_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500); 
 hist_vec_reco_tt_1D.push_back(h_tt_dPhi_EMissProj_reco_MET_gen_sqrtS_2500);//332
 hist_vec_reco_tt_1D.push_back(h_tt_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_dPhi_MHT_reco_MET_gen_sqrtS_2500); 
 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_750);//337
 hist_vec_reco_tt_1D.push_back(h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_rj1_sj_CosMinHelicityAngle_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_750);
 hist_vec_reco_tt_1D.push_back(h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500);//342
 hist_vec_reco_tt_1D.push_back(h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_750);
 hist_vec_reco_tt_1D.push_back(h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_rj2_sj_CosMinHelicityAngle_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_750);
 hist_vec_reco_tt_1D.push_back(h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500);//347
 hist_vec_reco_tt_1D.push_back(h_tt_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_MET_reco_over_E_tot_reco_sqrtS_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_MHT_reco_over_E_tot_reco_sqrtS_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_MET_reco_over_E_j1_j2_reco_sqrtS_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_MHT_reco_over_E_j1_j2_reco_sqrtS_2500);//352

 hist_vec_reco_tt_1D.push_back(h_tt_delta_mass_reco_rj1_rj2_sqrt_s_750);
 hist_vec_reco_tt_1D.push_back(h_tt_delta_mass_reco_rj1_rj2_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_delta_mass_reco_rj1_rj2_sqrt_s_2500);//355

 hist_vec_reco_tt_1D.push_back(h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750);
 hist_vec_reco_tt_1D.push_back(h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500);
 hist_vec_reco_tt_1D.push_back(h_tt_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500);//358

 for(unsigned int i=0;i<hist_vec_reco_tt_1D.size();i++){
   hist_vec_reco_tt_1D[i]->Sumw2();
   hist_vec_reco_tt_1D[i]->SetLineColor(kRed);
   hist_vec_reco_tt_1D[i]->SetLineWidth(2);
   hist_vec_reco_tt_1D[i]->GetYaxis()->SetTitle("Events");
 }
 std::vector<TH1F*>hist_vec_tt_parton;
 std::vector<TH1F*>hist_vec_gen_tt_1D;

 TH2F* h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_0_750 = new TH2F("h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_0_750->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500 = new TH2F("h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_2500 = new TH2F("h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_2500->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_tt_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500 = new TH2F("h_tt_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500","", n_bins_high,n_bins_low_METCorr,n_bins_high_METCorr , n_bins_high,n_bins_low_METCorr, n_bins_high_METCorr);
 h_tt_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetXaxis()->SetTitle("gj1 METCorr");
 h_tt_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetYaxis()->SetTitle("gj2 METCorr");

 TH2F* h_tt_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500 = new TH2F("h_tt_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500","", n_bins_high,n_bins_low_METCorr,n_bins_high_METCorr , n_bins_high,n_bins_low_METCorr, n_bins_high_METCorr);
 h_tt_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetXaxis()->SetTitle("rj1 METCorr");
 h_tt_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetYaxis()->SetTitle("rj2 METCorr");  

 std::vector<TH2F*>hist_vec_tt_2DHist;
 hist_vec_tt_2DHist.push_back(h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_0_750);
 hist_vec_tt_2DHist.push_back(h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500);
 hist_vec_tt_2DHist.push_back(h_tt_mass_j1_vs_mass_j2_reco_sqrt_s_2500);
 hist_vec_tt_2DHist.push_back(h_tt_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500);
 hist_vec_tt_2DHist.push_back(h_tt_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500);
 for(unsigned int i=0;i<hist_vec_tt_2DHist.size();i++){
   hist_vec_tt_2DHist[i]->Sumw2();
 }

 usePartonInfo=false;
 fillPartonInfo=false;
 fillGenInfo=false;

 std::vector<TH1F*> hist_vec_reco_tt_1D_reco_vs_gen_selection;  
 std::vector<TH2F*> hist_vec_reco_tt_2D_reco_vs_gen_selection;

 fill_HZ_histograms(file_CLIC_tt, hist_vec_tt_parton, hist_vec_gen_tt_1D, hist_vec_reco_tt_1D, hist_vec_tt_2DHist, hist_vec_reco_tt_1D_reco_vs_gen_selection, hist_vec_reco_tt_2D_reco_vs_gen_selection, usePartonInfo ,xsec_tt,fillPartonInfo,fillGenInfo);
  */


  //NOW _ee_qq

 TH1F* h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");
 h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");
 h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");

 TH1F* h_ee_qq_dPhi_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_dPhi_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_ee_qq_dPhi_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_dPhi_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_ee_qq_dPhi_j1_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_dPhi_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_ee_qq_dPhi_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");
 h_ee_qq_dPhi_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");
 h_ee_qq_dPhi_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");

 TH1F* h_ee_qq_dTheta_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_dTheta_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_ee_qq_dTheta_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_dTheta_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_ee_qq_dTheta_j1_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_dTheta_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_ee_qq_dTheta_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");
 h_ee_qq_dTheta_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");
 h_ee_qq_dTheta_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");

 TH1F* h_ee_qq_Theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_Theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_ee_qq_Theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_Theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_ee_qq_Theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_Theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_ee_qq_Theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#theta(j1)[#circ]");
 h_ee_qq_Theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#theta(j1)[#circ]");
 h_ee_qq_Theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("#theta(j1)[#circ]");

 TH1F* h_ee_qq_Theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_Theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_ee_qq_Theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_Theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_ee_qq_Theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_Theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_ee_qq_Theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 h_ee_qq_Theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 h_ee_qq_Theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 
 TH1F* h_ee_qq_mass_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_mass_j1_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_mass_j1_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_mass_j1_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 
 h_ee_qq_mass_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_ee_qq_mass_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_ee_qq_mass_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");

 TH1F* h_ee_qq_mass_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_mass_j2_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_mass_j2_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_mass_j2_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);

 h_ee_qq_mass_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_ee_qq_mass_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_ee_qq_mass_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_ee_qq_tau21_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_tau21_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau21_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_tau21_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau21_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_tau21_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_ee_qq_tau21_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_ee_qq_tau21_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_ee_qq_tau21_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j1)");

 TH1F* h_ee_qq_tau21_lrz_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_tau21_lrz_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau21_lrz_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_tau21_lrz_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau21_lrz_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_tau21_lrz_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_ee_qq_tau21_lrz_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_ee_qq_tau21_lrz_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_ee_qq_tau21_lrz_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j1)");

 TH1F* h_ee_qq_tau32_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_tau32_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau32_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_tau32_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau32_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_tau32_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_ee_qq_tau32_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_ee_qq_tau32_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_ee_qq_tau32_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j1)");

 TH1F* h_ee_qq_tau32_lrz_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_tau32_lrz_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau32_lrz_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_tau32_lrz_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau32_lrz_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_tau32_lrz_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_ee_qq_tau32_lrz_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_ee_qq_tau32_lrz_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_ee_qq_tau32_lrz_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j1)");

 TH1F* h_ee_qq_tau21_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_tau21_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau21_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_tau21_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau21_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_tau21_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_ee_qq_tau21_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_ee_qq_tau21_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_ee_qq_tau21_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j2)");

 TH1F* h_ee_qq_tau21_lrz_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_tau21_lrz_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau21_lrz_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_tau21_lrz_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau21_lrz_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_tau21_lrz_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_ee_qq_tau21_lrz_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_ee_qq_tau21_lrz_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_ee_qq_tau21_lrz_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j2)");

 TH1F* h_ee_qq_tau32_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_tau32_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau32_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_tau32_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau32_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_tau32_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_ee_qq_tau32_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_ee_qq_tau32_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_ee_qq_tau32_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j2)");

 TH1F* h_ee_qq_tau32_lrz_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_tau32_lrz_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau32_lrz_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_tau32_lrz_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_ee_qq_tau32_lrz_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_tau32_lrz_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_ee_qq_tau32_lrz_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_ee_qq_tau32_lrz_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_ee_qq_tau32_lrz_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 //N2 j1
 TH1F* h_ee_qq_N2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_ee_qq_N2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_ee_qq_N2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");

 TH1F* h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");

 TH1F* h_ee_qq_N2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_ee_qq_N2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_ee_qq_N2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");

 TH1F* h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");

 TH1F* h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");

 TH1F* h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 //N2 j2
 TH1F* h_ee_qq_N2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_ee_qq_N2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_ee_qq_N2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");

 TH1F* h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");

 TH1F* h_ee_qq_N2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_ee_qq_N2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_ee_qq_N2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");

 TH1F* h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");

 TH1F* h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");

 TH1F* h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");



 //N3 j1
 TH1F* h_ee_qq_N3_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_ee_qq_N3_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_ee_qq_N3_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");

 TH1F* h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");

 TH1F* h_ee_qq_N3_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_ee_qq_N3_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_ee_qq_N3_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");

 TH1F* h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");

 TH1F* h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");

 TH1F* h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 //N3 j2
 TH1F* h_ee_qq_N3_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_ee_qq_N3_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_ee_qq_N3_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");

 TH1F* h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");

 TH1F* h_ee_qq_N3_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_ee_qq_N3_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_ee_qq_N3_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");

 TH1F* h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");

 TH1F* h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");

 TH1F* h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");


 //now the C2 and C3 series
 //C2 j1
 TH1F* h_ee_qq_C2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_ee_qq_C2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_ee_qq_C2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");

 TH1F* h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");

 TH1F* h_ee_qq_C2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_ee_qq_C2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_ee_qq_C2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");

 TH1F* h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");

 TH1F* h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");

 TH1F* h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 //C2 j2
 TH1F* h_ee_qq_C2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_ee_qq_C2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_ee_qq_C2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");

 TH1F* h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");

 TH1F* h_ee_qq_C2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_ee_qq_C2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_ee_qq_C2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");

 TH1F* h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");

 TH1F* h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");

 TH1F* h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");



 //C3 j1
 TH1F* h_ee_qq_C3_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_ee_qq_C3_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_ee_qq_C3_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");

 TH1F* h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");

 TH1F* h_ee_qq_C3_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_ee_qq_C3_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_ee_qq_C3_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");

 TH1F* h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");

 TH1F* h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");

 TH1F* h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 //C3 j2
 TH1F* h_ee_qq_C3_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_ee_qq_C3_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_ee_qq_C3_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");

 TH1F* h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");

 TH1F* h_ee_qq_C3_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_ee_qq_C3_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_ee_qq_C3_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");

 TH1F* h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");

 TH1F* h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");

 TH1F* h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");

 //D2 series
 //D2 j1
 TH1F* h_ee_qq_D2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_ee_qq_D2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_ee_qq_D2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");

 TH1F* h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");

 TH1F* h_ee_qq_D2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_ee_qq_D2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_ee_qq_D2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");

 TH1F* h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");

 TH1F* h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");

 TH1F* h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 //D2 j2
 TH1F* h_ee_qq_D2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_ee_qq_D2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_ee_qq_D2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");

 TH1F* h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");

 TH1F* h_ee_qq_D2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_ee_qq_D2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_ee_qq_D2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");

 TH1F* h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");

 TH1F* h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");

 TH1F* h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");




 //D2 (1,2) series
 //D2 (1,2) j1
 TH1F* h_ee_qq_D2_1_2_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_1_2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_ee_qq_D2_1_2_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_1_2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_ee_qq_D2_1_2_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_1_2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_ee_qq_D2_1_2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_ee_qq_D2_1_2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_ee_qq_D2_1_2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");

 TH1F* h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");

 //D2_1_2 (1,2) j2
 TH1F* h_ee_qq_D2_1_2_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_1_2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_ee_qq_D2_1_2_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_1_2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_ee_qq_D2_1_2_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_1_2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_ee_qq_D2_1_2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_ee_qq_D2_1_2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_ee_qq_D2_1_2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");

 TH1F* h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");



 TH1F* h_ee_qq_subjet1_E_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_subjet1_E_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_ee_qq_subjet1_E_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_subjet1_E_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_ee_qq_subjet1_E_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_subjet1_E_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_ee_qq_subjet1_E_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");
 h_ee_qq_subjet1_E_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");
 h_ee_qq_subjet1_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");

 TH1F* h_ee_qq_subjet2_E_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_subjet2_E_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_ee_qq_subjet2_E_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_subjet2_E_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_ee_qq_subjet2_E_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_subjet2_E_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_ee_qq_subjet2_E_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");
 h_ee_qq_subjet2_E_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");
 h_ee_qq_subjet2_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");

 TH1F* h_ee_qq_subjet1_E_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_subjet1_E_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_ee_qq_subjet1_E_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_subjet1_E_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_ee_qq_subjet1_E_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_subjet1_E_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_ee_qq_subjet1_E_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");
 h_ee_qq_subjet1_E_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");
 h_ee_qq_subjet1_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");

 TH1F* h_ee_qq_subjet2_E_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_subjet2_E_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_ee_qq_subjet2_E_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_subjet2_E_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_ee_qq_subjet2_E_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_subjet2_E_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_ee_qq_subjet2_E_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");
 h_ee_qq_subjet2_E_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");
 h_ee_qq_subjet2_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");

 TH1F* h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}/E (j1)");
 h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}/E (j1)");
 h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}/E (j1)");

 TH1F* h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}/E (j2)");
 h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}/E (j2)");
 h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}/E (j2)");


 TH1F* h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1 [#circ]");

 TH1F* h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");


 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_ee_qq_sqrtS_reco_isoPh = new TH1F("h_ee_qq_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);//default
 TH1F* h_ee_qq_sqrtS_reco = new TH1F("h_ee_qq_sqrtS_reco","", n_bins_high, lim_energy_low,lim_energy_high);
 TH1F* h_ee_qq_sqrtS_reco_isoPh_inv = new TH1F("h_ee_qq_sqrtS_reco_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
 h_ee_qq_sqrtS_reco_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_isoPh_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_ee_qq_mass_j1_METProj_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_mass_j1_METProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j1_METProj_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_mass_j1_METProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j1_METProj_reco_sqrt_s_2500 = new TH1F("h_ee_qq_mass_j1_METProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_ee_qq_mass_j1_METProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_ee_qq_mass_j1_METProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_ee_qq_mass_j1_METProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 
 TH1F* h_ee_qq_mass_j2_METProj_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_mass_j2_METProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j2_METProj_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_mass_j2_METProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j2_METProj_reco_sqrt_s_2500 = new TH1F("h_ee_qq_mass_j2_METProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_ee_qq_mass_j2_METProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_ee_qq_mass_j2_METProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_ee_qq_mass_j2_METProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_2500 = new TH1F("h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 
 TH1F* h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_0_750 = new TH1F("h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_750_2500 = new TH1F("h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_2500 = new TH1F("h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_ee_qq_sqrtS_reco_isoPh_EMissCorr = new TH1F("h_ee_qq_sqrtS_reco_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high);
 h_ee_qq_sqrtS_reco_isoPh_EMissCorr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500 = new TH1F("h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_ee_qq_sqrtS_reco_j1_j2_isoPh = new TH1F("h_ee_qq_sqrtS_reco_j1_j2_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);
 h_ee_qq_sqrtS_reco_j1_j2_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 TH1F* h_ee_qq_sqrtS_reco_j1_j2_isoPh_EMissCorr = new TH1F("h_ee_qq_sqrtS_reco_j1_j2_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high);
 h_ee_qq_sqrtS_reco_j1_j2_isoPh_EMissCorr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 
 TH1F* h_ee_qq_E_tot_j1_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_E_tot_j1_j2_reco_sqrt_s_2500","", n_bins_high, lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_ee_qq_E_tot_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
 TH1F* h_ee_qq_E_tot_j1_j2_EMiss_reco_sqrt_s_2500 = new TH1F("h_ee_qq_E_tot_j1_j2_EMiss_reco_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_ee_qq_E_tot_j1_j2_EMiss_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");

 TH1F* h_ee_qq_E_tot_isoPh_reco_sqrt_s_2500 = new TH1F("h_ee_qq_E_tot_isoPh_reco_sqrt_s_2500","", n_bins_high, lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_ee_qq_E_tot_isoPh_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{tot} [GeV]");
 TH1F* h_ee_qq_E_tot_isoPh_EMiss_reco_sqrt_s_2500 = new TH1F("h_ee_qq_E_tot_isoPh_EMiss_reco_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_ee_qq_E_tot_isoPh_EMiss_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{tot}+E_{miss}^{jetProj}) [GeV]");
 
 TH1F* h_ee_qq_E_j1_reco_sqrt_s_2500 = new TH1F("h_ee_qq_E_j1_reco_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_ee_qq_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_ee_qq_E_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_ee_qq_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");
 TH1F* h_ee_qq_E_j1_reco_EMiss_sqrt_s_2500 = new TH1F("h_ee_qq_E_j1_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_ee_qq_E_j1_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_ee_qq_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_ee_qq_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_ee_qq_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");

 TH1F* h_ee_qq_E_j1_min_E_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_E_j1_min_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);
 h_ee_qq_E_j1_min_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2}) [GeV]");
 TH1F* h_ee_qq_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_ee_qq_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);
 h_ee_qq_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2}) [GeV]");

 TH1F* h_ee_qq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500 = new TH1F("h_ee_qq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_ee_qq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");
 TH1F* h_ee_qq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_ee_qq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_ee_qq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");

 TH1F* h_ee_qq_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_ee_qq_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_ee_qq_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MET_{jet-proj,reco}-MET_{gen})/MET_{gen}");
 TH1F* h_ee_qq_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_ee_qq_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_ee_qq_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MHT_{jet-proj}-MET_{reco})/MET_{gen}");
 TH1F* h_ee_qq_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_ee_qq_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_ee_qq_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MHT_{jet-proj,reco}-MET_{gen})/MET_{gen}");
 
 TH1F* h_ee_qq_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500 = new TH1F("h_ee_qq_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_ee_qq_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500->GetXaxis()->SetTitle("(E_{miss,jet-proj,reco}-E_{miss,gen}/E_{miss,gen}");
 TH1F* h_ee_qq_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500 = new TH1F("h_ee_qq_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_ee_qq_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500->GetXaxis()->SetTitle("(E^{MHT}_{miss,jet-proj,reco}-E_{miss,gen})/E_{miss,gen}");
 
 TH1F* h_ee_qq_dPhi_EMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_ee_qq_dPhi_EMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_ee_qq_dPhi_EMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_ee_qq_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_ee_qq_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_ee_qq_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_ee_qq_dPhi_MHT_reco_MET_gen_sqrtS_2500 = new TH1F("h_ee_qq_dPhi_MHT_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_ee_qq_dPhi_MHT_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MHT_{reco},MET_{gen}) [#circ]");
 
 TH1F* h_ee_qq_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_ee_qq_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_ee_qq_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#alpha(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_ee_qq_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_ee_qq_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_ee_qq_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#alpha(MET_{jet-proj,reco},MET_{gen}) [#circ]");


  TH1F* h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_750 = new TH1F("h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750 = new TH1F("h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");

  TH1F* h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_750 = new TH1F("h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750 = new TH1F("h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");

  TH1F* h_ee_qq_MET_reco_over_E_tot_reco_sqrtS_2500 = new TH1F("h_ee_qq_MET_reco_over_E_tot_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_ee_qq_MET_reco_over_E_tot_reco_sqrtS_2500->GetXaxis()->SetTitle("MET_{reco}/E_{tot}");
  TH1F* h_ee_qq_MHT_reco_over_E_tot_reco_sqrtS_2500 = new TH1F("h_ee_qq_MHT_reco_over_E_tot_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_ee_qq_MHT_reco_over_E_tot_reco_sqrtS_2500->GetXaxis()->SetTitle("MHT_{reco}/E_{tot}");
  TH1F* h_ee_qq_MET_reco_over_E_j1_j2_reco_sqrtS_2500 = new TH1F("h_ee_qq_MET_reco_over_E_j1_j2_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_ee_qq_MET_reco_over_E_j1_j2_reco_sqrtS_2500->GetXaxis()->SetTitle("MET_{reco}/(E_{j1}+E_{j2})");
  TH1F* h_ee_qq_MHT_reco_over_E_j1_j2_reco_sqrtS_2500 = new TH1F("h_ee_qq_MHT_reco_over_E_j1_j2_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_ee_qq_MHT_reco_over_E_j1_j2_reco_sqrtS_2500->GetXaxis()->SetTitle("MHT_{reco}/(E_{j1}+E_{j2})");

  TH1F* h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_750 = new TH1F("h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_750","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_750->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_750_2500 = new TH1F("h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_700_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_2500 = new TH1F("h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");

  TH1F* h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750 = new TH1F("h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500 = new TH1F("h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_700_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500 = new TH1F("h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");

 std::vector<TH1F*> hist_vec_reco_ee_qq_1D;  
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_j1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dPhi_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dPhi_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dPhi_j1_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dTheta_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dTheta_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dTheta_j1_j2_reco_sqrt_s_2500);//8
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_Theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_Theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_Theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_Theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_Theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_Theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j1_reco_sqrt_s_2500);//17
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_lrz_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_lrz_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_lrz_j1_reco_sqrt_s_2500);//26
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_lrz_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_lrz_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_lrz_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_j2_reco_sqrt_s_2500);//35
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_lrz_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_lrz_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau21_lrz_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_lrz_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_lrz_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_tau32_lrz_j2_reco_sqrt_s_2500);//44
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_j1_reco_sqrt_s_2500);//53
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//62

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta1_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_j2_reco_sqrt_s_2500);//71
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta2_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//80
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_j1_reco_sqrt_s_2500);//89
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500);//98
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta1_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_j2_reco_sqrt_s_2500);//107
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta2_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500);//116
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_j1_reco_sqrt_s_2500);//125

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//134

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_j2_reco_sqrt_s_2500);//143

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//152
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_j1_reco_sqrt_s_2500);//161

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500);//170

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_j2_reco_sqrt_s_2500);//179

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500);//188

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_j1_reco_sqrt_s_2500);//197

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//206

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_j2_reco_sqrt_s_2500);//215

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//224

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_E_theta_j1_reco_sqrt_s_2500);


 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_j2_reco_sqrt_s_2500);//233

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_D2_1_2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet2_E_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet2_E_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet2_E_j1_reco_sqrt_s_2500);//242

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet2_E_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet2_E_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet2_E_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_over_jetE_j1_reco_sqrt_s_2500);//251

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_subjet1_E_over_jetE_j2_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500);//260

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500);//269

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500);//278

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_isoPh);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_isoPh_inv);//281

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j1_METProj_reco_sqrt_s_0_750);//282
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j1_METProj_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j1_METProj_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j2_METProj_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j2_METProj_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j2_METProj_reco_sqrt_s_2500);//287
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j1_EMissProj_reco_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_mass_j2_EMissProj_reco_sqrt_s_2500); //293

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_isoPh_EMissCorr);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500);//297
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500);//303
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500);
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500);//309
 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500);//312
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_j1_j2_isoPh);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_sqrtS_reco_j1_j2_isoPh_EMissCorr);//314

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_tot_j1_j2_reco_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_tot_j1_j2_EMiss_reco_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_tot_isoPh_reco_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_tot_isoPh_EMiss_reco_sqrt_s_2500); 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_j1_reco_sqrt_s_2500);//319
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_j2_reco_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_j1_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_j2_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_j1_min_E_j2_reco_sqrt_s_2500);//323
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500);

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500);//327
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500); 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500); 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dPhi_EMissProj_reco_MET_gen_sqrtS_2500);//332
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dPhi_MHT_reco_MET_gen_sqrtS_2500); 
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_750);//337
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj1_sj_CosMinHelicityAngle_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500);//342
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj2_sj_CosMinHelicityAngle_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500);//347
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_MET_reco_over_E_tot_reco_sqrtS_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_MHT_reco_over_E_tot_reco_sqrtS_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_MET_reco_over_E_j1_j2_reco_sqrtS_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_MHT_reco_over_E_j1_j2_reco_sqrtS_2500);//352
 

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_2500);//355

 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500);
 hist_vec_reco_ee_qq_1D.push_back(h_ee_qq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500);//358

 for(unsigned int i=0;i<hist_vec_reco_ee_qq_1D.size();i++){
   hist_vec_reco_ee_qq_1D[i]->Sumw2();
   hist_vec_reco_ee_qq_1D[i]->SetLineColor(kBlue);
   hist_vec_reco_ee_qq_1D[i]->SetLineWidth(2);
   hist_vec_reco_ee_qq_1D[i]->GetYaxis()->SetTitle("Events");
 }
 std::vector<TH1F*>hist_vec_ee_qq_parton;
 std::vector<TH1F*>hist_vec_gen_ee_qq_1D;

 TH2F* h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_0_750 = new TH2F("h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_0_750->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500 = new TH2F("h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_2500 = new TH2F("h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_2500->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_ee_qq_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500 = new TH2F("h_ee_qq_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500","", n_bins_high,n_bins_low_METCorr,n_bins_high_METCorr , n_bins_high,n_bins_low_METCorr, n_bins_high_METCorr);
 h_ee_qq_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetXaxis()->SetTitle("gj1 METCorr");
 h_ee_qq_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetYaxis()->SetTitle("gj2 METCorr");

 TH2F* h_ee_qq_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500 = new TH2F("h_ee_qq_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500","", n_bins_high,n_bins_low_METCorr,n_bins_high_METCorr , n_bins_high,n_bins_low_METCorr, n_bins_high_METCorr);
 h_ee_qq_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetXaxis()->SetTitle("rj1 METCorr");
 h_ee_qq_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetYaxis()->SetTitle("rj2 METCorr");

 std::vector<TH2F*>hist_vec_ee_qq_2DHist;
 hist_vec_ee_qq_2DHist.push_back(h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_0_750);
 hist_vec_ee_qq_2DHist.push_back(h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500);
 hist_vec_ee_qq_2DHist.push_back(h_ee_qq_mass_j1_vs_mass_j2_reco_sqrt_s_2500);
 hist_vec_ee_qq_2DHist.push_back(h_ee_qq_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500);
 hist_vec_ee_qq_2DHist.push_back(h_ee_qq_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500);
 for(unsigned int i=0;i<hist_vec_ee_qq_2DHist.size();i++){
   hist_vec_ee_qq_2DHist[i]->Sumw2();
 }

 usePartonInfo=false;
 fillPartonInfo=false;
 fillGenInfo=false;

 std::vector<TH1F*> hist_vec_reco_ee_qq_1D_reco_vs_gen_selection;  
 std::vector<TH2F*> hist_vec_reco_ee_qq_2D_reco_vs_gen_selection;

 std::cout<<"before filling ee_qq"<<std::endl;
 fill_HZ_histograms(file_CLIC_ee_qq, hist_vec_ee_qq_parton, hist_vec_gen_ee_qq_1D, hist_vec_reco_ee_qq_1D, hist_vec_ee_qq_2DHist, hist_vec_reco_ee_qq_1D_reco_vs_gen_selection, hist_vec_reco_ee_qq_2D_reco_vs_gen_selection, usePartonInfo ,xsec_ee_qq,fillPartonInfo,fillGenInfo);
 std::cout<<"after filling ee_qq"<<std::endl;
  //NOW _qqqq

 TH1F* h_qqqq_dAlpha_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_dAlpha_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_qqqq_dAlpha_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_dAlpha_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_qqqq_dAlpha_j1_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_dAlpha_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_qqqq_dAlpha_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");
 h_qqqq_dAlpha_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");
 h_qqqq_dAlpha_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(j1,j2)[#circ]");

 TH1F* h_qqqq_dPhi_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_dPhi_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_qqqq_dPhi_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_dPhi_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_qqqq_dPhi_j1_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_dPhi_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_qqqq_dPhi_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");
 h_qqqq_dPhi_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");
 h_qqqq_dPhi_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#phi(j1,j2)[#circ]");

 TH1F* h_qqqq_dTheta_j1_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_dTheta_j1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_qqqq_dTheta_j1_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_dTheta_j1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_qqqq_dTheta_j1_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_dTheta_j1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_qqqq_dTheta_j1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");
 h_qqqq_dTheta_j1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");
 h_qqqq_dTheta_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#theta(j1,j2)[#circ]");

 TH1F* h_qqqq_Theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_Theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_qqqq_Theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_Theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_qqqq_Theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_Theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_qqqq_Theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#theta(j1)[#circ]");
 h_qqqq_Theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#theta(j1)[#circ]");
 h_qqqq_Theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("#theta(j1)[#circ]");

 TH1F* h_qqqq_Theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_Theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_qqqq_Theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_Theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);
 TH1F* h_qqqq_Theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_Theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle, n_bins_high_dangle);

 h_qqqq_Theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 h_qqqq_Theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 h_qqqq_Theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#theta(j2)[#circ]");
 
 TH1F* h_qqqq_mass_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_mass_j1_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_mass_j1_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_mass_j1_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 
 h_qqqq_mass_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_qqqq_mass_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_qqqq_mass_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");

 TH1F* h_qqqq_mass_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_mass_j2_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_mass_j2_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_mass_j2_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);

 h_qqqq_mass_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_qqqq_mass_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_qqqq_mass_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_qqqq_tau21_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_tau21_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau21_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_tau21_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau21_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_tau21_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_qqqq_tau21_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_qqqq_tau21_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_qqqq_tau21_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j1)");

 TH1F* h_qqqq_tau21_lrz_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_tau21_lrz_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau21_lrz_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_tau21_lrz_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau21_lrz_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_tau21_lrz_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_qqqq_tau21_lrz_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_qqqq_tau21_lrz_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j1)");
 h_qqqq_tau21_lrz_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j1)");

 TH1F* h_qqqq_tau32_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_tau32_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau32_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_tau32_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau32_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_tau32_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_qqqq_tau32_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_qqqq_tau32_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_qqqq_tau32_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j1)");

 TH1F* h_qqqq_tau32_lrz_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_tau32_lrz_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau32_lrz_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_tau32_lrz_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau32_lrz_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_tau32_lrz_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_qqqq_tau32_lrz_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_qqqq_tau32_lrz_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j1)");
 h_qqqq_tau32_lrz_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j1)");

 TH1F* h_qqqq_tau21_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_tau21_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau21_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_tau21_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau21_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_tau21_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_qqqq_tau21_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_qqqq_tau21_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_qqqq_tau21_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j2)");

 TH1F* h_qqqq_tau21_lrz_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_tau21_lrz_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau21_lrz_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_tau21_lrz_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau21_lrz_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_tau21_lrz_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_qqqq_tau21_lrz_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_qqqq_tau21_lrz_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{21}(j2)");
 h_qqqq_tau21_lrz_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{21}(j2)");

 TH1F* h_qqqq_tau32_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_tau32_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau32_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_tau32_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau32_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_tau32_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_qqqq_tau32_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_qqqq_tau32_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_qqqq_tau32_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j2)");

 TH1F* h_qqqq_tau32_lrz_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_tau32_lrz_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau32_lrz_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_tau32_lrz_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);
 TH1F* h_qqqq_tau32_lrz_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_tau32_lrz_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_tau_rat, n_bins_high_tau_rat);

 h_qqqq_tau32_lrz_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_qqqq_tau32_lrz_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 h_qqqq_tau32_lrz_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("tau_{32}(j2)");
 //N2 j1
 TH1F* h_qqqq_N2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_qqqq_N2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_qqqq_N2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");

 TH1F* h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");
 h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j1)");

 TH1F* h_qqqq_N2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_qqqq_N2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_qqqq_N2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");

 TH1F* h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");
 h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j1)");

 TH1F* h_qqqq_N2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_qqqq_N2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_qqqq_N2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");

 TH1F* h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j1)");
 //N2 j2
 TH1F* h_qqqq_N2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_qqqq_N2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_qqqq_N2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");

 TH1F* h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");
 h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(1)}(j2)");

 TH1F* h_qqqq_N2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_qqqq_N2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_qqqq_N2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");

 TH1F* h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");
 h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(2)}(j2)");

 TH1F* h_qqqq_N2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_qqqq_N2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_qqqq_N2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");

 TH1F* h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);
 TH1F* h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N2, n_bins_high_N2);

 h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");
 h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{2}^{(0.5)}(j2)");



 //N3 j1
 TH1F* h_qqqq_N3_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_qqqq_N3_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_qqqq_N3_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");

 TH1F* h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");
 h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j1)");

 TH1F* h_qqqq_N3_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_qqqq_N3_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_qqqq_N3_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");

 TH1F* h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");
 h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j1)");

 TH1F* h_qqqq_N3_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_qqqq_N3_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_qqqq_N3_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");

 TH1F* h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j1)");
 //N3 j2
 TH1F* h_qqqq_N3_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_qqqq_N3_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_qqqq_N3_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");

 TH1F* h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");
 h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(1)}(j2)");

 TH1F* h_qqqq_N3_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_qqqq_N3_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_qqqq_N3_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");

 TH1F* h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");
 h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(2)}(j2)");

 TH1F* h_qqqq_N3_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_qqqq_N3_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_qqqq_N3_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");

 TH1F* h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);
 TH1F* h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_N3, n_bins_high_N3);

 h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");
 h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("N_{3}^{(0.5)}(j2)");


 //now the C2 and C3 series
 //C2 j1
 TH1F* h_qqqq_C2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_qqqq_C2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_qqqq_C2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");

 TH1F* h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");
 h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j1)");

 TH1F* h_qqqq_C2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_qqqq_C2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_qqqq_C2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");

 TH1F* h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");
 h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j1)");

 TH1F* h_qqqq_C2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_qqqq_C2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_qqqq_C2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");

 TH1F* h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j1)");
 //C2 j2
 TH1F* h_qqqq_C2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_qqqq_C2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_qqqq_C2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");

 TH1F* h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");
 h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(1)}(j2)");

 TH1F* h_qqqq_C2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_qqqq_C2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_qqqq_C2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");

 TH1F* h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");
 h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(2)}(j2)");

 TH1F* h_qqqq_C2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_qqqq_C2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_qqqq_C2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");

 TH1F* h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);
 TH1F* h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C2, n_bins_high_C2);

 h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");
 h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{2}^{(0.5)}(j2)");



 //C3 j1
 TH1F* h_qqqq_C3_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_qqqq_C3_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_qqqq_C3_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");

 TH1F* h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");
 h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j1)");

 TH1F* h_qqqq_C3_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_qqqq_C3_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_qqqq_C3_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");

 TH1F* h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");
 h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j1)");

 TH1F* h_qqqq_C3_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_qqqq_C3_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_qqqq_C3_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");

 TH1F* h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j1)");
 //C3 j2
 TH1F* h_qqqq_C3_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_qqqq_C3_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_qqqq_C3_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");

 TH1F* h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");
 h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(1)}(j2)");

 TH1F* h_qqqq_C3_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_qqqq_C3_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_qqqq_C3_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");

 TH1F* h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");
 h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(2)}(j2)");

 TH1F* h_qqqq_C3_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_qqqq_C3_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_qqqq_C3_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");

 TH1F* h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);
 TH1F* h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_C3, n_bins_high_C3);

 h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");
 h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("C_{3}^{(0.5)}(j2)");

 //D2 series
 //D2 j1
 TH1F* h_qqqq_D2_beta1_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta1_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta1_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta1_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta1_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta1_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta1_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_qqqq_D2_beta1_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_qqqq_D2_beta1_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");

 TH1F* h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");
 h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j1)");

 TH1F* h_qqqq_D2_beta2_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta2_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta2_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_qqqq_D2_beta2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_qqqq_D2_beta2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");

 TH1F* h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");
 h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j1)");

 TH1F* h_qqqq_D2_beta0_5_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta0_5_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta0_5_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta0_5_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta0_5_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta0_5_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta0_5_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_qqqq_D2_beta0_5_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_qqqq_D2_beta0_5_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");

 TH1F* h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j1)");
 //D2 j2
 TH1F* h_qqqq_D2_beta1_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta1_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta1_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta1_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta1_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta1_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta1_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_qqqq_D2_beta1_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_qqqq_D2_beta1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");

 TH1F* h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");
 h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1)}(j2)");

 TH1F* h_qqqq_D2_beta2_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta2_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta2_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_qqqq_D2_beta2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_qqqq_D2_beta2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");

 TH1F* h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");
 h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(2)}(j2)");

 TH1F* h_qqqq_D2_beta0_5_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta0_5_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta0_5_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta0_5_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta0_5_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta0_5_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta0_5_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_qqqq_D2_beta0_5_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_qqqq_D2_beta0_5_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");

 TH1F* h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);
 TH1F* h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2, n_bins_high_D2);

 h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");
 h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(0.5)}(j2)");




 //D2 (1,2) series
 //D2 (1,2) j1
 TH1F* h_qqqq_D2_1_2_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_1_2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_qqqq_D2_1_2_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_1_2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_qqqq_D2_1_2_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_1_2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_qqqq_D2_1_2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_qqqq_D2_1_2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_qqqq_D2_1_2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");

 TH1F* h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");
 h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j1)");

 //D2_1_2 (1,2) j2
 TH1F* h_qqqq_D2_1_2_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_1_2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_qqqq_D2_1_2_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_1_2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_qqqq_D2_1_2_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_1_2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_qqqq_D2_1_2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_qqqq_D2_1_2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_qqqq_D2_1_2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");

 TH1F* h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);
 TH1F* h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_D2_1_2, n_bins_high_D2_1_2);

 h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");
 h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("D_{2}^{(1,2)}(j2)");



 TH1F* h_qqqq_subjet1_E_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_subjet1_E_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_qqqq_subjet1_E_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_subjet1_E_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_qqqq_subjet1_E_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_subjet1_E_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_qqqq_subjet1_E_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");
 h_qqqq_subjet1_E_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");
 h_qqqq_subjet1_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}(j1)[GeV]");

 TH1F* h_qqqq_subjet2_E_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_subjet2_E_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_qqqq_subjet2_E_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_subjet2_E_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_qqqq_subjet2_E_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_subjet2_E_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_qqqq_subjet2_E_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");
 h_qqqq_subjet2_E_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");
 h_qqqq_subjet2_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj2}(j1)[GeV]");

 TH1F* h_qqqq_subjet1_E_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_subjet1_E_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_qqqq_subjet1_E_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_subjet1_E_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_qqqq_subjet1_E_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_subjet1_E_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_qqqq_subjet1_E_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");
 h_qqqq_subjet1_E_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");
 h_qqqq_subjet1_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}(j2)[GeV]");

 TH1F* h_qqqq_subjet2_E_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_subjet2_E_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_qqqq_subjet2_E_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_subjet2_E_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 TH1F* h_qqqq_subjet2_E_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_subjet2_E_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE, n_bins_high_subjetE);
 h_qqqq_subjet2_E_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");
 h_qqqq_subjet2_E_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");
 h_qqqq_subjet2_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj2}(j2)[GeV]");

 TH1F* h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}/E (j1)");
 h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}/E (j1)");
 h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}/E (j1)");

 TH1F* h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 TH1F* h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_subjetE_rat, n_bins_high_subjetE_rat);
 h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("E_{sj1}/E (j2)");
 h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("E_{sj1}/E (j2)");
 h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{sj1}/E (j2)");


 TH1F* h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750 = new TH1F("h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1) [#circ]");
 h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j1 [#circ]");

 TH1F* h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750 = new TH1F("h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 TH1F* h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500","", n_bins_high_reco,n_bins_low_dangle_subjet, n_bins_high_dangle_subjet);
 h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");
 h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("#Delta#alpha(sj1,sj2) (j2) [#circ]");


 TH1F* h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 TH1F* h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 TH1F* h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 TH1F* h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_qqqq_sqrtS_reco_isoPh = new TH1F("h_qqqq_sqrtS_reco_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);//default
 TH1F* h_qqqq_sqrtS_reco = new TH1F("h_qqqq_sqrtS_reco","", n_bins_high, lim_energy_low,lim_energy_high);
 TH1F* h_qqqq_sqrtS_reco_isoPh_inv = new TH1F("h_qqqq_sqrtS_reco_isoPh_inv","", n_bins_high, lim_energy_low,lim_energy_high);
 h_qqqq_sqrtS_reco_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_isoPh_inv->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_qqqq_mass_j1_METProj_reco_sqrt_s_0_750 = new TH1F("h_qqqq_mass_j1_METProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j1_METProj_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_mass_j1_METProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j1_METProj_reco_sqrt_s_2500 = new TH1F("h_qqqq_mass_j1_METProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_qqqq_mass_j1_METProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_qqqq_mass_j1_METProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_qqqq_mass_j1_METProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 
 TH1F* h_qqqq_mass_j2_METProj_reco_sqrt_s_0_750 = new TH1F("h_qqqq_mass_j2_METProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j2_METProj_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_mass_j2_METProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j2_METProj_reco_sqrt_s_2500 = new TH1F("h_qqqq_mass_j2_METProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_qqqq_mass_j2_METProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_qqqq_mass_j2_METProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_qqqq_mass_j2_METProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_qqqq_mass_j1_EMissProj_reco_sqrt_s_0_750 = new TH1F("h_qqqq_mass_j1_EMissProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j1_EMissProj_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_mass_j1_EMissProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j1_EMissProj_reco_sqrt_s_2500 = new TH1F("h_qqqq_mass_j1_EMissProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_qqqq_mass_j1_EMissProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_qqqq_mass_j1_EMissProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 h_qqqq_mass_j1_EMissProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j1)[GeV]");
 
 TH1F* h_qqqq_mass_j2_EMissProj_reco_sqrt_s_0_750 = new TH1F("h_qqqq_mass_j2_EMissProj_reco_sqrt_s_0_750","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j2_EMissProj_reco_sqrt_s_750_2500 = new TH1F("h_qqqq_mass_j2_EMissProj_reco_sqrt_s_750_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 TH1F* h_qqqq_mass_j2_EMissProj_reco_sqrt_s_2500 = new TH1F("h_qqqq_mass_j2_EMissProj_reco_sqrt_s_2500","", n_bins_high_reco_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_qqqq_mass_j2_EMissProj_reco_sqrt_s_0_750->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_qqqq_mass_j2_EMissProj_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");
 h_qqqq_mass_j2_EMissProj_reco_sqrt_s_2500->GetXaxis()->SetTitle("mass(j2)[GeV]");

 TH1F* h_qqqq_sqrtS_reco_isoPh_EMissCorr = new TH1F("h_qqqq_sqrtS_reco_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high);
 h_qqqq_sqrtS_reco_isoPh_EMissCorr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_0_750, n_bins_high_sqrt_0_750);

 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_750_2500, n_bins_high_sqrt_750_2500);

 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);
 TH1F* h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500 = new TH1F("h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500","",n_bins_high_reco,n_bins_low_sqrt_2500, n_bins_high_sqrt_2500);

 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");

 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500->GetXaxis()->SetTitle("#sqrt{s} [GeV]");


 TH1F* h_qqqq_sqrtS_reco_j1_j2_isoPh = new TH1F("h_qqqq_sqrtS_reco_j1_j2_isoPh","", n_bins_high, lim_energy_low,lim_energy_high);
 h_qqqq_sqrtS_reco_j1_j2_isoPh->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 TH1F* h_qqqq_sqrtS_reco_j1_j2_isoPh_EMissCorr = new TH1F("h_qqqq_sqrtS_reco_j1_j2_isoPh_EMissCorr","", n_bins_high, lim_energy_low,lim_energy_high);
 h_qqqq_sqrtS_reco_j1_j2_isoPh_EMissCorr->GetXaxis()->SetTitle("#sqrt{s} [GeV]");
 
 TH1F* h_qqqq_E_tot_j1_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_E_tot_j1_j2_reco_sqrt_s_2500","", n_bins_high, lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_qqqq_E_tot_j1_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");
 TH1F* h_qqqq_E_tot_j1_j2_EMiss_reco_sqrt_s_2500 = new TH1F("h_qqqq_E_tot_j1_j2_EMiss_reco_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_qqqq_E_tot_j1_j2_EMiss_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}+E_{j2}) [GeV]");

 TH1F* h_qqqq_E_tot_isoPh_reco_sqrt_s_2500 = new TH1F("h_qqqq_E_tot_isoPh_reco_sqrt_s_2500","", n_bins_high, lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_qqqq_E_tot_isoPh_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{tot} [GeV]");
 TH1F* h_qqqq_E_tot_isoPh_EMiss_reco_sqrt_s_2500 = new TH1F("h_qqqq_E_tot_isoPh_EMiss_reco_sqrt_s_2500","", n_bins_high,lim_jet_energy_low_sqrt_s_2500,lim_energy_high);
 h_qqqq_E_tot_isoPh_EMiss_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{tot}+E_{miss}^{jetProj}) [GeV]");
 
 TH1F* h_qqqq_E_j1_reco_sqrt_s_2500 = new TH1F("h_qqqq_E_j1_reco_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_qqqq_E_j1_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_qqqq_E_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_qqqq_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");
 TH1F* h_qqqq_E_j1_reco_EMiss_sqrt_s_2500 = new TH1F("h_qqqq_E_j1_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_qqqq_E_j1_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j1} [GeV]");
 TH1F* h_qqqq_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_qqqq_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_energy_jet_low,lim_energy_jet_high);
 h_qqqq_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("E_{j2} [GeV]");

 TH1F* h_qqqq_E_j1_min_E_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_E_j1_min_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);
 h_qqqq_E_j1_min_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2}) [GeV]");
 TH1F* h_qqqq_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_qqqq_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_low,lim_delta_energy_jet_high);
 h_qqqq_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2}) [GeV]");

 TH1F* h_qqqq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500 = new TH1F("h_qqqq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_qqqq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");
 TH1F* h_qqqq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500 = new TH1F("h_qqqq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500","", n_bins_high, lim_delta_energy_jet_rel_low,lim_delta_energy_jet_rel_high);
 h_qqqq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500->GetXaxis()->SetTitle("(E_{j1}-E_{j2})/(E_{j1}+E_{j2})");

 TH1F* h_qqqq_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_qqqq_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_qqqq_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MET_{jet-proj,reco}-MET_{gen})/MET_{gen}");
 TH1F* h_qqqq_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_qqqq_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_qqqq_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MHT_{jet-proj}-MET_{reco})/MET_{gen}");
 TH1F* h_qqqq_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500 = new TH1F("h_qqqq_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_qqqq_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("(MHT_{jet-proj,reco}-MET_{gen})/MET_{gen}");
 
 TH1F* h_qqqq_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500 = new TH1F("h_qqqq_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_qqqq_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500->GetXaxis()->SetTitle("(E_{miss,jet-proj,reco}-E_{miss,gen}/E_{miss,gen}");
 TH1F* h_qqqq_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500 = new TH1F("h_qqqq_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500","", n_bins_high, lim_delta_energy_met_rel_low,lim_delta_energy_met_rel_high);
 h_qqqq_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500->GetXaxis()->SetTitle("(E^{MHT}_{miss,jet-proj,reco}-E_{miss,gen})/E_{miss,gen}");
 
 TH1F* h_qqqq_dPhi_EMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_qqqq_dPhi_EMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_qqqq_dPhi_EMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_qqqq_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_qqqq_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_qqqq_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_qqqq_dPhi_MHT_reco_MET_gen_sqrtS_2500 = new TH1F("h_qqqq_dPhi_MHT_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_qqqq_dPhi_MHT_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#Phi(MHT_{reco},MET_{gen}) [#circ]");
 
 TH1F* h_qqqq_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_qqqq_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_qqqq_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#alpha(MET_{jet-proj,reco},MET_{gen}) [#circ]");
 TH1F* h_qqqq_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500 = new TH1F("h_qqqq_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500","", n_bins_high_reco,lim_dalpha_low, lim_dalpha_high);
 h_qqqq_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500->GetXaxis()->SetTitle("#Delta#alpha(MET_{jet-proj,reco},MET_{gen}) [#circ]");


  TH1F* h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_750 = new TH1F("h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750 = new TH1F("h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");

  TH1F* h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_750 = new TH1F("h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500 = new TH1F("h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_2500 = new TH1F("h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("min(cos#theta_{H}(sj,j1))");
  TH1F* h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750 = new TH1F("h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500 = new TH1F("h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_700_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");
  TH1F* h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500 = new TH1F("h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500","", n_bins_high, lim_cosdalpha_low,lim_cosdalpha_high);
  h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500->GetXaxis()->SetTitle("max(cos#theta_{H}(sj,j1))");

  TH1F* h_qqqq_MET_reco_over_E_tot_reco_sqrtS_2500 = new TH1F("h_qqqq_MET_reco_over_E_tot_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_qqqq_MET_reco_over_E_tot_reco_sqrtS_2500->GetXaxis()->SetTitle("MET_{reco}/E_{tot}");
  TH1F* h_qqqq_MHT_reco_over_E_tot_reco_sqrtS_2500 = new TH1F("h_qqqq_MHT_reco_over_E_tot_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_qqqq_MHT_reco_over_E_tot_reco_sqrtS_2500->GetXaxis()->SetTitle("MHT_{reco}/E_{tot}");
  TH1F* h_qqqq_MET_reco_over_E_j1_j2_reco_sqrtS_2500 = new TH1F("h_qqqq_MET_reco_over_E_j1_j2_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_qqqq_MET_reco_over_E_j1_j2_reco_sqrtS_2500->GetXaxis()->SetTitle("MET_{reco}/(E_{j1}+E_{j2})");
  TH1F* h_qqqq_MHT_reco_over_E_j1_j2_reco_sqrtS_2500 = new TH1F("h_qqqq_MHT_reco_over_E_j1_j2_reco_sqrtS_2500","", n_bins_high_reco,0.,lim_delta_energy_jet_rel_high);
  h_qqqq_MHT_reco_over_E_j1_j2_reco_sqrtS_2500->GetXaxis()->SetTitle("MHT_{reco}/(E_{j1}+E_{j2})");

  TH1F* h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_750 = new TH1F("h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_750","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_750->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_750_2500 = new TH1F("h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_700_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_2500 = new TH1F("h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");

  TH1F* h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750 = new TH1F("h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500 = new TH1F("h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_700_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");
  TH1F* h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500 = new TH1F("h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500","", n_bins_low, lim_delta_mass_low,lim_delta_mass_high);
  h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500->GetXaxis()->SetTitle("#Delta mass(rj1,rj2)");

 std::vector<TH1F*> hist_vec_reco_qqqq_1D;  
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_j1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dPhi_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dPhi_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dPhi_j1_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dTheta_j1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dTheta_j1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dTheta_j1_j2_reco_sqrt_s_2500);//8
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_Theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_Theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_Theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_Theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_Theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_Theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j1_reco_sqrt_s_2500);//17
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_lrz_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_lrz_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_lrz_j1_reco_sqrt_s_2500);//26
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_lrz_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_lrz_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_lrz_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_j2_reco_sqrt_s_2500);//35
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_lrz_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_lrz_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau21_lrz_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_lrz_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_lrz_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_tau32_lrz_j2_reco_sqrt_s_2500);//44
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_j1_reco_sqrt_s_2500);//53
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//62

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta1_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_j2_reco_sqrt_s_2500);//71
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta2_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//80
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_j1_reco_sqrt_s_2500);//89
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_E_theta_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_E_theta_j1_reco_sqrt_s_2500);//98
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta1_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_j2_reco_sqrt_s_2500);//107
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta2_E_theta_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_N3_beta0_5_E_theta_j2_reco_sqrt_s_2500);//116
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_j1_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_j1_reco_sqrt_s_2500);//125

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//134

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_j2_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_j2_reco_sqrt_s_2500);//143

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//152
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_j1_reco_sqrt_s_2500);//161

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_E_theta_j1_reco_sqrt_s_2500);//170

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_j2_reco_sqrt_s_2500);//179

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_C3_beta0_5_E_theta_j2_reco_sqrt_s_2500);//188

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_j1_reco_sqrt_s_2500);//197

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_E_theta_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_E_theta_j1_reco_sqrt_s_2500);//206

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta1_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_j2_reco_sqrt_s_2500);//215

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_beta0_5_E_theta_j2_reco_sqrt_s_2500);//224

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_E_theta_j1_reco_sqrt_s_2500);


 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_j2_reco_sqrt_s_2500);//233

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_D2_1_2_E_theta_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet2_E_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet2_E_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet2_E_j1_reco_sqrt_s_2500);//242

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet2_E_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet2_E_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet2_E_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_over_jetE_j1_reco_sqrt_s_2500);//251

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_subjet1_E_over_jetE_j2_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_sj1_sj2_j1_reco_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_sj1_sj2_j2_reco_sqrt_s_2500);//260

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_0_750_gen_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_0_750_parton_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_750_2500_gen_sqrt_s_2500);//269

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_750_2500_parton_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_2500_gen_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_sqrt_s_2500_parton_sqrt_s_2500);//278

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_isoPh);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_isoPh_inv);//281

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j1_METProj_reco_sqrt_s_0_750);//282
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j1_METProj_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j1_METProj_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j2_METProj_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j2_METProj_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j2_METProj_reco_sqrt_s_2500);//287
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j1_EMissProj_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j1_EMissProj_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j1_EMissProj_reco_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j2_EMissProj_reco_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j2_EMissProj_reco_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_mass_j2_EMissProj_reco_sqrt_s_2500); //293

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_isoPh_EMissCorr);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_gen_EMissCorr_sqrt_s_2500);//297
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_0_750_parton_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_gen_EMissCorr_sqrt_s_2500);//303
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_750_2500_parton_sqrt_s_2500);
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_gen_EMissCorr_sqrt_s_2500);//309
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_0_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_EMissCorr_sqrt_s_2500_parton_sqrt_s_2500);//312
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_j1_j2_isoPh);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_sqrtS_reco_j1_j2_isoPh_EMissCorr);//314

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_tot_j1_j2_reco_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_tot_j1_j2_EMiss_reco_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_tot_isoPh_reco_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_tot_isoPh_EMiss_reco_sqrt_s_2500); 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_j1_reco_sqrt_s_2500);//319
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_j2_reco_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_j1_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_j2_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_j1_min_E_j2_reco_sqrt_s_2500);//323
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_j1_min_E_j2_reco_EMiss_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_E_j1_min_E_j2_over_E_j1_plus_E_j2_reco_EMiss_sqrt_s_2500);

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_EMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500);//327
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_MHMissProj_reco_MET_gen_over_MET_gen_sqrtS_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_MHT_reco_MET_gen_over_MET_gen_sqrtS_2500); 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_EMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_MHMiss_E_reco_MET_E_gen_over_MET_E_gen_sqrtS_2500); 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dPhi_EMissProj_reco_MET_gen_sqrtS_2500);//332
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dPhi_MHMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dPhi_MHT_reco_MET_gen_sqrtS_2500); 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_EMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_dAlpha_MHMissProj_reco_MET_gen_sqrtS_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_750);//337
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj1_sj_CosMinHelicityAngle_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj1_sj_CosMaxHelicityAngle_sqrt_s_2500);//342
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj2_sj_CosMinHelicityAngle_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_750_2500);//347
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_rj2_sj_CosMaxHelicityAngle_sqrt_s_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_MET_reco_over_E_tot_reco_sqrtS_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_MHT_reco_over_E_tot_reco_sqrtS_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_MET_reco_over_E_j1_j2_reco_sqrtS_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_MHT_reco_over_E_j1_j2_reco_sqrtS_2500);//352
 
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_2500);//355

 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_750_2500);
 hist_vec_reco_qqqq_1D.push_back(h_qqqq_delta_mass_EMiss_reco_rj1_rj2_sqrt_s_2500);//358

 for(unsigned int i=0;i<hist_vec_reco_qqqq_1D.size();i++){
   hist_vec_reco_qqqq_1D[i]->Sumw2();
   hist_vec_reco_qqqq_1D[i]->SetLineColor(kGreen-2);
   hist_vec_reco_qqqq_1D[i]->SetLineWidth(2);
   hist_vec_reco_qqqq_1D[i]->GetYaxis()->SetTitle("Events");
 }
 std::vector<TH1F*>hist_vec_qqqq_parton;
 std::vector<TH1F*>hist_vec_gen_qqqq_1D;

 TH2F* h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_0_750 = new TH2F("h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_0_750","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_0_750->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_0_750->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500 = new TH2F("h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_2500 = new TH2F("h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_2500","", n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass, n_bins_high_gen_mass,n_bins_low_jetmass, n_bins_high_jetmass);
 h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_2500->GetXaxis()->SetTitle("jet1 mass [GeV]");
 h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_2500->GetYaxis()->SetTitle("jet2 mass [GeV]");

 TH2F* h_qqqq_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500 = new TH2F("h_qqqq_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500","", n_bins_high,n_bins_low_METCorr,n_bins_high_METCorr , n_bins_high,n_bins_low_METCorr, n_bins_high_METCorr);
 h_qqqq_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetXaxis()->SetTitle("gj1 METCorr");
 h_qqqq_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetYaxis()->SetTitle("gj2 METCorr");

 TH2F* h_qqqq_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500 = new TH2F("h_qqqq_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500","", n_bins_high,n_bins_low_METCorr,n_bins_high_METCorr , n_bins_high,n_bins_low_METCorr, n_bins_high_METCorr);
 h_qqqq_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetXaxis()->SetTitle("rj1 METCorr");
 h_qqqq_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500->GetYaxis()->SetTitle("rj2 METCorr");

 std::vector<TH2F*>hist_vec_qqqq_2DHist;
 hist_vec_qqqq_2DHist.push_back(h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_0_750);
 hist_vec_qqqq_2DHist.push_back(h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_750_2500);
 hist_vec_qqqq_2DHist.push_back(h_qqqq_mass_j1_vs_mass_j2_reco_sqrt_s_2500);
 hist_vec_qqqq_2DHist.push_back(h_qqqq_gen_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500);
 hist_vec_qqqq_2DHist.push_back(h_qqqq_reco_j1_METCorr_vs_j2_METCorr_part_sqrt_s_2500);
 for(unsigned int i=0;i<hist_vec_qqqq_2DHist.size();i++){
   hist_vec_qqqq_2DHist[i]->Sumw2();
 }

 usePartonInfo=false;
 fillPartonInfo=false;
 fillGenInfo=false;

 std::vector<TH1F*> hist_vec_reco_qqqq_1D_reco_vs_gen_selection;  
 std::vector<TH2F*> hist_vec_reco_qqqq_2D_reco_vs_gen_selection;
 std::cout<<"before filling qqqq"<<std::endl;
 fill_HZ_histograms(file_CLIC_qqqq, hist_vec_qqqq_parton, hist_vec_gen_qqqq_1D, hist_vec_reco_qqqq_1D, hist_vec_qqqq_2DHist,hist_vec_reco_qqqq_1D_reco_vs_gen_selection, hist_vec_reco_qqqq_2D_reco_vs_gen_selection,  usePartonInfo ,xsec_qqqq,fillPartonInfo,fillGenInfo);
 std::cout<<"after filling qqqq"<<std::endl;
 
 std::cout<<"total numbers HZ/ee_qq/qqqq "<<h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_2500->Integral(0,h_HZ_delta_mass_reco_rj1_rj2_sqrt_s_2500->GetNbinsX()+1)<<"/"<<h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_2500->Integral(0,h_ee_qq_delta_mass_reco_rj1_rj2_sqrt_s_2500->GetNbinsX()+1)<<"/"<<h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_2500->Integral(0,h_qqqq_delta_mass_reco_rj1_rj2_sqrt_s_2500->GetNbinsX()+1)<<std::endl;
  
  file_histogram->Write();
  file_histogram->Close();

}
