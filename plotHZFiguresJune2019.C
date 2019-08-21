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
#include "TLegendEntry.h"
#include "TGraphErrors.h"
#include "TEfficiency.h"
#include "TGaxis.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TColor.h"
#include <vector>
#include "TGraphAsymmErrors.h"


Color_t color1 = kBlack;
Color_t color2 = kRed-7;
Color_t color3 = kBlue;
Color_t color4 = kGreen+2;
Color_t color5 = kCyan+1;
Color_t color6 = kRed+2;
Color_t color7 = kOrange;
Color_t color8 = kViolet+2;


void CLICdpLabel(std::string status);
void CLICdpLabel(Double_t x=0.20, Double_t y=0.86, std::string status="work in progress", Color_t color=kBlack);

void CLICdpStyle()
{
  gROOT->SetStyle("Plain"); /*Default white background for all plots*/

  TStyle *clicdpStyle = gROOT->GetStyle("Plain");

  /* set bkg color of all to kWhite: white, but not 0*/
  clicdpStyle->SetCanvasColor(kWhite);
  clicdpStyle->SetFrameFillColor(kWhite);
  clicdpStyle->SetStatColor(kWhite);
  clicdpStyle->SetPadColor(kWhite);
  clicdpStyle->SetFillColor(10);
  clicdpStyle->SetTitleFillColor(kWhite);
  
  
   /* SetPaperSize wants width & height in cm: A4 is 20,26 & US is 20,24*/
   clicdpStyle->SetPaperSize(20, 26); 
   /* No yellow border around histogram*/
   clicdpStyle->SetDrawBorder(0);
   /* remove border of canvas*/
   clicdpStyle->SetCanvasBorderMode(0);
   /* remove border of pads*/
   clicdpStyle->SetPadBorderMode(0);
   clicdpStyle->SetFrameBorderMode(0);
   clicdpStyle->SetLegendBorderSize(0);
  
   /* default text size*/
   clicdpStyle->SetTextSize(0.05);
   clicdpStyle->SetTitleSize(0.07,"xyz");
   clicdpStyle->SetLabelSize(0.06,"xyz");
   /* title offset: distance between given text and axis, here x,y,z*/
   clicdpStyle->SetLabelOffset(0.015,"xyz");
   clicdpStyle->SetTitleOffset(1.2,"yz"); //equivalent to: clicdpStyle->SetTitleYOffset(1.2);
   clicdpStyle->SetTitleOffset(1.0,"x");



   /* Use visible font for all text*/
   int font = 42; 
   clicdpStyle->SetTitleFont(font);
   clicdpStyle->SetTitleFontSize(0.06);
   clicdpStyle->SetStatFont(font);
   clicdpStyle->SetStatFontSize(0.07);
   clicdpStyle->SetTextFont(font);
   clicdpStyle->SetLabelFont(font,"xyz");
   clicdpStyle->SetTitleFont(font,"xyz");
   clicdpStyle->SetTitleBorderSize(0);
   clicdpStyle->SetStatBorderSize(1);
   //clicdpStyle->SetLegendFont(font);

   /* big marker points*/
   clicdpStyle->SetMarkerStyle(1);
   clicdpStyle->SetLineWidth(2);  
   clicdpStyle->SetMarkerSize(1.2);
   /*set palette in 2d histogram to nice and colorful one*/
   clicdpStyle->SetPalette(1,0); 

   /*No title for histograms*/
   clicdpStyle->SetOptTitle(0);
   /* show the errors on the stat box */
   clicdpStyle->SetOptStat(0); 
   /* show errors on fitted parameters*/
   clicdpStyle->SetOptFit(0); 
   /* number of decimals used for errors*/
   clicdpStyle->SetEndErrorSize(5);   

   /* set line width to 2 by default so that histograms are visible when printed small
      idea: emphasize the data, not the frame around*/
   clicdpStyle->SetHistLineWidth(2);
   clicdpStyle->SetFrameLineWidth(2);
   clicdpStyle->SetFuncWidth(2);
   clicdpStyle->SetHistLineColor(kBlack);
   clicdpStyle->SetFuncColor(kRed);
   clicdpStyle->SetLabelColor(kBlack,"xyz");

   //set the margins
   clicdpStyle->SetPadBottomMargin(0.18);
   clicdpStyle->SetPadTopMargin(0.08);
   clicdpStyle->SetPadRightMargin(0.08);
   clicdpStyle->SetPadLeftMargin(0.17);
   
   //set the number of divisions to show
   clicdpStyle->SetNdivisions(506, "xy");
   
   //turn off xy grids
   clicdpStyle->SetPadGridX(0);
   clicdpStyle->SetPadGridY(0);
   
   //set the tick mark style
   clicdpStyle->SetPadTickX(1);
   clicdpStyle->SetPadTickY(1);

   clicdpStyle->SetCanvasDefW(800);
   clicdpStyle->SetCanvasDefH(700);

   gROOT->ForceStyle();
}

/***********************************************************************/
/*                                                                     */
/*       draw legend for 2 histos on the same plot                     */
/*                                                                     */
/***********************************************************************/
void Draw2Legend(TH1 *histo1, 
                 TH1 *histo2,
                 const Char_t *label1, 
                 const Char_t *label2,
                 const Char_t *header="")
{

  Float_t max1 = histo1->GetMaximum();
  Float_t max2 = histo2->GetMaximum();
  if (max1 >= max2) histo1->SetMaximum(max1 * 1.3);
  else histo1->SetMaximum(max2 * 1.3);

  TLegend *legend = new TLegend(0.1436782,0.8072034,0.6408046,0.9851695,NULL,"brNDC");
  legend->SetTextAlign(22);
  legend->SetTextSize(0.1);
  legend->SetTextSize(0.06);
  legend->SetTextFont(42);
  
  TLegendEntry* entry1 = legend->AddEntry(histo1,label1,"LPF");  
  entry1->SetTextColor(histo1->GetLineColor());
  
  TLegendEntry* entry2 = legend->AddEntry(histo2,label2,"LPF");
  entry2->SetTextColor(histo2->GetLineColor());

  if (std::string(header) != "") legend->SetHeader(header); 

  legend->SetFillColor(kWhite); 
  legend->Draw();

  gPad->Update();
}

/***********************************************************************/
/*                                                                     */
/*       draw legend for 7 histos on the same plot                     */
/*                                                                     */
/***********************************************************************/
void Draw7Legend(TGraph *histo1, 
                 TGraph *histo2,
                 TGraph *histo3,
                 TGraph *histo4,
                 TGraph *histo5,
                 TGraph *histo6,
                 TGraph *histo7,
		 const Char_t *label1, 
                 const Char_t *label2,
		 const Char_t *label3, 
                 const Char_t *label4,
		 const Char_t *label5, 
                 const Char_t *label6,
		 const Char_t *label7, 
                 const Char_t *header="")
{

  TLegend *legend = new TLegend(0.006231824,0.3648567,0.8683008,0.9190962,NULL,"brNDC");
  legend->SetTextAlign(12);
  legend->SetTextSize(0.1);
  legend->SetTextSize(0.1);
  legend->SetTextFont(42);
  
  TLegendEntry* entry1 = legend->AddEntry(histo1,label1,"L");  
  entry1->SetTextColor(histo1->GetLineColor());
  
  TLegendEntry* entry2 = legend->AddEntry(histo2,label2,"L");
  entry2->SetTextColor(histo2->GetLineColor());

  TLegendEntry* entry3 = legend->AddEntry(histo3,label3,"L");
  entry3->SetTextColor(histo3->GetLineColor());

  TLegendEntry* entry4 = legend->AddEntry(histo4,label4,"L");
  entry4->SetTextColor(histo4->GetLineColor());

  TLegendEntry* entry5 = legend->AddEntry(histo5,label5,"L");
  entry5->SetTextColor(histo5->GetLineColor());

  TLegendEntry* entry6 = legend->AddEntry(histo6,label6,"L");
  entry6->SetTextColor(histo6->GetLineColor());

  TLegendEntry* entry7 = legend->AddEntry(histo7,label7,"L");
  entry7->SetTextColor(histo7->GetLineColor());

  if (std::string(header) != "") legend->SetHeader(header); 

  legend->SetFillColor(kWhite); 
  legend->Draw();

  gPad->Update();
}

/***********************************************************************/
/*                                                                     */
/*       draw statistics boxes for 2 histos                            */
/*       (needs h2->Draw("sames"))                                     */
/*                                                                     */
/***********************************************************************/
bool Draw2StatsBoxes(TH1 *histo1, TH1 *histo2)
{
  /*Get the stats box and set the right colors*/
  TPaveStats *statsh1 = (TPaveStats*)histo1->GetListOfFunctions()->FindObject("stats");
  if (statsh1 == NULL)
    {
      std::cout<<"\n Could not find statistics box for histo " << histo1->GetName()<<std::endl;
      return false;
    }
  statsh1->SetLineColor(histo1->GetLineColor());
  statsh1->SetTextColor(histo1->GetLineColor());
  statsh1->SetX1NDC(0.72);
  statsh1->SetY1NDC(0.68);
  statsh1->SetX2NDC(0.92);
  statsh1->SetY2NDC(0.92);

  TPaveStats *statsh2 = (TPaveStats*)histo2->GetListOfFunctions()->FindObject("stats");
  if (statsh2 == NULL)
    {
      std::cout<<"\n Could not find statistics box for histo " << histo2->GetName()<<std::endl;
      std::cout<<"You need to call: histo->Draw(\"sames\")"<<std::endl;
      return false;
    }
  statsh2->SetLineColor(histo2->GetLineColor());
  statsh2->SetTextColor(histo2->GetLineColor());

  /*need to move the stats box below the previous one, as by default they are stacked*/
  statsh2->SetX1NDC(0.72);
  statsh2->SetY1NDC(0.44);
  statsh2->SetX2NDC(0.92);
  statsh2->SetY2NDC(0.67);
  statsh2->Draw();
  
  return true;
}


void CLICdpLabel(std::string status){
  CLICdpLabel(0.20,0.86,status,kBlack);
}//end CLICdpLabel



// void CLICdpLabel(Double_t x=0.20, Double_t y=0.86, std::string status="work in progress", Color_t color=kBlack){

//   TLatex l; 
//   l.SetNDC();
//   l.SetTextFont(42);
//   l.SetTextColor(color);
//   l.SetTextSize(0.05);

//   std::string label = std::string("CLICdp ");
//   l.DrawLatex(x,y,label.c_str());

//   TLatex l2; 
//   l2.SetNDC();
//   l2.SetTextFont(42);
//   l2.SetTextColor(color);
//   l2.SetTextSize(0.035);

//   double dy = l2.GetTextSize()+0.005;

//   l2.DrawLatex(x,y-dy,status.c_str());

//   return;
// }//end CLICdpLabel


void CLICdpLabel(Double_t x, Double_t y, std::string status, Color_t color){

  TLatex l; 
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextColor(kBlack);
  l.SetTextSize(0.05);

  std::string label = std::string("CLICdp");
  l.DrawLatex(x,y,label.c_str());

  TLatex l2; 
  l2.SetNDC();
  l2.SetTextFont(42);
  l2.SetTextColor(color);
  l2.SetTextSize(0.035);

  //double dx = l2.GetTextSize()+0.005;
  double dx = 0.4*l.GetTextSize()*label.length();

  l2.DrawLatex(x+dx,y,status.c_str());

  return;
}//end CLICdpLabel


TCanvas* setUpperCanvas(const char* canvas_name) {
  TCanvas* c1= new TCanvas(canvas_name,canvas_name,0.,0.,800.,700.);
  c1->cd();
  gStyle->SetOptStat(0);
  //c1->Range(-186.894,-0.873515,1682.046,6.114605);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  //c1->SetGridx();
  //c1->SetGridy();
  //c1->SetRightMargin(0.0172);
  //c1->SetTopMargin(0.055);
  //c1->SetBottomMargin(0.138);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);
  
  return c1;
}

void drawFancyPlots(){


  //TLegend *leg_label_CLICdp = new TLegend(0.20,0.90,0.40,0.94, "CLICdp");
  //lewwlabel_CLICdp->SetBorderSize(0);
 
  CLICdpStyle();

  std::cout<<"do i get anywhere here"<<std::endl;

  TLatex* l = new TLatex();
  l->SetNDC();
  l->SetTextFont(42);
  l->SetTextColor(kBlack);
  l->SetTextSize(0.045);


  TLatex* lp = new TLatex();
  lp->SetNDC();
  lp->SetTextFont(42);
  lp->SetTextColor(kBlack);
  lp->SetTextSize(0.045);
  std::string label("CLICdp Work in Progress");
  double x=0.17, y=0.93;
  std::string label2("L=4 ab^{-1},e^{-} pol -80%");
  std::string label3("L=1 ab^{-1},e^{-} pol +80%");
  double x2=0.62, y2=0.93;
  
  TText* leg_label_CLICdp= new TText();
  leg_label_CLICdp->SetTextSize(0.045);
  leg_label_CLICdp->SetTextFont(42);

  bool draw_kappa_plots=true;

  std::cout<<"do i get anywhere here"<<std::endl;

  TFile* file_polm80_hzqq_SignalHistos_=TFile::Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/MVTrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth4NTrees400AdaBoostBeta020.root");  
  TFile* file_polm80_ee_qq_mqq_1TeV_BGHistos_=TFile::Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/MVTrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth4NTrees400AdaBoostBeta020.root"); 
  TFile* file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_=TFile::Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/MVTrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth4NTrees400AdaBoostBeta020.root"); 
  TFile* file_polm80_ee_qqqqqq_BGHistos_=TFile::Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/MVTrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth4NTrees400AdaBoostBeta020.root"); 

  TH1F* h_hzqq_BDT_polm80_jet1_mass=(TH1F*)file_polm80_hzqq_SignalHistos_->Get("0.35/h_jet1_mass");
  h_hzqq_BDT_polm80_jet1_mass->GetXaxis()->SetTitle("jet1 mass [GeV]");
  h_hzqq_BDT_polm80_jet1_mass->GetYaxis()->SetTitle("Events");
  h_hzqq_BDT_polm80_jet1_mass->SetLineWidth(3);
  TH1F* h_ee_qq_BDT_polm80_jet1_mass=(TH1F*)file_polm80_ee_qq_mqq_1TeV_BGHistos_->Get("0.35/h_jet1_mass");
  h_ee_qq_BDT_polm80_jet1_mass->GetXaxis()->SetTitle("jet1 mass [GeV]");
  h_ee_qq_BDT_polm80_jet1_mass->GetYaxis()->SetTitle("Events");
  h_ee_qq_BDT_polm80_jet1_mass->SetFillColor(kBlue);
  h_ee_qq_BDT_polm80_jet1_mass->SetLineColor(kBlue);
  TH1F* h_ee_qqqq_BDT_polm80_jet1_mass=(TH1F*)file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_->Get("0.35/h_jet1_mass");
  h_ee_qqqq_BDT_polm80_jet1_mass->GetXaxis()->SetTitle("jet1 mass [GeV]");
  h_ee_qqqq_BDT_polm80_jet1_mass->GetYaxis()->SetTitle("Events");
  h_ee_qqqq_BDT_polm80_jet1_mass->SetFillColor(kRed);
  h_ee_qqqq_BDT_polm80_jet1_mass->SetLineColor(kRed);
  TH1F* h_ee_qqqqqq_BDT_polm80_jet1_mass=(TH1F*)file_polm80_ee_qqqqqq_BGHistos_->Get("0.35/h_jet1_mass");
  h_ee_qqqqqq_BDT_polm80_jet1_mass->GetXaxis()->SetTitle("jet1 mass [GeV]");
  h_ee_qqqqqq_BDT_polm80_jet1_mass->GetYaxis()->SetTitle("Events");
  h_ee_qqqqqq_BDT_polm80_jet1_mass->SetFillColor(kGreen-2);
  h_ee_qqqqqq_BDT_polm80_jet1_mass->SetLineColor(kGreen-2);

  std::cout<<"do i get anywhere here"<<std::endl;

  THStack* hzqq_sig_BG_Stack_polm80= new THStack("hzqq_sig_BG_Stack_polm80", "");
  hzqq_sig_BG_Stack_polm80->Add(h_ee_qqqqqq_BDT_polm80_jet1_mass);
  hzqq_sig_BG_Stack_polm80->Add(h_ee_qqqq_BDT_polm80_jet1_mass);
  hzqq_sig_BG_Stack_polm80->Add(h_ee_qq_BDT_polm80_jet1_mass);
  hzqq_sig_BG_Stack_polm80->Add(h_hzqq_BDT_polm80_jet1_mass);

  TCanvas *canvas_h_BDT_Signal_polm80_thstack = setUpperCanvas("canvas_h_BDT_Signal_polm80_thstack");
  canvas_h_BDT_Signal_polm80_thstack->cd();
  //h_hzqq_BDT_polm80_jet1_mass->
  hzqq_sig_BG_Stack_polm80->Draw("hist");
  hzqq_sig_BG_Stack_polm80->GetXaxis()->SetTitle("jet1 mass [GeV]");
  hzqq_sig_BG_Stack_polm80->GetYaxis()->SetTitle("Events");
  canvas_h_BDT_Signal_polm80_thstack->Modified();

  TLegend* leg_h_BDT_Overview_polm80_test=new TLegend(0.25,0.65,0.65,0.90);
  //TLegend* leg_h_BDT_Overview_polm80_test=new TLegend(0.25,0.75,0.65,0.90);
  leg_h_BDT_Overview_polm80_test->SetBorderSize(0);
  leg_h_BDT_Overview_polm80_test->SetTextAlign(12);
  leg_h_BDT_Overview_polm80_test->SetTextSize(0.050);
  leg_h_BDT_Overview_polm80_test->SetTextFont(42);
  leg_h_BDT_Overview_polm80_test->SetMargin(0.15);
  leg_h_BDT_Overview_polm80_test->SetLineColor(1);
  leg_h_BDT_Overview_polm80_test->SetLineStyle(1);
  leg_h_BDT_Overview_polm80_test->SetLineWidth(1);
  leg_h_BDT_Overview_polm80_test->SetFillColor(0);
  //leg_h_BDT_Overview_polm80_test->SetFillStyle(1001);
  leg_h_BDT_Overview_polm80_test->SetFillStyle(0);
  leg_h_BDT_Overview_polm80_test->SetHeader("#sqrt{s}>2500 GeV");
  leg_h_BDT_Overview_polm80_test->AddEntry(h_hzqq_BDT_polm80_jet1_mass,"HZ");
  leg_h_BDT_Overview_polm80_test->AddEntry(h_ee_qq_BDT_polm80_jet1_mass,"ee#rightarrow qq");
  leg_h_BDT_Overview_polm80_test->AddEntry(h_ee_qqqq_BDT_polm80_jet1_mass,"ee#rightarrow qqqq");
  leg_h_BDT_Overview_polm80_test->AddEntry(h_ee_qqqqqq_BDT_polm80_jet1_mass,"ee#rightarrow qqqqqq");
  leg_h_BDT_Overview_polm80_test->Draw();

  l->DrawLatex(x,y,label.c_str());
  l->DrawLatex(x2,y2,label2.c_str());

  TFile* file_polp80_hzqq_SignalHistos_=TFile::Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/MVTrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth4NTrees400AdaBoostBeta020.root");  
  TFile* file_polp80_ee_qq_mqq_1TeV_BGHistos_=TFile::Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/MVTrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth4NTrees400AdaBoostBeta020.root"); 
  TFile* file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_=TFile::Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/MVTrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth4NTrees400AdaBoostBeta020.root"); 
  TFile* file_polp80_ee_qqqqqq_BGHistos_=TFile::Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/MVTrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth4NTrees400AdaBoostBeta020.root"); 


  TH1F* h_hzqq_BDT_polp80_jet1_mass=(TH1F*)file_polp80_hzqq_SignalHistos_->Get("0.2/h_jet1_mass");
  h_hzqq_BDT_polp80_jet1_mass->GetXaxis()->SetTitle("jet1 mass [GeV]");
  h_hzqq_BDT_polp80_jet1_mass->GetYaxis()->SetTitle("Events");
  h_hzqq_BDT_polp80_jet1_mass->SetLineWidth(3);
  TH1F* h_ee_qq_BDT_polp80_jet1_mass=(TH1F*)file_polp80_ee_qq_mqq_1TeV_BGHistos_->Get("0.2/h_jet1_mass");
  h_ee_qq_BDT_polp80_jet1_mass->GetXaxis()->SetTitle("jet1 mass [GeV]");
  h_ee_qq_BDT_polp80_jet1_mass->GetYaxis()->SetTitle("Events");
  h_ee_qq_BDT_polp80_jet1_mass->SetFillColor(kBlue);
  h_ee_qq_BDT_polp80_jet1_mass->SetLineColor(kBlue);
  TH1F* h_ee_qqqq_BDT_polp80_jet1_mass=(TH1F*)file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_->Get("0.2/h_jet1_mass");
  h_ee_qqqq_BDT_polp80_jet1_mass->GetXaxis()->SetTitle("jet1 mass [GeV]");
  h_ee_qqqq_BDT_polp80_jet1_mass->GetYaxis()->SetTitle("Events");
  h_ee_qqqq_BDT_polp80_jet1_mass->SetFillColor(kRed);
  h_ee_qqqq_BDT_polp80_jet1_mass->SetLineColor(kRed);
  TH1F* h_ee_qqqqqq_BDT_polp80_jet1_mass=(TH1F*)file_polp80_ee_qqqqqq_BGHistos_->Get("0.2/h_jet1_mass");
  h_ee_qqqqqq_BDT_polp80_jet1_mass->GetXaxis()->SetTitle("jet1 mass [GeV]");
  h_ee_qqqqqq_BDT_polp80_jet1_mass->GetYaxis()->SetTitle("Events");
  h_ee_qqqqqq_BDT_polp80_jet1_mass->SetFillColor(kGreen-2);
  h_ee_qqqqqq_BDT_polp80_jet1_mass->SetLineColor(kGreen-2);

  THStack* hzqq_sig_BG_Stack_polp80= new THStack("hzqq_sig_BG_Stack_polp80", "");
  hzqq_sig_BG_Stack_polp80->Add(h_ee_qqqqqq_BDT_polp80_jet1_mass);
  hzqq_sig_BG_Stack_polp80->Add(h_ee_qqqq_BDT_polp80_jet1_mass);
  hzqq_sig_BG_Stack_polp80->Add(h_ee_qq_BDT_polp80_jet1_mass);
  hzqq_sig_BG_Stack_polp80->Add(h_hzqq_BDT_polp80_jet1_mass);

  TCanvas *canvas_h_BDT_Signal_polp80_thstack = setUpperCanvas("canvas_h_BDT_Signal_polp80_thstack");
  canvas_h_BDT_Signal_polp80_thstack->cd();
  //h_hzqq_BDT_polp80_jet1_mass->
  hzqq_sig_BG_Stack_polp80->Draw("hist");
  hzqq_sig_BG_Stack_polp80->GetXaxis()->SetTitle("jet1 mass [GeV]");
  hzqq_sig_BG_Stack_polp80->GetYaxis()->SetTitle("Events");
  canvas_h_BDT_Signal_polp80_thstack->Modified();

  TLegend* leg_h_BDT_Overview_polp80_test=new TLegend(0.25,0.65,0.65,0.90);
  //TLegend* leg_h_BDT_Overview_polp80_test=new TLegend(0.25,0.75,0.65,0.90);
  leg_h_BDT_Overview_polp80_test->SetBorderSize(0);
  leg_h_BDT_Overview_polp80_test->SetTextAlign(12);
  leg_h_BDT_Overview_polp80_test->SetTextSize(0.050);
  leg_h_BDT_Overview_polp80_test->SetTextFont(42);
  leg_h_BDT_Overview_polp80_test->SetMargin(0.15);
  leg_h_BDT_Overview_polp80_test->SetLineColor(1);
  leg_h_BDT_Overview_polp80_test->SetLineStyle(1);
  leg_h_BDT_Overview_polp80_test->SetLineWidth(1);
  leg_h_BDT_Overview_polp80_test->SetFillColor(0);
  //leg_h_BDT_Overview_polp80_test->SetFillStyle(1001);
  leg_h_BDT_Overview_polp80_test->SetFillStyle(0);
  leg_h_BDT_Overview_polp80_test->SetHeader("#sqrt{s}>2500 GeV");
  leg_h_BDT_Overview_polp80_test->AddEntry(h_hzqq_BDT_polp80_jet1_mass,"HZ");
  leg_h_BDT_Overview_polp80_test->AddEntry(h_ee_qq_BDT_polp80_jet1_mass,"ee#rightarrow qq");
  leg_h_BDT_Overview_polp80_test->AddEntry(h_ee_qqqq_BDT_polp80_jet1_mass,"ee#rightarrow qqqq");
  leg_h_BDT_Overview_polp80_test->AddEntry(h_ee_qqqqqq_BDT_polp80_jet1_mass,"ee#rightarrow qqqqqq");
  leg_h_BDT_Overview_polp80_test->Draw();

  l->DrawLatex(x,y,label.c_str());
  l->DrawLatex(x2,y2,label3.c_str());



  std::cout<<"do i get anywhere here"<<std::endl;


  TFile* file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polm80=TFile::Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug5/MVATrainingWeights_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees300NCuts20__hzqq__ee_qq_mqq_1TeV__ee_qqqq_mqqqq_2TeV__ee_qqqqqq_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_allVar_no_Jet2_C2_D2.root");
  //eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/MVTrainingWeights_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NCuts20__hzqq__ee_qq_mqq_1TeV__ee_qqqq_mqqqq_2TeV__ee_qqqqqq_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_allVar.root");



  TFile* file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polp80=TFile::Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug5/MVATrainingWeights_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees250NCuts20__hzqq__ee_qq_mqq_1TeV__ee_qqqq_mqqqq_2TeV__ee_qqqqqq_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_allVar_no_Jet2_C2_D2.root");
  //eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/MVTrainingWeights_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NCuts20___hzqq__ee_qq_mqq_1TeV__ee_qqqq_mqqqq_2TeV__ee_qqqqqq_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_allVar.root");

  TH1F* h_BDT_Signal_polm80_test=(TH1F*)file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polm80->Get("dataset/Method_BDT/BDT/MVA_BDT_S");
  h_BDT_Signal_polm80_test->SetFillColor(kBlue-2);
  h_BDT_Signal_polm80_test->SetFillStyle(3001);
  h_BDT_Signal_polm80_test->GetXaxis()->SetTitle("BDT Score");
  h_BDT_Signal_polm80_test->GetYaxis()->SetTitle("Entries");

  TH1F* h_BDT_Signal_polm80_train=(TH1F*)file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polm80->Get("dataset/Method_BDT/BDT/MVA_BDT_Train_S");
  h_BDT_Signal_polm80_train->SetLineColor(kBlue);
  h_BDT_Signal_polm80_train->SetLineWidth(2);
  h_BDT_Signal_polm80_train->GetXaxis()->SetTitle("BDT Score");
  h_BDT_Signal_polm80_train->GetYaxis()->SetTitle("Entries");

  TH1F* h_BDT_BG_polm80_test=(TH1F*)file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polm80->Get("dataset/Method_BDT/BDT/MVA_BDT_B");
  h_BDT_BG_polm80_test->SetFillColor(kRed-2);
  h_BDT_BG_polm80_test->SetFillStyle(3002);
  h_BDT_BG_polm80_test->GetXaxis()->SetTitle("BDT Score");
  h_BDT_BG_polm80_test->GetYaxis()->SetTitle("Entries");

  TH1F* h_BDT_BG_polm80_train=(TH1F*)file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polm80->Get("dataset/Method_BDT/BDT/MVA_BDT_Train_B");
  h_BDT_BG_polm80_train->SetLineColor(kRed);
  h_BDT_BG_polm80_train->SetLineWidth(2);
  h_BDT_BG_polm80_train->GetXaxis()->SetTitle("BDT Score");
  h_BDT_BG_polm80_train->GetYaxis()->SetTitle("Entries");
 
  TCanvas *canvas_h_BDT_Signal_polm80_test = setUpperCanvas("canvas_h_BDT_Signal_polm80_test");
  canvas_h_BDT_Signal_polm80_test->cd();

  TLegend* leg_h_BDT_Signal_polm80_test=new TLegend(0.25,0.65,0.65,0.90);
  //TLegend* leg_h_BDT_Signal_polm80_test=new TLegend(0.25,0.75,0.65,0.90);
  leg_h_BDT_Signal_polm80_test->SetBorderSize(0);
  leg_h_BDT_Signal_polm80_test->SetTextAlign(12);
  leg_h_BDT_Signal_polm80_test->SetTextSize(0.050);
  leg_h_BDT_Signal_polm80_test->SetTextFont(42);
  leg_h_BDT_Signal_polm80_test->SetMargin(0.15);
  leg_h_BDT_Signal_polm80_test->SetLineColor(1);
  leg_h_BDT_Signal_polm80_test->SetLineStyle(1);
  leg_h_BDT_Signal_polm80_test->SetLineWidth(1);
  leg_h_BDT_Signal_polm80_test->SetFillColor(0);
  //leg_h_BDT_Signal_polm80_test->SetFillStyle(1001);
  leg_h_BDT_Signal_polm80_test->SetFillStyle(0);
  leg_h_BDT_Signal_polm80_test->SetHeader("#sqrt{s}>2500 GeV");
  leg_h_BDT_Signal_polm80_test->AddEntry(h_BDT_Signal_polm80_test->DrawCopy("hist,e"),"Signal,test");
  leg_h_BDT_Signal_polm80_test->AddEntry(h_BDT_BG_polm80_test->DrawCopy("hist,e,same"),"BG,test");
  leg_h_BDT_Signal_polm80_test->AddEntry(h_BDT_Signal_polm80_train->DrawCopy("hist,e,same"),"Signal,Train");
  leg_h_BDT_Signal_polm80_test->AddEntry(h_BDT_BG_polm80_train->DrawCopy("hist,e,same"),"BG,Train");
  leg_h_BDT_Signal_polm80_test->Draw();

  l->DrawLatex(x,y,label.c_str());
  l->DrawLatex(x2,y2,label2.c_str());


  TH1F* h_BDT_Signal_polp80_test=(TH1F*)file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polp80->Get("dataset/Method_BDT/BDT/MVA_BDT_S");
  h_BDT_Signal_polp80_test->SetFillColor(kBlue-2);
  h_BDT_Signal_polp80_test->SetFillStyle(3001);
  h_BDT_Signal_polp80_test->GetXaxis()->SetTitle("BDT Score");
  h_BDT_Signal_polp80_test->GetYaxis()->SetTitle("Entries");

  TH1F* h_BDT_Signal_polp80_train=(TH1F*)file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polp80->Get("dataset/Method_BDT/BDT/MVA_BDT_Train_S");
  h_BDT_Signal_polp80_train->SetLineColor(kBlue);
  h_BDT_Signal_polp80_train->SetLineWidth(2);
  h_BDT_Signal_polp80_train->GetXaxis()->SetTitle("BDT Score");
  h_BDT_Signal_polp80_train->GetYaxis()->SetTitle("Entries");

  TH1F* h_BDT_BG_polp80_test=(TH1F*)file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polp80->Get("dataset/Method_BDT/BDT/MVA_BDT_B");
  h_BDT_BG_polp80_test->SetFillColor(kRed-2);
  h_BDT_BG_polp80_test->SetFillStyle(3002);
  h_BDT_BG_polp80_test->GetXaxis()->SetTitle("BDT Score");
  h_BDT_BG_polp80_test->GetYaxis()->SetTitle("Entries");

  TH1F* h_BDT_BG_polp80_train=(TH1F*)file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polp80->Get("dataset/Method_BDT/BDT/MVA_BDT_Train_B");
  h_BDT_BG_polp80_train->SetLineColor(kRed);
  h_BDT_BG_polp80_train->SetLineWidth(2);
  h_BDT_BG_polp80_train->GetXaxis()->SetTitle("BDT Score");
  h_BDT_BG_polp80_train->GetYaxis()->SetTitle("Entries");
 
  TCanvas *canvas_h_BDT_Signal_polp80_test = setUpperCanvas("canvas_h_BDT_Signal_polp80_test");
  canvas_h_BDT_Signal_polp80_test->cd();

  TLegend* leg_h_BDT_Signal_polp80_test=new TLegend(0.25,0.65,0.65,0.90);
  //TLegend* leg_h_BDT_Signal_polp80_test=new TLegend(0.25,0.75,0.65,0.90);
  leg_h_BDT_Signal_polp80_test->SetBorderSize(0);
  leg_h_BDT_Signal_polp80_test->SetTextAlign(12);
  leg_h_BDT_Signal_polp80_test->SetTextSize(0.050);
  leg_h_BDT_Signal_polp80_test->SetTextFont(42);
  leg_h_BDT_Signal_polp80_test->SetMargin(0.15);
  leg_h_BDT_Signal_polp80_test->SetLineColor(1);
  leg_h_BDT_Signal_polp80_test->SetLineStyle(1);
  leg_h_BDT_Signal_polp80_test->SetLineWidth(1);
  leg_h_BDT_Signal_polp80_test->SetFillColor(0);
  //leg_h_BDT_Signal_polp80_test->SetFillStyle(1001);
  leg_h_BDT_Signal_polp80_test->SetFillStyle(0);
  leg_h_BDT_Signal_polp80_test->SetHeader("#sqrt{s}>2500 GeV");
  leg_h_BDT_Signal_polp80_test->AddEntry(h_BDT_Signal_polp80_test->DrawCopy("hist,e"),"Signal,test");
  leg_h_BDT_Signal_polp80_test->AddEntry(h_BDT_BG_polp80_test->DrawCopy("hist,e,same"),"BG,test");
  leg_h_BDT_Signal_polp80_test->AddEntry(h_BDT_Signal_polp80_train->DrawCopy("hist,e,same"),"Signal,Train");
  leg_h_BDT_Signal_polp80_test->AddEntry(h_BDT_BG_polp80_train->DrawCopy("hist,e,same"),"BG,Train");
  leg_h_BDT_Signal_polp80_test->Draw();

  l->DrawLatex(x,y,label.c_str());
  l->DrawLatex(x2,y2,label3.c_str());


}




