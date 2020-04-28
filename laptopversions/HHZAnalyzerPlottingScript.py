from ROOT import gROOT, TCanvas,THStack, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor,TMVA,TCut,TString,TDirectory,TLatex,TColor,kBlack, kBlue,kRed,kCyan,kGreen,kOrange,kYellow,kWhite
import ROOT as root
from math import cos, sin, pi, degrees, radians, pow, sqrt,acos
from array import array
#import numpy as np

#h_Zuds100_JER_wO_DR07->SetLineColor(kCyan+1);
#h_Zuds200_JER_wO_DR07->SetLineColor(kYellow+1);
#h_Zuds380_JER_wO_DR07->SetLineColor(kOrange);
#h_Zuds500_JER_wO_DR07->SetLineColor(kGreen-2)
#h_Zuds1500_JER_wO_DR07->SetLineColor(kBlue);
#h_Zuds3000_JER_wO_DR07->SetLineColor(kRed);


def DeltaPhi( Phi1, Phi2 ):
    deltaphi=abs(Phi1-Phi2)
    if (deltaphi>M_PI):
        deltaphi=2*M_PI-deltaphi        
    return deltaphi

def DeltaPhiDir( Phi1, Phi2 ):
    deltaphi=Phi1-Phi2
    if(deltaphi>M_PI):
        deltaphi=deltaphi-2*M_PI   
    if(deltaphi<(-M_PI)):
       deltaphi=2*M_PI+deltaphi
    return deltaphi

def setUpperCanvas(canvas_name) :
  c1= TCanvas(canvas_name,canvas_name,0,0,800,700);
  c1.cd();
  gStyle.SetOptStat(0);
  #c1.Range(-186.894,-0.873515,1682.046,6.114605);
  c1.SetFillColor(0);
  c1.SetBorderMode(0);
  c1.SetBorderSize(2);
  #c1.SetGridx();
  #c1.SetGridy();
  #c1.SetRightMargin(0.0172);
  c1.SetTopMargin(0.085);
  #c1.SetBottomMargin(0.138);
  c1.SetFrameBorderMode(0);
  c1.SetFrameBorderMode(0);
  return c1

def setUpperCanvasWColz(canvas_name) :
  c1= TCanvas(canvas_name,canvas_name,0,0,800,700);
  c1.cd();
  gStyle.SetOptStat(0);
  #c1.Range(-186.894,-0.873515,1682.046,6.114605);
  c1.SetFillColor(0);
  c1.SetBorderMode(0);
  c1.SetBorderSize(2);
  #c1.SetGridx();
  #c1.SetGridy();
  c1.SetRightMargin(0.160);
  c1.SetTopMargin(0.085);
  #c1.SetBottomMargin(0.138);
  c1.SetFrameBorderMode(0);
  c1.SetFrameBorderMode(0);
  return c1

    #c1= TCanvas(canvas_name,canvas_name,10,50,600,500)
    #c1.cd()
    #gPad.SetTopMargin(0.06)
    #return c1

def CLICdpStyle(): 
    gROOT.SetStyle("Plain") 
    gStyle.SetCanvasColor(root.kWhite)
    gStyle.SetFrameFillColor(root.kWhite)
    gStyle.SetStatColor(root.kWhite)
    gStyle.SetPadColor(root.kWhite)
    gStyle.SetFillColor(10)
    gStyle.SetTitleFillColor(root.kWhite)
  
    gStyle.SetPaperSize(20, 26) 
    
    gStyle.SetDrawBorder(0)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetFrameBorderMode(0)
    gStyle.SetLegendBorderSize(0)
   
    gStyle.SetTextSize(0.05)
    gStyle.SetTitleSize(0.06,"xyz")
    gStyle.SetLabelSize(0.06,"xyz")
    gStyle.SetLabelOffset(0.015,"xyz")
    gStyle.SetTitleOffset(1.2,"yz") 
    gStyle.SetTitleOffset(1.17,"x")
 
    font = 42 
    gStyle.SetTitleFont(font)
    gStyle.SetTitleFontSize(0.06)
    gStyle.SetStatFont(font)
    gStyle.SetStatFontSize(0.07)
    gStyle.SetTextFont(font)
    gStyle.SetLabelFont(font,"xyz")
    gStyle.SetTitleFont(font,"xyz")
    gStyle.SetTitleBorderSize(0)
    gStyle.SetStatBorderSize(1)
    gStyle.SetMarkerStyle(1)
    gStyle.SetLineWidth(2)  
    gStyle.SetMarkerSize(1.2)
    gStyle.SetPalette(1) 

    gStyle.SetOptTitle(0)
    gStyle.SetOptStat(0) 
    gStyle.SetOptFit(0) 
    gStyle.SetEndErrorSize(5)   

    gStyle.SetHistLineWidth(2)
    gStyle.SetFrameLineWidth(2)
    gStyle.SetFuncWidth(2)
    gStyle.SetHistLineColor(root.kBlack)
    gStyle.SetFuncColor(root.kBlack)
    gStyle.SetLabelColor(root.kBlack,"xyz")


    gStyle.SetPadBottomMargin(0.18)
    gStyle.SetPadTopMargin(0.11)
    gStyle.SetPadRightMargin(0.08)
    gStyle.SetPadLeftMargin(0.17)
    
    gStyle.SetNdivisions(506, "xy")
   
    gStyle.SetPadGridX(0)
    gStyle.SetPadGridY(0)
   
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)

    gStyle.SetCanvasDefW(800)
    gStyle.SetCanvasDefH(700)
    
    gROOT.ForceStyle()


def process_files():

    CLICdpStyle()

    l = TLatex();
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextColor(kBlack);
    l.SetTextSize(0.045);

    label="CLICdp"
    x=0.17
    y=0.93
    label2="L=4 ab^{-1},e^{-} pol -80%"
    label3="L=1 ab^{-1},e^{-} pol +80%"
    x2=0.62
    y2=0.93
    label2_noLumi="e^{-} pol -80%"
    label3_noLumi="e^{-} pol +80%"
    x2_noLumi=0.75
    x2_noLumi_normed=0.67
    y2_noLumi=0.93

    
    file_polm80_HHZ_SignalHistos_VLC10_njet3_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC10_NJets3_yij/polm80/test_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14_partonlevelOnly_noIsoPh_TrueMCJets_highCombPostProcess.root")

    h_sqrtS_parton_all_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_sqrtS_e1_e2_effective")
    h_sqrtS_parton_all_polm80.Rebin(4)
    h_sqrtS_parton_all_polm80.SetLineWidth(2)
    h_sqrtS_parton_all_polm80.SetLineColor(1)
    h_sqrtS_parton_all_polm80.GetXaxis().SetTitle('#sqrt{s}_{part}')
    h_sqrtS_parton_all_polm80.GetYaxis().SetTitle('Events')
    
    canvas_h_sqrtS_parton_all_polm80 = setUpperCanvas("canvas_hhz_h_sqrtS_parton_all_polm80");
    canvas_h_sqrtS_parton_all_polm80.cd()
    h_sqrtS_parton_all_polm80.Draw()

    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test=TLegend(0.20,0.61,0.60,0.87);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetBorderSize(0);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetTextAlign(12);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetTextSize(0.050);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetTextFont(42);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetMargin(0.15);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetLineColor(1);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetLineStyle(1);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetLineWidth(1);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetFillColor(0);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetFillStyle(0);
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.AddEntry(h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.DrawCopy("h"),"method1");
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.AddEntry(h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.DrawCopy("h,same"),"method2");
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.AddEntry(h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.DrawCopy("h,same"),"method3");
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.AddEntry(h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.DrawCopy("h,same"),"method4");
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.AddEntry(h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.DrawCopy("h,same"),"method5");
    #leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_sqrtS_parton_all_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_hhz_h_sqrtS_parton_all_polm80.eps")

    canvas_h_sqrtS_parton_vs_gj_all_polm80 = setUpperCanvas("canvas_hhz_h_sqrtS_parton_vs_genj_all_polm80");
    canvas_h_sqrtS_parton_vs_gj_all_polm80.cd()
    h_sqrtS_gj_sum_all_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_sqrtS_gj_sum_all")
    h_sqrtS_gj_sum_all_polm80.Rebin(4)
    h_sqrtS_gj_sum_all_polm80.SetLineWidth(2)
    h_sqrtS_gj_sum_all_polm80.SetLineColor(2)
    h_sqrtS_gj_sum_all_polm80.GetXaxis().SetTitle('#sqrt{s}')
    h_sqrtS_gj_sum_all_polm80.GetYaxis().SetTitle('Events')
    h_sqrtS_parton_all_polm80.GetXaxis().SetTitle('#sqrt{s}')

    h_sqrtS_gj_sum_plus_InvVec_all_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_sqrtS_gj_sum_plus_InvVec_all")
    h_sqrtS_gj_sum_plus_InvVec_all_polm80.Rebin(4)
    h_sqrtS_gj_sum_plus_InvVec_all_polm80.SetLineWidth(2)
    h_sqrtS_gj_sum_plus_InvVec_all_polm80.SetLineColor(4)
    h_sqrtS_gj_sum_plus_InvVec_all_polm80.GetXaxis().SetTitle('#sqrt{s}_{part}')
    h_sqrtS_gj_sum_plus_InvVec_all_polm80.GetYaxis().SetTitle('Events')
    
    
    canvas_h_sqrtS_parton_vs_gj_all_polm80 = setUpperCanvas("canvas_hhz_h_sqrtS_parton_vs_gj_all_polm80");
    canvas_h_sqrtS_parton_vs_gj_all_polm80.cd()

    leg_h_sqrtS_parton_vs_gj_all_polm80_test=TLegend(0.20,0.61,0.60,0.87);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetBorderSize(0);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetTextAlign(12);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetTextSize(0.050);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetTextFont(42);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetMargin(0.15);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetLineColor(1);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetLineStyle(1);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetLineWidth(1);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetFillColor(0);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetFillStyle(0);
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.SetHeader("all Events");
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.AddEntry(h_sqrtS_parton_all_polm80.DrawCopy("h"),"partons");
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.AddEntry(h_sqrtS_gj_sum_all_polm80.DrawCopy("h,same"),"genJets");
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.AddEntry(h_sqrtS_gj_sum_plus_InvVec_all_polm80.DrawCopy("h,same"),"genJets+Neutrinos");
    leg_h_sqrtS_parton_vs_gj_all_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);


    canvas_h_sqrtS_parton_vs_gj_all_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_hhz_h_sqrtS_parton_vs_gj_all_polm80.eps")

    h_sqrtS_parton_HHZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_sqrtS_e1_e2_effective_HHZ_bbbbqq")
    h_sqrtS_parton_HHZ_bbbbqq_polm80.Rebin(4)
    h_sqrtS_parton_HHZ_bbbbqq_polm80.SetLineWidth(2)
    h_sqrtS_parton_HHZ_bbbbqq_polm80.SetLineColor(1)
    h_sqrtS_parton_HHZ_bbbbqq_polm80.GetXaxis().SetTitle('#sqrt{s}_{part}')
    h_sqrtS_parton_HHZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_sqrtS_gj_sum_HHZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_sqrtS_gj_sum_HHZ_bbbbqq")
    h_sqrtS_gj_sum_HHZ_bbbbqq_polm80.Rebin(4)
    h_sqrtS_gj_sum_HHZ_bbbbqq_polm80.SetLineWidth(2)
    h_sqrtS_gj_sum_HHZ_bbbbqq_polm80.SetLineColor(2)
    h_sqrtS_gj_sum_HHZ_bbbbqq_polm80.GetXaxis().SetTitle('#sqrt{s}')
    h_sqrtS_gj_sum_HHZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_sqrtS_parton_HHZ_bbbbqq_polm80.GetXaxis().SetTitle('#sqrt{s}')

    h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq")
    h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq_polm80.Rebin(4)
    h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq_polm80.SetLineWidth(2)
    h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq_polm80.SetLineColor(4)
    h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq_polm80.GetXaxis().SetTitle('#sqrt{s}_{part}')
    h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')


    canvas_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80 = setUpperCanvas("canvas_hhz_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80");
    canvas_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80.cd()

    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test=TLegend(0.20,0.61,0.60,0.87);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.SetHeader("HHZ #rightarrow bbbbqq Events");
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.AddEntry(h_sqrtS_parton_HHZ_bbbbqq_polm80.DrawCopy("h"),"partons");
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.AddEntry(h_sqrtS_gj_sum_HHZ_bbbbqq_polm80.DrawCopy("h,same"),"genJets");
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.AddEntry(h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq_polm80.DrawCopy("h,same"),"genJets+Neutrinos");
    leg_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_hhz_h_sqrtS_parton_vs_gj_HHZ_bbbbqq_polm80.eps")
  



    h_dalpha_H1_H2_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_dalpha_H1_H2")
    h_dalpha_H1_H2_polm80.Rebin(4)
    h_dalpha_H1_H2_polm80.SetLineWidth(2)
    h_dalpha_H1_H2_polm80.SetLineColor(1)
    h_dalpha_H1_H2_polm80.GetXaxis().SetTitle('#Delta#alpha(B1,B2) [#circ]')
    h_dalpha_H1_H2_polm80.GetYaxis().SetTitle('Events')

    h_dalpha_H_Z_min_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_dalpha_H_Z_min")
    h_dalpha_H_Z_min_polm80.Rebin(4)
    h_dalpha_H_Z_min_polm80.SetLineWidth(2)
    h_dalpha_H_Z_min_polm80.SetLineColor(2)
    h_dalpha_H_Z_min_polm80.GetXaxis().SetTitle('#Delta#alpha(B1,B2)')
    h_dalpha_H_Z_min_polm80.GetYaxis().SetTitle('Events')

    h_dalpha_H_Z_max_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_dalpha_H_Z_max")
    h_dalpha_H_Z_max_polm80.Rebin(4)
    h_dalpha_H_Z_max_polm80.SetLineWidth(2)
    h_dalpha_H_Z_max_polm80.SetLineColor(4)
    h_dalpha_H_Z_max_polm80.GetXaxis().SetTitle('#Delta#alpha(B1,B2)')
    h_dalpha_H_Z_max_polm80.GetYaxis().SetTitle('Events')


    canvas_h_dalpha_B1_B2_polm80 = setUpperCanvas("canvas_hhz_h_dalpha_B1_B2_polm80");
    canvas_h_dalpha_B1_B2_polm80.cd()

    leg_h_dalpha_B1_B2_polm80_test=TLegend(0.20,0.61,0.60,0.87);
    leg_h_dalpha_B1_B2_polm80_test.SetBorderSize(0);
    leg_h_dalpha_B1_B2_polm80_test.SetTextAlign(12);
    leg_h_dalpha_B1_B2_polm80_test.SetTextSize(0.050);
    leg_h_dalpha_B1_B2_polm80_test.SetTextFont(42);
    leg_h_dalpha_B1_B2_polm80_test.SetMargin(0.15);
    leg_h_dalpha_B1_B2_polm80_test.SetLineColor(1);
    leg_h_dalpha_B1_B2_polm80_test.SetLineStyle(1);
    leg_h_dalpha_B1_B2_polm80_test.SetLineWidth(1);
    leg_h_dalpha_B1_B2_polm80_test.SetFillColor(0);
    leg_h_dalpha_B1_B2_polm80_test.SetFillStyle(0);
    leg_h_dalpha_B1_B2_polm80_test.SetHeader("HHZ Events");
    leg_h_dalpha_B1_B2_polm80_test.AddEntry(h_dalpha_H1_H2_polm80.DrawCopy("h"),"#Delta#alpha(H1,H2)");
    leg_h_dalpha_B1_B2_polm80_test.AddEntry(h_dalpha_H_Z_min_polm80.DrawCopy("h,same"),"min #Delta#alpha(H,Z)");
    leg_h_dalpha_B1_B2_polm80_test.AddEntry(h_dalpha_H_Z_max_polm80.DrawCopy("h,same"),"max #Delta#alpha(H,Z)");
    leg_h_dalpha_B1_B2_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_dalpha_B1_B2_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_hhz_h_dalpha_B1_B2_polm80.eps")

    h_theta_H1_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_theta_H1")
    h_theta_H1_polm80.Rebin(4)
    h_theta_H1_polm80.SetLineWidth(2)
    h_theta_H1_polm80.SetLineColor(1)
    h_theta_H1_polm80.GetXaxis().SetTitle('#theta [#circ]')
    h_theta_H1_polm80.GetYaxis().SetTitle('Events')

    h_theta_H2_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_theta_H2")
    h_theta_H2_polm80.Rebin(4)
    h_theta_H2_polm80.SetLineWidth(2)
    h_theta_H2_polm80.SetLineColor(2)
    h_theta_H2_polm80.GetXaxis().SetTitle('#theta [#circ]')
    h_theta_H2_polm80.GetYaxis().SetTitle('Events')

    h_theta_Z_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_theta_Z")
    h_theta_Z_polm80.Rebin(4)
    h_theta_Z_polm80.SetLineWidth(2)
    h_theta_Z_polm80.SetLineColor(4)
    h_theta_Z_polm80.GetXaxis().SetTitle('#theta [#circ]')
    h_theta_Z_polm80.GetYaxis().SetTitle('Events')


    canvas_theta_B_polm80 = setUpperCanvas("canvas_hhz_theta_B_polm80");
    canvas_theta_B_polm80.cd()

    leg_theta_B_polm80_test=TLegend(0.20,0.61,0.60,0.87);
    leg_theta_B_polm80_test.SetBorderSize(0);
    leg_theta_B_polm80_test.SetTextAlign(12);
    leg_theta_B_polm80_test.SetTextSize(0.050);
    leg_theta_B_polm80_test.SetTextFont(42);
    leg_theta_B_polm80_test.SetMargin(0.15);
    leg_theta_B_polm80_test.SetLineColor(1);
    leg_theta_B_polm80_test.SetLineStyle(1);
    leg_theta_B_polm80_test.SetLineWidth(1);
    leg_theta_B_polm80_test.SetFillColor(0);
    leg_theta_B_polm80_test.SetFillStyle(0);
    leg_theta_B_polm80_test.SetHeader("HHZ Events");
    leg_theta_B_polm80_test.AddEntry(h_theta_H2_polm80.DrawCopy("h"),"H2");
    leg_theta_B_polm80_test.AddEntry(h_theta_H1_polm80.DrawCopy("h,same"),"H1");
    leg_theta_B_polm80_test.AddEntry(h_theta_Z_polm80.DrawCopy("h,same"),"Z");
    leg_theta_B_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_theta_B_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_hhz_theta_B_polm80.eps")


    h_E_H1_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_H1_E")
    h_E_H1_polm80.Rebin(4)
    h_E_H1_polm80.SetLineWidth(2)
    h_E_H1_polm80.SetLineColor(1)
    h_E_H1_polm80.GetXaxis().SetTitle('E [GeV]')
    h_E_H1_polm80.GetYaxis().SetTitle('Events')

    h_E_H2_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_H2_E")
    h_E_H2_polm80.Rebin(4)
    h_E_H2_polm80.SetLineWidth(2)
    h_E_H2_polm80.SetLineColor(2)
    h_E_H2_polm80.GetXaxis().SetTitle('E [GeV]')
    h_E_H2_polm80.GetYaxis().SetTitle('Events')

    h_E_Z_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_Z_E")
    h_E_Z_polm80.Rebin(4)
    h_E_Z_polm80.SetLineWidth(2)
    h_E_Z_polm80.SetLineColor(4)
    h_E_Z_polm80.GetXaxis().SetTitle('E [GeV]')
    h_E_Z_polm80.GetYaxis().SetTitle('Events')


    canvas_E_B_polm80 = setUpperCanvas("canvas_hhz_E_B_polm80");
    canvas_E_B_polm80.cd()

    leg_E_B_polm80_test=TLegend(0.40,0.61,0.80,0.87);
    leg_E_B_polm80_test.SetBorderSize(0);
    leg_E_B_polm80_test.SetTextAlign(12);
    leg_E_B_polm80_test.SetTextSize(0.050);
    leg_E_B_polm80_test.SetTextFont(42);
    leg_E_B_polm80_test.SetMargin(0.15);
    leg_E_B_polm80_test.SetLineColor(1);
    leg_E_B_polm80_test.SetLineStyle(1);
    leg_E_B_polm80_test.SetLineWidth(1);
    leg_E_B_polm80_test.SetFillColor(0);
    leg_E_B_polm80_test.SetFillStyle(0);
    leg_E_B_polm80_test.SetHeader("HHZ Events");
    leg_E_B_polm80_test.AddEntry(h_E_Z_polm80.DrawCopy("h"),"Z");
    leg_E_B_polm80_test.AddEntry(h_E_H1_polm80.DrawCopy("h,same"),"H1");
    leg_E_B_polm80_test.AddEntry(h_E_H2_polm80.DrawCopy("h,same"),"H2");
    leg_E_B_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_E_B_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_hhz_E_B_polm80.eps")

    h_P_H1_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_H1_P")
    h_P_H1_polm80.Rebin(4)
    h_P_H1_polm80.SetLineWidth(2)
    h_P_H1_polm80.SetLineColor(1)
    h_P_H1_polm80.GetXaxis().SetTitle('p [GeV]')
    h_P_H1_polm80.GetYaxis().SetTitle('Events')

    h_P_H2_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_H2_P")
    h_P_H2_polm80.Rebin(4)
    h_P_H2_polm80.SetLineWidth(2)
    h_P_H2_polm80.SetLineColor(2)
    h_P_H2_polm80.GetXaxis().SetTitle('p [GeV]')
    h_P_H2_polm80.GetYaxis().SetTitle('Events')

    h_P_Z_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_Z_P")
    h_P_Z_polm80.Rebin(4)
    h_P_Z_polm80.SetLineWidth(2)
    h_P_Z_polm80.SetLineColor(4)
    h_P_Z_polm80.GetXaxis().SetTitle('p [GeV]')
    h_P_Z_polm80.GetYaxis().SetTitle('Events')


    canvas_P_B_polm80 = setUpperCanvas("canvas_hhz_P_B_polm80");
    canvas_P_B_polm80.cd()

    leg_P_B_polm80_test=TLegend(0.40,0.61,0.80,0.87);
    leg_P_B_polm80_test.SetBorderSize(0);
    leg_P_B_polm80_test.SetTextAlign(12);
    leg_P_B_polm80_test.SetTextSize(0.050);
    leg_P_B_polm80_test.SetTextFont(42);
    leg_P_B_polm80_test.SetMargin(0.15);
    leg_P_B_polm80_test.SetLineColor(1);
    leg_P_B_polm80_test.SetLineStyle(1);
    leg_P_B_polm80_test.SetLineWidth(1);
    leg_P_B_polm80_test.SetFillColor(0);
    leg_P_B_polm80_test.SetFillStyle(0);
    leg_P_B_polm80_test.SetHeader("HHZ Events");
    leg_P_B_polm80_test.AddEntry(h_P_Z_polm80.DrawCopy("h"),"Z");
    leg_P_B_polm80_test.AddEntry(h_P_H1_polm80.DrawCopy("h,same"),"H1");
    leg_P_B_polm80_test.AddEntry(h_P_H2_polm80.DrawCopy("h,same"),"H2");
    leg_P_B_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_P_B_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_hhz_P_B_polm80.eps")


    h_d_ij_H1_H2_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_d_ij_H1_H2")
    h_d_ij_H1_H2_polm80.Rebin(4)
    h_d_ij_H1_H2_polm80.SetLineWidth(2)
    h_d_ij_H1_H2_polm80.SetLineColor(1)
    h_d_ij_H1_H2_polm80.GetXaxis().SetTitle('d_{B1,B2}=1-cos#theta(B1,B2)')
    h_d_ij_H1_H2_polm80.GetYaxis().SetTitle('Events')

    h_d_ij_min_V1_V2_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_d_ij_min_V1_V2")
    h_d_ij_min_V1_V2_polm80.Rebin(4)
    h_d_ij_min_V1_V2_polm80.SetLineWidth(2)
    h_d_ij_min_V1_V2_polm80.SetLineColor(2)
    h_d_ij_min_V1_V2_polm80.GetXaxis().SetTitle('d_{B1,B2}=1-cos#theta(B1,B2)')
    h_d_ij_min_V1_V2_polm80.GetYaxis().SetTitle('Events')

    h_d_ij_max_V1_V2_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_d_ij_max_V1_V2")
    h_d_ij_max_V1_V2_polm80.Rebin(4)
    h_d_ij_max_V1_V2_polm80.SetLineWidth(2)
    h_d_ij_max_V1_V2_polm80.SetLineColor(4)
    h_d_ij_max_V1_V2_polm80.GetXaxis().SetTitle('d_{B1,B2}=1-cos#theta(B1,B2)')
    h_d_ij_max_V1_V2_polm80.GetYaxis().SetTitle('Events')

    canvas_d_ij_B_polm80 = setUpperCanvas("canvas_hhz_d_ij_B_polm80");
    canvas_d_ij_B_polm80.cd()

    leg_d_ij_B_polm80_test=TLegend(0.70,0.61,0.90,0.87);
    leg_d_ij_B_polm80_test.SetBorderSize(0);
    leg_d_ij_B_polm80_test.SetTextAlign(12);
    leg_d_ij_B_polm80_test.SetTextSize(0.050);
    leg_d_ij_B_polm80_test.SetTextFont(42);
    leg_d_ij_B_polm80_test.SetMargin(0.15);
    leg_d_ij_B_polm80_test.SetLineColor(1);
    leg_d_ij_B_polm80_test.SetLineStyle(1);
    leg_d_ij_B_polm80_test.SetLineWidth(1);
    leg_d_ij_B_polm80_test.SetFillColor(0);
    leg_d_ij_B_polm80_test.SetFillStyle(0);
    leg_d_ij_B_polm80_test.SetHeader("HHZ Events");
    leg_d_ij_B_polm80_test.AddEntry(h_d_ij_max_V1_V2_polm80.DrawCopy("h"),"d_{ij} max");
    leg_d_ij_B_polm80_test.AddEntry(h_d_ij_min_V1_V2_polm80.DrawCopy("h,same"),"d_{ij} min");
    leg_d_ij_B_polm80_test.AddEntry(h_d_ij_H1_H2_polm80.DrawCopy("h,same"),"d_{ij} (H1,H2)");
    leg_d_ij_B_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_d_ij_B_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_hhz_d_ij_B_polm80.eps")


    h_d_ij_H1_qqbar_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_d_ij_H1_qqbar")
    h_d_ij_H1_qqbar_polm80.Rebin(2)
    h_d_ij_H1_qqbar_polm80.SetLineWidth(2)
    h_d_ij_H1_qqbar_polm80.SetLineColor(1)
    h_d_ij_H1_qqbar_polm80.GetXaxis().SetTitle('d_{q1,q2}=1-cos#theta(q1,q2)')
    h_d_ij_H1_qqbar_polm80.GetYaxis().SetTitle('Events')
    h_d_ij_H1_qqbar_polm80.SetMinimum(0)
    h_d_ij_H1_qqbar_polm80.GetXaxis().SetRangeUser(0,2.0)

    h_d_ij_H2_qqbar_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_d_ij_H2_qqbar")
    h_d_ij_H2_qqbar_polm80.Rebin(2)
    h_d_ij_H2_qqbar_polm80.SetLineWidth(2)
    h_d_ij_H2_qqbar_polm80.SetLineColor(2)
    h_d_ij_H2_qqbar_polm80.GetXaxis().SetTitle('d_{B1,B2}=1-cos#theta(B1,B2)')
    h_d_ij_H2_qqbar_polm80.GetYaxis().SetTitle('Events')

    h_d_ij_Z_qqbar_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_d_ij_Z_qqbar")
    h_d_ij_Z_qqbar_polm80.Rebin(2)
    h_d_ij_Z_qqbar_polm80.SetLineWidth(2)
    h_d_ij_Z_qqbar_polm80.SetLineColor(4)
    h_d_ij_Z_qqbar_polm80.GetXaxis().SetTitle('d_{B1,B2}=1-cos#theta(B1,B2)')
    h_d_ij_Z_qqbar_polm80.GetYaxis().SetTitle('Events')

    canvas_d_ij_qqbar_fromB_polm80 = setUpperCanvas("canvas_hhz_d_ij_qqbar_fromB_polm80");
    canvas_d_ij_qqbar_fromB_polm80.cd()

    leg_d_ij_qqbar_fromB_polm80_test=TLegend(0.70,0.61,0.90,0.87);
    leg_d_ij_qqbar_fromB_polm80_test.SetBorderSize(0);
    leg_d_ij_qqbar_fromB_polm80_test.SetTextAlign(12);
    leg_d_ij_qqbar_fromB_polm80_test.SetTextSize(0.050);
    leg_d_ij_qqbar_fromB_polm80_test.SetTextFont(42);
    leg_d_ij_qqbar_fromB_polm80_test.SetMargin(0.15);
    leg_d_ij_qqbar_fromB_polm80_test.SetLineColor(1);
    leg_d_ij_qqbar_fromB_polm80_test.SetLineStyle(1);
    leg_d_ij_qqbar_fromB_polm80_test.SetLineWidth(1);
    leg_d_ij_qqbar_fromB_polm80_test.SetFillColor(0);
    leg_d_ij_qqbar_fromB_polm80_test.SetFillStyle(0);
    leg_d_ij_qqbar_fromB_polm80_test.SetHeader("HHZ Events");
    leg_d_ij_qqbar_fromB_polm80_test.AddEntry(h_d_ij_H1_qqbar_polm80.DrawCopy("h"),"H1 d_{ij}");
    leg_d_ij_qqbar_fromB_polm80_test.AddEntry(h_d_ij_H2_qqbar_polm80.DrawCopy("h,same"),"H2 d_{ij}");
    leg_d_ij_qqbar_fromB_polm80_test.AddEntry(h_d_ij_Z_qqbar_polm80.DrawCopy("h,same"),"Z d_{ij}");
    leg_d_ij_qqbar_fromB_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_d_ij_qqbar_fromB_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_hhz_d_ij_qqbar_fromB_polm80.eps")



    h_njet3_VLC10_mass_gj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_mass_comb_gj1_HHZ_bbbbqq")
    h_njet3_VLC10_mass_gj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet3_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet3_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet3_VLC10_mass_gj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet mass [GeV]')
    h_njet3_VLC10_mass_gj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet3_VLC10_mass_gj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_mass_comb_gj2_HHZ_bbbbqq")
    h_njet3_VLC10_mass_gj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet3_VLC10_mass_gj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet3_VLC10_mass_gj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet3_VLC10_mass_gj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet mass [GeV]')
    h_njet3_VLC10_mass_gj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet3_VLC10_mass_gj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_mass_comb_gj3_HHZ_bbbbqq")
    h_njet3_VLC10_mass_gj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet3_VLC10_mass_gj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet3_VLC10_mass_gj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet3_VLC10_mass_gj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet Mass [GeV]')
    h_njet3_VLC10_mass_gj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    print h_njet3_VLC10_mass_gj3_HZZ_bbbbqq_polm80.GetBinContent(10)

    canvas_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet3_VLC10_hhz_gj_mass_HZZ_bbbbqq_polm80");
    canvas_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=3");
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet3_VLC10_mass_gj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"gj1");
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet3_VLC10_mass_gj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"gj2");
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet3_VLC10_mass_gj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"gj3");
    leg_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet3_VLC10_gj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet3_VLC10_hhz_gj_mass_B_HZZ_bbbbqq_polm80.eps")


    h_njet3_VLC10_mass_rj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_mass_comb_rfj_rj1_HHZ_bbbbqq")
    h_njet3_VLC10_mass_rj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet3_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet3_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet3_VLC10_mass_rj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet3_VLC10_mass_rj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet3_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetMaximum(15)

    h_njet3_VLC10_mass_rj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_mass_comb_rfj_rj2_HHZ_bbbbqq")
    h_njet3_VLC10_mass_rj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet3_VLC10_mass_rj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet3_VLC10_mass_rj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet3_VLC10_mass_rj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet3_VLC10_mass_rj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet3_VLC10_mass_rj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet3_.Get("h_mass_comb_rfj_rj3_HHZ_bbbbqq")
    h_njet3_VLC10_mass_rj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet3_VLC10_mass_rj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet3_VLC10_mass_rj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet3_VLC10_mass_rj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet3_VLC10_mass_rj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')


    canvas_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet3_VLC10_hhz_rfj_rj_mass_HZZ_bbbbqq_polm80");
    canvas_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=3");
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet3_VLC10_mass_rj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"j1");
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet3_VLC10_mass_rj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j2");
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet3_VLC10_mass_rj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j3");
    leg_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet3_VLC10_rj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet3_VLC10_hhz_rfj_rj_mass_HZZ_bbbbqq_polm80.eps")


    file_polm80_HHZ_SignalHistos_VLC10_njet4_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC10_NJets4_yij/polm80/test_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14_partonlevelOnly_noIsoPh_TrueMCJets_highCombPostProcess.root")


    h_njet4_VLC10_mass_gj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet4_.Get("h_mass_comb_gj1_HHZ_bbbbqq")
    h_njet4_VLC10_mass_gj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet4_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet4_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet4_VLC10_mass_gj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet mass [GeV]')
    h_njet4_VLC10_mass_gj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet4_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet4_VLC10_mass_gj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet4_.Get("h_mass_comb_gj2_HHZ_bbbbqq")
    h_njet4_VLC10_mass_gj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet4_VLC10_mass_gj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet4_VLC10_mass_gj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet4_VLC10_mass_gj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet mass [GeV]')
    h_njet4_VLC10_mass_gj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet4_VLC10_mass_gj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet4_.Get("h_mass_comb_gj3_HHZ_bbbbqq")
    h_njet4_VLC10_mass_gj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet4_VLC10_mass_gj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet4_VLC10_mass_gj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet4_VLC10_mass_gj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet Mass [GeV]')
    h_njet4_VLC10_mass_gj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    canvas_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet4_VLC10_hhz_gj_mass_HZZ_bbbbqq_polm80");
    canvas_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=4");
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet4_VLC10_mass_gj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"gj1");
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet4_VLC10_mass_gj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"gj2");
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet4_VLC10_mass_gj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"gj3");
    leg_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet4_VLC10_gj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet4_VLC10_hhz_gj_mass_B_HZZ_bbbbqq_polm80.eps")


    h_njet4_VLC10_mass_rj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet4_.Get("h_mass_comb_rfj_rj1_HHZ_bbbbqq")
    h_njet4_VLC10_mass_rj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet4_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet4_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet4_VLC10_mass_rj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet4_VLC10_mass_rj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet4_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetMaximum(15)
    h_njet4_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet4_VLC10_mass_rj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet4_.Get("h_mass_comb_rfj_rj2_HHZ_bbbbqq")
    h_njet4_VLC10_mass_rj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet4_VLC10_mass_rj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet4_VLC10_mass_rj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet4_VLC10_mass_rj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet4_VLC10_mass_rj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet4_VLC10_mass_rj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet4_.Get("h_mass_comb_rfj_rj3_HHZ_bbbbqq")
    h_njet4_VLC10_mass_rj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet4_VLC10_mass_rj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet4_VLC10_mass_rj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet4_VLC10_mass_rj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet4_VLC10_mass_rj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')


    canvas_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet4_VLC10_hhz_comb_rfj_rj_mass_HZZ_bbbbqq_polm80");
    canvas_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=4");
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet4_VLC10_mass_rj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"j1");
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet4_VLC10_mass_rj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j2");
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet4_VLC10_mass_rj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j3");
    leg_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet4_VLC10_rj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet4_VLC10_hhz_rfj_rj_mass_HZZ_bbbbqq_polm80.eps")



    file_polm80_HHZ_SignalHistos_VLC10_njet5_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC10_NJets5_yij/polm80/test_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14_partonlevelOnly_noIsoPh_TrueMCJets_highCombPostProcess.root")


    h_njet5_VLC10_mass_gj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet5_.Get("h_mass_comb_gj1_HHZ_bbbbqq")
    h_njet5_VLC10_mass_gj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet5_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet5_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet5_VLC10_mass_gj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet mass [GeV]')
    h_njet5_VLC10_mass_gj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet5_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet5_VLC10_mass_gj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet5_.Get("h_mass_comb_gj2_HHZ_bbbbqq")
    h_njet5_VLC10_mass_gj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet5_VLC10_mass_gj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet5_VLC10_mass_gj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet5_VLC10_mass_gj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet mass [GeV]')
    h_njet5_VLC10_mass_gj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet5_VLC10_mass_gj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet5_.Get("h_mass_comb_gj3_HHZ_bbbbqq")
    h_njet5_VLC10_mass_gj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet5_VLC10_mass_gj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet5_VLC10_mass_gj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet5_VLC10_mass_gj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet Mass [GeV]')
    h_njet5_VLC10_mass_gj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    canvas_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet5_VLC10_hhz_gj_mass_HZZ_bbbbqq_polm80");
    canvas_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=5");
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet5_VLC10_mass_gj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"gj1");
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet5_VLC10_mass_gj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"gj2");
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet5_VLC10_mass_gj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"gj3");
    leg_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet5_VLC10_gj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet5_VLC10_hhz_gj_mass_B_HZZ_bbbbqq_polm80.eps")


    h_njet5_VLC10_mass_rj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet5_.Get("h_mass_comb_rfj_rj1_HHZ_bbbbqq")
    h_njet5_VLC10_mass_rj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet5_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet5_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet5_VLC10_mass_rj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet5_VLC10_mass_rj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet5_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetMaximum(15)
    h_njet5_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet5_VLC10_mass_rj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet5_.Get("h_mass_comb_rfj_rj2_HHZ_bbbbqq")
    h_njet5_VLC10_mass_rj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet5_VLC10_mass_rj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet5_VLC10_mass_rj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet5_VLC10_mass_rj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet5_VLC10_mass_rj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet5_VLC10_mass_rj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet5_.Get("h_mass_comb_rfj_rj3_HHZ_bbbbqq")
    h_njet5_VLC10_mass_rj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet5_VLC10_mass_rj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet5_VLC10_mass_rj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet5_VLC10_mass_rj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet5_VLC10_mass_rj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')


    canvas_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet5_VLC10_hhz_comb_rfj_rj_mass_HZZ_bbbbqq_polm80");
    canvas_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=5");
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet5_VLC10_mass_rj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"j1");
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet5_VLC10_mass_rj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j2");
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet5_VLC10_mass_rj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j3");
    leg_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet5_VLC10_rj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet5_VLC10_hhz_rfj_rj_mass_HZZ_bbbbqq_polm80.eps")


    file_polm80_HHZ_SignalHistos_VLC10_njet6_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC10_NJets6_yij/polm80/test_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14_partonlevelOnly_noIsoPh_TrueMCJets_highCombPostProcess.root")


    h_njet6_VLC10_mass_gj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_mass_comb_gj1_HHZ_bbbbqq")
    h_njet6_VLC10_mass_gj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet6_VLC10_mass_gj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet mass [GeV]')
    h_njet6_VLC10_mass_gj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet6_VLC10_mass_gj1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet6_VLC10_mass_gj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_mass_comb_gj2_HHZ_bbbbqq")
    h_njet6_VLC10_mass_gj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_mass_gj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_mass_gj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet6_VLC10_mass_gj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet mass [GeV]')
    h_njet6_VLC10_mass_gj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet6_VLC10_mass_gj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_mass_comb_gj3_HHZ_bbbbqq")
    h_njet6_VLC10_mass_gj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_mass_gj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_mass_gj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet6_VLC10_mass_gj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('genjet Mass [GeV]')
    h_njet6_VLC10_mass_gj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    canvas_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet6_VLC10_hhz_gj_mass_HZZ_bbbbqq_polm80");
    canvas_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=6"); 
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_mass_gj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"gj1");
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_mass_gj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"gj2");
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_mass_gj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"gj3");
    leg_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet6_VLC10_gj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet6_VLC10_hhz_gj_mass_B_HZZ_bbbbqq_polm80.eps")


    h_njet6_VLC10_mass_rj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_mass_comb_rj1_HHZ_bbbbqq")
    h_njet6_VLC10_mass_rj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet6_VLC10_mass_rj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_mass_rj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet6_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetMaximum(15)
    h_njet6_VLC10_mass_rj1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet6_VLC10_mass_rj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_mass_comb_rj2_HHZ_bbbbqq")
    h_njet6_VLC10_mass_rj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_mass_rj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_mass_rj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet6_VLC10_mass_rj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_mass_rj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet6_VLC10_mass_rj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_mass_comb_rj3_HHZ_bbbbqq")
    h_njet6_VLC10_mass_rj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_mass_rj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_mass_rj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet6_VLC10_mass_rj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_mass_rj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')


    canvas_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet6_VLC10_hhz_comb_rj_mass_HZZ_bbbbqq_polm80");
    canvas_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=6");
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_mass_rj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"j1");
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_mass_rj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j2");
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_mass_rj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j3");
    leg_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet6_VLC10_rj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet6_VLC10_hhz_rj_mass_B_HZZ_bbbbqq_polm80.eps")

    h_njet6_VLC10_mass_rfj_rj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_mass_comb_rfj_rj1_HHZ_bbbbqq")
    h_njet6_VLC10_mass_rfj_rj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet6_VLC10_mass_rfj_rj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_mass_rfj_rj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet6_VLC10_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetMaximum(15)
    h_njet6_VLC10_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet6_VLC10_mass_rfj_rj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_mass_comb_rfj_rj2_HHZ_bbbbqq")
    h_njet6_VLC10_mass_rfj_rj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_mass_rfj_rj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_mass_rfj_rj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet6_VLC10_mass_rfj_rj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_mass_rfj_rj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet6_VLC10_mass_rfj_rj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_mass_comb_rfj_rj3_HHZ_bbbbqq")
    h_njet6_VLC10_mass_rfj_rj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_mass_rfj_rj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_mass_rfj_rj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet6_VLC10_mass_rfj_rj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_mass_rfj_rj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')


    canvas_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet6_VLC10_hhz_comb_rfj_rj_mass_HZZ_bbbbqq_polm80");
    canvas_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=6");
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_mass_rfj_rj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"j1");
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_mass_rfj_rj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j2");
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_mass_rfj_rj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j3");
    leg_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet6_VLC10_rfj_rj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet6_VLC10_hhz_rfj_rj_mass_B_HZZ_bbbbqq_polm80.eps")





    h_njet6_VLC10_dalpha_H1_comb_rfj_rj_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_dalpha_H1_comb_rfj_rj_match_HHZ_bbbbqq")
    h_njet6_VLC10_dalpha_H1_comb_rfj_rj_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_dalpha_H1_comb_rfj_rj_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_dalpha_H1_comb_rfj_rj_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet6_VLC10_dalpha_H1_comb_rfj_rj_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('minimum #Delta#alpha(part,jet) [#circ]')
    h_njet6_VLC10_dalpha_H1_comb_rfj_rj_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet6_VLC10_dalpha_H1_comb_rfj_rj_HZZ_bbbbqq_polm80.SetMaximum(60)
    h_njet6_VLC10_dalpha_H1_comb_rfj_rj_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet6_VLC10_dalpha_H2_comb_rfj_rj_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_dalpha_H2_comb_rfj_rj_match_HHZ_bbbbqq")
    h_njet6_VLC10_dalpha_H2_comb_rfj_rj_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_dalpha_H2_comb_rfj_rj_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_dalpha_H2_comb_rfj_rj_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet6_VLC10_dalpha_H2_comb_rfj_rj_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_dalpha_H2_comb_rfj_rj_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet6_VLC10_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq")
    h_njet6_VLC10_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq.Rebin(4)
    h_njet6_VLC10_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq.SetLineWidth(2)
    h_njet6_VLC10_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq.SetLineColor(4)
    h_njet6_VLC10_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq.GetYaxis().SetTitle('Events')


    canvas_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet6_VLC10_hhz_rfj_rj_dalpha_H_Z_vs_rfj_jets_HZZ_bbbbqq_polm80.eps");
    canvas_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80.cd()

    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=6");
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_dalpha_H1_comb_rfj_rj_HZZ_bbbbqq_polm80.DrawCopy("h"),"H1-matched jet");
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_dalpha_H2_comb_rfj_rj_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"H2-matched jet");
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_dalpha_Z_comb_rfj_rj_match_HHZ_bbbbqq.DrawCopy("h,same"),"Z-matched jet");
    leg_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet6_VLC10_rfj_rj_dalpha_parton_comb_jet_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet6_VLC10_hhz_rfj_rj_dalpha_H_Z_vs_rfj_jets_HZZ_bbbbqq_polm80.eps")







    h_njet6_VLC10_P_reco_over_P_parton_H1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_H1_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq")
    h_njet6_VLC10_P_reco_over_P_parton_H1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_P_reco_over_P_parton_H1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_P_reco_over_P_parton_H1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet6_VLC10_P_reco_over_P_parton_H1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('p_{reco}/p_{parton}')
    h_njet6_VLC10_P_reco_over_P_parton_H1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet6_VLC10_P_reco_over_P_parton_H1_HZZ_bbbbqq_polm80.SetMaximum(25)
    h_njet6_VLC10_P_reco_over_P_parton_H1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet6_VLC10_P_reco_over_P_parton_H2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_H2_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq")
    h_njet6_VLC10_P_reco_over_P_parton_H2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_P_reco_over_P_parton_H2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_P_reco_over_P_parton_H2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet6_VLC10_P_reco_over_P_parton_H2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_P_reco_over_P_parton_H2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet6_VLC10_P_reco_over_P_parton_Z_HHZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_Z_comb_rfj_rj_match_P_reco_over_P_part_HHZ_bbbbqq")
    h_njet6_VLC10_P_reco_over_P_parton_Z_HHZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_P_reco_over_P_parton_Z_HHZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_P_reco_over_P_parton_Z_HHZ_bbbbqq_polm80.SetLineColor(4)
    h_njet6_VLC10_P_reco_over_P_parton_Z_HHZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_P_reco_over_P_parton_Z_HHZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')


    canvas_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet6_VLC10_hhz_rfj_rj_P_reco_over_P_parton_H_Z_vs_rfj_jets_HZZ_bbbbqq_polm80.eps");
    canvas_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80.cd()

    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test=TLegend(0.50,0.61,0.90,0.87);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=6");
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_P_reco_over_P_parton_H1_HZZ_bbbbqq_polm80.DrawCopy("h"),"H1-matched jet");
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_P_reco_over_P_parton_H2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"H2-matched jet");
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_P_reco_over_P_parton_Z_HHZ_bbbbqq_polm80.DrawCopy("h,same"),"Z-matched jet");
    leg_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet6_VLC10_rfj_rj_P_reco_over_P_partonHZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet6_VLC10_hhz_rfj_rj_reco_P_over_parton_P_rfj_jets_HZZ_bbbbqq_polm80.eps")


    h_njet6_VLC10_BTag_rjf_rj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_BTagMax_comb_rfj_rj1_HHZ_bbbbqq")
    h_njet6_VLC10_BTag_rjf_rj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_BTag_rjf_rj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_BTag_rjf_rj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet6_VLC10_BTag_rjf_rj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('BTag')
    h_njet6_VLC10_BTag_rjf_rj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet6_VLC10_BTag_rjf_rj1_HZZ_bbbbqq_polm80.SetMaximum(20)
    h_njet6_VLC10_BTag_rjf_rj1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet6_VLC10_BTag_rjf_rj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_BTagMax_comb_rfj_rj2_HHZ_bbbbqq")
    h_njet6_VLC10_BTag_rjf_rj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_BTag_rjf_rj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_BTag_rjf_rj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet6_VLC10_BTag_rjf_rj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_BTag_rjf_rj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet6_VLC10_BTag_rjf_rj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_BTagMax_comb_rfj_rj3_HHZ_bbbbqq")
    h_njet6_VLC10_BTag_rjf_rj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_BTag_rjf_rj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_BTag_rjf_rj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet6_VLC10_BTag_rjf_rj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_BTag_rjf_rj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')


    canvas_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet6_VLC10_hhz_BTagMax_rfj_rj_HZZ_bbbbqq_polm80");
    canvas_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80.cd()

    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.SetHeader("HHZ Events, N_{jet}=6");
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_BTag_rjf_rj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"j1");
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_BTag_rjf_rj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j2");
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC10_BTag_rjf_rj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"j3");
    leg_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet6_VLC10_BTagMax_rjf_rj_B_HZZ_bbbbqq_polm80.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/canvas_njet6_VLC10_hhz_BTagMax_rfj_rj_B_HZZ_bbbbqq_polm80.eps")






    h_njet6_VLC10_dalpha_comb_rfj_rj1_Hs_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_dalpha_comb_rfj_rj1_Hs_HHZ_bbbbqq")
    h_njet6_VLC10_dalpha_comb_rfj_rj1_Hs_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_dalpha_comb_rfj_rj1_Hs_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_dalpha_comb_rfj_rj1_Hs_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet6_VLC10_dalpha_comb_rfj_rj1_Hs_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('#Delta#alpha(part,jet) [#circ]')
    h_njet6_VLC10_dalpha_comb_rfj_rj1_Hs_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet6_VLC10_dalpha_comb_rfj_rj1_Hs_HZZ_bbbbqq_polm80.SetMaximum(60)
    h_njet6_VLC10_dalpha_comb_rfj_rj1_Hs_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet6_VLC10_alpha_comb_rfj_rj2_Hs_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_dalpha_comb_rfj_rj2_Hs_HHZ_bbbbqq")
    h_njet6_VLC10_alpha_comb_rfj_rj2_Hs_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC10_alpha_comb_rfj_rj2_Hs_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC10_alpha_comb_rfj_rj2_Hs_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet6_VLC10_alpha_comb_rfj_rj2_Hs_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_alpha_comb_rfj_rj2_Hs_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet6_VLC10_dalpha_comb_rfj_rj3_Z_HHZ_bbbbqq=file_polm80_HHZ_SignalHistos_VLC10_njet6_.Get("h_dalpha_comb_rfj_rj3_Z_HHZ_bbbbqq")
    h_njet6_VLC10_dalpha_comb_rfj_rj3_Z_HHZ_bbbbqq.Rebin(4)
    h_njet6_VLC10_dalpha_comb_rfj_rj3_Z_HHZ_bbbbqq.SetLineWidth(2)
    h_njet6_VLC10_dalpha_comb_rfj_rj3_Z_HHZ_bbbbqq.SetLineColor(4)
    h_njet6_VLC10_dalpha_comb_rfj_rj3_Z_HHZ_bbbbqq.GetXaxis().SetTitle('jet mass [GeV]')
    h_njet6_VLC10_dalpha_comb_rfj_rj3_Z_HHZ_bbbbqq.GetYaxis().SetTitle('Events')




    return None


process_files()
#root.gApplication.Run()

#
