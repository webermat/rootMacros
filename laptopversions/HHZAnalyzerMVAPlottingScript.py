from ROOT import gROOT, TCanvas,THStack, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TGraph2D, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor,TMVA,TCut,TString,TDirectory,TLatex,TColor,kBlack, kBlue,kRed,kCyan,kGreen,kOrange,kYellow,kWhite,kViolet,TLine,TArrow,TMath
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



#dm35
#BDT_GiniIndexNormNumEventsMaxDepth3NTrees400AdaBoostBeta020NCuts20.root --> go for this one
#0.35 polm significance,purity  20.9102981094 55.7510280984 784.273549631 229.694612712 324.408336341 68.366258144
#0.2 polp significance,purity  10.4412052721 63.7540377843 170.999000728 45.9298094958 42.9342992902 8.3536288621
#dm40
#BDT_GiniIndexNormNumEventsMaxDepth3NTrees400AdaBoostBeta020NCuts20.root
#0.35 polm significance,purity  20.725849016 52.0769530315 824.857816041 289.927689523 389.484276652 79.6512577534
#0.2 polm significance,purity  10.487648454 60.2347424738 182.603536725 57.2681302428 53.409231782 9.87228033846
#dm45
#BDT_GiniIndexNormNumEventsMaxDepth3NTrees400AdaBoostBeta020NCuts20.root
#0.35 polm significance,purity  20.4781264333 48.3653113549 867.054610983 339.372769147 491.468963504 94.823658973
#0.25 polm significance,purity  10.308988609 66.2359806654 160.449418992 37.2819378823 37.6321724653 6.87553374973

#dm35
#BDT_GiniIndexNormNumEventsMaxDepth2NTrees400AdaBoostBeta020NCuts20.root
# 0.35 polm significance,purity  20.5903012301 50.6065111975 837.75880754 314.6502226 431.249436677 71.7783291936
# 0.25 polp significance,purity  10.4104523399 68.1162298383 159.106747657 35.7445385903 32.5886862874 6.14128707821
#dm40
#BDT_GiniIndexNormNumEventsMaxDepth2NTrees400AdaBoostBeta020NCuts20.root
# 0.4 polm significance,purity  20.4845410358 56.0589276914 748.527377754 227.896589249 304.982672751 53.8443337828
# 0.3 polp significance,purity  10.302970444 71.4547798783 148.557171613 31.1323400289 22.501712501 5.71253242326
#dm45
#BDT_GiniIndexNormNumEventsMaxDepth2NTrees400AdaBoostBeta020NCuts20.root
# 0.4 polm significance,purity  20.3958685867 53.8349295933 772.716633141 258.912124842 346.747836173 56.9677798003
# 0.25 polp significance,purity  10.3856654545 64.5621957704 167.066881239 42.8550107479 40.8651757836 7.98182824627

#dm35
#BDT_GiniIndexNormNumEventsMaxDepth4NTrees400AdaBoostBeta020NCuts20.root
# 0.3 polm significance,purity  21.1480786181 52.0658424541 858.991630897 276.44265312 428.335584819 86.0481124148
# 0.2 polp significance,purity  10.5091721536 67.1084955895 164.573350042 37.4741127789 35.9510105848 7.23627458912
#dm40
#BDT_GiniIndexNormNumEventsMaxDepth4NTrees400AdaBoostBeta020NCuts20.root
# 0.35 polm significance,purity  20.9473768707 57.7294527197 760.084457919 200.027578145 291.384716392 65.1352648437
# 0.2 polp significance,purity  10.5884401794 65.8602518287 170.231759399 42.4706608951 37.3735321164 8.39824814667



#dm35 BDT_GiniIndexNormNumEventsMaxDepth3NTrees400AdaBoostBeta015
# 0.35 polm significance,purity  20.8748205593 49.6874065099 876.999151275 322.741238564 477.871001542 87.4216680527
# 0.25 polp significance,purity  10.4542843877 64.6756031637 168.984990805 43.4315354824 41.2531368732 7.61122567378
#dm40 BDT_GiniIndexNormNumEventsMaxDepth3NTrees400AdaBoostBeta015
# 0.4 polm significance,purity  20.5486061101 52.1155302006 810.209953621 274.644658625 395.311976969 74.4754733369
# 0.25 polm significance,purity  10.5173378211 62.1094240928 178.095991805 52.0794066936 47.4605035782 9.10961908419
#dm45 BDT_GiniIndexNormNumEventsMaxDepth3NTrees400AdaBoostBeta015
# 0.4 polm significance,purity  20.3709062917 49.8215205263 832.920831725 317.796741396 441.933556557 79.1582067832
# 0.3 polm significance,purity  10.4470464859 67.6579526386 161.312567145 34.014964208 36.4682905674 6.62783834973


#dm35 BDT_GiniIndexNormNumEventsMaxDepth3NTrees400AdaBoostBeta020NCuts25.root
# 0.3 polm significance,purity  20.7427773628 47.6657449346 902.666712359 347.913264602 547.803367138 95.3596648425
# 0.2 polm significance,purity  10.4936469096 63.8233368175 172.533482194 45.161109969 44.2275008559 8.40768100621
#dm40
#BDT_GiniIndexNormNumEventsMaxDepth3NTrees400AdaBoostBeta020NCuts25.root
#0.3 polm significance,purity  20.5449648868 44.401101572 950.642140076 414.439372033 658.529595315 117.421674237
#0.2 polm significance,purity  10.5212789865 60.9095820697 181.740389198 53.6168064773 52.8919504881 10.1281912316
#dm45
#BDT_GiniIndexNormNumEventsMaxDepth3NTrees400AdaBoostBeta020NCuts25.root
#0.4 polm significance,purity  20.757054158 61.1387116603 704.717658609 172.608067304 224.366221964 50.9618118629
#0.25 polm significance,purity  10.3150896444 66.6330312014 159.682176322 34.2071390897 38.6667333841 7.08812383562

#dm35 BDT_GiniIndexNormNumEventsMaxDepth3NTrees400AdaBoostBeta020NCuts15.root
# 0.35 polm significance,purity  20.6808073384 55.0528304982 776.882475063 239.134112686 327.3221789 67.8193300962
# 0.2 polm significance,purity  10.4462714467 61.6714291902 176.945124462 50.5420080721 49.9175882936 9.5111602949
#dm40
# 0.35 polm significance,purity  20.6849263859 52.8443536635 809.672462478 260.710146785 387.541696846 74.2591843605
# 0.25 polm significance,purity  10.493020931 67.8513508038 162.27162312 34.5914888829 36.3389704227 5.95546005148
#dm45
# 0.4 polm significance,purity  20.4875660688 56.1057180873 748.124037504 227.896595597 298.183704078 59.2141859382
# 0.25 polm significance,purity  10.572066337 68.7555639741 162.559333637 30.3636408895 37.1148920655 6.3929387728


def process_files():

    CLICdpStyle()

    l = TLatex();
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextColor(kBlack);
    l.SetTextSize(0.045);

    label="CLICdp"
    x=0.17
    x0=0.35
    y=0.93
    label2="L=4 ab^{-1},e^{-} pol -80%"
    label3="L=1 ab^{-1},e^{-} pol +80%"
    label4="L=4 ab^{-1},e^{-} pol -80% + L=1 ab^{-1},e^{-} pol +80%"
    x2=0.62
    x3=0.30
    y2=0.93
    label2_noLumi="e^{-} pol -80%"
    label3_noLumi="e^{-} pol +80%"
    x2_noLumi=0.75
    x2_noLumi_normed=0.67
    y2_noLumi=0.93

    directory=[];

    
    #NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTag_noLTagMins.root
    
    directory.append("-0.2")
    #directory.append("-0.15")
    directory.append("-0.1")
    #directory.append("-0.05")
    directory.append("0.0")
    #directory.append("0.05")
    directory.append("0.1")
    directory.append("0.15")
    #directory.append("0.175")
    directory.append("0.2")
    #directory.append("0.225")
    directory.append("0.25")
    #directory.append("0.275")
    directory.append("0.3")
    directory.append("0.325")
    directory.append("0.35")
    directory.append("0.375")
    directory.append("0.385")
    directory.append("0.4")
    directory.append("0.41")
    directory.append("0.425")
    directory.append("0.435")
    directory.append("0.44")
    directory.append("0.45")
    directory.append("0.46")
    directory.append("0.47")
    directory.append("0.48")
    directory.append("0.49")
    directory.append("0.5")
    directory.append("0.51")
    directory.append("0.525")
    directory.append("0.55")
    directory.append("0.6")
    directory.append("0.625")
    directory.append("0.63")
    directory.append("0.64")
    directory.append("0.65")
    directory.append("0.675")
    directory.append("0.7")
    directory.append("0.725")
    directory.append("0.735")
    directory.append("0.74")
    directory.append("0.745")
    directory.append("0.75")
    directory.append("0.76")
    directory.append("0.77")
    directory.append("0.775")
    directory.append("0.8")
    directory.append("0.81")
    directory.append("0.825")
    directory.append("0.83")
    directory.append("0.85")
    directory.append("0.875")
    directory.append("0.9")
    directory.append("0.925")
    directory.append("0.95")
    directory.append("0.975")
    directory.append("1.0")
  

    n_graphs=int(len(directory))
    
    file_polm80_preselection_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/testfile_preselectionhistos_polm80.root")
 
 
 
    h_comb_jet1_mass_HHZ_polm80_AllEvents_massCuts_=file_polm80_preselection_.Get("hhqq_AllEvents/h_comb_jet1_mass");
    h_comb_jet1_mass_HHZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("H1 candidate mass [GeV]");
    h_comb_jet1_mass_HHZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_comb_jet1_mass_HHZ_polm80_massCuts_=file_polm80_preselection_.Get("hhqq/h_comb_jet1_mass");
    h_comb_jet1_mass_HHZ_polm80_massCuts_.GetXaxis().SetTitle("H1 candidate mass [GeV]");
    h_comb_jet1_mass_HHZ_polm80_massCuts_.SetLineWidth(3);
    h_comb_jet1_mass_HHZ_polm80_massCuts_.SetFillColor(kWhite)
    h_comb_jet1_mass_HHZ_polm80_massCuts_.SetLineColor(kBlack)
 
 
    h_comb_jet1_mass_hzqq_polm80_massCuts_=file_polm80_preselection_.Get("hzqq/h_comb_jet1_mass");
    h_comb_jet1_mass_hzqq_polm80_massCuts_.GetXaxis().SetTitle("comb jet1 mass [GeV]");
    h_comb_jet1_mass_hzqq_polm80_massCuts_.SetLineWidth(3);
    h_comb_jet1_mass_hzqq_polm80_massCuts_.SetFillColor(kCyan+1)
    h_comb_jet1_mass_hzqq_polm80_massCuts_.SetLineColor(kCyan+1)
    h_comb_jet1_mass_hzqq_polm80_massCuts_.SetFillStyle(3002)
    
    h_comb_jet1_mass_ee_qq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qq/h_comb_jet1_mass")
    h_comb_jet1_mass_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_comb_jet1_mass_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_comb_jet1_mass_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_comb_jet1_mass_ee_qqqq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qqqq/h_comb_jet1_mass");
    h_comb_jet1_mass_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_comb_jet1_mass_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_comb_jet1_mass_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qqqqqq/h_comb_jet1_mass");
    h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);
    
    h_comb_jet1_mass_WWH_qqqqH_polm80_massCuts_=file_polm80_preselection_.Get("WWH_qqqqH/h_comb_jet1_mass");
    h_comb_jet1_mass_WWH_qqqqH_polm80_massCuts_.SetFillColor(kOrange);
    h_comb_jet1_mass_WWH_qqqqH_polm80_massCuts_.SetLineColor(kOrange);
    h_comb_jet1_mass_WWH_qqqqH_polm80_massCuts_.SetFillStyle(3002);
    
    h_comb_jet1_mass_ZZH_qqqqH_polm80_massCuts_=file_polm80_preselection_.Get("ZZH_qqqqH/h_comb_jet1_mass");
    h_comb_jet1_mass_ZZH_qqqqH_polm80_massCuts_.SetFillColor(kViolet+2);
    h_comb_jet1_mass_ZZH_qqqqH_polm80_massCuts_.SetLineColor(kViolet+2);
    h_comb_jet1_mass_ZZH_qqqqH_polm80_massCuts_.SetFillStyle(3002);
    
    h_tot_norm_comb_jet1_mass_BG = h_comb_jet1_mass_ee_qq_polm80_massCuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_comb_jet1_mass_BG.Add(h_comb_jet1_mass_ee_qqqq_polm80_massCuts_);
    h_tot_norm_comb_jet1_mass_BG.Add(h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_);
    h_tot_norm_comb_jet1_mass_BG.Add(h_comb_jet1_mass_hzqq_polm80_massCuts_);
    norm_tot_BG_to_SIG = h_comb_jet1_mass_HHZ_polm80_massCuts_.Integral(0,h_comb_jet1_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)/(h_comb_jet1_mass_hzqq_polm80_massCuts_.Integral(0,h_comb_jet1_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_comb_jet1_mass_ee_qq_polm80_massCuts_.Integral(0,h_comb_jet1_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_comb_jet1_mass_ee_qqqq_polm80_massCuts_.Integral(0,h_comb_jet1_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_.Integral(0,h_comb_jet1_mass_HHZ_polm80_massCuts_.GetNbinsX()+1))
    h_tot_norm_comb_jet1_mass_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_comb_jet1_mass_BG.SetLineColor(kBlack)
    h_tot_norm_comb_jet1_mass_BG.SetFillColor(0)
    
    print 'scale or range, norm here ',h_tot_norm_comb_jet1_mass_BG.Integral(),h_comb_jet1_mass_HHZ_polm80_massCuts_.Integral(),norm_tot_BG_to_SIG
    #0-1650 with 55 bins, cut to 150 GeV
    #h_comb_jet1_mass_hzqq_polm80_massCuts_.Rebin()
    #h_comb_jet1_mass_hzqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    #h_comb_jet1_mass_ee_qq_polm80_massCuts_.Scale(nor_tot_BG_to_SIG)
    #h_comb_jet1_mass_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    #h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    
    #h_comb_jet1_mass_ee_qq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_comb_jet1_mass_ee_qqqq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_comb_jet1_mass_HZ_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    
    hhqq_BG_comb_jet1_mass_polm80_massCuts_= THStack("hhqq_BG_comb_jet1_mass_polm80_massCuts_", "");
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.Add(h_comb_jet1_mass_ZZH_qqqqH_polm80_massCuts_);
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.Add(h_comb_jet1_mass_WWH_qqqqH_polm80_massCuts_);
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.Add(h_comb_jet1_mass_hzqq_polm80_massCuts_);
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.Add(h_comb_jet1_mass_ee_qq_polm80_massCuts_);
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.Add(h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_);
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.Add(h_comb_jet1_mass_ee_qqqq_polm80_massCuts_);
    
    canvas_h_SIG_BG_comb_jet1_mass_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_comb_jet1_mass_polm80_massCuts_thstack");
    canvas_h_SIG_BG_comb_jet1_mass_polm80_massCuts_thstack.cd();
    #h_comb_jet1_mass_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_comb_jet1_mass_BG.Draw("hist,e")
    #h_tot_norm_comb_jet1_mass_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.Draw("hist");
    #hhqq_BG_comb_jet1_mass_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.GetXaxis().SetTitle("H1 candidate mass [GeV]");
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.GetYaxis().SetTitle("Events");
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.SetMaximum(180000.)
    h_comb_jet1_mass_HHZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_comb_jet1_mass_polm80_massCuts_thstack.Modified();

    h_comb_jet1_mass_HHZ_polm80_massCuts_.Scale(50000.)
    
    line = TLine(75,0,75,90000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(79,80000,100,80000,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_=TLegend(0.20,0.60,0.590,0.88);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_HHZ_polm80_massCuts_,"HHZ#rightarrowbbbbqq x 50000");
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_hzqq_polm80_massCuts_,"HZ#rightarrowHqq");
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_WWH_qqqqH_polm80_massCuts_,"WWH#rightarrowqqqqH");
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_ZZH_qqqqH_polm80_massCuts_,"ZZH#rightarrowqqqqH");
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_comb_jet1_mass_polm80_massCuts_.Draw();
    
    
    l.DrawLatex(x0,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_comb_jet1_mass_polm80_massCuts_thstack.Print("h_comb_jet1_mass_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000.eps")
    canvas_h_SIG_BG_comb_jet1_mass_polm80_massCuts_thstack.cd()
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.SetMaximum(10.e10)
    hhqq_BG_comb_jet1_mass_polm80_massCuts_.SetMinimum(0.1)
    line = TLine(75,0,75,750000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    arrow = TArrow(79,120000,100,120000,0.025,"|>")
    arrow.Draw();
    canvas_h_SIG_BG_comb_jet1_mass_polm80_massCuts_thstack.SetLogy()
    canvas_h_SIG_BG_comb_jet1_mass_polm80_massCuts_thstack.Print("h_comb_jet1_mass_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000_logy.eps")
    
    h_comb_jet1_mass_HHZ_polm80_massCuts_.Scale(1./50000.)
    
    h_comb_jet1_mass_hzqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_comb_jet1_mass_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_comb_jet1_mass_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    
    hhqq_norm_BG_comb_jet1_mass_polm80_massCuts_= THStack("hhqq_norm_BG_comb_jet1_mass_polm80_massCuts_", "");
    hhqq_norm_BG_comb_jet1_mass_polm80_massCuts_.Add(h_comb_jet1_mass_hzqq_polm80_massCuts_);
    hhqq_norm_BG_comb_jet1_mass_polm80_massCuts_.Add(h_comb_jet1_mass_ee_qq_polm80_massCuts_);
    hhqq_norm_BG_comb_jet1_mass_polm80_massCuts_.Add(h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_);
    hhqq_norm_BG_comb_jet1_mass_polm80_massCuts_.Add(h_comb_jet1_mass_ee_qqqq_polm80_massCuts_);
    
    canvas_h_SIG_norm_BG_comb_jet1_mass_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_norm_BG_comb_jet1_mass_polm80_massCuts_thstack");
    canvas_h_SIG_norm_BG_comb_jet1_mass_polm80_massCuts_thstack.cd();
    #h_comb_jet1_mass_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_comb_jet1_mass_BG.Draw("hist,e")
    #h_tot_norm_comb_jet1_mass_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_comb_jet1_mass_polm80_massCuts_.Draw("hist");
    #hhqq_BG_comb_jet1_mass_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_comb_jet1_mass_polm80_massCuts_.GetXaxis().SetTitle("H1 candidate mass [GeV]");
    hhqq_norm_BG_comb_jet1_mass_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hhqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetMaximum(2.5)
    h_comb_jet1_mass_HHZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_norm_BG_comb_jet1_mass_polm80_massCuts_thstack.Modified();

    line = TLine(75,0,75,1.5)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(79,1.2,100,1.2,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();

    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_=TLegend(0.20,0.63,0.590,0.88);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_HHZ_polm80_massCuts_,"HHZ#rightarrowbbbbqq");
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_hzqq_polm80_massCuts_,"HZ#rightarrowHqq");
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.AddEntry(h_comb_jet1_mass_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_norm_BG_comb_jet1_mass_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_norm_BG_comb_jet1_mass_polm80_massCuts_thstack.Print("h_comb_jet1_mass_polm80_normed_hzqq_ee_qqqqqq_qqqq_qq_and_hhqq.eps")
    
    #comb jet2
    h_comb_jet2_mass_HHZ_polm80_AllEvents_massCuts_=file_polm80_preselection_.Get("hhqq_AllEvents/h_comb_jet2_mass");
    h_comb_jet2_mass_HHZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("H2 candidate mass [GeV]");
    h_comb_jet2_mass_HHZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_comb_jet2_mass_HHZ_polm80_massCuts_=file_polm80_preselection_.Get("hhqq/h_comb_jet2_mass");
    h_comb_jet2_mass_HHZ_polm80_massCuts_.GetXaxis().SetTitle("H2 candidate mass [GeV]");
    h_comb_jet2_mass_HHZ_polm80_massCuts_.SetLineWidth(3);
    h_comb_jet2_mass_HHZ_polm80_massCuts_.SetFillColor(kWhite)
    h_comb_jet2_mass_HHZ_polm80_massCuts_.SetLineColor(kBlack)
    

    h_comb_jet2_mass_hzqq_polm80_massCuts_=file_polm80_preselection_.Get("hzqq/h_comb_jet2_mass");
    h_comb_jet2_mass_hzqq_polm80_massCuts_.GetXaxis().SetTitle("H2 candidate mass [GeV]");
    h_comb_jet2_mass_hzqq_polm80_massCuts_.SetLineWidth(3);
    h_comb_jet2_mass_hzqq_polm80_massCuts_.SetFillColor(kCyan+1)
    h_comb_jet2_mass_hzqq_polm80_massCuts_.SetLineColor(kCyan+1)
    h_comb_jet2_mass_hzqq_polm80_massCuts_.SetFillStyle(3002)
    
    h_comb_jet2_mass_ee_qq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qq/h_comb_jet2_mass")
    h_comb_jet2_mass_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_comb_jet2_mass_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_comb_jet2_mass_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_comb_jet2_mass_ee_qqqq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qqqq/h_comb_jet2_mass");
    h_comb_jet2_mass_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_comb_jet2_mass_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_comb_jet2_mass_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qqqqqq/h_comb_jet2_mass");
    h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);
    
    h_comb_jet2_mass_WWH_qqqqH_polm80_massCuts_=file_polm80_preselection_.Get("WWH_qqqqH/h_comb_jet2_mass");
    h_comb_jet2_mass_WWH_qqqqH_polm80_massCuts_.SetFillColor(kOrange);
    h_comb_jet2_mass_WWH_qqqqH_polm80_massCuts_.SetLineColor(kOrange);
    h_comb_jet2_mass_WWH_qqqqH_polm80_massCuts_.SetFillStyle(3002);
    
    h_comb_jet2_mass_ZZH_qqqqH_polm80_massCuts_=file_polm80_preselection_.Get("ZZH_qqqqH/h_comb_jet2_mass");
    h_comb_jet2_mass_ZZH_qqqqH_polm80_massCuts_.SetFillColor(kViolet+2);
    h_comb_jet2_mass_ZZH_qqqqH_polm80_massCuts_.SetLineColor(kViolet+2);
    h_comb_jet2_mass_ZZH_qqqqH_polm80_massCuts_.SetFillStyle(3002);
    
    h_tot_norm_comb_jet2_mass_BG = h_comb_jet2_mass_ee_qq_polm80_massCuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_comb_jet2_mass_BG.Add(h_comb_jet2_mass_ee_qqqq_polm80_massCuts_);
    h_tot_norm_comb_jet2_mass_BG.Add(h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_);
    h_tot_norm_comb_jet2_mass_BG.Add(h_comb_jet2_mass_hzqq_polm80_massCuts_);
    norm_tot_BG_to_SIG = h_comb_jet2_mass_HHZ_polm80_massCuts_.Integral(0,h_comb_jet2_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)/(h_comb_jet2_mass_hzqq_polm80_massCuts_.Integral(0,h_comb_jet2_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_comb_jet2_mass_ee_qq_polm80_massCuts_.Integral(0,h_comb_jet2_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_comb_jet2_mass_ee_qqqq_polm80_massCuts_.Integral(0,h_comb_jet2_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_.Integral(0,h_comb_jet2_mass_HHZ_polm80_massCuts_.GetNbinsX()+1))
    h_tot_norm_comb_jet2_mass_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_comb_jet2_mass_BG.SetLineColor(kBlack)
    h_tot_norm_comb_jet2_mass_BG.SetFillColor(0)
    
    print 'scale or range, norm here ',h_tot_norm_comb_jet2_mass_BG.Integral(),h_comb_jet2_mass_HHZ_polm80_massCuts_.Integral(),norm_tot_BG_to_SIG
    #0-1650 with 55 bins, cut to 150 GeV
    #h_comb_jet2_mass_hzqq_polm80_massCuts_.Rebin()
    #h_comb_jet2_mass_hzqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    #h_comb_jet2_mass_ee_qq_polm80_massCuts_.Scale(nor_tot_BG_to_SIG)
    #h_comb_jet2_mass_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    #h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    
    #h_comb_jet2_mass_ee_qq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_comb_jet2_mass_ee_qqqq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_comb_jet2_mass_HZ_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    
    hhqq_BG_comb_jet2_mass_polm80_massCuts_= THStack("hhqq_BG_comb_jet2_mass_polm80_massCuts_", "");
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.Add(h_comb_jet2_mass_ZZH_qqqqH_polm80_massCuts_);
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.Add(h_comb_jet2_mass_WWH_qqqqH_polm80_massCuts_);
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.Add(h_comb_jet2_mass_hzqq_polm80_massCuts_);
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.Add(h_comb_jet2_mass_ee_qq_polm80_massCuts_);
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.Add(h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_);
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.Add(h_comb_jet2_mass_ee_qqqq_polm80_massCuts_);
    
    
    canvas_h_SIG_BG_comb_jet2_mass_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_comb_jet2_mass_polm80_massCuts_thstack");
    canvas_h_SIG_BG_comb_jet2_mass_polm80_massCuts_thstack.cd();
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.Draw("hist");
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.GetXaxis().SetTitle("H2 candidate mass [GeV]");
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.GetYaxis().SetTitle("Events");
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.SetMaximum(210000.)
    h_comb_jet2_mass_HHZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_comb_jet2_mass_polm80_massCuts_thstack.Modified();
    
    h_comb_jet2_mass_HHZ_polm80_massCuts_.Scale(50000.)
    canvas_h_SIG_BG_comb_jet2_mass_polm80_massCuts_thstack.Update()
    
    line = TLine(75,0,75,110000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(79,90000,100,90000,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_=TLegend(0.30,0.60,0.690,0.88);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_HHZ_polm80_massCuts_,"HHZ#rightarrowbbbbqq x 50000");
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_hzqq_polm80_massCuts_,"HZ#rightarrowHqq");
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_WWH_qqqqH_polm80_massCuts_,"WWH#rightarrowqqqqH");
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_ZZH_qqqqH_polm80_massCuts_,"ZZH#rightarrowqqqqH");
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_comb_jet2_mass_polm80_massCuts_.Draw();
    
    
    l.DrawLatex(x0,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_comb_jet2_mass_polm80_massCuts_thstack.Print("h_comb_jet2_mass_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000.eps")
    canvas_h_SIG_BG_comb_jet2_mass_polm80_massCuts_thstack.cd()
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.SetMaximum(10.e10)
    hhqq_BG_comb_jet2_mass_polm80_massCuts_.SetMinimum(0.1)
    canvas_h_SIG_BG_comb_jet2_mass_polm80_massCuts_thstack.SetLogy()
    canvas_h_SIG_BG_comb_jet2_mass_polm80_massCuts_thstack.Print("h_comb_jet2_mass_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000_logy.eps")
    
    
    h_comb_jet2_mass_HHZ_polm80_massCuts_.Scale(1./50000.)
    
    h_comb_jet2_mass_hzqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_comb_jet2_mass_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_comb_jet2_mass_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    
    hhqq_norm_BG_comb_jet2_mass_polm80_massCuts_= THStack("hhqq_norm_BG_comb_jet2_mass_polm80_massCuts_", "");
    hhqq_norm_BG_comb_jet2_mass_polm80_massCuts_.Add(h_comb_jet2_mass_hzqq_polm80_massCuts_);
    hhqq_norm_BG_comb_jet2_mass_polm80_massCuts_.Add(h_comb_jet2_mass_ee_qq_polm80_massCuts_);
    hhqq_norm_BG_comb_jet2_mass_polm80_massCuts_.Add(h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_);
    hhqq_norm_BG_comb_jet2_mass_polm80_massCuts_.Add(h_comb_jet2_mass_ee_qqqq_polm80_massCuts_);
    
    canvas_h_SIG_norm_BG_comb_jet2_mass_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_norm_BG_comb_jet2_mass_polm80_massCuts_thstack");
    canvas_h_SIG_norm_BG_comb_jet2_mass_polm80_massCuts_thstack.cd();
    #h_comb_jet2_mass_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_comb_jet2_mass_BG.Draw("hist,e")
    #h_tot_norm_comb_jet2_mass_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_comb_jet2_mass_polm80_massCuts_.Draw("hist");
    #hhqq_BG_comb_jet2_mass_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_comb_jet2_mass_polm80_massCuts_.GetXaxis().SetTitle("H2 candidate mass [GeV]");
    hhqq_norm_BG_comb_jet2_mass_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hhqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetMaximum(2.5)
    h_comb_jet2_mass_HHZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_norm_BG_comb_jet2_mass_polm80_massCuts_thstack.Modified();
    
    line = TLine(75,0,75,1.5)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(79,1.4,100,1.4,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_=TLegend(0.20,0.63,0.590,0.88);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_HHZ_polm80_massCuts_,"HHZ#rightarrowbbbbqq");
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_hzqq_polm80_massCuts_,"HZ#rightarrowHqq");
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.AddEntry(h_comb_jet2_mass_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_norm_BG_comb_jet2_mass_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_norm_BG_comb_jet2_mass_polm80_massCuts_thstack.Print("h_comb_jet2_mass_polm80_normed_hzqq_ee_qqqqqq_qqqq_qq_and_hhqq.eps")
    
    #comb jet3
    h_comb_jet3_mass_HHZ_polm80_AllEvents_massCuts_=file_polm80_preselection_.Get("hhqq_AllEvents/h_comb_jet3_mass");
    h_comb_jet3_mass_HHZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("comb jet3 mass [GeV]");
    h_comb_jet3_mass_HHZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_comb_jet3_mass_HHZ_polm80_massCuts_=file_polm80_preselection_.Get("hhqq/h_comb_jet3_mass");
    h_comb_jet3_mass_HHZ_polm80_massCuts_.GetXaxis().SetTitle("Z candidate mass [GeV]");
    h_comb_jet3_mass_HHZ_polm80_massCuts_.SetLineWidth(3);
    h_comb_jet3_mass_HHZ_polm80_massCuts_.SetFillColor(kWhite)
    h_comb_jet3_mass_HHZ_polm80_massCuts_.SetLineColor(kBlack)
    
    
    h_comb_jet3_mass_hzqq_polm80_massCuts_=file_polm80_preselection_.Get("hzqq/h_comb_jet3_mass");
    h_comb_jet3_mass_hzqq_polm80_massCuts_.GetXaxis().SetTitle("Z candidate mass [GeV]");
    h_comb_jet3_mass_hzqq_polm80_massCuts_.SetLineWidth(3);
    h_comb_jet3_mass_hzqq_polm80_massCuts_.SetFillColor(kCyan+1)
    h_comb_jet3_mass_hzqq_polm80_massCuts_.SetLineColor(kCyan+1)
    h_comb_jet3_mass_hzqq_polm80_massCuts_.SetFillStyle(3002)
    
    h_comb_jet3_mass_ee_qq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qq/h_comb_jet3_mass")
    h_comb_jet3_mass_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_comb_jet3_mass_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_comb_jet3_mass_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_comb_jet3_mass_ee_qqqq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qqqq/h_comb_jet3_mass");
    h_comb_jet3_mass_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_comb_jet3_mass_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_comb_jet3_mass_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qqqqqq/h_comb_jet3_mass");
    h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);
    
    h_comb_jet3_mass_WWH_qqqqH_polm80_massCuts_=file_polm80_preselection_.Get("WWH_qqqqH/h_comb_jet3_mass");
    h_comb_jet3_mass_WWH_qqqqH_polm80_massCuts_.SetFillColor(kOrange);
    h_comb_jet3_mass_WWH_qqqqH_polm80_massCuts_.SetLineColor(kOrange);
    h_comb_jet3_mass_WWH_qqqqH_polm80_massCuts_.SetFillStyle(3002);
    
    h_comb_jet3_mass_ZZH_qqqqH_polm80_massCuts_=file_polm80_preselection_.Get("ZZH_qqqqH/h_comb_jet3_mass");
    h_comb_jet3_mass_ZZH_qqqqH_polm80_massCuts_.SetFillColor(kViolet+2);
    h_comb_jet3_mass_ZZH_qqqqH_polm80_massCuts_.SetLineColor(kViolet+2);
    h_comb_jet3_mass_ZZH_qqqqH_polm80_massCuts_.SetFillStyle(3002);
    
    h_tot_norm_comb_jet3_mass_BG = h_comb_jet3_mass_ee_qq_polm80_massCuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_comb_jet3_mass_BG.Add(h_comb_jet3_mass_ee_qqqq_polm80_massCuts_);
    h_tot_norm_comb_jet3_mass_BG.Add(h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_);
    h_tot_norm_comb_jet3_mass_BG.Add(h_comb_jet3_mass_hzqq_polm80_massCuts_);
    norm_tot_BG_to_SIG = h_comb_jet3_mass_HHZ_polm80_massCuts_.Integral(0,h_comb_jet3_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)/(h_comb_jet3_mass_hzqq_polm80_massCuts_.Integral(0,h_comb_jet3_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_comb_jet3_mass_ee_qq_polm80_massCuts_.Integral(0,h_comb_jet3_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_comb_jet3_mass_ee_qqqq_polm80_massCuts_.Integral(0,h_comb_jet3_mass_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_.Integral(0,h_comb_jet3_mass_HHZ_polm80_massCuts_.GetNbinsX()+1))
    h_tot_norm_comb_jet3_mass_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_comb_jet3_mass_BG.SetLineColor(kBlack)
    h_tot_norm_comb_jet3_mass_BG.SetFillColor(0)
    
    print 'scale or range, norm here ',h_tot_norm_comb_jet3_mass_BG.Integral(),h_comb_jet3_mass_HHZ_polm80_massCuts_.Integral(),norm_tot_BG_to_SIG
    #0-1650 with 55 bins, cut to 150 GeV
    #h_comb_jet3_mass_hzqq_polm80_massCuts_.Rebin()
    #h_comb_jet3_mass_hzqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    #h_comb_jet3_mass_ee_qq_polm80_massCuts_.Scale(nor_tot_BG_to_SIG)
    #h_comb_jet3_mass_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    #h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    
    #h_comb_jet3_mass_ee_qq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_comb_jet3_mass_ee_qqqq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_comb_jet3_mass_HZ_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    
    hhqq_BG_comb_jet3_mass_polm80_massCuts_= THStack("hhqq_BG_comb_jet3_mass_polm80_massCuts_", "");
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.Add(h_comb_jet3_mass_ZZH_qqqqH_polm80_massCuts_);
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.Add(h_comb_jet3_mass_WWH_qqqqH_polm80_massCuts_);
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.Add(h_comb_jet3_mass_hzqq_polm80_massCuts_);
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.Add(h_comb_jet3_mass_ee_qq_polm80_massCuts_);
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.Add(h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_);
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.Add(h_comb_jet3_mass_ee_qqqq_polm80_massCuts_);
    
    
    canvas_h_SIG_BG_comb_jet3_mass_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_comb_jet3_mass_polm80_massCuts_thstack");
    canvas_h_SIG_BG_comb_jet3_mass_polm80_massCuts_thstack.cd();
    #h_comb_jet3_mass_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_comb_jet3_mass_BG.Draw("hist,e")
    #h_tot_norm_comb_jet3_mass_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.Draw("hist");
    #hhqq_BG_comb_jet3_mass_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.GetXaxis().SetTitle("Z candidate mass [GeV]");
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.GetYaxis().SetTitle("Events");
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.SetMaximum(530000.)
    h_comb_jet3_mass_HHZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_comb_jet3_mass_polm80_massCuts_thstack.Modified();
    
    h_comb_jet3_mass_HHZ_polm80_massCuts_.Scale(50000.)
    
    line = TLine(50,0,50,275000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(53,225000,75,225000,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    line2 = TLine(150,0,150,275000)
    line2.SetLineColor(kBlack);
    line2.SetLineWidth(2);
    line2.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow2 = TArrow(130,225000,147,225000,0.025,"<|")
    #arrow.SetAngle(0);
    arrow2.SetLineWidth(2);
    arrow2.Draw();
    
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_=TLegend(0.35,0.60,0.740,0.88);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_HHZ_polm80_massCuts_,"HHZ#rightarrowbbbbqq x 50000");
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_hzqq_polm80_massCuts_,"HZ#rightarrowHqq");
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_WWH_qqqqH_polm80_massCuts_,"WWH#rightarrowqqqqH");
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_ZZH_qqqqH_polm80_massCuts_,"ZZH#rightarrowqqqqH");
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_comb_jet3_mass_polm80_massCuts_.Draw();
    
    l.DrawLatex(x0,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_comb_jet3_mass_polm80_massCuts_thstack.Print("h_comb_jet3_mass_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000.eps")
    canvas_h_SIG_BG_comb_jet3_mass_polm80_massCuts_thstack.cd()
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.SetMaximum(10.e10)
    hhqq_BG_comb_jet3_mass_polm80_massCuts_.SetMinimum(0.1)
    canvas_h_SIG_BG_comb_jet3_mass_polm80_massCuts_thstack.SetLogy()
    canvas_h_SIG_BG_comb_jet3_mass_polm80_massCuts_thstack.Print("h_comb_jet3_mass_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000_logy.eps")
    
    h_comb_jet3_mass_HHZ_polm80_massCuts_.Scale(1./50000.)
    
    h_comb_jet3_mass_hzqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_comb_jet3_mass_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_comb_jet3_mass_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    
    hhqq_norm_BG_comb_jet3_mass_polm80_massCuts_= THStack("hhqq_norm_BG_comb_jet3_mass_polm80_massCuts_", "");
    hhqq_norm_BG_comb_jet3_mass_polm80_massCuts_.Add(h_comb_jet3_mass_hzqq_polm80_massCuts_);
    hhqq_norm_BG_comb_jet3_mass_polm80_massCuts_.Add(h_comb_jet3_mass_ee_qq_polm80_massCuts_);
    hhqq_norm_BG_comb_jet3_mass_polm80_massCuts_.Add(h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_);
    hhqq_norm_BG_comb_jet3_mass_polm80_massCuts_.Add(h_comb_jet3_mass_ee_qqqq_polm80_massCuts_);
    
    canvas_h_SIG_norm_BG_comb_jet3_mass_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_norm_BG_comb_jet3_mass_polm80_massCuts_thstack");
    canvas_h_SIG_norm_BG_comb_jet3_mass_polm80_massCuts_thstack.cd();
    #h_comb_jet3_mass_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_comb_jet3_mass_BG.Draw("hist,e")
    #h_tot_norm_comb_jet3_mass_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_comb_jet3_mass_polm80_massCuts_.Draw("hist");
    #hhqq_BG_comb_jet3_mass_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_comb_jet3_mass_polm80_massCuts_.GetXaxis().SetTitle("Z candidate mass [GeV]");
    hhqq_norm_BG_comb_jet3_mass_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hhqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetMaximum(3.5)
    h_comb_jet3_mass_HHZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_norm_BG_comb_jet3_mass_polm80_massCuts_thstack.Modified();
    
    line = TLine(50,0,50,2)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(53,1.8,75,1.8,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    line2 = TLine(150,0,150,2)
    line2.SetLineColor(kBlack);
    line2.SetLineWidth(2);
    line2.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow2 = TArrow(130,1.8,147,1.8,0.025,"<|")
    #arrow.SetAngle(0);
    arrow2.SetLineWidth(2);
    arrow2.Draw();
    
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_=TLegend(0.50,0.63,0.890,0.88);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_HHZ_polm80_massCuts_,"HHZ#rightarrowbbbbqq");
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_hzqq_polm80_massCuts_,"HZ#rightarrowHqq");
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.AddEntry(h_comb_jet3_mass_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_norm_BG_comb_jet3_mass_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_norm_BG_comb_jet3_mass_polm80_massCuts_thstack.Print("h_comb_jet3_mass_polm80_normed_hzqq_ee_qqqqqq_qqqq_qq_and_hhqq.eps")
  
    #jet1 theta
    h_jet1_theta_HHZ_polm80_AllEvents_ECuts_=file_polm80_preselection_.Get("hhqq_AllEvents/h_jet1_theta");
    h_jet1_theta_HHZ_polm80_AllEvents_ECuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    h_jet1_theta_HHZ_polm80_AllEvents_ECuts_.SetLineWidth(3);
    h_jet1_theta_HHZ_polm80_ECuts_=file_polm80_preselection_.Get("hhqq/h_jet1_theta");
    h_jet1_theta_HHZ_polm80_ECuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    h_jet1_theta_HHZ_polm80_ECuts_.SetLineWidth(3);
    h_jet1_theta_HHZ_polm80_ECuts_.SetFillColor(kWhite)
    h_jet1_theta_HHZ_polm80_ECuts_.SetLineColor(kBlack)
  
  
    h_jet1_theta_hzqq_polm80_ECuts_=file_polm80_preselection_.Get("hzqq/h_jet1_theta");
    h_jet1_theta_hzqq_polm80_ECuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    h_jet1_theta_hzqq_polm80_ECuts_.SetLineWidth(3);
    h_jet1_theta_hzqq_polm80_ECuts_.SetFillColor(kCyan+1)
    h_jet1_theta_hzqq_polm80_ECuts_.SetLineColor(kCyan+1)
    h_jet1_theta_hzqq_polm80_ECuts_.SetFillStyle(3002)
    
    h_jet1_theta_ee_qq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qq/h_jet1_theta")
    h_jet1_theta_ee_qq_polm80_ECuts_.SetFillColor(kBlue);
    h_jet1_theta_ee_qq_polm80_ECuts_.SetLineColor(kBlue);
    h_jet1_theta_ee_qq_polm80_ECuts_.SetFillStyle(3002);
    h_jet1_theta_ee_qqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqq/h_jet1_theta");
    h_jet1_theta_ee_qqqq_polm80_ECuts_.SetFillColor(kRed);
    h_jet1_theta_ee_qqqq_polm80_ECuts_.SetLineColor(kRed);
    h_jet1_theta_ee_qqqq_polm80_ECuts_.SetFillStyle(3002);
    h_jet1_theta_ee_qqqqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqqqq/h_jet1_theta");
    h_jet1_theta_ee_qqqqqq_polm80_ECuts_.SetFillColor(kGreen-2);
    h_jet1_theta_ee_qqqqqq_polm80_ECuts_.SetLineColor(kGreen-2);
    h_jet1_theta_ee_qqqqqq_polm80_ECuts_.SetFillStyle(3002);
    
    h_jet1_theta_WWH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("WWH_qqqqH/h_jet1_theta");
    h_jet1_theta_WWH_qqqqH_polm80_ECuts_.SetFillColor(kOrange);
    h_jet1_theta_WWH_qqqqH_polm80_ECuts_.SetLineColor(kOrange);
    h_jet1_theta_WWH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_jet1_theta_ZZH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("ZZH_qqqqH/h_jet1_theta");
    h_jet1_theta_ZZH_qqqqH_polm80_ECuts_.SetFillColor(kViolet+2);
    h_jet1_theta_ZZH_qqqqH_polm80_ECuts_.SetLineColor(kViolet+2);
    h_jet1_theta_ZZH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
      
    h_tot_norm_jet1_theta_BG = h_jet1_theta_ee_qq_polm80_ECuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_jet1_theta_BG.Add(h_jet1_theta_ee_qqqq_polm80_ECuts_);
    h_tot_norm_jet1_theta_BG.Add(h_jet1_theta_ee_qqqqqq_polm80_ECuts_);
    h_tot_norm_jet1_theta_BG.Add(h_jet1_theta_hzqq_polm80_ECuts_);
    norm_tot_BG_to_SIG= h_jet1_theta_HHZ_polm80_ECuts_.Integral(0,h_jet1_theta_HHZ_polm80_ECuts_.GetNbinsX()+1)/(h_jet1_theta_hzqq_polm80_ECuts_.Integral(0,h_jet1_theta_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet1_theta_ee_qq_polm80_ECuts_.Integral(0,h_jet1_theta_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet1_theta_ee_qqqq_polm80_ECuts_.Integral(0,h_jet1_theta_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet1_theta_ee_qqqqqq_polm80_ECuts_.Integral(0,h_jet1_theta_HHZ_polm80_ECuts_.GetNbinsX()+1))
    h_tot_norm_jet1_theta_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_jet1_theta_BG.SetLineColor(kBlack)
    h_tot_norm_jet1_theta_BG.SetFillColor(0)
      
    print 'scale or range, norm here ',h_tot_norm_jet1_theta_BG.Integral(),h_jet1_theta_HHZ_polm80_ECuts_.Integral(),norm_tot_BG_to_SIG
    #0-1650 with 55 bins, cut to 150 GeV
    #h_jet1_theta_hzqq_polm80_ECuts_.Rebin()
    #h_jet1_theta_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet1_theta_ee_qq_polm80_ECuts_.Scale(nor_tot_BG_to_SIG)
    #h_jet1_theta_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet1_theta_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
      
    #h_jet1_theta_ee_qq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet1_theta_ee_qqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet1_theta_ee_qqqqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet1_theta_HZ_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    
    hhqq_BG_jet1_theta_polm80_ECuts_= THStack("hhqq_BG_jet1_theta_polm80_ECuts_", "");
    hhqq_BG_jet1_theta_polm80_ECuts_.Add(h_jet1_theta_ZZH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet1_theta_polm80_ECuts_.Add(h_jet1_theta_WWH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet1_theta_polm80_ECuts_.Add(h_jet1_theta_hzqq_polm80_ECuts_);
    hhqq_BG_jet1_theta_polm80_ECuts_.Add(h_jet1_theta_ee_qq_polm80_ECuts_);
    hhqq_BG_jet1_theta_polm80_ECuts_.Add(h_jet1_theta_ee_qqqqqq_polm80_ECuts_);
    hhqq_BG_jet1_theta_polm80_ECuts_.Add(h_jet1_theta_ee_qqqq_polm80_ECuts_);
      
    canvas_h_SIG_BG_jet1_theta_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_theta_polm80_ECuts_thstack");
    canvas_h_SIG_BG_jet1_theta_polm80_ECuts_thstack.cd();
    #h_jet1_theta_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet1_theta_BG.Draw("hist,e")
    #h_tot_norm_jet1_theta_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet1_theta_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet1_theta_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet1_theta_polm80_ECuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    hhqq_BG_jet1_theta_polm80_ECuts_.GetYaxis().SetTitle("Events");
    hhqq_BG_jet1_theta_polm80_ECuts_.GetYaxis().SetTitleOffset(1.5);
    hhqq_BG_jet1_theta_polm80_ECuts_.SetMaximum(750000.)
    h_jet1_theta_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_theta_polm80_ECuts_thstack.Modified();
    
    h_jet1_theta_HHZ_polm80_ECuts_.Scale(50000.)
    
    line = TLine(10,0,10,725000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(13,650000,25,650000,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    line2 = TLine(170,0,170,725000)
    line2.SetLineColor(kBlack);
    line2.SetLineWidth(2);
    line2.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow2 = TArrow(155,650000,167,650000,0.025,"<|")
    #arrow.SetAngle(0);
    arrow2.SetLineWidth(2);
    arrow2.Draw();
    
    leg_hzqq_BG_jet1_theta_polm80_ECuts_=TLegend(0.30,0.60,0.690,0.88);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq x 50000");
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_WWH_qqqqH_polm80_ECuts_,"WWH#rightarrowqqqqH");
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_ZZH_qqqqH_polm80_ECuts_,"ZZH#rightarrowqqqqH");
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_theta_polm80_ECuts_.Draw();
    
    l.DrawLatex(x0,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet1_theta_polm80_ECuts_thstack.Print("h_jet1_theta_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_50000.eps")
    canvas_h_SIG_BG_jet1_theta_polm80_ECuts_thstack.cd()
    hhqq_BG_jet1_theta_polm80_ECuts_.SetMaximum(10.e9)
    hhqq_BG_jet1_theta_polm80_ECuts_.SetMinimum(0.1)
    canvas_h_SIG_BG_jet1_theta_polm80_ECuts_thstack.SetLogy()
    canvas_h_SIG_BG_jet1_theta_polm80_ECuts_thstack.Print("h_jet1_theta_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000_logy.eps")
    
    h_jet1_theta_HHZ_polm80_ECuts_.Scale(1./50000.)
    
    h_jet1_theta_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_theta_ee_qq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_theta_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_theta_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    
    hhqq_norm_BG_jet1_theta_polm80_ECuts_= THStack("hhqq_norm_BG_jet1_theta_polm80_ECuts_", "");
    hhqq_norm_BG_jet1_theta_polm80_ECuts_.Add(h_jet1_theta_hzqq_polm80_ECuts_);
    hhqq_norm_BG_jet1_theta_polm80_ECuts_.Add(h_jet1_theta_ee_qq_polm80_ECuts_);
    hhqq_norm_BG_jet1_theta_polm80_ECuts_.Add(h_jet1_theta_ee_qqqqqq_polm80_ECuts_);
    hhqq_norm_BG_jet1_theta_polm80_ECuts_.Add(h_jet1_theta_ee_qqqq_polm80_ECuts_);
      
    canvas_h_SIG_norm_BG_jet1_theta_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_norm_BG_jet1_theta_polm80_ECuts_thstack");
    canvas_h_SIG_norm_BG_jet1_theta_polm80_ECuts_thstack.cd();
    #h_jet1_theta_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet1_theta_BG.Draw("hist,e")
    #h_tot_norm_jet1_theta_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet1_theta_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet1_theta_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet1_theta_polm80_ECuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    hhqq_norm_BG_jet1_theta_polm80_ECuts_.GetYaxis().SetTitle("A.U.");
    hhqq_norm_BG_jet1_theta_polm80_ECuts_.SetMaximum(5.5)
    h_jet1_theta_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_norm_BG_jet1_theta_polm80_ECuts_thstack.Modified();
      
    line = TLine(10,0,10,5.40)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(13,4.5,25,4.5,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    line2 = TLine(170,0,170,5.40)
    line2.SetLineColor(kBlack);
    line2.SetLineWidth(2);
    line2.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow2 = TArrow(155,4.5,167,4.5,0.025,"<|")
    #arrow.SetAngle(0);
    arrow2.SetLineWidth(2);
    arrow2.Draw();
      
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_=TLegend(0.30,0.63,0.690,0.88);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq");
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.AddEntry(h_jet1_theta_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_norm_BG_jet1_theta_polm80_ECuts_.Draw();
      
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_norm_BG_jet1_theta_polm80_ECuts_thstack.Print("h_jet1_theta_polm80_normed_hzqq_ee_qqqqqq_qqqq_qq_and_hhqq.eps")
    
    #comb jet2
    h_jet2_theta_HHZ_polm80_AllEvents_ECuts_=file_polm80_preselection_.Get("hhqq_AllEvents/h_jet2_theta");
    h_jet2_theta_HHZ_polm80_AllEvents_ECuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    h_jet2_theta_HHZ_polm80_AllEvents_ECuts_.SetLineWidth(3);
    h_jet2_theta_HHZ_polm80_ECuts_=file_polm80_preselection_.Get("hhqq/h_jet2_theta");
    h_jet2_theta_HHZ_polm80_ECuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    h_jet2_theta_HHZ_polm80_ECuts_.SetLineWidth(3);
    h_jet2_theta_HHZ_polm80_ECuts_.SetFillColor(kWhite)
    h_jet2_theta_HHZ_polm80_ECuts_.SetLineColor(kBlack)
    
      
    h_jet2_theta_hzqq_polm80_ECuts_=file_polm80_preselection_.Get("hzqq/h_jet2_theta");
    h_jet2_theta_hzqq_polm80_ECuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    h_jet2_theta_hzqq_polm80_ECuts_.SetLineWidth(3);
    h_jet2_theta_hzqq_polm80_ECuts_.SetFillColor(kCyan+1)
    h_jet2_theta_hzqq_polm80_ECuts_.SetLineColor(kCyan+1)
    h_jet2_theta_hzqq_polm80_ECuts_.SetFillStyle(3002)
    
    h_jet2_theta_ee_qq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qq/h_jet2_theta")
    h_jet2_theta_ee_qq_polm80_ECuts_.SetFillColor(kBlue);
    h_jet2_theta_ee_qq_polm80_ECuts_.SetLineColor(kBlue);
    h_jet2_theta_ee_qq_polm80_ECuts_.SetFillStyle(3002);
    h_jet2_theta_ee_qqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqq/h_jet2_theta");
    h_jet2_theta_ee_qqqq_polm80_ECuts_.SetFillColor(kRed);
    h_jet2_theta_ee_qqqq_polm80_ECuts_.SetLineColor(kRed);
    h_jet2_theta_ee_qqqq_polm80_ECuts_.SetFillStyle(3002);
    h_jet2_theta_ee_qqqqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqqqq/h_jet2_theta");
    h_jet2_theta_ee_qqqqqq_polm80_ECuts_.SetFillColor(kGreen-2);
    h_jet2_theta_ee_qqqqqq_polm80_ECuts_.SetLineColor(kGreen-2);
    h_jet2_theta_ee_qqqqqq_polm80_ECuts_.SetFillStyle(3002);

    h_jet2_theta_WWH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("WWH_qqqqH/h_jet2_theta");
    h_jet2_theta_WWH_qqqqH_polm80_ECuts_.SetFillColor(kOrange);
    h_jet2_theta_WWH_qqqqH_polm80_ECuts_.SetLineColor(kOrange);
    h_jet2_theta_WWH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_jet2_theta_ZZH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("ZZH_qqqqH/h_jet2_theta");
    h_jet2_theta_ZZH_qqqqH_polm80_ECuts_.SetFillColor(kViolet+2);
    h_jet2_theta_ZZH_qqqqH_polm80_ECuts_.SetLineColor(kViolet+2);
    h_jet2_theta_ZZH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_tot_norm_jet2_theta_BG = h_jet2_theta_ee_qq_polm80_ECuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_jet2_theta_BG.Add(h_jet2_theta_ee_qqqq_polm80_ECuts_);
    h_tot_norm_jet2_theta_BG.Add(h_jet2_theta_ee_qqqqqq_polm80_ECuts_);
    h_tot_norm_jet2_theta_BG.Add(h_jet2_theta_hzqq_polm80_ECuts_);
    norm_tot_BG_to_SIG= h_jet2_theta_HHZ_polm80_ECuts_.Integral(0,h_jet2_theta_HHZ_polm80_ECuts_.GetNbinsX()+1)/(h_jet2_theta_hzqq_polm80_ECuts_.Integral(0,h_jet2_theta_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet2_theta_ee_qq_polm80_ECuts_.Integral(0,h_jet2_theta_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet2_theta_ee_qqqq_polm80_ECuts_.Integral(0,h_jet2_theta_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet2_theta_ee_qqqqqq_polm80_ECuts_.Integral(0,h_jet2_theta_HHZ_polm80_ECuts_.GetNbinsX()+1))
    h_tot_norm_jet2_theta_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_jet2_theta_BG.SetLineColor(kBlack)
    h_tot_norm_jet2_theta_BG.SetFillColor(0)
    
    print 'scale or range, norm here ',h_tot_norm_jet2_theta_BG.Integral(),h_jet2_theta_HHZ_polm80_ECuts_.Integral(),norm_tot_BG_to_SIG
    #0-1650 with 55 bins, cut to 150 GeV
    #h_jet2_theta_hzqq_polm80_ECuts_.Rebin()
    #h_jet2_theta_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet2_theta_ee_qq_polm80_ECuts_.Scale(nor_tot_BG_to_SIG)
    #h_jet2_theta_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet2_theta_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    
    #h_jet2_theta_ee_qq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet2_theta_ee_qqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet2_theta_ee_qqqqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet2_theta_HZ_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    
    hhqq_BG_jet2_theta_polm80_ECuts_= THStack("hhqq_BG_jet2_theta_polm80_ECuts_", "");
    hhqq_BG_jet2_theta_polm80_ECuts_.Add(h_jet2_theta_ZZH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet2_theta_polm80_ECuts_.Add(h_jet2_theta_WWH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet2_theta_polm80_ECuts_.Add(h_jet2_theta_hzqq_polm80_ECuts_);
    hhqq_BG_jet2_theta_polm80_ECuts_.Add(h_jet2_theta_ee_qq_polm80_ECuts_);
    hhqq_BG_jet2_theta_polm80_ECuts_.Add(h_jet2_theta_ee_qqqqqq_polm80_ECuts_);
    hhqq_BG_jet2_theta_polm80_ECuts_.Add(h_jet2_theta_ee_qqqq_polm80_ECuts_);
    
      
    canvas_h_SIG_BG_jet2_theta_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_theta_polm80_ECuts_thstack");
    canvas_h_SIG_BG_jet2_theta_polm80_ECuts_thstack.cd();
    #h_jet2_theta_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet2_theta_BG.Draw("hist,e")
    #h_tot_norm_jet2_theta_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet2_theta_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet2_theta_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet2_theta_polm80_ECuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    hhqq_BG_jet2_theta_polm80_ECuts_.GetYaxis().SetTitle("Events");
    hhqq_BG_jet2_theta_polm80_ECuts_.GetYaxis().SetTitleOffset(1.5)
    hhqq_BG_jet2_theta_polm80_ECuts_.SetMaximum(750000.)
    h_jet2_theta_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_theta_polm80_ECuts_thstack.Modified();
      
    h_jet2_theta_HHZ_polm80_ECuts_.Scale(50000.)
    
    line = TLine(10,0,10,725000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(13,650000,25,650000,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    line2 = TLine(170,0,170,725000)
    line2.SetLineColor(kBlack);
    line2.SetLineWidth(2);
    line2.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow2 = TArrow(155,650000,167,650000,0.025,"<|")
    #arrow.SetAngle(0);
    arrow2.SetLineWidth(2);
    arrow2.Draw();
    
    leg_hzqq_BG_jet2_theta_polm80_ECuts_=TLegend(0.30,0.60,0.690,0.88);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq x 50000");
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_WWH_qqqqH_polm80_ECuts_,"WWH#rightarrowqqqqH");
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_ZZH_qqqqH_polm80_ECuts_,"ZZH#rightarrowqqqqH");
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_theta_polm80_ECuts_.Draw();
      
    l.DrawLatex(x0,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet2_theta_polm80_ECuts_thstack.Print("h_jet2_theta_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000.eps")
    canvas_h_SIG_BG_jet2_theta_polm80_ECuts_thstack.cd()
    hhqq_BG_jet2_theta_polm80_ECuts_.SetMaximum(10.e9)
    hhqq_BG_jet2_theta_polm80_ECuts_.SetMinimum(0.1)
    canvas_h_SIG_BG_jet2_theta_polm80_ECuts_thstack.SetLogy()
    canvas_h_SIG_BG_jet2_theta_polm80_ECuts_thstack.Print("h_jet2_theta_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000_logy.eps")
    h_jet2_theta_HHZ_polm80_ECuts_.Scale(1./50000.)
      
    h_jet2_theta_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_theta_ee_qq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_theta_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_theta_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
      
    hhqq_norm_BG_jet2_theta_polm80_ECuts_= THStack("hhqq_norm_BG_jet2_theta_polm80_ECuts_", "");
    hhqq_norm_BG_jet2_theta_polm80_ECuts_.Add(h_jet2_theta_hzqq_polm80_ECuts_);
    hhqq_norm_BG_jet2_theta_polm80_ECuts_.Add(h_jet2_theta_ee_qq_polm80_ECuts_);
    hhqq_norm_BG_jet2_theta_polm80_ECuts_.Add(h_jet2_theta_ee_qqqqqq_polm80_ECuts_);
    hhqq_norm_BG_jet2_theta_polm80_ECuts_.Add(h_jet2_theta_ee_qqqq_polm80_ECuts_);
      
    canvas_h_SIG_norm_BG_jet2_theta_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_norm_BG_jet2_theta_polm80_ECuts_thstack");
    canvas_h_SIG_norm_BG_jet2_theta_polm80_ECuts_thstack.cd();
    #h_jet2_theta_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet2_theta_BG.Draw("hist,e")
    #h_tot_norm_jet2_theta_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet2_theta_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet2_theta_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet2_theta_polm80_ECuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    hhqq_norm_BG_jet2_theta_polm80_ECuts_.GetYaxis().SetTitle("A.U.");
    hhqq_norm_BG_jet2_theta_polm80_ECuts_.SetMaximum(5.5)
    h_jet2_theta_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_norm_BG_jet2_theta_polm80_ECuts_thstack.Modified();
    
    line = TLine(10,0,10,5.40)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(13,4.5,25,4.5,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    line2 = TLine(170,0,170,5.40)
    line2.SetLineColor(kBlack);
    line2.SetLineWidth(2);
    line2.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow2 = TArrow(155,4.5,167,4.5,0.025,"<|")
    #arrow.SetAngle(0);
    arrow2.SetLineWidth(2);
    arrow2.Draw();
      
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_=TLegend(0.30,0.63,0.690,0.88);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq");
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.AddEntry(h_jet2_theta_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_norm_BG_jet2_theta_polm80_ECuts_.Draw();
      
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_norm_BG_jet2_theta_polm80_ECuts_thstack.Print("h_jet2_theta_polm80_normed_hzqq_ee_qqqqqq_qqqq_qq_and_hhqq.eps")
  
  
    #jet1 E
    h_jet1_E_HHZ_polm80_AllEvents_ECuts_=file_polm80_preselection_.Get("hhqq_AllEvents/h_jet1_E");
    h_jet1_E_HHZ_polm80_AllEvents_ECuts_.GetXaxis().SetTitle("jet1 E [GeV]");
    h_jet1_E_HHZ_polm80_AllEvents_ECuts_.SetLineWidth(3);
    h_jet1_E_HHZ_polm80_ECuts_=file_polm80_preselection_.Get("hhqq/h_jet1_E");
    h_jet1_E_HHZ_polm80_ECuts_.GetXaxis().SetTitle("jet1 E [GeV]");
    h_jet1_E_HHZ_polm80_ECuts_.SetLineWidth(3);
    h_jet1_E_HHZ_polm80_ECuts_.SetFillColor(kWhite)
    h_jet1_E_HHZ_polm80_ECuts_.SetLineColor(kBlack)
    

    h_jet1_E_hzqq_polm80_ECuts_=file_polm80_preselection_.Get("hzqq/h_jet1_E");
    h_jet1_E_hzqq_polm80_ECuts_.GetXaxis().SetTitle("jet1 E [GeV]");
    h_jet1_E_hzqq_polm80_ECuts_.SetLineWidth(3);
    h_jet1_E_hzqq_polm80_ECuts_.SetFillColor(kCyan+1)
    h_jet1_E_hzqq_polm80_ECuts_.SetLineColor(kCyan+1)
    h_jet1_E_hzqq_polm80_ECuts_.SetFillStyle(3002)
    
    h_jet1_E_ee_qq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qq/h_jet1_E")
    h_jet1_E_ee_qq_polm80_ECuts_.SetFillColor(kBlue);
    h_jet1_E_ee_qq_polm80_ECuts_.SetLineColor(kBlue);
    h_jet1_E_ee_qq_polm80_ECuts_.SetFillStyle(3002);
    h_jet1_E_ee_qqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqq/h_jet1_E");
    h_jet1_E_ee_qqqq_polm80_ECuts_.SetFillColor(kRed);
    h_jet1_E_ee_qqqq_polm80_ECuts_.SetLineColor(kRed);
    h_jet1_E_ee_qqqq_polm80_ECuts_.SetFillStyle(3002);
    h_jet1_E_ee_qqqqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqqqq/h_jet1_E");
    h_jet1_E_ee_qqqqqq_polm80_ECuts_.SetFillColor(kGreen-2);
    h_jet1_E_ee_qqqqqq_polm80_ECuts_.SetLineColor(kGreen-2);
    h_jet1_E_ee_qqqqqq_polm80_ECuts_.SetFillStyle(3002);
    
    h_jet1_E_WWH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("WWH_qqqqH/h_jet1_E");
    h_jet1_E_WWH_qqqqH_polm80_ECuts_.SetFillColor(kOrange);
    h_jet1_E_WWH_qqqqH_polm80_ECuts_.SetLineColor(kOrange);
    h_jet1_E_WWH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_jet1_E_ZZH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("ZZH_qqqqH/h_jet1_E");
    h_jet1_E_ZZH_qqqqH_polm80_ECuts_.SetFillColor(kViolet+2);
    h_jet1_E_ZZH_qqqqH_polm80_ECuts_.SetLineColor(kViolet+2);
    h_jet1_E_ZZH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_tot_norm_jet1_E_BG = h_jet1_E_ee_qq_polm80_ECuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_jet1_E_BG.Add(h_jet1_E_ee_qqqq_polm80_ECuts_);
    h_tot_norm_jet1_E_BG.Add(h_jet1_E_ee_qqqqqq_polm80_ECuts_);
    h_tot_norm_jet1_E_BG.Add(h_jet1_E_hzqq_polm80_ECuts_);
    norm_tot_BG_to_SIG = h_jet1_E_HHZ_polm80_ECuts_.Integral(0,h_jet1_E_HHZ_polm80_ECuts_.GetNbinsX()+1)/(h_jet1_E_hzqq_polm80_ECuts_.Integral(0,h_jet1_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet1_E_ee_qq_polm80_ECuts_.Integral(0,h_jet1_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet1_E_ee_qqqq_polm80_ECuts_.Integral(0,h_jet1_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet1_E_ee_qqqqqq_polm80_ECuts_.Integral(0,h_jet1_E_HHZ_polm80_ECuts_.GetNbinsX()+1))
    h_tot_norm_jet1_E_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_jet1_E_BG.SetLineColor(kBlack)
    h_tot_norm_jet1_E_BG.SetFillColor(0)
    
    print 'scale or range, norm here ',h_tot_norm_jet1_E_BG.Integral(),h_jet1_E_HHZ_polm80_ECuts_.Integral(),norm_tot_BG_to_SIG
    #0-1650 with 55 bins, cut to 150 GeV
    #h_jet1_E_hzqq_polm80_ECuts_.Rebin()
    #h_jet1_E_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet1_E_ee_qq_polm80_ECuts_.Scale(nor_tot_BG_to_SIG)
    #h_jet1_E_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet1_E_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    
    #h_jet1_E_ee_qq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet1_E_ee_qqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet1_E_ee_qqqqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet1_E_HZ_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    
    hhqq_BG_jet1_E_polm80_ECuts_= THStack("hhqq_BG_jet1_E_polm80_ECuts_", "");
    hhqq_BG_jet1_E_polm80_ECuts_.Add(h_jet1_E_ZZH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet1_E_polm80_ECuts_.Add(h_jet1_E_WWH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet1_E_polm80_ECuts_.Add(h_jet1_E_hzqq_polm80_ECuts_);
    hhqq_BG_jet1_E_polm80_ECuts_.Add(h_jet1_E_ee_qq_polm80_ECuts_);
    hhqq_BG_jet1_E_polm80_ECuts_.Add(h_jet1_E_ee_qqqqqq_polm80_ECuts_);
    hhqq_BG_jet1_E_polm80_ECuts_.Add(h_jet1_E_ee_qqqq_polm80_ECuts_);
    
    
    canvas_h_SIG_BG_jet1_E_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_E_polm80_ECuts_thstack");
    canvas_h_SIG_BG_jet1_E_polm80_ECuts_thstack.cd();
    #h_jet1_E_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet1_E_BG.Draw("hist,e")
    #h_tot_norm_jet1_E_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet1_E_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet1_E_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet1_E_polm80_ECuts_.GetXaxis().SetTitle("jet1 E [GeV]");
    hhqq_BG_jet1_E_polm80_ECuts_.GetYaxis().SetTitle("Events");
    hhqq_BG_jet1_E_polm80_ECuts_.GetYaxis().SetTitleOffset(1.5);
    hhqq_BG_jet1_E_polm80_ECuts_.SetMaximum(1400000.)
    h_jet1_E_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_E_polm80_ECuts_thstack.Modified();
    
    h_jet1_E_HHZ_polm80_ECuts_.Scale(50000.)
    
    line = TLine(150,0,150,1400000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(165,750000,320,750000,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_BG_jet1_E_polm80_ECuts_=TLegend(0.40,0.60,0.790,0.88);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq x 50000");
    leg_hzqq_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_WWH_qqqqH_polm80_ECuts_,"WWH#rightarrowqqqqH");
    leg_hzqq_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_ZZH_qqqqH_polm80_ECuts_,"ZZH#rightarrowqqqqH");
    leg_hzqq_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_E_polm80_ECuts_.Draw();
    
    l.DrawLatex(x0,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet1_E_polm80_ECuts_thstack.Print("h_jet1_E_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000.eps")
    canvas_h_SIG_BG_jet1_E_polm80_ECuts_thstack.cd()
    hhqq_BG_jet1_E_polm80_ECuts_.SetMaximum(10.e9)
    hhqq_BG_jet1_E_polm80_ECuts_.SetMinimum(0.1)
    canvas_h_SIG_BG_jet1_E_polm80_ECuts_thstack.SetLogy()
    canvas_h_SIG_BG_jet1_E_polm80_ECuts_thstack.Print("h_jet1_E_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000_logy.eps")
    
    h_jet1_E_HHZ_polm80_ECuts_.Scale(1./50000.)
    
    h_jet1_E_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_E_ee_qq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_E_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_E_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    
    hhqq_norm_BG_jet1_E_polm80_ECuts_= THStack("hhqq_norm_BG_jet1_E_polm80_ECuts_", "");
    hhqq_norm_BG_jet1_E_polm80_ECuts_.Add(h_jet1_E_hzqq_polm80_ECuts_);
    hhqq_norm_BG_jet1_E_polm80_ECuts_.Add(h_jet1_E_ee_qq_polm80_ECuts_);
    hhqq_norm_BG_jet1_E_polm80_ECuts_.Add(h_jet1_E_ee_qqqqqq_polm80_ECuts_);
    hhqq_norm_BG_jet1_E_polm80_ECuts_.Add(h_jet1_E_ee_qqqq_polm80_ECuts_);
    
    canvas_h_SIG_norm_BG_jet1_E_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_norm_BG_jet1_E_polm80_ECuts_thstack");
    canvas_h_SIG_norm_BG_jet1_E_polm80_ECuts_thstack.cd();
    #h_jet1_E_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet1_E_BG.Draw("hist,e")
    #h_tot_norm_jet1_E_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet1_E_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet1_E_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet1_E_polm80_ECuts_.GetXaxis().SetTitle("jet1 E [GeV]");
    hhqq_norm_BG_jet1_E_polm80_ECuts_.GetYaxis().SetTitle("A.U.");
    hhqq_norm_BG_jet1_E_polm80_ECuts_.SetMaximum(10.5)
    h_jet1_E_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_norm_BG_jet1_E_polm80_ECuts_thstack.Modified();
    
    line = TLine(150,0,150,10.5)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(165,5.5,320,5.5,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_=TLegend(0.40,0.63,0.790,0.88);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq");
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.AddEntry(h_jet1_E_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_norm_BG_jet1_E_polm80_ECuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_norm_BG_jet1_E_polm80_ECuts_thstack.Print("h_jet1_E_polm80_normed_hzqq_ee_qqqqqq_qqqq_qq_and_hhqq.eps")
    
    #comb jet2
    h_jet2_E_HHZ_polm80_AllEvents_ECuts_=file_polm80_preselection_.Get("hhqq_AllEvents/h_jet2_E");
    h_jet2_E_HHZ_polm80_AllEvents_ECuts_.GetXaxis().SetTitle("jet2 E [GeV]");
    h_jet2_E_HHZ_polm80_AllEvents_ECuts_.SetLineWidth(3);
    h_jet2_E_HHZ_polm80_ECuts_=file_polm80_preselection_.Get("hhqq/h_jet2_E");
    h_jet2_E_HHZ_polm80_ECuts_.GetXaxis().SetTitle("jet2 E [GeV]");
    h_jet2_E_HHZ_polm80_ECuts_.SetLineWidth(3);
    h_jet2_E_HHZ_polm80_ECuts_.SetFillColor(kWhite)
    h_jet2_E_HHZ_polm80_ECuts_.SetLineColor(kBlack)
    
    
    h_jet2_E_hzqq_polm80_ECuts_=file_polm80_preselection_.Get("hzqq/h_jet2_E");
    h_jet2_E_hzqq_polm80_ECuts_.GetXaxis().SetTitle("jet2 E [GeV]");
    h_jet2_E_hzqq_polm80_ECuts_.SetLineWidth(3);
    h_jet2_E_hzqq_polm80_ECuts_.SetFillColor(kCyan+1)
    h_jet2_E_hzqq_polm80_ECuts_.SetLineColor(kCyan+1)
    h_jet2_E_hzqq_polm80_ECuts_.SetFillStyle(3002)
    
    h_jet2_E_ee_qq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qq/h_jet2_E")
    h_jet2_E_ee_qq_polm80_ECuts_.SetFillColor(kBlue);
    h_jet2_E_ee_qq_polm80_ECuts_.SetLineColor(kBlue);
    h_jet2_E_ee_qq_polm80_ECuts_.SetFillStyle(3002);
    h_jet2_E_ee_qqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqq/h_jet2_E");
    h_jet2_E_ee_qqqq_polm80_ECuts_.SetFillColor(kRed);
    h_jet2_E_ee_qqqq_polm80_ECuts_.SetLineColor(kRed);
    h_jet2_E_ee_qqqq_polm80_ECuts_.SetFillStyle(3002);
    h_jet2_E_ee_qqqqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqqqq/h_jet2_E");
    h_jet2_E_ee_qqqqqq_polm80_ECuts_.SetFillColor(kGreen-2);
    h_jet2_E_ee_qqqqqq_polm80_ECuts_.SetLineColor(kGreen-2);
    h_jet2_E_ee_qqqqqq_polm80_ECuts_.SetFillStyle(3002);
    h_jet2_E_WWH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("WWH_qqqqH/h_jet2_E");
    h_jet2_E_WWH_qqqqH_polm80_ECuts_.SetFillColor(kOrange);
    h_jet2_E_WWH_qqqqH_polm80_ECuts_.SetLineColor(kOrange);
    h_jet2_E_WWH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_jet2_E_ZZH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("ZZH_qqqqH/h_jet2_E");
    h_jet2_E_ZZH_qqqqH_polm80_ECuts_.SetFillColor(kViolet+2);
    h_jet2_E_ZZH_qqqqH_polm80_ECuts_.SetLineColor(kViolet+2);
    h_jet2_E_ZZH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_tot_norm_jet2_E_BG = h_jet2_E_ee_qq_polm80_ECuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_jet2_E_BG.Add(h_jet2_E_ee_qqqq_polm80_ECuts_);
    h_tot_norm_jet2_E_BG.Add(h_jet2_E_ee_qqqqqq_polm80_ECuts_);
    h_tot_norm_jet2_E_BG.Add(h_jet2_E_hzqq_polm80_ECuts_);
    norm_tot_BG_to_SIG = h_jet2_E_HHZ_polm80_ECuts_.Integral(0,h_jet2_E_HHZ_polm80_ECuts_.GetNbinsX()+1)/(h_jet2_E_hzqq_polm80_ECuts_.Integral(0,h_jet2_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet2_E_ee_qq_polm80_ECuts_.Integral(0,h_jet2_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet2_E_ee_qqqq_polm80_ECuts_.Integral(0,h_jet2_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet2_E_ee_qqqqqq_polm80_ECuts_.Integral(0,h_jet2_E_HHZ_polm80_ECuts_.GetNbinsX()+1))
    h_tot_norm_jet2_E_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_jet2_E_BG.SetLineColor(kBlack)
    h_tot_norm_jet2_E_BG.SetFillColor(0)
    
    print 'scale or range, norm here ',h_tot_norm_jet2_E_BG.Integral(),h_jet2_E_HHZ_polm80_ECuts_.Integral(),norm_tot_BG_to_SIG
    #0-1650 with 55 bins, cut to 150 GeV
    #h_jet2_E_hzqq_polm80_ECuts_.Rebin()
    #h_jet2_E_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet2_E_ee_qq_polm80_ECuts_.Scale(nor_tot_BG_to_SIG)
    #h_jet2_E_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet2_E_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    
    #h_jet2_E_ee_qq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet2_E_ee_qqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet2_E_ee_qqqqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet2_E_HZ_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    
    hhqq_BG_jet2_E_polm80_ECuts_= THStack("hhqq_BG_jet2_E_polm80_ECuts_", "");
    hhqq_BG_jet2_E_polm80_ECuts_.Add(h_jet2_E_ZZH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet2_E_polm80_ECuts_.Add(h_jet2_E_WWH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet2_E_polm80_ECuts_.Add(h_jet2_E_hzqq_polm80_ECuts_);
    hhqq_BG_jet2_E_polm80_ECuts_.Add(h_jet2_E_ee_qq_polm80_ECuts_);
    hhqq_BG_jet2_E_polm80_ECuts_.Add(h_jet2_E_ee_qqqqqq_polm80_ECuts_);
    hhqq_BG_jet2_E_polm80_ECuts_.Add(h_jet2_E_ee_qqqq_polm80_ECuts_);
    
    
    canvas_h_SIG_BG_jet2_E_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_E_polm80_ECuts_thstack");
    canvas_h_SIG_BG_jet2_E_polm80_ECuts_thstack.cd();
    #h_jet2_E_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet2_E_BG.Draw("hist,e")
    #h_tot_norm_jet2_E_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet2_E_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet2_E_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet2_E_polm80_ECuts_.GetXaxis().SetTitle("jet2 E [GeV]");
    hhqq_BG_jet2_E_polm80_ECuts_.GetYaxis().SetTitle("Events");
    hhqq_BG_jet2_E_polm80_ECuts_.GetYaxis().SetTitleOffset(1.5);
    hhqq_BG_jet2_E_polm80_ECuts_.SetMaximum(2200000.)
    h_jet2_E_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_E_polm80_ECuts_thstack.Modified();
    
    h_jet2_E_HHZ_polm80_ECuts_.Scale(50000.)
    
    line = TLine(100,0,100,2200000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(115,1100000,270,1100000,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_BG_jet2_E_polm80_ECuts_=TLegend(0.40,0.60,0.790,0.88);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq x 50000");
    leg_hzqq_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_WWH_qqqqH_polm80_ECuts_,"WWH#rightarrowqqqqH");
    leg_hzqq_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_ZZH_qqqqH_polm80_ECuts_,"ZZH#rightarrowqqqqH");
    leg_hzqq_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_E_polm80_ECuts_.Draw();
    
    l.DrawLatex(x0,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet2_E_polm80_ECuts_thstack.Print("h_jet2_E_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000.eps")
    canvas_h_SIG_BG_jet2_E_polm80_ECuts_thstack.cd()
    hhqq_BG_jet2_E_polm80_ECuts_.SetMaximum(10.e9)
    hhqq_BG_jet2_E_polm80_ECuts_.SetMinimum(0.1)
    canvas_h_SIG_BG_jet2_E_polm80_ECuts_thstack.SetLogy()
    canvas_h_SIG_BG_jet2_E_polm80_ECuts_thstack.Print("h_jet2_E_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000_logy.eps")
    
    h_jet2_E_HHZ_polm80_ECuts_.Scale(1./50000.)
    
    h_jet2_E_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_E_ee_qq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_E_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_E_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    
    hhqq_norm_BG_jet2_E_polm80_ECuts_= THStack("hhqq_norm_BG_jet2_E_polm80_ECuts_", "");
    hhqq_norm_BG_jet2_E_polm80_ECuts_.Add(h_jet2_E_hzqq_polm80_ECuts_);
    hhqq_norm_BG_jet2_E_polm80_ECuts_.Add(h_jet2_E_ee_qq_polm80_ECuts_);
    hhqq_norm_BG_jet2_E_polm80_ECuts_.Add(h_jet2_E_ee_qqqqqq_polm80_ECuts_);
    hhqq_norm_BG_jet2_E_polm80_ECuts_.Add(h_jet2_E_ee_qqqq_polm80_ECuts_);
    
    canvas_h_SIG_norm_BG_jet2_E_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_norm_BG_jet2_E_polm80_ECuts_thstack");
    canvas_h_SIG_norm_BG_jet2_E_polm80_ECuts_thstack.cd();
    #h_jet2_E_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet2_E_BG.Draw("hist,e")
    #h_tot_norm_jet2_E_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet2_E_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet2_E_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet2_E_polm80_ECuts_.GetXaxis().SetTitle("jet2 E [GeV]");
    hhqq_norm_BG_jet2_E_polm80_ECuts_.GetYaxis().SetTitle("A.U.");
    hhqq_norm_BG_jet2_E_polm80_ECuts_.SetMaximum(18.5)
    h_jet2_E_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_norm_BG_jet2_E_polm80_ECuts_thstack.Modified();
    
    line = TLine(100,0,100,17.75)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(115,8,270,8,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_=TLegend(0.40,0.63,0.790,0.88);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq");
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.AddEntry(h_jet2_E_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_norm_BG_jet2_E_polm80_ECuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_norm_BG_jet2_E_polm80_ECuts_thstack.Print("h_jet2_E_polm80_normed_hzqq_ee_qqqqqq_qqqq_qq_and_hhqq.eps")
    
    #comb jet3
    h_jet3_E_HHZ_polm80_AllEvents_ECuts_=file_polm80_preselection_.Get("hhqq_AllEvents/h_jet3_E");
    h_jet3_E_HHZ_polm80_AllEvents_ECuts_.GetXaxis().SetTitle("jet3 E [GeV]");
    h_jet3_E_HHZ_polm80_AllEvents_ECuts_.SetLineWidth(3);
    h_jet3_E_HHZ_polm80_ECuts_=file_polm80_preselection_.Get("hhqq/h_jet3_E");
    h_jet3_E_HHZ_polm80_ECuts_.GetXaxis().SetTitle("jet3 E [GeV]");
    h_jet3_E_HHZ_polm80_ECuts_.SetLineWidth(3);
    h_jet3_E_HHZ_polm80_ECuts_.SetFillColor(kWhite)
    h_jet3_E_HHZ_polm80_ECuts_.SetLineColor(kBlack)
    
    
    h_jet3_E_hzqq_polm80_ECuts_=file_polm80_preselection_.Get("hzqq/h_jet3_E");
    h_jet3_E_hzqq_polm80_ECuts_.GetXaxis().SetTitle("jet3 E [GeV]");
    h_jet3_E_hzqq_polm80_ECuts_.SetLineWidth(3);
    h_jet3_E_hzqq_polm80_ECuts_.SetFillColor(kCyan+1)
    h_jet3_E_hzqq_polm80_ECuts_.SetLineColor(kCyan+1)
    h_jet3_E_hzqq_polm80_ECuts_.SetFillStyle(3002)
    
    h_jet3_E_ee_qq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qq/h_jet3_E")
    h_jet3_E_ee_qq_polm80_ECuts_.SetFillColor(kBlue);
    h_jet3_E_ee_qq_polm80_ECuts_.SetLineColor(kBlue);
    h_jet3_E_ee_qq_polm80_ECuts_.SetFillStyle(3002);
    h_jet3_E_ee_qqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqq/h_jet3_E");
    h_jet3_E_ee_qqqq_polm80_ECuts_.SetFillColor(kRed);
    h_jet3_E_ee_qqqq_polm80_ECuts_.SetLineColor(kRed);
    h_jet3_E_ee_qqqq_polm80_ECuts_.SetFillStyle(3002);
    h_jet3_E_ee_qqqqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqqqq/h_jet3_E");
    h_jet3_E_ee_qqqqqq_polm80_ECuts_.SetFillColor(kGreen-2);
    h_jet3_E_ee_qqqqqq_polm80_ECuts_.SetLineColor(kGreen-2);
    h_jet3_E_ee_qqqqqq_polm80_ECuts_.SetFillStyle(3002);
    h_jet3_E_WWH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("WWH_qqqqH/h_jet3_E");
    h_jet3_E_WWH_qqqqH_polm80_ECuts_.SetFillColor(kOrange);
    h_jet3_E_WWH_qqqqH_polm80_ECuts_.SetLineColor(kOrange);
    h_jet3_E_WWH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_jet3_E_ZZH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("ZZH_qqqqH/h_jet3_E");
    h_jet3_E_ZZH_qqqqH_polm80_ECuts_.SetFillColor(kViolet+2);
    h_jet3_E_ZZH_qqqqH_polm80_ECuts_.SetLineColor(kViolet+2);
    h_jet3_E_ZZH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_tot_norm_jet3_E_BG = h_jet3_E_ee_qq_polm80_ECuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_jet3_E_BG.Add(h_jet3_E_ee_qqqq_polm80_ECuts_);
    h_tot_norm_jet3_E_BG.Add(h_jet3_E_ee_qqqqqq_polm80_ECuts_);
    h_tot_norm_jet3_E_BG.Add(h_jet3_E_hzqq_polm80_ECuts_);
    norm_tot_BG_to_SIG = h_jet3_E_HHZ_polm80_ECuts_.Integral(0,h_jet3_E_HHZ_polm80_ECuts_.GetNbinsX()+1)/(h_jet3_E_hzqq_polm80_ECuts_.Integral(0,h_jet3_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet3_E_ee_qq_polm80_ECuts_.Integral(0,h_jet3_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet3_E_ee_qqqq_polm80_ECuts_.Integral(0,h_jet3_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet3_E_ee_qqqqqq_polm80_ECuts_.Integral(0,h_jet3_E_HHZ_polm80_ECuts_.GetNbinsX()+1))
    h_tot_norm_jet3_E_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_jet3_E_BG.SetLineColor(kBlack)
    h_tot_norm_jet3_E_BG.SetFillColor(0)
    
    print 'scale or range, norm here ',h_tot_norm_jet3_E_BG.Integral(),h_jet3_E_HHZ_polm80_ECuts_.Integral(),norm_tot_BG_to_SIG
    #0-1650 with 55 bins, cut to 150 GeV
    #h_jet3_E_hzqq_polm80_ECuts_.Rebin()
    #h_jet3_E_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet3_E_ee_qq_polm80_ECuts_.Scale(nor_tot_BG_to_SIG)
    #h_jet3_E_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet3_E_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    
    #h_jet3_E_ee_qq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet3_E_ee_qqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet3_E_ee_qqqqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet3_E_HZ_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    
    hhqq_BG_jet3_E_polm80_ECuts_= THStack("hhqq_BG_jet3_E_polm80_ECuts_", "");
    hhqq_BG_jet3_E_polm80_ECuts_.Add(h_jet3_E_ZZH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet3_E_polm80_ECuts_.Add(h_jet3_E_WWH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet3_E_polm80_ECuts_.Add(h_jet3_E_hzqq_polm80_ECuts_);
    hhqq_BG_jet3_E_polm80_ECuts_.Add(h_jet3_E_ee_qq_polm80_ECuts_);
    hhqq_BG_jet3_E_polm80_ECuts_.Add(h_jet3_E_ee_qqqqqq_polm80_ECuts_);
    hhqq_BG_jet3_E_polm80_ECuts_.Add(h_jet3_E_ee_qqqq_polm80_ECuts_);
    
    
    canvas_h_SIG_BG_jet3_E_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet3_E_polm80_ECuts_thstack");
    canvas_h_SIG_BG_jet3_E_polm80_ECuts_thstack.cd();
    #h_jet3_E_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet3_E_BG.Draw("hist,e")
    #h_tot_norm_jet3_E_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet3_E_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet3_E_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet3_E_polm80_ECuts_.GetXaxis().SetTitle("jet3 E [GeV]");
    hhqq_BG_jet3_E_polm80_ECuts_.GetYaxis().SetTitle("Events");
    hhqq_BG_jet3_E_polm80_ECuts_.GetYaxis().SetTitleOffset(1.5);
    hhqq_BG_jet3_E_polm80_ECuts_.SetMaximum(2200000.)
    h_jet3_E_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet3_E_polm80_ECuts_thstack.Modified();
    
    h_jet3_E_HHZ_polm80_ECuts_.Scale(50000.)
    
    line = TLine(50,0,50,2200000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(65,1100000,140,1100000,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_BG_jet3_E_polm80_ECuts_=TLegend(0.40,0.60,0.790,0.88);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq x 50000");
    leg_hzqq_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_WWH_qqqqH_polm80_ECuts_,"WWH#rightarrowqqqqH");
    leg_hzqq_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_ZZH_qqqqH_polm80_ECuts_,"ZZH#rightarrowqqqqH");
    leg_hzqq_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet3_E_polm80_ECuts_.Draw();
    
    l.DrawLatex(x0,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet3_E_polm80_ECuts_thstack.Print("h_jet3_E_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000.eps")
    canvas_h_SIG_BG_jet3_E_polm80_ECuts_thstack.cd()
    hhqq_BG_jet3_E_polm80_ECuts_.SetMaximum(10.e9)
    hhqq_BG_jet3_E_polm80_ECuts_.SetMinimum(0.1)
    canvas_h_SIG_BG_jet3_E_polm80_ECuts_thstack.SetLogy()
    canvas_h_SIG_BG_jet3_E_polm80_ECuts_thstack.Print("h_jet3_E_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000_logy.eps")
    
    h_jet3_E_HHZ_polm80_ECuts_.Scale(1./50000.)
    
    h_jet3_E_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet3_E_ee_qq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet3_E_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet3_E_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    
    hhqq_norm_BG_jet3_E_polm80_ECuts_= THStack("hhqq_norm_BG_jet3_E_polm80_ECuts_", "");
    hhqq_norm_BG_jet3_E_polm80_ECuts_.Add(h_jet3_E_hzqq_polm80_ECuts_);
    hhqq_norm_BG_jet3_E_polm80_ECuts_.Add(h_jet3_E_ee_qq_polm80_ECuts_);
    hhqq_norm_BG_jet3_E_polm80_ECuts_.Add(h_jet3_E_ee_qqqqqq_polm80_ECuts_);
    hhqq_norm_BG_jet3_E_polm80_ECuts_.Add(h_jet3_E_ee_qqqq_polm80_ECuts_);
    
    canvas_h_SIG_norm_BG_jet3_E_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_norm_BG_jet3_E_polm80_ECuts_thstack");
    canvas_h_SIG_norm_BG_jet3_E_polm80_ECuts_thstack.cd();
    #h_jet3_E_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet3_E_BG.Draw("hist,e")
    #h_tot_norm_jet3_E_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet3_E_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet3_E_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet3_E_polm80_ECuts_.GetXaxis().SetTitle("jet3 E [GeV]");
    hhqq_norm_BG_jet3_E_polm80_ECuts_.GetYaxis().SetTitle("A.U.");
    hhqq_norm_BG_jet3_E_polm80_ECuts_.SetMaximum(15.5)
    h_jet3_E_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_norm_BG_jet3_E_polm80_ECuts_thstack.Modified();
    
    line = TLine(50,0,50,15.25)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(65,6,140,6,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_=TLegend(0.40,0.60,0.790,0.88);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq");
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.AddEntry(h_jet3_E_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_norm_BG_jet3_E_polm80_ECuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_norm_BG_jet3_E_polm80_ECuts_thstack.Print("h_jet3_E_polm80_normed_hzqq_ee_qqqqqq_qqqq_qq_and_hhqq.eps")
    
    #comb jet4
    h_jet4_E_HHZ_polm80_AllEvents_ECuts_=file_polm80_preselection_.Get("hhqq_AllEvents/h_jet4_E");
    h_jet4_E_HHZ_polm80_AllEvents_ECuts_.GetXaxis().SetTitle("jet4 E [GeV]");
    h_jet4_E_HHZ_polm80_AllEvents_ECuts_.SetLineWidth(3);
    h_jet4_E_HHZ_polm80_ECuts_=file_polm80_preselection_.Get("hhqq/h_jet4_E");
    h_jet4_E_HHZ_polm80_ECuts_.GetXaxis().SetTitle("jet4 E [GeV]");
    h_jet4_E_HHZ_polm80_ECuts_.SetLineWidth(3);
    h_jet4_E_HHZ_polm80_ECuts_.SetFillColor(kWhite)
    h_jet4_E_HHZ_polm80_ECuts_.SetLineColor(kBlack)
    
    
    h_jet4_E_hzqq_polm80_ECuts_=file_polm80_preselection_.Get("hzqq/h_jet4_E");
    h_jet4_E_hzqq_polm80_ECuts_.GetXaxis().SetTitle("jet4 E [GeV]");
    h_jet4_E_hzqq_polm80_ECuts_.SetLineWidth(3);
    h_jet4_E_hzqq_polm80_ECuts_.SetFillColor(kCyan+1)
    h_jet4_E_hzqq_polm80_ECuts_.SetLineColor(kCyan+1)
    h_jet4_E_hzqq_polm80_ECuts_.SetFillStyle(3002)
    
    h_jet4_E_ee_qq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qq/h_jet4_E")
    h_jet4_E_ee_qq_polm80_ECuts_.SetFillColor(kBlue);
    h_jet4_E_ee_qq_polm80_ECuts_.SetLineColor(kBlue);
    h_jet4_E_ee_qq_polm80_ECuts_.SetFillStyle(3002);
    h_jet4_E_ee_qqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqq/h_jet4_E");
    h_jet4_E_ee_qqqq_polm80_ECuts_.SetFillColor(kRed);
    h_jet4_E_ee_qqqq_polm80_ECuts_.SetLineColor(kRed);
    h_jet4_E_ee_qqqq_polm80_ECuts_.SetFillStyle(3002);
    h_jet4_E_ee_qqqqqq_polm80_ECuts_=file_polm80_preselection_.Get("ee_qqqqqq/h_jet4_E");
    h_jet4_E_ee_qqqqqq_polm80_ECuts_.SetFillColor(kGreen-2);
    h_jet4_E_ee_qqqqqq_polm80_ECuts_.SetLineColor(kGreen-2);
    h_jet4_E_ee_qqqqqq_polm80_ECuts_.SetFillStyle(3002);
    h_jet4_E_WWH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("WWH_qqqqH/h_jet4_E");
    h_jet4_E_WWH_qqqqH_polm80_ECuts_.SetFillColor(kOrange);
    h_jet4_E_WWH_qqqqH_polm80_ECuts_.SetLineColor(kOrange);
    h_jet4_E_WWH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_jet4_E_ZZH_qqqqH_polm80_ECuts_=file_polm80_preselection_.Get("ZZH_qqqqH/h_jet4_E");
    h_jet4_E_ZZH_qqqqH_polm80_ECuts_.SetFillColor(kViolet+2);
    h_jet4_E_ZZH_qqqqH_polm80_ECuts_.SetLineColor(kViolet+2);
    h_jet4_E_ZZH_qqqqH_polm80_ECuts_.SetFillStyle(3002);
    
    h_tot_norm_jet4_E_BG = h_jet4_E_ee_qq_polm80_ECuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_jet4_E_BG.Add(h_jet4_E_ee_qqqq_polm80_ECuts_);
    h_tot_norm_jet4_E_BG.Add(h_jet4_E_ee_qqqqqq_polm80_ECuts_);
    h_tot_norm_jet4_E_BG.Add(h_jet4_E_hzqq_polm80_ECuts_);
    norm_tot_BG_to_SIG = h_jet4_E_HHZ_polm80_ECuts_.Integral(0,h_jet4_E_HHZ_polm80_ECuts_.GetNbinsX()+1)/(h_jet4_E_hzqq_polm80_ECuts_.Integral(0,h_jet4_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet4_E_ee_qq_polm80_ECuts_.Integral(0,h_jet4_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet4_E_ee_qqqq_polm80_ECuts_.Integral(0,h_jet4_E_HHZ_polm80_ECuts_.GetNbinsX()+1)+h_jet4_E_ee_qqqqqq_polm80_ECuts_.Integral(0,h_jet4_E_HHZ_polm80_ECuts_.GetNbinsX()+1))
    h_tot_norm_jet4_E_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_jet4_E_BG.SetLineColor(kBlack)
    h_tot_norm_jet4_E_BG.SetFillColor(0)
    
    print 'scale or range, norm here ',h_tot_norm_jet4_E_BG.Integral(),h_jet4_E_HHZ_polm80_ECuts_.Integral(),norm_tot_BG_to_SIG
    #0-1650 with 55 bins, cut to 150 GeV
    #h_jet4_E_hzqq_polm80_ECuts_.Rebin()
    #h_jet4_E_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet4_E_ee_qq_polm80_ECuts_.Scale(nor_tot_BG_to_SIG)
    #h_jet4_E_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    #h_jet4_E_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    
    #h_jet4_E_ee_qq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet4_E_ee_qqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet4_E_ee_qqqqqq_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    #h_jet4_E_HZ_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    
    hhqq_BG_jet4_E_polm80_ECuts_= THStack("hhqq_BG_jet4_E_polm80_ECuts_", "");
    hhqq_BG_jet4_E_polm80_ECuts_.Add(h_jet4_E_ZZH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet4_E_polm80_ECuts_.Add(h_jet4_E_WWH_qqqqH_polm80_ECuts_);
    hhqq_BG_jet4_E_polm80_ECuts_.Add(h_jet4_E_hzqq_polm80_ECuts_);
    hhqq_BG_jet4_E_polm80_ECuts_.Add(h_jet4_E_ee_qq_polm80_ECuts_);
    hhqq_BG_jet4_E_polm80_ECuts_.Add(h_jet4_E_ee_qqqqqq_polm80_ECuts_);
    hhqq_BG_jet4_E_polm80_ECuts_.Add(h_jet4_E_ee_qqqq_polm80_ECuts_);
    
    
    canvas_h_SIG_BG_jet4_E_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet4_E_polm80_ECuts_thstack");
    canvas_h_SIG_BG_jet4_E_polm80_ECuts_thstack.cd();
    #h_jet4_E_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet4_E_BG.Draw("hist,e")
    #h_tot_norm_jet4_E_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet4_E_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet4_E_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_jet4_E_polm80_ECuts_.GetXaxis().SetTitle("jet4 E [GeV]");
    hhqq_BG_jet4_E_polm80_ECuts_.GetYaxis().SetTitle("Events");
    hhqq_BG_jet4_E_polm80_ECuts_.GetYaxis().SetTitleOffset(1.5);
    hhqq_BG_jet4_E_polm80_ECuts_.SetMaximum(2500000.)
    h_jet4_E_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet4_E_polm80_ECuts_thstack.Modified();
    
    h_jet4_E_HHZ_polm80_ECuts_.Scale(50000.)
    
    line = TLine(50,0,50,2400000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(65,1100000,140,1100000,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_BG_jet4_E_polm80_ECuts_=TLegend(0.40,0.60,0.790,0.88);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq x 50000");
    leg_hzqq_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_WWH_qqqqH_polm80_ECuts_,"WWH#rightarrowqqqqH");
    leg_hzqq_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_ZZH_qqqqH_polm80_ECuts_,"ZZH#rightarrowqqqqH");
    leg_hzqq_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet4_E_polm80_ECuts_.Draw();
    
    l.DrawLatex(x0,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet4_E_polm80_ECuts_thstack.Print("h_jet4_E_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000.eps")
    canvas_h_SIG_BG_jet4_E_polm80_ECuts_thstack.cd()
    hhqq_BG_jet4_E_polm80_ECuts_.SetMaximum(10.e9)
    hhqq_BG_jet4_E_polm80_ECuts_.SetMinimum(0.1)
    canvas_h_SIG_BG_jet4_E_polm80_ECuts_thstack.SetLogy()
    canvas_h_SIG_BG_jet4_E_polm80_ECuts_thstack.Print("h_jet4_E_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000_logy.eps")
    
    h_jet4_E_HHZ_polm80_ECuts_.Scale(1./50000.)
    
    h_jet4_E_hzqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet4_E_ee_qq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet4_E_ee_qqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    h_jet4_E_ee_qqqqqq_polm80_ECuts_.Scale(norm_tot_BG_to_SIG)
    
    hhqq_norm_BG_jet4_E_polm80_ECuts_= THStack("hhqq_norm_BG_jet4_E_polm80_ECuts_", "");
    hhqq_norm_BG_jet4_E_polm80_ECuts_.Add(h_jet4_E_hzqq_polm80_ECuts_);
    hhqq_norm_BG_jet4_E_polm80_ECuts_.Add(h_jet4_E_ee_qq_polm80_ECuts_);
    hhqq_norm_BG_jet4_E_polm80_ECuts_.Add(h_jet4_E_ee_qqqqqq_polm80_ECuts_);
    hhqq_norm_BG_jet4_E_polm80_ECuts_.Add(h_jet4_E_ee_qqqq_polm80_ECuts_);
    
    canvas_h_SIG_norm_BG_jet4_E_polm80_ECuts_thstack = setUpperCanvas("canvas_h_SIG_norm_BG_jet4_E_polm80_ECuts_thstack");
    canvas_h_SIG_norm_BG_jet4_E_polm80_ECuts_thstack.cd();
    #h_jet4_E_HZ_polm80_ECuts_.Draw("hist,e")
    #h_tot_norm_jet4_E_BG.Draw("hist,e")
    #h_tot_norm_jet4_E_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet4_E_polm80_ECuts_.Draw("hist");
    #hhqq_BG_jet4_E_polm80_ECuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_jet4_E_polm80_ECuts_.GetXaxis().SetTitle("jet4 E [GeV]");
    hhqq_norm_BG_jet4_E_polm80_ECuts_.GetYaxis().SetTitle("A.U.");
    hhqq_norm_BG_jet4_E_polm80_ECuts_.SetMaximum(18.5)
    h_jet4_E_HHZ_polm80_ECuts_.Draw("hist,e,same")
    canvas_h_SIG_norm_BG_jet4_E_polm80_ECuts_thstack.Modified();
    
    line = TLine(50,0,50,17.75)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(65,7.5,140,7.5,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_=TLegend(0.40,0.63,0.790,0.88);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.SetBorderSize(0);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.SetTextAlign(12);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.SetTextSize(0.050);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.SetTextFont(42);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.SetMargin(0.15);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.SetLineColor(1);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.SetLineStyle(1);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.SetLineWidth(1);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.SetFillColor(0);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.SetFillStyle(0);
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_HHZ_polm80_ECuts_,"HHZ#rightarrowbbbbqq");
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_ee_qqqqqq_polm80_ECuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_ee_qqqq_polm80_ECuts_,"ee#rightarrow qqqq");
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_hzqq_polm80_ECuts_,"HZ#rightarrowHqq");
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.AddEntry(h_jet4_E_ee_qq_polm80_ECuts_,"ee#rightarrow qq");
    leg_hzqq_norm_BG_jet4_E_polm80_ECuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_norm_BG_jet4_E_polm80_ECuts_thstack.Print("h_jet4_E_polm80_normed_hzqq_ee_qqqqqq_qqqq_qq_and_hhqq.eps")
    
    #BTag_sum_max3
    h_BTag_sum_max3_HHZ_polm80_AllEvents_massCuts_=file_polm80_preselection_.Get("hhqq_AllEvents/h_BTag_sum_max3");
    h_BTag_sum_max3_HHZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("#sum BTag (max 3 jets)");
    h_BTag_sum_max3_HHZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_BTag_sum_max3_HHZ_polm80_massCuts_=file_polm80_preselection_.Get("hhqq/h_BTag_sum_max3");
    h_BTag_sum_max3_HHZ_polm80_massCuts_.GetXaxis().SetTitle("#sum BTag (max 3 jets)");
    h_BTag_sum_max3_HHZ_polm80_massCuts_.SetLineWidth(3);
    h_BTag_sum_max3_HHZ_polm80_massCuts_.SetFillColor(kWhite)
    h_BTag_sum_max3_HHZ_polm80_massCuts_.SetLineColor(kBlack)
    
    
    h_BTag_sum_max3_hzqq_polm80_massCuts_=file_polm80_preselection_.Get("hzqq/h_BTag_sum_max3");
    h_BTag_sum_max3_hzqq_polm80_massCuts_.GetXaxis().SetTitle("#sum BTag (max 3 jets)");
    h_BTag_sum_max3_hzqq_polm80_massCuts_.SetLineWidth(3);
    h_BTag_sum_max3_hzqq_polm80_massCuts_.SetFillColor(kCyan+1)
    h_BTag_sum_max3_hzqq_polm80_massCuts_.SetLineColor(kCyan+1)
    h_BTag_sum_max3_hzqq_polm80_massCuts_.SetFillStyle(3002)
    
    h_BTag_sum_max3_ee_qq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qq/h_BTag_sum_max3")
    h_BTag_sum_max3_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_BTag_sum_max3_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_BTag_sum_max3_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_BTag_sum_max3_ee_qqqq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qqqq/h_BTag_sum_max3");
    h_BTag_sum_max3_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_BTag_sum_max3_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_BTag_sum_max3_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_=file_polm80_preselection_.Get("ee_qqqqqq/h_BTag_sum_max3");
    h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);
    h_BTag_sum_max3_WWH_qqqqH_polm80_massCuts_=file_polm80_preselection_.Get("WWH_qqqqH/h_BTag_sum_max3");
    h_BTag_sum_max3_WWH_qqqqH_polm80_massCuts_.SetFillColor(kOrange);
    h_BTag_sum_max3_WWH_qqqqH_polm80_massCuts_.SetLineColor(kOrange);
    h_BTag_sum_max3_WWH_qqqqH_polm80_massCuts_.SetFillStyle(3002);
    
    h_BTag_sum_max3_ZZH_qqqqH_polm80_massCuts_=file_polm80_preselection_.Get("ZZH_qqqqH/h_BTag_sum_max3");
    h_BTag_sum_max3_ZZH_qqqqH_polm80_massCuts_.SetFillColor(kViolet+2);
    h_BTag_sum_max3_ZZH_qqqqH_polm80_massCuts_.SetLineColor(kViolet+2);
    h_BTag_sum_max3_ZZH_qqqqH_polm80_massCuts_.SetFillStyle(3002);
    
    h_tot_norm_h_BTag_sum_max3_BG = h_BTag_sum_max3_ee_qq_polm80_massCuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_h_BTag_sum_max3_BG.Add(h_BTag_sum_max3_ee_qqqq_polm80_massCuts_);
    h_tot_norm_h_BTag_sum_max3_BG.Add(h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_);
    h_tot_norm_h_BTag_sum_max3_BG.Add(h_BTag_sum_max3_hzqq_polm80_massCuts_);
    norm_tot_BG_to_SIG = h_BTag_sum_max3_HHZ_polm80_massCuts_.Integral(0,h_BTag_sum_max3_HHZ_polm80_massCuts_.GetNbinsX()+1)/(h_BTag_sum_max3_hzqq_polm80_massCuts_.Integral(0,h_BTag_sum_max3_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_BTag_sum_max3_ee_qq_polm80_massCuts_.Integral(0,h_BTag_sum_max3_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_BTag_sum_max3_ee_qqqq_polm80_massCuts_.Integral(0,h_BTag_sum_max3_HHZ_polm80_massCuts_.GetNbinsX()+1)+h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_.Integral(0,h_BTag_sum_max3_HHZ_polm80_massCuts_.GetNbinsX()+1))
    h_tot_norm_h_BTag_sum_max3_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_h_BTag_sum_max3_BG.SetLineColor(kBlack)
    h_tot_norm_h_BTag_sum_max3_BG.SetFillColor(0)
    
    print 'scale or range, norm here ',h_tot_norm_h_BTag_sum_max3_BG.Integral(),h_BTag_sum_max3_HHZ_polm80_massCuts_.Integral(),norm_tot_BG_to_SIG
    #0-1650 with 55 bins, cut to 150 GeV
    #h_BTag_sum_max3_hzqq_polm80_massCuts_.Rebin()
    #h_BTag_sum_max3_hzqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    #h_BTag_sum_max3_ee_qq_polm80_massCuts_.Scale(nor_tot_BG_to_SIG)
    #h_BTag_sum_max3_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    #h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    
    #h_BTag_sum_max3_ee_qq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_BTag_sum_max3_ee_qqqq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    #h_BTag_sum_max3_HZ_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    
    hhqq_BG_BTag_sum_max3_polm80_massCuts_= THStack("hhqq_BG_BTag_sum_max3_polm80_massCuts_", "");
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.Add(h_BTag_sum_max3_ZZH_qqqqH_polm80_massCuts_);
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.Add(h_BTag_sum_max3_WWH_qqqqH_polm80_massCuts_);
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.Add(h_BTag_sum_max3_hzqq_polm80_massCuts_);
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.Add(h_BTag_sum_max3_ee_qq_polm80_massCuts_);
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.Add(h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_);
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.Add(h_BTag_sum_max3_ee_qqqq_polm80_massCuts_);
    
    
    canvas_h_SIG_BG_BTag_sum_max3_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_BTag_sum_max3_polm80_massCuts_thstack");
    canvas_h_SIG_BG_BTag_sum_max3_polm80_massCuts_thstack.cd();
    #h_BTag_sum_max3_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_h_BTag_sum_max3_BG.Draw("hist,e")
    #h_tot_norm_h_BTag_sum_max3_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.Draw("hist");
    #hhqq_BG_BTag_sum_max3_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.GetXaxis().SetTitle("#sum BTag (max 3 jets)");
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.GetYaxis().SetTitle("Events");
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.SetMaximum(450000.)
    h_BTag_sum_max3_HHZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_BTag_sum_max3_polm80_massCuts_thstack.Modified();
    
    h_BTag_sum_max3_HHZ_polm80_massCuts_.Scale(50000.)
    
    line = TLine(2.2,0,2.2,280000)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(2.3,150000,2.9,150000,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_=TLegend(0.35,0.60,0.740,0.88);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_HHZ_polm80_massCuts_,"HHZ#rightarrowbbbbqq x 50000");
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_hzqq_polm80_massCuts_,"HZ#rightarrowHqq");
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_WWH_qqqqH_polm80_massCuts_,"WWH#rightarrowqqqqH");
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_ZZH_qqqqH_polm80_massCuts_,"ZZH#rightarrowqqqqH");
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_BTag_sum_max3_polm80_massCuts_.Draw();
    
    l.DrawLatex(x0,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_BTag_sum_max3_polm80_massCuts_thstack.Print("h_BTag_sum_max3_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000.eps")
    canvas_h_SIG_BG_BTag_sum_max3_polm80_massCuts_thstack.cd()
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.SetMaximum(10.e9)
    hhqq_BG_BTag_sum_max3_polm80_massCuts_.SetMinimum(0.1)
    canvas_h_SIG_BG_BTag_sum_max3_polm80_massCuts_thstack.SetLogy()
    canvas_h_SIG_BG_BTag_sum_max3_polm80_massCuts_thstack.Print("h_BTag_sum_max3_polm80_hzqq_ee_qqqqqq_qqqq_qq_WWH_ZZH_and_hhqq_50000_logy.eps")
    
    h_BTag_sum_max3_HHZ_polm80_massCuts_.Scale(1./50000.)
    
    h_BTag_sum_max3_hzqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_BTag_sum_max3_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_BTag_sum_max3_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    
    hhqq_norm_BG_BTag_sum_max3_polm80_massCuts_= THStack("hhqq_norm_BG_BTag_sum_max3_polm80_massCuts_", "");
    hhqq_norm_BG_BTag_sum_max3_polm80_massCuts_.Add(h_BTag_sum_max3_hzqq_polm80_massCuts_);
    hhqq_norm_BG_BTag_sum_max3_polm80_massCuts_.Add(h_BTag_sum_max3_ee_qq_polm80_massCuts_);
    hhqq_norm_BG_BTag_sum_max3_polm80_massCuts_.Add(h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_);
    hhqq_norm_BG_BTag_sum_max3_polm80_massCuts_.Add(h_BTag_sum_max3_ee_qqqq_polm80_massCuts_);
    
    canvas_h_SIG_norm_BG_BTag_sum_max3_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_norm_BG_BTag_sum_max3_polm80_massCuts_thstack");
    canvas_h_SIG_norm_BG_BTag_sum_max3_polm80_massCuts_thstack.cd();
    #h_BTag_sum_max3_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_h_BTag_sum_max3_BG.Draw("hist,e")
    #h_tot_norm_h_BTag_sum_max3_BG.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_BTag_sum_max3_polm80_massCuts_.Draw("hist");
    #hhqq_BG_BTag_sum_max3_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    hhqq_norm_BG_BTag_sum_max3_polm80_massCuts_.GetXaxis().SetTitle("#sum BTag (max 3 jets)");
    hhqq_norm_BG_BTag_sum_max3_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hhqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetMaximum(3.5)
    h_BTag_sum_max3_HHZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_norm_BG_BTag_sum_max3_polm80_massCuts_thstack.Modified();
    
    line = TLine(1,0,1,2.1)
    line.SetLineColor(kBlack);
    line.SetLineWidth(2);
    line.Draw()
    #x1,y1,x2,y2,arrow size, direction
    arrow = TArrow(1.1,1.5,1.6,1.5,0.025,"|>")
    #arrow.SetAngle(0);
    arrow.SetLineWidth(2);
    arrow.Draw();
    
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_=TLegend(0.50,0.63,0.890,0.88);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_HHZ_polm80_massCuts_,"HHZ#rightarrowbbbbqq");
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_hzqq_polm80_massCuts_,"HZ#rightarrowHqq");
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.AddEntry(h_BTag_sum_max3_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_norm_BG_BTag_sum_max3_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_norm_BG_BTag_sum_max3_polm80_massCuts_thstack.Print("h_BTag_sum_max3_polm80_normed_hzqq_ee_qqqqqq_qqqq_qq_and_hhqq.eps")
    '''
    
    file_BDT_Hhqq_sig_vs_BG_polm80 = root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/MVATrainingWeights_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_hzqq__ee_qq__ee_qqqq__ee_qqqqqq.root")
    
    h_BDT_Signal_polm80_test=file_BDT_Hhqq_sig_vs_BG_polm80.Get("dataset/Method_BDT/BDT/MVA_BDT_S");
    h_BDT_Signal_polm80_test.SetFillColor(kBlue-2);
    h_BDT_Signal_polm80_test.SetFillStyle(3001);
    h_BDT_Signal_polm80_test.GetXaxis().SetTitle("BDT Score");
    h_BDT_Signal_polm80_test.GetYaxis().SetTitle("a.u.");
    h_BDT_Signal_polm80_test.SetMinimum(0)
    h_BDT_Signal_polm80_test.SetMaximum(2.75)
    
    h_BDT_Signal_polm80_train=file_BDT_Hhqq_sig_vs_BG_polm80.Get("dataset/Method_BDT/BDT/MVA_BDT_Train_S");
    h_BDT_Signal_polm80_train.SetLineColor(kBlue);
    h_BDT_Signal_polm80_train.SetLineWidth(2);
    h_BDT_Signal_polm80_train.GetXaxis().SetTitle("BDT Score");
    h_BDT_Signal_polm80_train.GetYaxis().SetTitle("Entries");
    
    h_BDT_BG_polm80_test=file_BDT_Hhqq_sig_vs_BG_polm80.Get("dataset/Method_BDT/BDT/MVA_BDT_B");
    h_BDT_BG_polm80_test.SetFillColor(kRed-2);
    h_BDT_BG_polm80_test.SetFillStyle(3002);
    h_BDT_BG_polm80_test.GetXaxis().SetTitle("BDT Score");
    h_BDT_BG_polm80_test.GetYaxis().SetTitle("Entries");
    
    h_BDT_BG_polm80_train=file_BDT_Hhqq_sig_vs_BG_polm80.Get("dataset/Method_BDT/BDT/MVA_BDT_Train_B");
    h_BDT_BG_polm80_train.SetLineColor(kRed);
    h_BDT_BG_polm80_train.SetLineWidth(2);
    h_BDT_BG_polm80_train.GetXaxis().SetTitle("BDT Score");
    h_BDT_BG_polm80_train.GetYaxis().SetTitle("Entries");
    
    canvas_h_BDT_Signal_polm80_test = setUpperCanvas("canvas_h_BDT_Signal_polm80_test");
    canvas_h_BDT_Signal_polm80_test.cd();

    leg_h_BDT_Signal_polm80_test=TLegend(0.475,0.65,0.875,0.90);
    leg_h_BDT_Signal_polm80_test.SetBorderSize(0);
    leg_h_BDT_Signal_polm80_test.SetTextAlign(12);
    leg_h_BDT_Signal_polm80_test.SetTextSize(0.050);
    leg_h_BDT_Signal_polm80_test.SetTextFont(42);
    leg_h_BDT_Signal_polm80_test.SetMargin(0.15);
    leg_h_BDT_Signal_polm80_test.SetLineColor(1);
    leg_h_BDT_Signal_polm80_test.SetLineStyle(1);
    leg_h_BDT_Signal_polm80_test.SetLineWidth(1);
    leg_h_BDT_Signal_polm80_test.SetFillColor(0);
    leg_h_BDT_Signal_polm80_test.SetFillStyle(0);
    leg_h_BDT_Signal_polm80_test.AddEntry(h_BDT_Signal_polm80_test.DrawCopy("hist,e"),"Signal,test");
    leg_h_BDT_Signal_polm80_test.AddEntry(h_BDT_BG_polm80_test.DrawCopy("hist,e,same"),"BG,test");
    leg_h_BDT_Signal_polm80_test.AddEntry(h_BDT_Signal_polm80_train.DrawCopy("hist,e,same"),"Signal,train");
    leg_h_BDT_Signal_polm80_test.AddEntry(h_BDT_BG_polm80_train.DrawCopy("hist,e,same"),"BG,train");
    leg_h_BDT_Signal_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_BDT_Signal_polm80_test.Update()
    canvas_h_BDT_Signal_polm80_test.Modify()
    
    canvas_h_BDT_Signal_polm80_test.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/h_BDT_OverTrainingChecks_polm80.eps")

    file_BDT_Hhqq_sig_vs_BG_polp80 = root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/MVATrainingWeights_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_hzqq__ee_qq__ee_qqqq__ee_qqqqqq.root")
    
    h_BDT_Signal_polp80_test=file_BDT_Hhqq_sig_vs_BG_polp80.Get("dataset/Method_BDT/BDT/MVA_BDT_S");
    h_BDT_Signal_polp80_test.SetFillColor(kBlue-2);
    h_BDT_Signal_polp80_test.SetFillStyle(3001);
    h_BDT_Signal_polp80_test.GetXaxis().SetTitle("BDT Score");
    h_BDT_Signal_polp80_test.GetYaxis().SetTitle("a.u.");
    h_BDT_Signal_polp80_test.SetMinimum(0)
    h_BDT_Signal_polp80_test.SetMaximum(2.75)
    
    h_BDT_Signal_polp80_train=file_BDT_Hhqq_sig_vs_BG_polp80.Get("dataset/Method_BDT/BDT/MVA_BDT_Train_S");
    h_BDT_Signal_polp80_train.SetLineColor(kBlue);
    h_BDT_Signal_polp80_train.SetLineWidth(2);
    h_BDT_Signal_polp80_train.GetXaxis().SetTitle("BDT Score");
    h_BDT_Signal_polp80_train.GetYaxis().SetTitle("Entries");
    
    h_BDT_BG_polp80_test=file_BDT_Hhqq_sig_vs_BG_polp80.Get("dataset/Method_BDT/BDT/MVA_BDT_B");
    h_BDT_BG_polp80_test.SetFillColor(kRed-2);
    h_BDT_BG_polp80_test.SetFillStyle(3002);
    h_BDT_BG_polp80_test.GetXaxis().SetTitle("BDT Score");
    h_BDT_BG_polp80_test.GetYaxis().SetTitle("Entries");

    h_BDT_BG_polp80_train=file_BDT_Hhqq_sig_vs_BG_polp80.Get("dataset/Method_BDT/BDT/MVA_BDT_Train_B");
    h_BDT_BG_polp80_train.SetLineColor(kRed);
    h_BDT_BG_polp80_train.SetLineWidth(2);
    h_BDT_BG_polp80_train.GetXaxis().SetTitle("BDT Score");
    h_BDT_BG_polp80_train.GetYaxis().SetTitle("Entries");
    
    canvas_h_BDT_Signal_polp80_test = setUpperCanvas("canvas_h_BDT_Signal_polp80_test");
    canvas_h_BDT_Signal_polp80_test.cd();

    leg_h_BDT_Signal_polp80_test=TLegend(0.475,0.65,0.875,0.90);
    leg_h_BDT_Signal_polp80_test.SetBorderSize(0);
    leg_h_BDT_Signal_polp80_test.SetTextAlign(12);
    leg_h_BDT_Signal_polp80_test.SetTextSize(0.050);
    leg_h_BDT_Signal_polp80_test.SetTextFont(42);
    leg_h_BDT_Signal_polp80_test.SetMargin(0.15);
    leg_h_BDT_Signal_polp80_test.SetLineColor(1);
    leg_h_BDT_Signal_polp80_test.SetLineStyle(1);
    leg_h_BDT_Signal_polp80_test.SetLineWidth(1);
    leg_h_BDT_Signal_polp80_test.SetFillColor(0);
    leg_h_BDT_Signal_polp80_test.SetFillStyle(0);
    leg_h_BDT_Signal_polp80_test.AddEntry(h_BDT_Signal_polp80_test.DrawCopy("hist,e"),"Signal,test");
    leg_h_BDT_Signal_polp80_test.AddEntry(h_BDT_BG_polp80_test.DrawCopy("hist,e,same"),"BG,test");
    leg_h_BDT_Signal_polp80_test.AddEntry(h_BDT_Signal_polp80_train.DrawCopy("hist,e,same"),"Signal,train");
    leg_h_BDT_Signal_polp80_test.AddEntry(h_BDT_BG_polp80_train.DrawCopy("hist,e,same"),"BG,train");
    leg_h_BDT_Signal_polp80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_BDT_Signal_polp80_test.Update()
    canvas_h_BDT_Signal_polp80_test.Modify()

    canvas_h_BDT_Signal_polp80_test.Print("~/plotsHHZ_VLC_Jets_rfJets_BTag_NJet3_to_NJet6_noIsoP_191213/h_BDT_OverTrainingChecks_polp80.eps")


    #canvas_h_BDT_Signal_polm80_test.Print("h_BDT_OverTrainingChecks_polm80.pdf")
    
    
    file_polm80_HHZ_SignalHistos_VLC11_njet6_ =root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis/ntuple_HHZ_hhqq_14364_polm80_3TeV_wO_CLIC_o3_v14.root")


    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC11_njet6_.Get("h_mass_comb_rfj_rj1_HHZ_bbbbqq")
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('comb jet mass [GeV]')
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetMaximum(15)
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC11_njet6_.Get("h_mass_comb_rfj_rj2_HHZ_bbbbqq")
    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('comb jet mass [GeV]')
    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC11_njet6_.Get("h_mass_comb_rfj_rj3_HHZ_bbbbqq")
    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('comb jet mass [GeV]')
    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')


    canvas_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet6_VLC11_hhz_rfj_rj_mass_HZZ_bbbbqq_polm80");
    canvas_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHqq#rightarrowb#bar{b} b#bar{b} qq");
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"rj1");
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"rj2");
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"rj3");
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/canvas_njet6_VLC11_hhz_rfj_rj_mass_B_HZZ_bbbbqq_polm80.eps")



    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC11_njet6_.Get("h_mass_comb_rfj_rj1_HHZ_bbbbqq")
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetLineColor(1)
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('comb jet BTagMax')
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetMaximum(15)
    h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.SetMinimum(0)

    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC11_njet6_.Get("h_mass_comb_rfj_rj2_HHZ_bbbbqq")
    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.SetLineColor(2)
    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('comb jet BTagMax')
    h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')

    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80=file_polm80_HHZ_SignalHistos_VLC11_njet6_.Get("h_mass_comb_rfj_rj3_HHZ_bbbbqq")
    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.Rebin(4)
    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.SetLineWidth(2)
    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.SetLineColor(4)
    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.GetXaxis().SetTitle('comb jet BTagMax')
    h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.GetYaxis().SetTitle('Events')


    canvas_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80 = setUpperCanvas("canvas_njet6_VLC11_hhz_rfj_rj_mass_HZZ_bbbbqq_polm80");
    canvas_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80.cd()

    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test=TLegend(0.30,0.61,0.70,0.87);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetBorderSize(0);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextAlign(12);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextSize(0.050);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetTextFont(42);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetMargin(0.15);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineColor(1);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineStyle(1);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetLineWidth(1);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillColor(0);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetFillStyle(0);
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.SetHeader("HHqq#rightarrowb#bar{b} b#bar{b} qq");
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC11_mass_rfj_rj1_HZZ_bbbbqq_polm80.DrawCopy("h"),"rj1");
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC11_mass_rfj_rj2_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"rj2");
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.AddEntry(h_njet6_VLC11_mass_rfj_rj3_HZZ_bbbbqq_polm80.DrawCopy("h,same"),"rj3");
    leg_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_njet6_VLC11_rfj_rj_mass_B_HZZ_bbbbqq_polm80.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/canvas_njet6_VLC11_comb_BTag_HZZ_bbbbqq_polm80.eps")

   '''
    file_polm80_hhqq_SignalHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")
    file_polm80_hhqq_SignalHistos_AllEvents_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root")
    file_polm80_hzqq_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__hzqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")

    file_polm80_ee_qq_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__ee_qq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")
    file_polm80_ee_qqqq_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__ee_qqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")
    file_polm80_ee_qqqqqq_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__ee_qqqqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")

    file_polm80_WWH_qqqqH_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__WWH_qqqqH__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")

    file_polm80_ZZH_qqqqH_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__ZZH_qqqqH__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")

    x_polm80_BDTScore = array( 'f' )
    y_polm80_significance = array( 'f' )
    y_polm80_purity = array( 'f' )
    y_polm80_efficiency = array( 'f' )
    norm_polm=file_polm80_hhqq_SignalHistos_.Get("-0.2/h_BDT_output").Integral()

    y_polm80_AllEvents_significance = array( 'f' )
    y_polm80_AllEvents_purity = array( 'f' )
    y_polm80_AllEvents_efficiency = array( 'f' )
    norm_polm_AllEvents=file_polm80_hhqq_SignalHistos_AllEvents_.Get("-0.2/h_BDT_output").Integral()
    #norm_polm=h_mass_sig_hzqq_norm_polm.Integral()
    #norm_polm = 10.
    print 'norm_polm',norm_polm

    for dir_ind in directory:
        x_polm80_BDTScore.append(float(dir_ind))
        h_BTag_sum_all_polm80_sig_hhqq=file_polm80_hhqq_SignalHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_sig_hhqq_AllEvents=file_polm80_hhqq_SignalHistos_AllEvents_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_sig_hzqq_BGHistos=file_polm80_hzqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_sig_ee_qq_BGHistos=file_polm80_ee_qq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos=file_polm80_ee_qqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos=file_polm80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos=file_polm80_WWH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos=file_polm80_ZZH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        #print dir_ind,"integral of signal mass",h_BTag_sum_all_polm80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_BTag_sum_all_polm80_sig_hhqq.Integral()>0:
            #print '#',dir_ind,"polm significance,purity ",h_BTag_sum_all_polm80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()),h_BTag_sum_all_polm80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polm80_efficiency.append(100.*h_BTag_sum_all_polm80_sig_hhqq.Integral()/norm_polm)
            y_polm80_purity.append(100.*h_BTag_sum_all_polm80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos.Integral()))
            y_polm80_significance.append(h_BTag_sum_all_polm80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos.Integral()))
            print 'at polm ',dir_ind,'signif/pur/eff/events',h_BTag_sum_all_polm80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_sig_hhqq.Integral()/norm_polm,h_BTag_sum_all_polm80_sig_hhqq.Integral()
        else:
            y_polm80_efficiency.append(0)
            y_polm80_purity.append(0)
            y_polm80_significance.append(0)
      #print dir_ind,"integral of signal mass",h_BTag_sum_all_polm80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()>0:
            #print '#',dir_ind,"polm significance,purity ",h_BTag_sum_all_polm80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()),h_BTag_sum_all_polm80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polm80_AllEvents_efficiency.append(100.*h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()/norm_polm_AllEvents)
            y_polm80_AllEvents_purity.append(100.*h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()/(h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos.Integral()))
            y_polm80_AllEvents_significance.append(h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()/sqrt(h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos.Integral()))
            print 'at polm ',dir_ind,'signif/pur/eff/events all ',h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()/sqrt(h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()/(h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()/norm_polm,h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral()
            print 'at polm ',dir_ind,'evt sig/sig all/hzqq/qq/qqqq/qqqqqq/WWH/ZZH',h_BTag_sum_all_polm80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral(),h_BTag_sum_all_polm80_sig_hzqq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos.Integral()
        else:
            y_polm80_AllEvents_efficiency.append(0)
            y_polm80_AllEvents_purity.append(0)
            y_polm80_AllEvents_significance.append(0)
    canvas_polm80_BDT_significance = setUpperCanvas("canvas_polm80_BDT_significance");
    canvas_polm80_BDT_significance.cd()
    graph_polm80_significance = TGraph( n_graphs, x_polm80_BDTScore, y_polm80_significance )
    graph_polm80_significance.SetMarkerStyle(20)
    graph_polm80_significance.SetMarkerColor(1)
    graph_polm80_significance.GetXaxis().SetTitle('BDT Score')
    graph_polm80_significance.GetYaxis().SetTitle('significance [#sigma]')
    graph_polm80_significance.Draw('AP')
    #canvas_polm80_BDT_significance.Update()
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polm80_BDT_significance.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polm80_sigificance_vs_BDT.eps")
    
    canvas_polm80_BDT_efficiency = setUpperCanvas("canvas_polm80_BDT_efficiency");
    canvas_polm80_BDT_efficiency.cd()
    graph_polm80_efficiency = TGraph( n_graphs, x_polm80_BDTScore, y_polm80_efficiency )
    graph_polm80_efficiency.SetMarkerStyle(21)
    graph_polm80_efficiency.SetMarkerColor(1)
    graph_polm80_efficiency.GetXaxis().SetTitle('BDT Score')
    graph_polm80_efficiency.GetYaxis().SetTitle('signal efficiency [%]')
    graph_polm80_efficiency.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polm80_BDT_efficiency.Print("hhqq_bbar_vs_totBG_polm80_signal_efficiency_vs_BDT.eps")
    canvas_polm80_BDT_purity = setUpperCanvas("canvas_polm80_BDT_purity");
    canvas_polm80_BDT_purity.cd()
    graph_polm80_purity = TGraph( n_graphs, x_polm80_BDTScore, y_polm80_purity )
    graph_polm80_purity.SetMarkerStyle(22)
    graph_polm80_purity.SetMarkerColor(2)
    graph_polm80_purity.GetXaxis().SetTitle('BDT Score')
    graph_polm80_purity.GetYaxis().SetTitle('purity [%]')
    graph_polm80_purity.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polm80_BDT_purity.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polm80_signal_purity_vs_BDT.eps")


    canvas_polm80_AllEvents_BDT_significance = setUpperCanvas("canvas_polm80_AllEvents_BDT_significance");
    canvas_polm80_AllEvents_BDT_significance.cd()
    graph_polm80_AllEvents_significance = TGraph( n_graphs, x_polm80_BDTScore, y_polm80_AllEvents_significance )
    graph_polm80_AllEvents_significance.SetMarkerStyle(20)
    graph_polm80_AllEvents_significance.SetMarkerColor(1)
    graph_polm80_AllEvents_significance.GetXaxis().SetTitle('BDT Score')
    graph_polm80_AllEvents_significance.GetYaxis().SetTitle('significance [#sigma]')
    graph_polm80_AllEvents_significance.Draw('AP')
    #canvas_polm80_AllEvents_BDT_significance.Update()
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polm80_AllEvents_BDT_significance.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polm80_AllEvents_sigificance_vs_BDT.eps")
    
    canvas_polm80_AllEvents_BDT_efficiency = setUpperCanvas("canvas_polm80_AllEvents_BDT_efficiency");
    canvas_polm80_AllEvents_BDT_efficiency.cd()
    graph_polm80_AllEvents_efficiency = TGraph( n_graphs, x_polm80_BDTScore, y_polm80_AllEvents_efficiency )
    graph_polm80_AllEvents_efficiency.SetMarkerStyle(21)
    graph_polm80_AllEvents_efficiency.SetMarkerColor(1)
    graph_polm80_AllEvents_efficiency.GetXaxis().SetTitle('BDT Score')
    graph_polm80_AllEvents_efficiency.GetYaxis().SetTitle('signal efficiency [%]')
    graph_polm80_AllEvents_efficiency.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polm80_AllEvents_BDT_efficiency.Print("hhqq_bbar_vs_totBG_polm80_AllEvents_signal_efficiency_vs_BDT.eps")
    canvas_polm80_AllEvents_BDT_purity = setUpperCanvas("canvas_polm80_AllEvents_BDT_purity");
    canvas_polm80_AllEvents_BDT_purity.cd()
    graph_polm80_AllEvents_purity = TGraph( n_graphs, x_polm80_BDTScore, y_polm80_AllEvents_purity )
    graph_polm80_AllEvents_purity.SetMarkerStyle(22)
    graph_polm80_AllEvents_purity.SetMarkerColor(2)
    graph_polm80_AllEvents_purity.GetXaxis().SetTitle('BDT Score')
    graph_polm80_AllEvents_purity.GetYaxis().SetTitle('purity [%]')
    graph_polm80_AllEvents_purity.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polm80_AllEvents_BDT_purity.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polm80_AllEvents_signal_purity_vs_BDT.eps")


    file_polp80_hhqq_SignalHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")
    file_polp80_hhqq_SignalHistos_AllEvents_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root")
    file_polp80_hzqq_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__hzqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")
    
    file_polp80_ee_qq_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__ee_qq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")
    file_polp80_ee_qqqq_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__ee_qqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")
    file_polp80_ee_qqqqqq_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__ee_qqqqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")
    
    file_polp80_WWH_qqqqH_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__WWH_qqqqH__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")
    
    file_polp80_ZZH_qqqqH_BGHistos_=root.TFile.Open("/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_2_20_tight_Mass_Cuts/MVATrainingReader_BDT__hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ__HHZ__ZZH_qqqqH__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root")
    
    x_polp80_BDTScore = array( 'f' )
    y_polp80_significance = array( 'f' )
    y_polp80_purity = array( 'f' )
    y_polp80_efficiency = array( 'f' )
    norm_polp=file_polp80_hhqq_SignalHistos_.Get("-0.2/h_BDT_output").Integral()
    
    y_polp80_AllEvents_significance = array( 'f' )
    y_polp80_AllEvents_purity = array( 'f' )
    y_polp80_AllEvents_efficiency = array( 'f' )
    norm_polp_AllEvents=file_polp80_hhqq_SignalHistos_AllEvents_.Get("-0.2/h_BDT_output").Integral()
    #norm_polp=h_mass_sig_hzqq_norm_polp.Integral()
    #norm_polp = 10.
    print 'norm_polp',norm_polp
    
    for dir_ind in directory:
        x_polp80_BDTScore.append(float(dir_ind))
        h_BTag_sum_all_polp80_sig_hhqq=file_polp80_hhqq_SignalHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_hhqq_AllEvents=file_polp80_hhqq_SignalHistos_AllEvents_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_hzqq_BGHistos=file_polp80_hzqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_ee_qq_BGHistos=file_polp80_ee_qq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos=file_polp80_ee_qqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos=file_polp80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos=file_polp80_WWH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos=file_polp80_ZZH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        #print dir_ind,"integral of signal mass",h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_BTag_sum_all_polp80_sig_hhqq.Integral()>0:
            #print '#',dir_ind,"polm significance,purity ",h_BTag_sum_all_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polp80_efficiency.append(100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/norm_polp)
            y_polp80_purity.append(100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
            y_polp80_significance.append(h_BTag_sum_all_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
            print 'at polp ',dir_ind,'signif/pur/eff/events',h_BTag_sum_all_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/norm_polp,h_BTag_sum_all_polp80_sig_hhqq.Integral()
        else:
            y_polp80_efficiency.append(0)
            y_polp80_purity.append(0)
            y_polp80_significance.append(0)
        #print dir_ind,"integral of signal mass",h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()>0:
            #print '#',dir_ind,"polm significance,purity ",h_BTag_sum_all_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polp80_AllEvents_efficiency.append(100.*h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/norm_polp_AllEvents)
            y_polp80_AllEvents_purity.append(100.*h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/(h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
            y_polp80_AllEvents_significance.append(h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
            print 'at polp ',dir_ind,'signif/pur/eff/events all ',h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/(h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/norm_polp,h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()
            print 'at polp ',dir_ind,'evt sig/sig all/hzqq/qq/qqqq/qqqqqq/WWH/ZZH',h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral(),h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos.Integral()
        else:
            y_polp80_AllEvents_efficiency.append(0)
            y_polp80_AllEvents_purity.append(0)
            y_polp80_AllEvents_significance.append(0)
    canvas_polp80_BDT_significance = setUpperCanvas("canvas_polp80_BDT_significance");
    canvas_polp80_BDT_significance.cd()
    graph_polp80_significance = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_significance )
    graph_polp80_significance.SetMarkerStyle(20)
    graph_polp80_significance.SetMarkerColor(1)
    graph_polp80_significance.GetXaxis().SetTitle('BDT Score')
    graph_polp80_significance.GetYaxis().SetTitle('significance [#sigma]')
    graph_polp80_significance.Draw('AP')
    #canvas_polp80_BDT_significance.Update()
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_BDT_significance.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polp80_sigificance_vs_BDT.eps")

    canvas_polp80_BDT_efficiency = setUpperCanvas("canvas_polp80_BDT_efficiency");
    canvas_polp80_BDT_efficiency.cd()
    graph_polp80_efficiency = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_efficiency )
    graph_polp80_efficiency.SetMarkerStyle(21)
    graph_polp80_efficiency.SetMarkerColor(1)
    graph_polp80_efficiency.GetXaxis().SetTitle('BDT Score')
    graph_polp80_efficiency.GetYaxis().SetTitle('signal efficiency [%]')
    graph_polp80_efficiency.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_BDT_efficiency.Print("hhqq_bbar_vs_totBG_polp80_signal_efficiency_vs_BDT.eps")
    canvas_polp80_BDT_purity = setUpperCanvas("canvas_polp80_BDT_purity");
    canvas_polp80_BDT_purity.cd()
    graph_polp80_purity = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_purity )
    graph_polp80_purity.SetMarkerStyle(22)
    graph_polp80_purity.SetMarkerColor(2)
    graph_polp80_purity.GetXaxis().SetTitle('BDT Score')
    graph_polp80_purity.GetYaxis().SetTitle('purity [%]')
    graph_polp80_purity.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_BDT_purity.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polp80_signal_purity_vs_BDT.eps")


    canvas_polp80_AllEvents_BDT_significance = setUpperCanvas("canvas_polp80_AllEvents_BDT_significance");
    canvas_polp80_AllEvents_BDT_significance.cd()
    graph_polp80_AllEvents_significance = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_AllEvents_significance )
    graph_polp80_AllEvents_significance.SetMarkerStyle(20)
    graph_polp80_AllEvents_significance.SetMarkerColor(1)
    graph_polp80_AllEvents_significance.GetXaxis().SetTitle('BDT Score')
    graph_polp80_AllEvents_significance.GetYaxis().SetTitle('significance [#sigma]')
    graph_polp80_AllEvents_significance.Draw('AP')
    #canvas_polp80_AllEvents_BDT_significance.Update()
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_AllEvents_BDT_significance.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polp80_AllEvents_sigificance_vs_BDT.eps")
    
    canvas_polp80_AllEvents_BDT_efficiency = setUpperCanvas("canvas_polp80_AllEvents_BDT_efficiency");
    canvas_polp80_AllEvents_BDT_efficiency.cd()
    graph_polp80_AllEvents_efficiency = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_AllEvents_efficiency )
    graph_polp80_AllEvents_efficiency.SetMarkerStyle(21)
    graph_polp80_AllEvents_efficiency.SetMarkerColor(1)
    graph_polp80_AllEvents_efficiency.GetXaxis().SetTitle('BDT Score')
    graph_polp80_AllEvents_efficiency.GetYaxis().SetTitle('signal efficiency [%]')
    graph_polp80_AllEvents_efficiency.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_AllEvents_BDT_efficiency.Print("hhqq_bbar_vs_totBG_polp80_AllEvents_signal_efficiency_vs_BDT.eps")
    canvas_polp80_AllEvents_BDT_purity = setUpperCanvas("canvas_polp80_AllEvents_BDT_purity");
    canvas_polp80_AllEvents_BDT_purity.cd()
    graph_polp80_AllEvents_purity = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_AllEvents_purity )
    graph_polp80_AllEvents_purity.SetMarkerStyle(22)
    graph_polp80_AllEvents_purity.SetMarkerColor(2)
    graph_polp80_AllEvents_purity.GetXaxis().SetTitle('BDT Score')
    graph_polp80_AllEvents_purity.GetYaxis().SetTitle('purity [%]')
    graph_polp80_AllEvents_purity.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_AllEvents_BDT_purity.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polp80_AllEvents_signal_purity_vs_BDT.eps")
    

    x_polm80_and_polp80_BDTScore = array( 'f' )
    y_polm80_and_polp80_significance = array( 'f' )
    y_polm80_and_polp80_purity = array( 'f' )
    y_polm80_and_polp80_efficiency = array( 'f' )
    norm_polm80_and_polp80=file_polm80_hhqq_SignalHistos_.Get("-0.2/h_BDT_output").Integral()+file_polp80_hhqq_SignalHistos_.Get("-0.2/h_BDT_output").Integral()

    y_polm80_and_polp80_AllEvents_significance = array( 'f' )
    y_polm80_and_polp80_AllEvents_purity = array( 'f' )
    y_polm80_and_polp80_AllEvents_efficiency = array( 'f' )
    norm_polm80_and_polp80_AllEvents=file_polp80_hhqq_SignalHistos_AllEvents_.Get("-0.2/h_BDT_output").Integral()+file_polm80_hhqq_SignalHistos_AllEvents_.Get("-0.2/h_BDT_output").Integral()
    #norm_polm=h_mass_sig_hzqq_norm_polm.Integral()
    #norm_polm = 10.
    print 'norm_polm80_and_polp80',norm_polm80_and_polp80
    
    for dir_ind in directory:
        x_polm80_and_polp80_BDTScore.append(float(dir_ind))
        h_BTag_sum_all_polm80_and_polp80_sig_hhqq=file_polm80_hhqq_SignalHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_hhqq_SignalHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents=file_polm80_hhqq_SignalHistos_AllEvents_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_hhqq_SignalHistos_AllEvents_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos=file_polm80_hzqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_hzqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos=file_polm80_ee_qq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_ee_qq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos=file_polm80_ee_qqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_ee_qqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos=file_polm80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos=file_polm80_WWH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_WWH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos=file_polm80_ZZH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_ZZH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        #print dir_ind,"integral of signal mass",h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()>0:
            #print '#',dir_ind,"polm significance,purity ",h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()),h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polm80_and_polp80_efficiency.append(100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/norm_polm80_and_polp80)
            y_polm80_and_polp80_purity.append(100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
            y_polm80_and_polp80_significance.append(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
            print 'at polm plus polp ',dir_ind,'signif/pur/eff/events',h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/norm_polm80_and_polp80,h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()
        else:
            y_polm80_and_polp80_efficiency.append(0)
            y_polm80_and_polp80_purity.append(0)
            y_polm80_and_polp80_significance.append(0)
        #print dir_ind,"integral of signal mass",h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()>0:
            #print '#',dir_ind,"polm significance,purity ",h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()),h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polm80_and_polp80_AllEvents_efficiency.append(100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/norm_polm80_and_polp80_AllEvents)
            y_polm80_and_polp80_AllEvents_purity.append(100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
            y_polm80_and_polp80_AllEvents_significance.append(h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
            print 'at polm plus polp ',dir_ind,'signif/pur/eff/events all ',h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/norm_polm80_and_polp80,h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()
            print 'at polm plus polp ',dir_ind,'evt sig/sig all/hzqq/qq/qqqq/qqqqqq/WWH/ZZH',h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()
        else:
            y_polm80_and_polp80_AllEvents_efficiency.append(0)
            y_polm80_and_polp80_AllEvents_purity.append(0)
            y_polm80_and_polp80_AllEvents_significance.append(0)

    for dir_ind2 in directory:
        for dir_ind in directory:
            x_polm80_and_polp80_BDTScore.append(float(dir_ind))
            h_BTag_sum_all_polm80_and_polp80_sig_hhqq=file_polm80_hhqq_SignalHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_hhqq_SignalHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents=file_polm80_hhqq_SignalHistos_AllEvents_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_hhqq_SignalHistos_AllEvents_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos=file_polm80_hzqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_hzqq_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos=file_polm80_ee_qq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_ee_qq_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos=file_polm80_ee_qqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_ee_qqqq_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos=file_polm80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_ee_qqqqqq_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos=file_polm80_WWH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_WWH_qqqqH_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos=file_polm80_ZZH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")+file_polp80_ZZH_qqqqH_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_sig_hhqq=file_polm80_hhqq_SignalHistos_.Get(dir_ind+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_sig_hhqq_AllEvents=file_polm80_hhqq_SignalHistos_AllEvents_.Get(dir_ind+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_sig_hzqq_BGHistos=file_polm80_hzqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_sig_ee_qq_BGHistos=file_polm80_ee_qq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos=file_polm80_ee_qqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos=file_polm80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos=file_polm80_WWH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
            h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos=file_polm80_ZZH_qqqqH_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
            h_BTag_sum_all_polp80_sig_hhqq=file_polp80_hhqq_SignalHistos_.Get(dir_ind+"/h_BTag_sum_all")
            h_BTag_sum_all_polp80_sig_hhqq_AllEvents=file_polp80_hhqq_SignalHistos_AllEvents_.Get(dir_ind+"/h_BTag_sum_all")
            h_BTag_sum_all_polp80_sig_hzqq_BGHistos=file_polp80_hzqq_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polp80_sig_ee_qq_BGHistos=file_polp80_ee_qq_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos=file_polp80_ee_qqqq_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos=file_polp80_ee_qqqqqq_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos=file_polp80_WWH_qqqqH_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
            h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos=file_polp80_ZZH_qqqqH_BGHistos_.Get(dir_ind2+"/h_BTag_sum_all")
                    #print dir_ind,"integral of signal mass",h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()
            if h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()>0:
                #print '#',dir_ind,"polm significance,purity ",h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()),h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()
                y_polm80_and_polp80_efficiency.append(100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/norm_polm80_and_polp80)
                y_polm80_and_polp80_purity.append(100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
                y_polm80_and_polp80_significance.append(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
                print 'at polm plus polp different BDT' ,dir_ind,dir_ind2,'signif/pur/eff/events',h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/norm_polm80_and_polp80,h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()
            else:
                y_polm80_and_polp80_efficiency.append(0)
                y_polm80_and_polp80_purity.append(0)
                y_polm80_and_polp80_significance.append(0)
                #print dir_ind,"integral of signal mass",h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()
            if h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()>0:
                #print '#',dir_ind,"polm significance,purity ",h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()),h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()
                y_polm80_and_polp80_AllEvents_efficiency.append(100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/norm_polm80_and_polp80_AllEvents)
                y_polm80_and_polp80_AllEvents_purity.append(100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
                y_polm80_and_polp80_AllEvents_significance.append(h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()))
                print 'at polm plus polp both BDT all',dir_ind,dir_ind2,'signif/pur/eff/events all ',h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/sqrt(h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/(h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral()+h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()),100.*h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()/norm_polm80_and_polp80,h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral()
                print 'at polm plus polp both BDT ',dir_ind,dir_ind2,'evt sig/sig all/hzqq/qq/qqqq/qqqqqq/WWH/ZZH',h_BTag_sum_all_polm80_and_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_hhqq_AllEvents.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_hzqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ee_qqqqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_WWH_qqqqH_BGHistos.Integral(),h_BTag_sum_all_polm80_and_polp80_sig_ZZH_qqqqH_BGHistos.Integral()
                print 'at polm both BDT ',dir_ind,dir_ind2,'evt sig/sig all/hzqq/qq/qqqq/qqqqqq/WWH/ZZH',h_BTag_sum_all_polm80_sig_hhqq.Integral(),h_BTag_sum_all_polm80_sig_hhqq_AllEvents.Integral(),h_BTag_sum_all_polm80_sig_hzqq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ee_qqqqqq_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_WWH_qqqqH_BGHistos.Integral(),h_BTag_sum_all_polm80_sig_ZZH_qqqqH_BGHistos.Integral()
                print 'at polp both BDT ',dir_ind,dir_ind2,'evt sig/sig all/hzqq/qq/qqqq/qqqqqq/WWH/ZZH',h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral(),h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_WWH_qqqqH_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ZZH_qqqqH_BGHistos.Integral()
            else:
                y_polm80_and_polp80_AllEvents_efficiency.append(0)
                y_polm80_and_polp80_AllEvents_purity.append(0)
                y_polm80_and_polp80_AllEvents_significance.append(0)

    """
    canvas_polm80_vs_polp80_BDT_significance = setUpperCanvas("canvas_polm80_vs_polp80_BDT_significance");
    canvas_polm80_vs_polp80_BDT_significance.cd()
    graph_polm80_vs_polp80_significance = TGraph2D( n_graphs, x1_polm80_vs_polp80_BDTScore,x2_polm80_vs_polp80_BDTScore, y_polm80_vs_polp80_significance )
    graph_polm80_vs_polp80_significance.SetMarkerStyle(20)
    graph_polm80_vs_polp80_significance.SetMarkerColor(1)
    graph_polm80_vs_polp80_significance.GetXaxis().SetTitle('BDT Score (-80%)')
    graph_polm80_vs_polp80_significance.GetXaxis().SetTitle('BDT Score (+80%)')
    graph_polm80_vs_polp80_significance.GetZaxis().SetTitle('significance [#sigma]')
    graph_polm80_vs_polp80_significance.Draw('COL')
    #canvas_polm80_vs_polp80_BDT_significance.Update()
    l.DrawLatex(x,y,label);
    l.DrawLatex(x3,y2,label4);
    canvas_polm80_vs_polp80_BDT_significance.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polm80_vs_polp80_signal_purity_vs_BDT_V1.eps")
    print 'do i ever get here 2'
   
    canvas_polm80_vs_polp80_BDT_efficiency = setUpperCanvas("canvas_polm80_vs_polp80_BDT_efficiency");
    canvas_polm80_vs_polp80_BDT_efficiency.cd()
    graph_polm80_vs_polp80_efficiency = TGraph2D( n_graphs, x1_polm80_vs_polp80_BDTScore,x2_polm80_vs_polp80_BDTScore, y_polm80_vs_polp80_efficiency )
    graph_polm80_vs_polp80_efficiency.SetMarkerStyle(21)
    graph_polm80_vs_polp80_efficiency.SetMarkerColor(1)
    graph_polm80_vs_polp80_efficiency.GetXaxis().SetTitle('BDT Score (-80%)')
    graph_polm80_vs_polp80_efficiency.GetXaxis().SetTitle('BDT Score (+80%)')
    graph_polm80_vs_polp80_efficiency.GetZaxis().SetTitle('signal efficiency [%]')
    graph_polm80_vs_polp80_efficiency.Draw('COL')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x3,y2,label2);
    canvas_polm80_vs_polp80_BDT_efficiency.Print("hhqq_bbar_vs_totBG_polm80_vs_polp80_signal_efficiency_vs_BDT.eps")
    print 'do i ever get here 2 b'

    canvas_polm80_vs_polp80_BDT_purity = setUpperCanvas("canvas_polm80_vs_polp80_BDT_purity");
    canvas_polm80_vs_polp80_BDT_purity.cd()
    graph_polm80_vs_polp80_purity = TGraph2D( n_graphs, x1_polm80_vs_polp80_BDTScore,x2_polm80_vs_polp80_BDTScore, y_polm80_vs_polp80_purity )
    graph_polm80_vs_polp80_purity.SetMarkerStyle(22)
    graph_polm80_vs_polp80_purity.SetMarkerColor(2)
    graph_polm80_vs_polp80_purity.GetXaxis().SetTitle('BDT Score (-80%)')
    graph_polm80_vs_polp80_purity.GetXaxis().SetTitle('BDT Score (+80%)')
    graph_polm80_vs_polp80_purity.GetZaxis().SetTitle('purity [%]')
    graph_polm80_vs_polp80_purity.Draw('COL')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x3,y2,label4);
    canvas_polm80_vs_polp80_BDT_purity.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polm80_vs_polp80_signal_purity_vs_BDT.eps")

    print 'do i ever get here 3'

    canvas_polm80_vs_polp80_AllEvents_BDT_significance = setUpperCanvas("canvas_polm80_vs_polp80_AllEvents_BDT_significance");
    canvas_polm80_vs_polp80_AllEvents_BDT_significance.cd()
    graph_polm80_vs_polp80_AllEvents_significance = TGraph2D( n_graphs, x1_polm80_vs_polp80_BDTScore,x2_polm80_vs_polp80_BDTScore, y_polm80_vs_polp80_AllEvents_significance )
    graph_polm80_vs_polp80_AllEvents_significance.SetMarkerStyle(20)
    graph_polm80_vs_polp80_AllEvents_significance.SetMarkerColor(1)
    graph_polm80_vs_polp80_AllEvents_significance.GetXaxis().SetTitle('BDT Score (-80%)')
    graph_polm80_vs_polp80_AllEvents_significance.GetXaxis().SetTitle('BDT Score (+80%)')
    graph_polm80_vs_polp80_AllEvents_significance.GetZaxis().SetTitle('significance [#sigma]')

    graph_polm80_vs_polp80_AllEvents_significance.Draw('COL')
    #canvas_polm80_vs_polp80_AllEvents_BDT_significance.Update()
    l.DrawLatex(x,y,label);
    l.DrawLatex(x3,y2,label2);
    canvas_polm80_vs_polp80_AllEvents_BDT_significance.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polm80_vs_polp80_AllEvents_sigificance_vs_BDT.eps")
    
    canvas_polm80_vs_polp80_AllEvents_BDT_efficiency = setUpperCanvas("canvas_polm80_vs_polp80_AllEvents_BDT_efficiency");
    canvas_polm80_vs_polp80_AllEvents_BDT_efficiency.cd()
    graph_polm80_vs_polp80_AllEvents_efficiency = TGraph2D( n_graphs, x1_polm80_vs_polp80_BDTScore,x2_polm80_vs_polp80_BDTScore, y_polm80_vs_polp80_AllEvents_efficiency )
    graph_polm80_vs_polp80_AllEvents_efficiency.SetMarkerStyle(21)
    graph_polm80_vs_polp80_AllEvents_efficiency.SetMarkerColor(1)
    graph_polm80_vs_polp80_AllEvents_efficiency.GetXaxis().SetTitle('BDT Score (-80%)')
    graph_polm80_vs_polp80_AllEvents_efficiency.GetXaxis().SetTitle('BDT Score (+80%)')
    graph_polm80_vs_polp80_AllEvents_efficiency.GetZaxis().SetTitle('signal efficiency [%]')
    graph_polm80_vs_polp80_AllEvents_efficiency.Draw('COL')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x3,y2,label4);
    canvas_polm80_vs_polp80_AllEvents_BDT_efficiency.Print("hhqq_bbar_vs_totBG_polm80_vs_polp80_AllEvents_signal_efficiency_vs_BDT.eps")
    canvas_polm80_vs_polp80_AllEvents_BDT_purity = setUpperCanvas("canvas_polm80_vs_polp80_AllEvents_BDT_purity");
    canvas_polm80_vs_polp80_AllEvents_BDT_purity.cd()
    graph_polm80_vs_polp80_AllEvents_purity = TGraph2D( n_graphs, x1_polm80_vs_polp80_BDTScore,x2_polm80_vs_polp80_BDTScore, y_polm80_vs_polp80_AllEvents_purity )
    graph_polm80_vs_polp80_AllEvents_purity.SetMarkerStyle(22)
    graph_polm80_vs_polp80_AllEvents_purity.SetMarkerColor(2)
    graph_polm80_vs_polp80_AllEvents_purity.GetXaxis().SetTitle('BDT Score (-80%)')
    graph_polm80_vs_polp80_AllEvents_purity.GetXaxis().SetTitle('BDT Score (+80%)')
    graph_polm80_vs_polp80_AllEvents_purity.GetZaxis().SetTitle('purity [%]')
    graph_polm80_vs_polp80_AllEvents_purity.Draw('COL')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x3,y2,label4);
    canvas_polm80_vs_polp80_AllEvents_BDT_purity.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polm80_vs_polp80_AllEvents_signal_purity_vs_BDT.eps")


    print 'do i ever get here 4'

    x_polp80_BDTScore = array( 'f' )
    y_polp80_significance = array( 'f' )
    y_polp80_purity = array( 'f' )
    y_polp80_efficiency = array( 'f' )
    norm_polp=file_polp80_hhqq_SignalHistos_.Get("-0.2/h_BDT_output").Integral()

    y_polp80_AllEvents_significance = array( 'f' )
    y_polp80_AllEvents_purity = array( 'f' )
    y_polp80_AllEvents_efficiency = array( 'f' )
    norm_polp_AllEvents=file_polp80_hhqq_SignalHistos_AllEvents_.Get("-0.2/h_BDT_output").Integral()
    #norm_polp=h_mass_sig_hzqq_norm_polp.Integral()
    #norm_polp = 10.
    print 'norm_polp',norm_polp

    for dir_ind in directory:
        x_polp80_BDTScore.append(float(dir_ind))
        h_BTag_sum_all_polp80_sig_hhqq=file_polp80_hhqq_SignalHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_hhqq_AllEvents=file_polp80_hhqq_SignalHistos_AllEvents_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_hzqq_BGHistos=file_polp80_hzqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_ee_qq_BGHistos=file_polp80_ee_qq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos=file_polp80_ee_qqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all")
        h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos=file_polp80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_BTag_sum_all") 
        #print dir_ind,"integral of signal mass",h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_BTag_sum_all_polp80_sig_hhqq.Integral()>0:
            #print '#',dir_ind,"polp significance,purity ",h_BTag_sum_all_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polp80_efficiency.append(100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/norm_polp)
            y_polp80_purity.append(100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()))
            y_polp80_significance.append(h_BTag_sum_all_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()))
            print 'at polp ',dir_ind,'signif/pur/eff/events',h_BTag_sum_all_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/norm_polp,h_BTag_sum_all_polp80_sig_hhqq.Integral()
        else:
            y_polp80_efficiency.append(0)
            y_polp80_purity.append(0)
            y_polp80_significance.append(0)
      #print dir_ind,"integral of signal mass",h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()>0:
            #print '#',dir_ind,"polp significance,purity ",h_BTag_sum_all_polp80_sig_hhqq.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq.Integral()/(h_BTag_sum_all_polp80_sig_hhqq.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polp80_AllEvents_efficiency.append(100.*h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/norm_polp_AllEvents)
            y_polp80_AllEvents_purity.append(100.*h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/(h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()))
            y_polp80_AllEvents_significance.append(h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()))
            print 'at polp ',dir_ind,'signif/pur/eff/events all ',h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/sqrt(h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/(h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()+h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral()+h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()/norm_polp,h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral()
            print 'at polp ',dir_ind,'evt sig/sig all/hzqq/qq/qqqq/qqqqqq',h_BTag_sum_all_polp80_sig_hhqq.Integral(),h_BTag_sum_all_polp80_sig_hhqq_AllEvents.Integral(),h_BTag_sum_all_polp80_sig_hzqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqq_BGHistos.Integral(),h_BTag_sum_all_polp80_sig_ee_qqqqqq_BGHistos.Integral()
        else:
            y_polp80_AllEvents_efficiency.append(0)
            y_polp80_AllEvents_purity.append(0)
            y_polp80_AllEvents_significance.append(0)


    canvas_polp80_BDT_significance = setUpperCanvas("canvas_polp80_BDT_significance");
    canvas_polp80_BDT_significance.cd()
    graph_polp80_significance = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_significance )
    graph_polp80_significance.SetMarkerStyle(20)
    graph_polp80_significance.SetMarkerColor(1)
    graph_polp80_significance.GetXaxis().SetTitle('BDT Score')
    graph_polp80_significance.GetYaxis().SetTitle('significance [#sigma]')
    graph_polp80_significance.Draw('AP')
    #canvas_polp80_BDT_significance.Update()
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_BDT_significance.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polp80_sigificance_vs_BDT.eps")
    
    canvas_polp80_BDT_efficiency = setUpperCanvas("canvas_polp80_BDT_efficiency");
    canvas_polp80_BDT_efficiency.cd()
    graph_polp80_efficiency = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_efficiency )
    graph_polp80_efficiency.SetMarkerStyle(21)
    graph_polp80_efficiency.SetMarkerColor(1)
    graph_polp80_efficiency.GetXaxis().SetTitle('BDT Score')
    graph_polp80_efficiency.GetYaxis().SetTitle('signal efficiency [%]')
    graph_polp80_efficiency.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_BDT_efficiency.Print("hhqq_bbar_vs_totBG_polp80_signal_efficiency_vs_BDT.eps")
    canvas_polp80_BDT_purity = setUpperCanvas("canvas_polp80_BDT_purity");
    canvas_polp80_BDT_purity.cd()
    graph_polp80_purity = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_purity )
    graph_polp80_purity.SetMarkerStyle(22)
    graph_polp80_purity.SetMarkerColor(2)
    graph_polp80_purity.GetXaxis().SetTitle('BDT Score')
    graph_polp80_purity.GetYaxis().SetTitle('purity [%]')
    graph_polp80_purity.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_BDT_purity.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polp80_signal_purity_vs_BDT.eps")


    canvas_polp80_AllEvents_BDT_significance = setUpperCanvas("canvas_polp80_AllEvents_BDT_significance");
    canvas_polp80_AllEvents_BDT_significance.cd()
    graph_polp80_AllEvents_significance = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_AllEvents_significance )
    graph_polp80_AllEvents_significance.SetMarkerStyle(20)
    graph_polp80_AllEvents_significance.SetMarkerColor(1)
    graph_polp80_AllEvents_significance.GetXaxis().SetTitle('BDT Score')
    graph_polp80_AllEvents_significance.GetYaxis().SetTitle('significance [#sigma]')
    graph_polp80_AllEvents_significance.Draw('AP')
    #canvas_polp80_AllEvents_BDT_significance.Update()
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_AllEvents_BDT_significance.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polp80_AllEvents_sigificance_vs_BDT.eps")
    
    canvas_polp80_AllEvents_BDT_efficiency = setUpperCanvas("canvas_polp80_AllEvents_BDT_efficiency");
    canvas_polp80_AllEvents_BDT_efficiency.cd()
    graph_polp80_AllEvents_efficiency = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_AllEvents_efficiency )
    graph_polp80_AllEvents_efficiency.SetMarkerStyle(21)
    graph_polp80_AllEvents_efficiency.SetMarkerColor(1)
    graph_polp80_AllEvents_efficiency.GetXaxis().SetTitle('BDT Score')
    graph_polp80_AllEvents_efficiency.GetYaxis().SetTitle('signal efficiency [%]')
    graph_polp80_AllEvents_efficiency.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_AllEvents_BDT_efficiency.Print("hhqq_bbar_vs_totBG_polp80_AllEvents_signal_efficiency_vs_BDT.eps")
    canvas_polp80_AllEvents_BDT_purity = setUpperCanvas("canvas_polp80_AllEvents_BDT_purity");
    canvas_polp80_AllEvents_BDT_purity.cd()
    graph_polp80_AllEvents_purity = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_AllEvents_purity )
    graph_polp80_AllEvents_purity.SetMarkerStyle(22)
    graph_polp80_AllEvents_purity.SetMarkerColor(2)
    graph_polp80_AllEvents_purity.GetXaxis().SetTitle('BDT Score')
    graph_polp80_AllEvents_purity.GetYaxis().SetTitle('purity [%]')
    graph_polp80_AllEvents_purity.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_polp80_AllEvents_BDT_purity.Print("~/plotsHHqq_VLC11_rfJets_BTag_NJet6/hhqq_bbar_vs_totBG_polp80_AllEvents_signal_purity_vs_BDT.eps")

    """
    

    return None


process_files()
#root.gApplication.Run()

#
