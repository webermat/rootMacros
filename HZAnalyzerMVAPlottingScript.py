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

    directory=[];
    directory.append("-0.2")
    directory.append("-0.15")
    directory.append("-0.1")
    directory.append("-0.05")
    directory.append("0.0")
    directory.append("0.05")
    directory.append("0.1")
    directory.append("0.15")
    directory.append("0.175")
    directory.append("0.2")
    directory.append("0.225")
    directory.append("0.25")
    directory.append("0.275")
    directory.append("0.3")
    directory.append("0.325")
    directory.append("0.35")
    directory.append("0.375")
    directory.append("0.4")
    directory.append("0.425")
    directory.append("0.45")
    directory.append("0.5")
    directory.append("0.55")
    directory.append("0.6")
    directory.append("0.65")
    directory.append("0.7")
    directory.append("0.75")
    directory.append("0.8")
    directory.append("0.85")
    directory.append("0.875")
    directory.append("0.9")
    directory.append("0.925")
    directory.append("0.95")
    directory.append("0.975")
    directory.append("1.0")

    n_graphs=int(len(directory))
    #MVATrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root
    #jet1_D2_C2_C3_tau21_d21_N2_N3_jet2_D2_C2_C3_tau21_d21_N2_N3
    #jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21
    #jet1_D2_C2_C3_jet2_D2_C2_C3
    #jet1_D2_C2_C3_tau21_d21_N2_N3_jet2_no_subStruct #jet1_D2_C2_C3_tau21_jet2_no_subStruct
    #jet1_D2_C2_C3_jet2_no_subStruct

    #AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root

    
    file_polm80_HZ_SignalHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root") 
    file_polm80_HZ_SignalHistos_AllEvents_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21_AllEvents.root") 
    file_polm80_HZ_SignalHistos_AllEventsParton_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_hzqq_noMass_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_withPartonHistos_AllEvents.root") 
    #file_polm80_HZ_SignalHistos_AllEvents_noCorr_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_AllEvents.root") 
    file_polm80_ee_qq_mqq_1TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root") 
    file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root") 
    file_polm80_ee_qqqqqq_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root") 


    #file_polm80_HZ_SignalHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21_BTagUp0_01_origBDT.root") 
    #file_polm80_HZ_SignalHistos_AllEvents_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21_BTagUp0_01_origBDT_AllEvents.root") 
    #file_polm80_HZ_SignalHistos_AllEventsParton_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_hzqq_noMass_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_withPartonHistos_AllEvents.root") 
    #file_polm80_HZ_SignalHistos_AllEvents_noCorr_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_AllEvents_origBDT.root") 
    #file_polm80_ee_qq_mqq_1TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21_BTagUp0_01_origBDT.root") 
    #file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21_BTagUp0_01_origBDT.root") 
    #file_polm80_ee_qqqqqq_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21_BTagUp0_01_origBDT.root") 



    h_2D_polm80_HZ_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polm80_HZ_SignalHistos_.Get("0.375/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polm80_HZ_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polm80_HZ_SignalHistos_.Get("0.375/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    print h_2D_polm80_HZ_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.GetBinContent(40,40)

    h_2D_polm80_HZ_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polm80_HZ_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")

    h_2D_polm80_HZ_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polm80_HZ_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
  
    h_2D_polm80_HZ_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polm80_HZ_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")

    h_2D_polm80_HZ_AllEvents_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_2D_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM")
    h_2D_polm80_HZ_AllEvents_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_2D_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM")
  

    h_2D_polm80_ee_qq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polm80_ee_qq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")

    h_2D_polm80_ee_qqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polm80_ee_qqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")

    h_2D_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
 
    A_pos_polm80_theta1_vs_theta2_signal=h_2D_polm80_HZ_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    B_neg_polm80_theta1_vs_theta2_signal=h_2D_polm80_HZ_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    alpha_polm80_theta1_vs_theta2_HZ=(A_pos_polm80_theta1_vs_theta2_signal+B_neg_polm80_theta1_vs_theta2_signal)/(abs(A_pos_polm80_theta1_vs_theta2_signal)+abs(B_neg_polm80_theta1_vs_theta2_signal))
    alpha_polm80_theta1_vs_theta2_HZ_error=sqrt(4.*abs(A_pos_polm80_theta1_vs_theta2_signal)*abs(B_neg_polm80_theta1_vs_theta2_signal)/pow((abs(A_pos_polm80_theta1_vs_theta2_signal)+abs(B_neg_polm80_theta1_vs_theta2_signal)),3))

    A_pos_polm80_theta1_vs_theta2_signal_all=h_2D_polm80_HZ_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    B_neg_polm80_theta1_vs_theta2_signal_all=h_2D_polm80_HZ_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    alpha_polm80_theta1_vs_theta2_HZ_all=(A_pos_polm80_theta1_vs_theta2_signal_all+B_neg_polm80_theta1_vs_theta2_signal_all)/(abs(A_pos_polm80_theta1_vs_theta2_signal_all)+abs(B_neg_polm80_theta1_vs_theta2_signal_all))
    alpha_polm80_theta1_vs_theta2_HZ_all_error=sqrt(4.*abs(A_pos_polm80_theta1_vs_theta2_signal_all)*abs(B_neg_polm80_theta1_vs_theta2_signal_all)/pow((abs(A_pos_polm80_theta1_vs_theta2_signal_all)+abs(B_neg_polm80_theta1_vs_theta2_signal_all)),3))
    
    A_pos_parton_polm80_theta1_vs_theta2_signal_all=h_2D_polm80_HZ_AllEvents_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Integral()
    B_neg_parton_polm80_theta1_vs_theta2_signal_all=h_2D_polm80_HZ_AllEvents_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Integral()
    alpha_parton_polm80_theta1_vs_theta2_HZ_all=(A_pos_parton_polm80_theta1_vs_theta2_signal_all+B_neg_parton_polm80_theta1_vs_theta2_signal_all)/(abs(A_pos_parton_polm80_theta1_vs_theta2_signal_all)+abs(B_neg_parton_polm80_theta1_vs_theta2_signal_all))
    alpha_parton_polm80_theta1_vs_theta2_HZ_all_error=sqrt(4.*abs(A_pos_parton_polm80_theta1_vs_theta2_signal_all)*abs(B_neg_parton_polm80_theta1_vs_theta2_signal_all)/pow((abs(A_pos_parton_polm80_theta1_vs_theta2_signal_all)+abs(B_neg_parton_polm80_theta1_vs_theta2_signal_all)),3))


    A_pos_polm80_theta1_vs_theta2_signal_BDTm020=h_2D_polm80_HZ_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020.Integral()
    B_neg_polm80_theta1_vs_theta2_signal_BDTm020=h_2D_polm80_HZ_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020.Integral()
    alpha_polm80_theta1_vs_theta2_HZ_BDTm020=(A_pos_polm80_theta1_vs_theta2_signal_BDTm020+B_neg_polm80_theta1_vs_theta2_signal_BDTm020)/(abs(A_pos_polm80_theta1_vs_theta2_signal_BDTm020)+abs(B_neg_polm80_theta1_vs_theta2_signal_BDTm020))
    alpha_polm80_theta1_vs_theta2_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polm80_theta1_vs_theta2_signal_BDTm020)*abs(B_neg_polm80_theta1_vs_theta2_signal_BDTm020)/pow((abs(A_pos_polm80_theta1_vs_theta2_signal_BDTm020)+abs(B_neg_polm80_theta1_vs_theta2_signal_BDTm020)),3))
    A_pos_polm80_theta1_vs_theta2_signal_all_BDTm020=h_2D_polm80_HZ_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020.Integral()
    B_neg_polm80_theta1_vs_theta2_signal_all_BDTm020=h_2D_polm80_HZ_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020.Integral()
    alpha_polm80_theta1_vs_theta2_HZ_all_BDTm020=(A_pos_polm80_theta1_vs_theta2_signal_all_BDTm020+B_neg_polm80_theta1_vs_theta2_signal_all_BDTm020)/(abs(A_pos_polm80_theta1_vs_theta2_signal_all_BDTm020)+abs(B_neg_polm80_theta1_vs_theta2_signal_all_BDTm020))
    alpha_polm80_theta1_vs_theta2_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polm80_theta1_vs_theta2_signal_all_BDTm020)*abs(B_neg_polm80_theta1_vs_theta2_signal_all_BDTm020)/pow((abs(A_pos_polm80_theta1_vs_theta2_signal_all_BDTm020)+abs(B_neg_polm80_theta1_vs_theta2_signal_all_BDTm020)),3))

    A_pos_polm80_theta1_vs_theta2_signal_all_withBG=h_2D_polm80_HZ_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polm80_ee_qq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polm80_ee_qqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    B_neg_polm80_theta1_vs_theta2_signal_all_withBG=h_2D_polm80_HZ_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polm80_ee_qq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polm80_ee_qqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    alpha_polm80_theta1_vs_theta2_HZ_all_withBG=(A_pos_polm80_theta1_vs_theta2_signal_all_withBG+B_neg_polm80_theta1_vs_theta2_signal_all_withBG)/(abs(A_pos_polm80_theta1_vs_theta2_signal_all_withBG)+abs(B_neg_polm80_theta1_vs_theta2_signal_all_withBG))
    alpha_polm80_theta1_vs_theta2_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polm80_theta1_vs_theta2_signal_all_withBG)*abs(B_neg_polm80_theta1_vs_theta2_signal_all_withBG)/pow((abs(A_pos_polm80_theta1_vs_theta2_signal_all_withBG)+abs(B_neg_polm80_theta1_vs_theta2_signal_all_withBG)),3))

    A_pos_polm80_theta1_vs_theta2_allBG=h_2D_polm80_ee_qq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polm80_ee_qqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    B_neg_polm80_theta1_vs_theta2_allBG=h_2D_polm80_ee_qq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polm80_ee_qqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    alpha_polm80_theta1_vs_theta2_allBG=(A_pos_polm80_theta1_vs_theta2_allBG+B_neg_polm80_theta1_vs_theta2_allBG)/(abs(A_pos_polm80_theta1_vs_theta2_allBG)+abs(B_neg_polm80_theta1_vs_theta2_allBG))
    alpha_polm80_theta1_vs_theta2_allBG_error=sqrt(4.*abs(A_pos_polm80_theta1_vs_theta2_allBG)*abs(B_neg_polm80_theta1_vs_theta2_allBG)/pow((abs(A_pos_polm80_theta1_vs_theta2_allBG)+abs(B_neg_polm80_theta1_vs_theta2_allBG)),3))

    print 'alpha_polm80_theta1_vs_theta2_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_theta1_vs_theta2_HZ,alpha_polm80_theta1_vs_theta2_HZ_all,alpha_polm80_theta1_vs_theta2_HZ_BDTm020,alpha_polm80_theta1_vs_theta2_HZ_all_BDTm020,alpha_polm80_theta1_vs_theta2_HZ_all_withBG
    print 'errors of alpha_polm80_theta1_vs_theta2_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_theta1_vs_theta2_HZ_error,alpha_polm80_theta1_vs_theta2_HZ_all_error,alpha_polm80_theta1_vs_theta2_HZ_BDTm020_error,alpha_polm80_theta1_vs_theta2_HZ_all_BDTm020_error,alpha_polm80_theta1_vs_theta2_HZ_all_withBG_error


    print 'alpha_polm80_theta1_vs_theta2_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polm80_theta1_vs_theta2_HZ,alpha_polm80_theta1_vs_theta2_HZ_all,alpha_polm80_theta1_vs_theta2_allBG,alpha_polm80_theta1_vs_theta2_HZ_all_withBG,alpha_parton_polm80_theta1_vs_theta2_HZ_all
    print 'errors of alpha_polm80_theta1_vs_theta2_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polm80_theta1_vs_theta2_HZ_error,alpha_polm80_theta1_vs_theta2_HZ_all_error,alpha_polm80_theta1_vs_theta2_allBG_error,alpha_polm80_theta1_vs_theta2_HZ_all_withBG_error,alpha_parton_polm80_theta1_vs_theta2_HZ_all_error



    h_polm80_HZ_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_HZ_SignalHistos_.Get("0.375/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polm80_HZ_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_HZ_SignalHistos_.Get("0.375/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polm80_HZ_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polm80_HZ_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")

    h_polm80_HZ_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polm80_HZ_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polm80_HZ_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polm80_HZ_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")

    h_polm80_HZ_AllEvents_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom")
    h_polm80_HZ_AllEvents_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom")
  

    h_polm80_ee_qq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polm80_ee_qq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")

    h_polm80_ee_qqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polm80_ee_qqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")

    h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")

    A_pos_polm80_theta1_signal=h_polm80_HZ_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    B_neg_polm80_theta1_signal=h_polm80_HZ_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    alpha_polm80_theta1_HZ=(A_pos_polm80_theta1_signal+B_neg_polm80_theta1_signal)/(abs(A_pos_polm80_theta1_signal)+abs(B_neg_polm80_theta1_signal))
    alpha_polm80_theta1_HZ_error=sqrt(4.*abs(A_pos_polm80_theta1_signal)*abs(B_neg_polm80_theta1_signal)/pow((abs(A_pos_polm80_theta1_signal)+abs(B_neg_polm80_theta1_signal)),3))

    A_pos_polm80_theta1_signal_all=h_polm80_HZ_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    B_neg_polm80_theta1_signal_all=h_polm80_HZ_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    alpha_polm80_theta1_HZ_all=(A_pos_polm80_theta1_signal_all+B_neg_polm80_theta1_signal_all)/(abs(A_pos_polm80_theta1_signal_all)+abs(B_neg_polm80_theta1_signal_all))
    alpha_polm80_theta1_HZ_all_error=sqrt(4.*abs(A_pos_polm80_theta1_signal_all)*abs(B_neg_polm80_theta1_signal_all)/pow((abs(A_pos_polm80_theta1_signal_all)+abs(B_neg_polm80_theta1_signal_all)),3))
 
    A_pos_parton_polm80_theta1_signal_all=h_polm80_HZ_AllEvents_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral()
    B_neg_parton_polm80_theta1_signal_all=h_polm80_HZ_AllEvents_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral()
    alpha_parton_polm80_theta1_HZ_all=(A_pos_parton_polm80_theta1_signal_all+B_neg_parton_polm80_theta1_signal_all)/(abs(A_pos_parton_polm80_theta1_signal_all)+abs(B_neg_parton_polm80_theta1_signal_all))
    alpha_parton_polm80_theta1_HZ_all_error=sqrt(4.*abs(A_pos_parton_polm80_theta1_signal_all)*abs(B_neg_parton_polm80_theta1_signal_all)/pow((abs(A_pos_parton_polm80_theta1_signal_all)+abs(B_neg_parton_polm80_theta1_signal_all)),3))



    A_pos_polm80_theta1_signal_BDTm020=h_polm80_HZ_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020.Integral()
    B_neg_polm80_theta1_signal_BDTm020=h_polm80_HZ_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020.Integral()
    alpha_polm80_theta1_HZ_BDTm020=(A_pos_polm80_theta1_signal_BDTm020+B_neg_polm80_theta1_signal_BDTm020)/(abs(A_pos_polm80_theta1_signal_BDTm020)+abs(B_neg_polm80_theta1_signal_BDTm020))
    alpha_polm80_theta1_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polm80_theta1_signal_BDTm020)*abs(B_neg_polm80_theta1_signal_BDTm020)/pow((abs(A_pos_polm80_theta1_signal_BDTm020)+abs(B_neg_polm80_theta1_signal_BDTm020)),3))
    A_pos_polm80_theta1_signal_all_BDTm020=h_polm80_HZ_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020.Integral()
    B_neg_polm80_theta1_signal_all_BDTm020=h_polm80_HZ_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020.Integral()
    alpha_polm80_theta1_HZ_all_BDTm020=(A_pos_polm80_theta1_signal_all_BDTm020+B_neg_polm80_theta1_signal_all_BDTm020)/(abs(A_pos_polm80_theta1_signal_all_BDTm020)+abs(B_neg_polm80_theta1_signal_all_BDTm020))
    alpha_polm80_theta1_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polm80_theta1_signal_all_BDTm020)*abs(B_neg_polm80_theta1_signal_all_BDTm020)/pow((abs(A_pos_polm80_theta1_signal_all_BDTm020)+abs(B_neg_polm80_theta1_signal_all_BDTm020)),3))

    A_pos_polm80_theta1_signal_all_withBG=h_polm80_HZ_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polm80_ee_qq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    B_neg_polm80_theta1_signal_all_withBG=h_polm80_HZ_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polm80_ee_qq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    alpha_polm80_theta1_HZ_all_withBG=(A_pos_polm80_theta1_signal_all_withBG+B_neg_polm80_theta1_signal_all_withBG)/(abs(A_pos_polm80_theta1_signal_all_withBG)+abs(B_neg_polm80_theta1_signal_all_withBG))
    alpha_polm80_theta1_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polm80_theta1_signal_all_withBG)*abs(B_neg_polm80_theta1_signal_all_withBG)/pow((abs(A_pos_polm80_theta1_signal_all_withBG)+abs(B_neg_polm80_theta1_signal_all_withBG)),3))

    A_pos_polm80_theta1_allBG=h_polm80_ee_qq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    B_neg_polm80_theta1_allBG=h_polm80_ee_qq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    alpha_polm80_theta1_allBG=(A_pos_polm80_theta1_allBG+B_neg_polm80_theta1_allBG)/(abs(A_pos_polm80_theta1_allBG)+abs(B_neg_polm80_theta1_allBG))
    alpha_polm80_theta1_allBG_error=sqrt(4.*abs(A_pos_polm80_theta1_allBG)*abs(B_neg_polm80_theta1_allBG)/pow((abs(A_pos_polm80_theta1_allBG)+abs(B_neg_polm80_theta1_allBG)),3))
  
    print 'alpha_polm80_theta1_HZ/HZallevents/HZ BDTm020/HZ all BDTm020/ BG',alpha_polm80_theta1_HZ,alpha_polm80_theta1_HZ_all,alpha_polm80_theta1_HZ_BDTm020,alpha_polm80_theta1_HZ_all_BDTm020,alpha_polm80_theta1_HZ_all_withBG
    print 'errors of alpha_polm80_theta1_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_theta1_HZ_error,alpha_polm80_theta1_HZ_all_error,alpha_polm80_theta1_HZ_BDTm020_error,alpha_polm80_theta1_HZ_all_BDTm020_error,alpha_polm80_theta1_HZ_all_withBG_error

    print 'alpha_polm80_theta1_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polm80_theta1_HZ,alpha_polm80_theta1_HZ_all,alpha_polm80_theta1_allBG,alpha_polm80_theta1_HZ_all_withBG,alpha_parton_polm80_theta1_HZ_all
    print 'errors of alpha_polm80_theta1_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polm80_theta1_HZ_error,alpha_polm80_theta1_HZ_all_error,alpha_polm80_theta1_allBG_error,alpha_polm80_theta1_HZ_all_withBG_error,alpha_parton_polm80_theta1_HZ_all_error



    h_polm80_HZ_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_.Get("0.375/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_.Get("0.375/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_HZ_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_HZ_AllEvents_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep")
    h_polm80_HZ_AllEvents_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep")

    h_polm80_HZ_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")


    A_pos_polm80_phi1_signal=h_polm80_HZ_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi1_signal=h_polm80_HZ_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi1_HZ=(A_pos_polm80_phi1_signal+B_neg_polm80_phi1_signal)/(abs(A_pos_polm80_phi1_signal)+abs(B_neg_polm80_phi1_signal))
    alpha_polm80_phi1_HZ_error=sqrt(4.*abs(A_pos_polm80_phi1_signal)*abs(B_neg_polm80_phi1_signal)/pow((abs(A_pos_polm80_phi1_signal)+abs(B_neg_polm80_phi1_signal)),3))

    A_pos_polm80_phi1_signal_all=h_polm80_HZ_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi1_signal_all=h_polm80_HZ_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi1_HZ_all=(A_pos_polm80_phi1_signal_all+B_neg_polm80_phi1_signal_all)/(abs(A_pos_polm80_phi1_signal_all)+abs(B_neg_polm80_phi1_signal_all))
    alpha_polm80_phi1_HZ_all_error=sqrt(4.*abs(A_pos_polm80_phi1_signal_all)*abs(B_neg_polm80_phi1_signal_all)/pow((abs(A_pos_polm80_phi1_signal_all)+abs(B_neg_polm80_phi1_signal_all)),3))
    
    A_pos_parton_polm80_phi1_signal_all=h_polm80_HZ_AllEvents_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    B_neg_parton_polm80_phi1_signal_all=h_polm80_HZ_AllEvents_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    alpha_parton_polm80_phi1_HZ_all=(A_pos_parton_polm80_phi1_signal_all+B_neg_parton_polm80_phi1_signal_all)/(abs(A_pos_parton_polm80_phi1_signal_all)+abs(B_neg_parton_polm80_phi1_signal_all))
    alpha_parton_polm80_phi1_HZ_all_error=sqrt(4.*abs(A_pos_parton_polm80_phi1_signal_all)*abs(B_neg_parton_polm80_phi1_signal_all)/pow((abs(A_pos_parton_polm80_phi1_signal_all)+abs(B_neg_parton_polm80_phi1_signal_all)),3))



    A_pos_polm80_phi1_signal_BDTm020=h_polm80_HZ_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polm80_phi1_signal_BDTm020=h_polm80_HZ_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polm80_phi1_HZ_BDTm020=(A_pos_polm80_phi1_signal_BDTm020+B_neg_polm80_phi1_signal_BDTm020)/(abs(A_pos_polm80_phi1_signal_BDTm020)+abs(B_neg_polm80_phi1_signal_BDTm020))
    alpha_polm80_phi1_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polm80_phi1_signal_BDTm020)*abs(B_neg_polm80_phi1_signal_BDTm020)/pow((abs(A_pos_polm80_phi1_signal_BDTm020)+abs(B_neg_polm80_phi1_signal_BDTm020)),3))
    A_pos_polm80_phi1_signal_all_BDTm020=h_polm80_HZ_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polm80_phi1_signal_all_BDTm020=h_polm80_HZ_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polm80_phi1_HZ_all_BDTm020=(A_pos_polm80_phi1_signal_all_BDTm020+B_neg_polm80_phi1_signal_all_BDTm020)/(abs(A_pos_polm80_phi1_signal_all_BDTm020)+abs(B_neg_polm80_phi1_signal_all_BDTm020))
    alpha_polm80_phi1_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polm80_phi1_signal_all_BDTm020)*abs(B_neg_polm80_phi1_signal_all_BDTm020)/pow((abs(A_pos_polm80_phi1_signal_all_BDTm020)+abs(B_neg_polm80_phi1_signal_all_BDTm020)),3))

    A_pos_polm80_phi1_signal_all_withBG=h_polm80_HZ_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi1_signal_all_withBG=h_polm80_HZ_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi1_HZ_all_withBG=(A_pos_polm80_phi1_signal_all_withBG+B_neg_polm80_phi1_signal_all_withBG)/(abs(A_pos_polm80_phi1_signal_all_withBG)+abs(B_neg_polm80_phi1_signal_all_withBG))
    alpha_polm80_phi1_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polm80_phi1_signal_all_withBG)*abs(B_neg_polm80_phi1_signal_all_withBG)/pow((abs(A_pos_polm80_phi1_signal_all_withBG)+abs(B_neg_polm80_phi1_signal_all_withBG)),3))

  
    A_pos_polm80_phi1_allBG=h_polm80_ee_qq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi1_allBG=h_polm80_ee_qq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi1_allBG=(A_pos_polm80_phi1_allBG+B_neg_polm80_phi1_allBG)/(abs(A_pos_polm80_phi1_allBG)+abs(B_neg_polm80_phi1_allBG))
    alpha_polm80_phi1_allBG_error=sqrt(4.*abs(A_pos_polm80_phi1_allBG)*abs(B_neg_polm80_phi1_allBG)/pow((abs(A_pos_polm80_phi1_allBG)+abs(B_neg_polm80_phi1_allBG)),3))

    print 'alpha_polm80_phi1_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_phi1_HZ,alpha_polm80_phi1_HZ_all,alpha_polm80_phi1_HZ_BDTm020,alpha_polm80_phi1_HZ_all_BDTm020,alpha_polm80_phi1_HZ_all_withBG
    print 'errors of alpha_polm80_phi1_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_phi1_HZ_error,alpha_polm80_phi1_HZ_all_error,alpha_polm80_phi1_HZ_BDTm020_error,alpha_polm80_phi1_HZ_all_BDTm020_error,alpha_polm80_phi1_HZ_all_withBG_error

    print 'alpha_polm80_phi1_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polm80_phi1_HZ,alpha_polm80_phi1_HZ_all,alpha_polm80_phi1_allBG,alpha_polm80_phi1_HZ_all_withBG,alpha_parton_polm80_phi1_HZ_all
    print 'errors of alpha_polm80_phi1_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polm80_phi1_HZ_error,alpha_polm80_phi1_HZ_all_error,alpha_polm80_phi1_allBG_error,alpha_polm80_phi1_HZ_all_withBG_error,alpha_parton_polm80_phi1_HZ_all_error


    h_polm80_HZ_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_.Get("0.375/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_.Get("0.375/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_HZ_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_HZ_AllEvents_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep")
    h_polm80_HZ_AllEvents_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep")

    h_polm80_HZ_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qqqq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    A_pos_polm80_phi2_signal=h_polm80_HZ_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi2_signal=h_polm80_HZ_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi2_HZ=(A_pos_polm80_phi2_signal+B_neg_polm80_phi2_signal)/(abs(A_pos_polm80_phi2_signal)+abs(B_neg_polm80_phi2_signal))
    alpha_polm80_phi2_HZ_error=sqrt(4.*abs(A_pos_polm80_phi2_signal)*abs(B_neg_polm80_phi2_signal)/pow((abs(A_pos_polm80_phi2_signal)+abs(B_neg_polm80_phi2_signal)),3))

    A_pos_polm80_phi2_signal_all=h_polm80_HZ_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi2_signal_all=h_polm80_HZ_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi2_HZ_all=(A_pos_polm80_phi2_signal_all+B_neg_polm80_phi2_signal_all)/(abs(A_pos_polm80_phi2_signal_all)+abs(B_neg_polm80_phi2_signal_all))
    alpha_polm80_phi2_HZ_all_error=sqrt(4.*abs(A_pos_polm80_phi2_signal_all)*abs(B_neg_polm80_phi2_signal_all)/pow((abs(A_pos_polm80_phi2_signal_all)+abs(B_neg_polm80_phi2_signal_all)),3))
    
    A_pos_parton_polm80_phi2_signal_all=h_polm80_HZ_AllEvents_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    B_neg_parton_polm80_phi2_signal_all=h_polm80_HZ_AllEvents_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    alpha_parton_polm80_phi2_HZ_all=(A_pos_parton_polm80_phi2_signal_all+B_neg_parton_polm80_phi2_signal_all)/(abs(A_pos_parton_polm80_phi2_signal_all)+abs(B_neg_parton_polm80_phi2_signal_all))
    alpha_parton_polm80_phi2_HZ_all_error=sqrt(4.*abs(A_pos_parton_polm80_phi2_signal_all)*abs(B_neg_parton_polm80_phi2_signal_all)/pow((abs(A_pos_parton_polm80_phi2_signal_all)+abs(B_neg_parton_polm80_phi2_signal_all)),3))


    A_pos_polm80_phi2_signal_BDTm020=h_polm80_HZ_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polm80_phi2_signal_BDTm020=h_polm80_HZ_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polm80_phi2_HZ_BDTm020=(A_pos_polm80_phi2_signal_BDTm020+B_neg_polm80_phi2_signal_BDTm020)/(abs(A_pos_polm80_phi2_signal_BDTm020)+abs(B_neg_polm80_phi2_signal_BDTm020))
    alpha_polm80_phi2_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polm80_phi2_signal_BDTm020)*abs(B_neg_polm80_phi2_signal_BDTm020)/pow((abs(A_pos_polm80_phi2_signal_BDTm020)+abs(B_neg_polm80_phi2_signal_BDTm020)),3))
    A_pos_polm80_phi2_signal_all_BDTm020=h_polm80_HZ_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polm80_phi2_signal_all_BDTm020=h_polm80_HZ_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polm80_phi2_HZ_all_BDTm020=(A_pos_polm80_phi2_signal_all_BDTm020+B_neg_polm80_phi2_signal_all_BDTm020)/(abs(A_pos_polm80_phi2_signal_all_BDTm020)+abs(B_neg_polm80_phi2_signal_all_BDTm020))
    alpha_polm80_phi2_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polm80_phi2_signal_all_BDTm020)*abs(B_neg_polm80_phi2_signal_all_BDTm020)/pow((abs(A_pos_polm80_phi2_signal_all_BDTm020)+abs(B_neg_polm80_phi2_signal_all_BDTm020)),3))

    A_pos_polm80_phi2_signal_all_withBG=h_polm80_HZ_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi2_signal_all_withBG=h_polm80_HZ_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi2_HZ_all_withBG=(A_pos_polm80_phi2_signal_all_withBG+B_neg_polm80_phi2_signal_all_withBG)/(abs(A_pos_polm80_phi2_signal_all_withBG)+abs(B_neg_polm80_phi2_signal_all_withBG))
    alpha_polm80_phi2_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polm80_phi2_signal_all_withBG)*abs(B_neg_polm80_phi2_signal_all_withBG)/pow((abs(A_pos_polm80_phi2_signal_all_withBG)+abs(B_neg_polm80_phi2_signal_all_withBG)),3))

    A_pos_polm80_phi2_allBG=h_polm80_ee_qq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi2_allBG=h_polm80_ee_qq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi2_allBG=(A_pos_polm80_phi2_allBG+B_neg_polm80_phi2_allBG)/(abs(A_pos_polm80_phi2_allBG)+abs(B_neg_polm80_phi2_allBG))
    alpha_polm80_phi2_allBG_error=sqrt(4.*abs(A_pos_polm80_phi2_allBG)*abs(B_neg_polm80_phi2_allBG)/pow((abs(A_pos_polm80_phi2_allBG)+abs(B_neg_polm80_phi2_allBG)),3))
  
    print 'alpha_polm80_phi2_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_phi2_HZ,alpha_polm80_phi2_HZ_all,alpha_polm80_phi2_HZ_BDTm020,alpha_polm80_phi2_HZ_all_BDTm020,alpha_polm80_phi2_HZ_all_withBG
    print 'errors of alpha_polm80_phi2_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_phi2_HZ_error,alpha_polm80_phi2_HZ_all_error,alpha_polm80_phi2_HZ_BDTm020_error,alpha_polm80_phi2_HZ_all_BDTm020_error,alpha_polm80_phi2_HZ_all_withBG_error

    print 'alpha_polm80_phi2_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polm80_phi2_HZ,alpha_polm80_phi2_HZ_all,alpha_polm80_phi2_allBG,alpha_polm80_phi2_HZ_all_withBG,alpha_parton_polm80_phi2_HZ_all
    print 'errors of alpha_polm80_phi2_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polm80_phi2_HZ_error,alpha_polm80_phi2_HZ_all_error,alpha_polm80_phi2_allBG_error,alpha_polm80_phi2_HZ_all_withBG_error,alpha_parton_polm80_phi2_HZ_all_error


    h_polm80_HZ_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_.Get("0.375/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_.Get("0.375/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_HZ_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_HZ_AllEvents_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep")
    h_polm80_HZ_AllEvents_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep")

    h_polm80_HZ_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qqqq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    A_pos_polm80_phi3_signal=h_polm80_HZ_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi3_signal=h_polm80_HZ_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi3_HZ=(A_pos_polm80_phi3_signal+B_neg_polm80_phi3_signal)/(abs(A_pos_polm80_phi3_signal)+abs(B_neg_polm80_phi3_signal))
    alpha_polm80_phi3_HZ_error=sqrt(4.*abs(A_pos_polm80_phi3_signal)*abs(B_neg_polm80_phi3_signal)/pow((abs(A_pos_polm80_phi3_signal)+abs(B_neg_polm80_phi3_signal)),3))

    A_pos_polm80_phi3_signal_all=h_polm80_HZ_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi3_signal_all=h_polm80_HZ_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi3_HZ_all=(A_pos_polm80_phi3_signal_all+B_neg_polm80_phi3_signal_all)/(abs(A_pos_polm80_phi3_signal_all)+abs(B_neg_polm80_phi3_signal_all))
    alpha_polm80_phi3_HZ_all_error=sqrt(4.*abs(A_pos_polm80_phi3_signal_all)*abs(B_neg_polm80_phi3_signal_all)/pow((abs(A_pos_polm80_phi3_signal_all)+abs(B_neg_polm80_phi3_signal_all)),3))
    
    A_pos_parton_polm80_phi3_signal_all=h_polm80_HZ_AllEvents_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    B_neg_parton_polm80_phi3_signal_all=h_polm80_HZ_AllEvents_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    alpha_parton_polm80_phi3_HZ_all=(A_pos_parton_polm80_phi3_signal_all+B_neg_parton_polm80_phi3_signal_all)/(abs(A_pos_parton_polm80_phi3_signal_all)+abs(B_neg_parton_polm80_phi3_signal_all))
    alpha_parton_polm80_phi3_HZ_all_error=sqrt(4.*abs(A_pos_parton_polm80_phi3_signal_all)*abs(B_neg_parton_polm80_phi3_signal_all)/pow((abs(A_pos_parton_polm80_phi3_signal_all)+abs(B_neg_parton_polm80_phi3_signal_all)),3))


    A_pos_polm80_phi3_signal_BDTm020=h_polm80_HZ_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polm80_phi3_signal_BDTm020=h_polm80_HZ_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polm80_phi3_HZ_BDTm020=(A_pos_polm80_phi3_signal_BDTm020+B_neg_polm80_phi3_signal_BDTm020)/(abs(A_pos_polm80_phi3_signal_BDTm020)+abs(B_neg_polm80_phi3_signal_BDTm020))
    alpha_polm80_phi3_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polm80_phi3_signal_BDTm020)*abs(B_neg_polm80_phi3_signal_BDTm020)/pow((abs(A_pos_polm80_phi3_signal_BDTm020)+abs(B_neg_polm80_phi3_signal_BDTm020)),3))
    A_pos_polm80_phi3_signal_all_BDTm020=h_polm80_HZ_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polm80_phi3_signal_all_BDTm020=h_polm80_HZ_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polm80_phi3_HZ_all_BDTm020=(A_pos_polm80_phi3_signal_all_BDTm020+B_neg_polm80_phi3_signal_all_BDTm020)/(abs(A_pos_polm80_phi3_signal_all_BDTm020)+abs(B_neg_polm80_phi3_signal_all_BDTm020))
    alpha_polm80_phi3_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polm80_phi3_signal_all_BDTm020)*abs(B_neg_polm80_phi3_signal_all_BDTm020)/pow((abs(A_pos_polm80_phi3_signal_all_BDTm020)+abs(B_neg_polm80_phi3_signal_all_BDTm020)),3))

    A_pos_polm80_phi3_signal_all_withBG=h_polm80_HZ_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi3_signal_all_withBG=h_polm80_HZ_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi3_HZ_all_withBG=(A_pos_polm80_phi3_signal_all_withBG+B_neg_polm80_phi3_signal_all_withBG)/(abs(A_pos_polm80_phi3_signal_all_withBG)+abs(B_neg_polm80_phi3_signal_all_withBG))
    alpha_polm80_phi3_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polm80_phi3_signal_all_withBG)*abs(B_neg_polm80_phi3_signal_all_withBG)/pow((abs(A_pos_polm80_phi3_signal_all_withBG)+abs(B_neg_polm80_phi3_signal_all_withBG)),3))

    A_pos_polm80_phi3_allBG=h_polm80_ee_qq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi3_allBG=h_polm80_ee_qq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi3_allBG=(A_pos_polm80_phi3_allBG+B_neg_polm80_phi3_allBG)/(abs(A_pos_polm80_phi3_allBG)+abs(B_neg_polm80_phi3_allBG))
    alpha_polm80_phi3_allBG_error=sqrt(4.*abs(A_pos_polm80_phi3_allBG)*abs(B_neg_polm80_phi3_allBG)/pow((abs(A_pos_polm80_phi3_allBG)+abs(B_neg_polm80_phi3_allBG)),3))
  
    print 'alpha_polm80_phi3_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_phi3_HZ,alpha_polm80_phi3_HZ_all,alpha_polm80_phi3_HZ_BDTm020,alpha_polm80_phi3_HZ_all_BDTm020,alpha_polm80_phi3_HZ_all_withBG
    print 'errors of alpha_polm80_phi3_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_phi3_HZ_error,alpha_polm80_phi3_HZ_all_error,alpha_polm80_phi3_HZ_BDTm020_error,alpha_polm80_phi3_HZ_all_BDTm020_error,alpha_polm80_phi3_HZ_all_withBG_error


    print 'alpha_polm80_phi3_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polm80_phi3_HZ,alpha_polm80_phi3_HZ_all,alpha_polm80_phi3_allBG,alpha_polm80_phi3_HZ_all_withBG,alpha_parton_polm80_phi3_HZ_all
    print 'errors of alpha_polm80_phi3_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polm80_phi3_HZ_error,alpha_polm80_phi3_HZ_all_error,alpha_polm80_phi3_allBG_error,alpha_polm80_phi3_HZ_all_withBG_error,alpha_parton_polm80_phi3_HZ_all_error




    h_polm80_HZ_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_.Get("0.375/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_.Get("0.375/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_HZ_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_.Get("-0.2/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_HZ_AllEvents_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep")
    h_polm80_HZ_AllEvents_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep=file_polm80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep")

    h_polm80_HZ_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_HZ_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qqqq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    A_pos_polm80_phi4_signal=h_polm80_HZ_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi4_signal=h_polm80_HZ_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi4_HZ=(A_pos_polm80_phi4_signal+B_neg_polm80_phi4_signal)/(abs(A_pos_polm80_phi4_signal)+abs(B_neg_polm80_phi4_signal))
    alpha_polm80_phi4_HZ_error=sqrt(4.*abs(A_pos_polm80_phi4_signal)*abs(B_neg_polm80_phi4_signal)/pow((abs(A_pos_polm80_phi4_signal)+abs(B_neg_polm80_phi4_signal)),3))

    A_pos_polm80_phi4_signal_all=h_polm80_HZ_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi4_signal_all=h_polm80_HZ_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi4_HZ_all=(A_pos_polm80_phi4_signal_all+B_neg_polm80_phi4_signal_all)/(abs(A_pos_polm80_phi4_signal_all)+abs(B_neg_polm80_phi4_signal_all))
    alpha_polm80_phi4_HZ_all_error=sqrt(4.*abs(A_pos_polm80_phi4_signal_all)*abs(B_neg_polm80_phi4_signal_all)/pow((abs(A_pos_polm80_phi4_signal_all)+abs(B_neg_polm80_phi4_signal_all)),3))
    
    A_pos_parton_polm80_phi4_signal_all=h_polm80_HZ_AllEvents_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    B_neg_parton_polm80_phi4_signal_all=h_polm80_HZ_AllEvents_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    alpha_parton_polm80_phi4_HZ_all=(A_pos_parton_polm80_phi4_signal_all+B_neg_parton_polm80_phi4_signal_all)/(abs(A_pos_parton_polm80_phi4_signal_all)+abs(B_neg_parton_polm80_phi4_signal_all))
    alpha_parton_polm80_phi4_HZ_all_error=sqrt(4.*abs(A_pos_parton_polm80_phi4_signal_all)*abs(B_neg_parton_polm80_phi4_signal_all)/pow((abs(A_pos_parton_polm80_phi4_signal_all)+abs(B_neg_parton_polm80_phi4_signal_all)),3))


    A_pos_polm80_phi4_signal_BDTm020=h_polm80_HZ_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polm80_phi4_signal_BDTm020=h_polm80_HZ_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polm80_phi4_HZ_BDTm020=(A_pos_polm80_phi4_signal_BDTm020+B_neg_polm80_phi4_signal_BDTm020)/(abs(A_pos_polm80_phi4_signal_BDTm020)+abs(B_neg_polm80_phi4_signal_BDTm020))
    alpha_polm80_phi4_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polm80_phi4_signal_BDTm020)*abs(B_neg_polm80_phi4_signal_BDTm020)/pow((abs(A_pos_polm80_phi4_signal_BDTm020)+abs(B_neg_polm80_phi4_signal_BDTm020)),3))
    A_pos_polm80_phi4_signal_all_BDTm020=h_polm80_HZ_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polm80_phi4_signal_all_BDTm020=h_polm80_HZ_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polm80_phi4_HZ_all_BDTm020=(A_pos_polm80_phi4_signal_all_BDTm020+B_neg_polm80_phi4_signal_all_BDTm020)/(abs(A_pos_polm80_phi4_signal_all_BDTm020)+abs(B_neg_polm80_phi4_signal_all_BDTm020))
    alpha_polm80_phi4_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polm80_phi4_signal_all_BDTm020)*abs(B_neg_polm80_phi4_signal_all_BDTm020)/pow((abs(A_pos_polm80_phi4_signal_all_BDTm020)+abs(B_neg_polm80_phi4_signal_all_BDTm020)),3))

    A_pos_polm80_phi4_signal_all_withBG=h_polm80_HZ_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi4_signal_all_withBG=h_polm80_HZ_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi4_HZ_all_withBG=(A_pos_polm80_phi4_signal_all_withBG+B_neg_polm80_phi4_signal_all_withBG)/(abs(A_pos_polm80_phi4_signal_all_withBG)+abs(B_neg_polm80_phi4_signal_all_withBG))
    alpha_polm80_phi4_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polm80_phi4_signal_all_withBG)*abs(B_neg_polm80_phi4_signal_all_withBG)/pow((abs(A_pos_polm80_phi4_signal_all_withBG)+abs(B_neg_polm80_phi4_signal_all_withBG)),3))

    A_pos_polm80_phi4_allBG=h_polm80_ee_qq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polm80_phi4_allBG=h_polm80_ee_qq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polm80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polm80_phi4_allBG=(A_pos_polm80_phi4_allBG+B_neg_polm80_phi4_allBG)/(abs(A_pos_polm80_phi4_allBG)+abs(B_neg_polm80_phi4_allBG))
    alpha_polm80_phi4_allBG_error=sqrt(4.*abs(A_pos_polm80_phi4_allBG)*abs(B_neg_polm80_phi4_allBG)/pow((abs(A_pos_polm80_phi4_allBG)+abs(B_neg_polm80_phi4_allBG)),3))
  
    print 'alpha_polm80_phi4_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_phi4_HZ,alpha_polm80_phi4_HZ_all,alpha_polm80_phi4_HZ_BDTm020,alpha_polm80_phi4_HZ_all_BDTm020,alpha_polm80_phi4_HZ_all_withBG
    print 'errors of alpha_polm80_phi4_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polm80_phi4_HZ_error,alpha_polm80_phi4_HZ_all_error,alpha_polm80_phi4_HZ_BDTm020_error,alpha_polm80_phi4_HZ_all_BDTm020_error,alpha_polm80_phi4_HZ_all_withBG_error

    print 'alpha_polm80_phi4_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polm80_phi4_HZ,alpha_polm80_phi4_HZ_all,alpha_polm80_phi4_allBG,alpha_polm80_phi4_HZ_all_withBG,alpha_parton_polm80_phi4_HZ_all
    print 'errors of alpha_polm80_phi4_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polm80_phi4_HZ_error,alpha_polm80_phi4_HZ_all_error,alpha_polm80_phi4_allBG_error,alpha_polm80_phi4_HZ_all_withBG_error,alpha_parton_polm80_phi4_HZ_all_error



    h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom")
    h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.Rebin(4)
    h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.SetLineWidth(2)
    h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.SetLineColor(1)
    h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.GetXaxis().SetTitle('cos(#theta_{1}(reco))*sgn(cos(#theta_{1}(part))')
    h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.GetYaxis().SetTitle('Events')

    h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom")
    h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.Rebin(4)
    h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.SetLineWidth(2)
    h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.SetLineColor(kRed)
    h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.GetXaxis().SetTitle('cos(#theta2(reco))*sgn(cos(#theta(part))')
    h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.GetYaxis().SetTitle('Events')
 
    h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom")
    h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.Rebin(4)
    h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.SetLineWidth(2)
    h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.SetLineColor(kBlue)
    h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.GetXaxis().SetTitle('cos(#theta2(reco))*sgn(cos(#theta(part))')
    h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.GetYaxis().SetTitle('Events')

    h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom")
    h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.Rebin(4)
    h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.SetLineWidth(2)
    h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.SetLineColor(kGreen-2)
    h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.GetXaxis().SetTitle('cos(#theta2(reco))*sgn(cos(#theta(part))')
    h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.GetYaxis().SetTitle('Events')

    h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom")
    h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.Rebin(4)
    h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.SetLineWidth(2)
    h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.SetLineColor(kCyan+1)
    h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.GetXaxis().SetTitle('cos(#theta2(reco))*sgn(cos(#theta(part))')
    h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.GetYaxis().SetTitle('Events')
    
    canvas_cosAngle_sgnCosPart_methods_quark_Charge_rel_polm80 = setUpperCanvas("canvas_cosAngle_sgnCosPart_methods_jcE_rel_polm80");
    canvas_cosAngle_sgnCosPart_methods_quark_Charge_rel_polm80.cd()

    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test=TLegend(0.20,0.61,0.60,0.87);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetBorderSize(0);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetTextAlign(12);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetTextSize(0.050);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetTextFont(42);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetMargin(0.15);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetLineColor(1);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetLineStyle(1);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetLineWidth(1);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetFillColor(0);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.SetFillStyle(0);
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.AddEntry(h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.DrawCopy("h"),"method1");
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.AddEntry(h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.DrawCopy("h,same"),"method2");
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.AddEntry(h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.DrawCopy("h,same"),"method3");
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.AddEntry(h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.DrawCopy("h,same"),"method4");
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.AddEntry(h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom_polm80.DrawCopy("h,same"),"method5");
    leg_h_cosAngle_sgnCosPart_methods_jcE_rel_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_cosAngle_sgnCosPart_methods_quark_Charge_rel_polm80.Print("hzqq_bbar_cosTheta1_reco_times_sgnCosTheta1_part_methods_quark_Charge.eps")

    h_jetE_reco_over_parton_H_matched_orig_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_jetE_reco_over_parton_H_matched_orig")
    h_jetE_reco_over_parton_H_matched_orig_polm80.SetLineWidth(2)
    h_jetE_reco_over_parton_H_matched_orig_polm80.SetLineColor(2)
    h_jetE_reco_over_parton_H_matched_orig_polm80.GetXaxis().SetTitle('H-jet E/H-parton E')
    h_jetE_reco_over_parton_H_matched_orig_polm80.GetYaxis().SetTitle('Events')
 
    h_jetE_reco_over_parton_H_matched_corr_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_jetE_reco_over_parton_H_matched_corr")
    h_jetE_reco_over_parton_H_matched_corr_polm80.SetLineWidth(2)
    h_jetE_reco_over_parton_H_matched_corr_polm80.SetLineColor(4)
    h_jetE_reco_over_parton_H_matched_corr_polm80.GetXaxis().SetTitle('H-jet E/H-parton E')
    h_jetE_reco_over_parton_H_matched_corr_polm80.GetYaxis().SetTitle('Events')

    canvas_H_matched_jetE_rel_polm80 = setUpperCanvas("canvas_H_matched_jetE_rel_polm80");
    canvas_H_matched_jetE_rel_polm80.cd()
    h_jetE_reco_over_parton_H_matched_corr_polm80.SetMinimum(0)

    leg_h_H_matched_jetE_rel_polm80_test=TLegend(0.20,0.72,0.60,0.87);
    leg_h_H_matched_jetE_rel_polm80_test.SetBorderSize(0);
    leg_h_H_matched_jetE_rel_polm80_test.SetTextAlign(12);
    leg_h_H_matched_jetE_rel_polm80_test.SetTextSize(0.050);
    leg_h_H_matched_jetE_rel_polm80_test.SetTextFont(42);
    leg_h_H_matched_jetE_rel_polm80_test.SetMargin(0.15);
    leg_h_H_matched_jetE_rel_polm80_test.SetLineColor(1);
    leg_h_H_matched_jetE_rel_polm80_test.SetLineStyle(1);
    leg_h_H_matched_jetE_rel_polm80_test.SetLineWidth(1);
    leg_h_H_matched_jetE_rel_polm80_test.SetFillColor(0);
    leg_h_H_matched_jetE_rel_polm80_test.SetFillStyle(0);
    leg_h_H_matched_jetE_rel_polm80_test.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_H_matched_jetE_rel_polm80_test.AddEntry(h_jetE_reco_over_parton_H_matched_corr_polm80.DrawCopy("h,e"),"corrected jet");
    leg_h_H_matched_jetE_rel_polm80_test.AddEntry(h_jetE_reco_over_parton_H_matched_orig_polm80.DrawCopy("h,e,same"),"original jet");
    leg_h_H_matched_jetE_rel_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_H_matched_jetE_rel_polm80.Print("hzqq_bbar_polm80_H_jet_E_corr_vs_E_orig.eps")


    h_jetP_reco_over_parton_H_matched_orig_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_jetP_reco_over_parton_H_matched_orig")
    h_jetP_reco_over_parton_H_matched_orig_polm80.SetLineWidth(2)
    h_jetP_reco_over_parton_H_matched_orig_polm80.SetLineColor(2)
    h_jetP_reco_over_parton_H_matched_orig_polm80.GetXaxis().SetTitle('H-jet P/H-parton P')
    h_jetP_reco_over_parton_H_matched_orig_polm80.GetYaxis().SetTitle('Events')
 
    h_jetP_reco_over_parton_H_matched_corr_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_jetP_reco_over_parton_H_matched_corr")
    h_jetP_reco_over_parton_H_matched_corr_polm80.SetLineWidth(2)
    h_jetP_reco_over_parton_H_matched_corr_polm80.SetLineColor(4)
    h_jetP_reco_over_parton_H_matched_corr_polm80.GetXaxis().SetTitle('H-jet P/H-parton P')
    h_jetP_reco_over_parton_H_matched_corr_polm80.GetYaxis().SetTitle('Events')
    h_jetP_reco_over_parton_H_matched_corr_polm80.SetMinimum(0)

    canvas_H_matched_jetP_rel_polm80 = setUpperCanvas("canvas_H_matched_jetP_rel_polm80");
    canvas_H_matched_jetP_rel_polm80.cd()

    leg_h_H_matched_jetP_rel_polm80_test=TLegend(0.20,0.72,0.60,0.87);
    leg_h_H_matched_jetP_rel_polm80_test.SetBorderSize(0);
    leg_h_H_matched_jetP_rel_polm80_test.SetTextAlign(12);
    leg_h_H_matched_jetP_rel_polm80_test.SetTextSize(0.050);
    leg_h_H_matched_jetP_rel_polm80_test.SetTextFont(42);
    leg_h_H_matched_jetP_rel_polm80_test.SetMargin(0.15);
    leg_h_H_matched_jetP_rel_polm80_test.SetLineColor(1);
    leg_h_H_matched_jetP_rel_polm80_test.SetLineStyle(1);
    leg_h_H_matched_jetP_rel_polm80_test.SetLineWidth(1);
    leg_h_H_matched_jetP_rel_polm80_test.SetFillColor(0);
    leg_h_H_matched_jetP_rel_polm80_test.SetFillStyle(0);
    leg_h_H_matched_jetP_rel_polm80_test.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_H_matched_jetP_rel_polm80_test.AddEntry(h_jetP_reco_over_parton_H_matched_corr_polm80.DrawCopy("h,e"),"corrected jet");
    leg_h_H_matched_jetP_rel_polm80_test.AddEntry(h_jetP_reco_over_parton_H_matched_orig_polm80.DrawCopy("h,e,same"),"original jet");
    leg_h_H_matched_jetP_rel_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_H_matched_jetP_rel_polm80.Print("hzqq_bbar_polm80_H_jet_P_corr_vs_P_orig.eps")

    h_jetE_reco_over_parton_Z_matched_orig_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_jetE_reco_over_parton_Z_matched_orig")
    h_jetE_reco_over_parton_Z_matched_orig_polm80.SetLineWidth(2)
    h_jetE_reco_over_parton_Z_matched_orig_polm80.SetLineColor(2)
    h_jetE_reco_over_parton_Z_matched_orig_polm80.GetXaxis().SetTitle('Z-jet E/Z-parton E')
    h_jetE_reco_over_parton_Z_matched_orig_polm80.GetYaxis().SetTitle('Events')
 
    h_jetE_reco_over_parton_Z_matched_corr_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_jetE_reco_over_parton_Z_matched_corr")
    h_jetE_reco_over_parton_Z_matched_corr_polm80.SetLineWidth(2)
    h_jetE_reco_over_parton_Z_matched_corr_polm80.SetLineColor(4)
    h_jetE_reco_over_parton_Z_matched_corr_polm80.GetXaxis().SetTitle('Z-jet E/Z-parton E')
    h_jetE_reco_over_parton_Z_matched_corr_polm80.GetYaxis().SetTitle('Events')
    h_jetE_reco_over_parton_Z_matched_corr_polm80.SetMinimum(0)

    canvas_Z_matched_jetE_rel_polm80 = setUpperCanvas("canvas_Z_matched_jetE_rel_polm80");
    canvas_Z_matched_jetE_rel_polm80.cd()

    leg_h_Z_matched_jetE_rel_polm80_test=TLegend(0.20,0.72,0.60,0.87);
    leg_h_Z_matched_jetE_rel_polm80_test.SetBorderSize(0);
    leg_h_Z_matched_jetE_rel_polm80_test.SetTextAlign(12);
    leg_h_Z_matched_jetE_rel_polm80_test.SetTextSize(0.050);
    leg_h_Z_matched_jetE_rel_polm80_test.SetTextFont(42);
    leg_h_Z_matched_jetE_rel_polm80_test.SetMargin(0.15);
    leg_h_Z_matched_jetE_rel_polm80_test.SetLineColor(1);
    leg_h_Z_matched_jetE_rel_polm80_test.SetLineStyle(1);
    leg_h_Z_matched_jetE_rel_polm80_test.SetLineWidth(1);
    leg_h_Z_matched_jetE_rel_polm80_test.SetFillColor(0);
    leg_h_Z_matched_jetE_rel_polm80_test.SetFillStyle(0);
    leg_h_Z_matched_jetE_rel_polm80_test.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_Z_matched_jetE_rel_polm80_test.AddEntry(h_jetE_reco_over_parton_Z_matched_corr_polm80.DrawCopy("h,e"),"corrected jet");
    leg_h_Z_matched_jetE_rel_polm80_test.AddEntry(h_jetE_reco_over_parton_Z_matched_orig_polm80.DrawCopy("h,e,same"),"original jet");
    leg_h_Z_matched_jetE_rel_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_Z_matched_jetE_rel_polm80.Print("hzqq_bbar_polm80_Z_jet_E_corr_vs_E_orig.eps")


    h_jetP_reco_over_parton_Z_matched_orig_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_jetP_reco_over_parton_Z_matched_orig")
    h_jetP_reco_over_parton_Z_matched_orig_polm80.SetLineWidth(2)
    h_jetP_reco_over_parton_Z_matched_orig_polm80.SetLineColor(2)
    h_jetP_reco_over_parton_Z_matched_orig_polm80.GetXaxis().SetTitle('Z-jet P/Z-parton P')
    h_jetP_reco_over_parton_Z_matched_orig_polm80.GetYaxis().SetTitle('Events')
 
    h_jetP_reco_over_parton_Z_matched_corr_polm80=file_polm80_HZ_SignalHistos_.Get("-0.2/h_jetP_reco_over_parton_Z_matched_corr")
    h_jetP_reco_over_parton_Z_matched_corr_polm80.SetLineWidth(2)
    h_jetP_reco_over_parton_Z_matched_corr_polm80.SetLineColor(4)
    h_jetP_reco_over_parton_Z_matched_corr_polm80.GetXaxis().SetTitle('Z-jet P/Z-parton P')
    h_jetP_reco_over_parton_Z_matched_corr_polm80.GetYaxis().SetTitle('Events')
    h_jetP_reco_over_parton_Z_matched_corr_polm80.SetMinimum(0)

    canvas_Z_matched_jetP_rel_polm80 = setUpperCanvas("canvas_Z_matched_jetP_rel_polm80");
    canvas_Z_matched_jetP_rel_polm80.cd()

    leg_h_Z_matched_jetP_rel_polm80_test=TLegend(0.20,0.72,0.60,0.87);
    leg_h_Z_matched_jetP_rel_polm80_test.SetBorderSize(0);
    leg_h_Z_matched_jetP_rel_polm80_test.SetTextAlign(12);
    leg_h_Z_matched_jetP_rel_polm80_test.SetTextSize(0.050);
    leg_h_Z_matched_jetP_rel_polm80_test.SetTextFont(42);
    leg_h_Z_matched_jetP_rel_polm80_test.SetMargin(0.15);
    leg_h_Z_matched_jetP_rel_polm80_test.SetLineColor(1);
    leg_h_Z_matched_jetP_rel_polm80_test.SetLineStyle(1);
    leg_h_Z_matched_jetP_rel_polm80_test.SetLineWidth(1);
    leg_h_Z_matched_jetP_rel_polm80_test.SetFillColor(0);
    leg_h_Z_matched_jetP_rel_polm80_test.SetFillStyle(0);
    leg_h_Z_matched_jetP_rel_polm80_test.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_Z_matched_jetP_rel_polm80_test.AddEntry(h_jetP_reco_over_parton_Z_matched_corr_polm80.DrawCopy("h,e"),"corrected jet");
    leg_h_Z_matched_jetP_rel_polm80_test.AddEntry(h_jetP_reco_over_parton_Z_matched_orig_polm80.DrawCopy("h,e,same"),"original jet");
    leg_h_Z_matched_jetP_rel_polm80_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_Z_matched_jetP_rel_polm80.Print("hzqq_bbar_polm80_Z_jet_P_corr_vs_P_orig.eps")





    h_jetE_reco_over_parton_H_matched_orig_polm80_AllEvents=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_jetE_reco_over_parton_H_matched_orig")
    h_jetE_reco_over_parton_H_matched_orig_polm80_AllEvents.SetLineWidth(2)
    h_jetE_reco_over_parton_H_matched_orig_polm80_AllEvents.SetLineColor(2)
    h_jetE_reco_over_parton_H_matched_orig_polm80_AllEvents.GetXaxis().SetTitle('H-jet E/H-parton E')
    h_jetE_reco_over_parton_H_matched_orig_polm80_AllEvents.GetYaxis().SetTitle('Events')
 
    h_jetE_reco_over_parton_H_matched_corr_polm80_AllEvents=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_jetE_reco_over_parton_H_matched_corr")
    h_jetE_reco_over_parton_H_matched_corr_polm80_AllEvents.SetLineWidth(2)
    h_jetE_reco_over_parton_H_matched_corr_polm80_AllEvents.SetLineColor(4)
    h_jetE_reco_over_parton_H_matched_corr_polm80_AllEvents.GetXaxis().SetTitle('H-jet E/H-parton E')
    h_jetE_reco_over_parton_H_matched_corr_polm80_AllEvents.GetYaxis().SetTitle('Events')

    canvas_H_matched_jetE_rel_polm80_AllEvents = setUpperCanvas("canvas_H_matched_jetE_rel_polm80_AllEvents");
    canvas_H_matched_jetE_rel_polm80_AllEvents.cd()
    h_jetE_reco_over_parton_H_matched_corr_polm80_AllEvents.SetMinimum(0)

    leg_h_H_matched_jetE_rel_polm80_AllEvents_test=TLegend(0.20,0.72,0.60,0.87);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetBorderSize(0);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetTextAlign(12);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetTextSize(0.050);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetTextFont(42);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetMargin(0.15);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetLineColor(1);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetLineStyle(1);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetLineWidth(1);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetFillColor(0);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetFillStyle(0);
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.AddEntry(h_jetE_reco_over_parton_H_matched_corr_polm80_AllEvents.DrawCopy("h,e"),"corrected jet");
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.AddEntry(h_jetE_reco_over_parton_H_matched_orig_polm80_AllEvents.DrawCopy("h,e,same"),"original jet");
    leg_h_H_matched_jetE_rel_polm80_AllEvents_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_H_matched_jetE_rel_polm80_AllEvents.Print("hzqq_AllEvents_polm80_H_jet_E_corr_vs_E_orig.eps")


    h_jetP_reco_over_parton_H_matched_orig_polm80_AllEvents=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_jetP_reco_over_parton_H_matched_orig")
    h_jetP_reco_over_parton_H_matched_orig_polm80_AllEvents.SetLineWidth(2)
    h_jetP_reco_over_parton_H_matched_orig_polm80_AllEvents.SetLineColor(2)
    h_jetP_reco_over_parton_H_matched_orig_polm80_AllEvents.GetXaxis().SetTitle('H-jet P/H-parton P')
    h_jetP_reco_over_parton_H_matched_orig_polm80_AllEvents.GetYaxis().SetTitle('Events')
 
    h_jetP_reco_over_parton_H_matched_corr_polm80_AllEvents=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_jetP_reco_over_parton_H_matched_corr")
    h_jetP_reco_over_parton_H_matched_corr_polm80_AllEvents.SetLineWidth(2)
    h_jetP_reco_over_parton_H_matched_corr_polm80_AllEvents.SetLineColor(4)
    h_jetP_reco_over_parton_H_matched_corr_polm80_AllEvents.GetXaxis().SetTitle('H-jet P/H-parton P')
    h_jetP_reco_over_parton_H_matched_corr_polm80_AllEvents.GetYaxis().SetTitle('Events')
    h_jetP_reco_over_parton_H_matched_corr_polm80_AllEvents.SetMinimum(0)

    canvas_H_matched_jetP_rel_polm80_AllEvents = setUpperCanvas("canvas_H_matched_jetP_rel_polm80_AllEvents");
    canvas_H_matched_jetP_rel_polm80_AllEvents.cd()

    leg_h_H_matched_jetP_rel_polm80_AllEvents_test=TLegend(0.20,0.72,0.60,0.87);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetBorderSize(0);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetTextAlign(12);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetTextSize(0.050);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetTextFont(42);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetMargin(0.15);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetLineColor(1);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetLineStyle(1);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetLineWidth(1);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetFillColor(0);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetFillStyle(0);
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.AddEntry(h_jetP_reco_over_parton_H_matched_corr_polm80_AllEvents.DrawCopy("h,e"),"corrected jet");
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.AddEntry(h_jetP_reco_over_parton_H_matched_orig_polm80_AllEvents.DrawCopy("h,e,same"),"original jet");
    leg_h_H_matched_jetP_rel_polm80_AllEvents_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_H_matched_jetP_rel_polm80_AllEvents.Print("hzqq_AllEvents_polm80_H_jet_P_corr_vs_P_orig.eps")

    h_jetE_reco_over_parton_Z_matched_orig_polm80_AllEvents=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_jetE_reco_over_parton_Z_matched_orig")
    h_jetE_reco_over_parton_Z_matched_orig_polm80_AllEvents.SetLineWidth(2)
    h_jetE_reco_over_parton_Z_matched_orig_polm80_AllEvents.SetLineColor(2)
    h_jetE_reco_over_parton_Z_matched_orig_polm80_AllEvents.GetXaxis().SetTitle('Z-jet E/Z-parton E')
    h_jetE_reco_over_parton_Z_matched_orig_polm80_AllEvents.GetYaxis().SetTitle('Events')
 
    h_jetE_reco_over_parton_Z_matched_corr_polm80_AllEvents=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_jetE_reco_over_parton_Z_matched_corr")
    h_jetE_reco_over_parton_Z_matched_corr_polm80_AllEvents.SetLineWidth(2)
    h_jetE_reco_over_parton_Z_matched_corr_polm80_AllEvents.SetLineColor(4)
    h_jetE_reco_over_parton_Z_matched_corr_polm80_AllEvents.GetXaxis().SetTitle('Z-jet E/Z-parton E')
    h_jetE_reco_over_parton_Z_matched_corr_polm80_AllEvents.GetYaxis().SetTitle('Events')
    h_jetE_reco_over_parton_Z_matched_corr_polm80_AllEvents.SetMinimum(0)

    canvas_Z_matched_jetE_rel_polm80_AllEvents = setUpperCanvas("canvas_Z_matched_jetE_rel_polm80_AllEvents");
    canvas_Z_matched_jetE_rel_polm80_AllEvents.cd()

    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test=TLegend(0.20,0.72,0.60,0.87);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetBorderSize(0);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetTextAlign(12);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetTextSize(0.050);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetTextFont(42);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetMargin(0.15);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetLineColor(1);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetLineStyle(1);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetLineWidth(1);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetFillColor(0);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetFillStyle(0);
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.AddEntry(h_jetE_reco_over_parton_Z_matched_corr_polm80_AllEvents.DrawCopy("h,e"),"corrected jet");
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.AddEntry(h_jetE_reco_over_parton_Z_matched_orig_polm80_AllEvents.DrawCopy("h,e,same"),"original jet");
    leg_h_Z_matched_jetE_rel_polm80_AllEvents_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_Z_matched_jetE_rel_polm80_AllEvents.Print("hzqq_AllEvents_polm80_Z_jet_E_corr_vs_E_orig.eps")


    h_jetP_reco_over_parton_Z_matched_orig_polm80_AllEvents=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_jetP_reco_over_parton_Z_matched_orig")
    h_jetP_reco_over_parton_Z_matched_orig_polm80_AllEvents.SetLineWidth(2)
    h_jetP_reco_over_parton_Z_matched_orig_polm80_AllEvents.SetLineColor(2)
    h_jetP_reco_over_parton_Z_matched_orig_polm80_AllEvents.GetXaxis().SetTitle('Z-jet P/Z-parton P')
    h_jetP_reco_over_parton_Z_matched_orig_polm80_AllEvents.GetYaxis().SetTitle('Events')
 
    h_jetP_reco_over_parton_Z_matched_corr_polm80_AllEvents=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_jetP_reco_over_parton_Z_matched_corr")
    h_jetP_reco_over_parton_Z_matched_corr_polm80_AllEvents.SetLineWidth(2)
    h_jetP_reco_over_parton_Z_matched_corr_polm80_AllEvents.SetLineColor(4)
    h_jetP_reco_over_parton_Z_matched_corr_polm80_AllEvents.GetXaxis().SetTitle('Z-jet P/Z-parton P')
    h_jetP_reco_over_parton_Z_matched_corr_polm80_AllEvents.GetYaxis().SetTitle('Events')
    h_jetP_reco_over_parton_Z_matched_corr_polm80_AllEvents.SetMinimum(0)

    canvas_Z_matched_jetP_rel_polm80_AllEvents = setUpperCanvas("canvas_Z_matched_jetP_rel_polm80_AllEvents");
    canvas_Z_matched_jetP_rel_polm80_AllEvents.cd()

    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test=TLegend(0.20,0.72,0.60,0.87);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetBorderSize(0);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetTextAlign(12);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetTextSize(0.050);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetTextFont(42);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetMargin(0.15);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetLineColor(1);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetLineStyle(1);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetLineWidth(1);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetFillColor(0);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetFillStyle(0);
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.AddEntry(h_jetP_reco_over_parton_Z_matched_corr_polm80_AllEvents.DrawCopy("h,e"),"corrected jet");
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.AddEntry(h_jetP_reco_over_parton_Z_matched_orig_polm80_AllEvents.DrawCopy("h,e,same"),"original jet");
    leg_h_Z_matched_jetP_rel_polm80_AllEvents_test.Draw();

    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_Z_matched_jetP_rel_polm80_AllEvents.Print("hzqq_AllEvents_polm80_Z_jet_P_corr_vs_P_orig.eps")






    x_polm80_BDTScore = array( 'f' )
    y_polm80_significance = array( 'f' )
    y_polm80_purity = array( 'f' )
    y_polm80_efficiency = array( 'f' )
    norm_polm=file_polm80_HZ_SignalHistos_.Get("-0.2/h_BDT_output").Integral()

    y_polm80_AllEvents_significance = array( 'f' )
    y_polm80_AllEvents_purity = array( 'f' )
    y_polm80_AllEvents_efficiency = array( 'f' )
    norm_polm_AllEvents=file_polm80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_BDT_output").Integral()
    #norm_polm=h_mass_sig_hzqq_norm_polm.Integral()
    #norm_polm = 10.
    print 'norm_polm',norm_polm

    for dir_ind in directory:
        x_polm80_BDTScore.append(float(dir_ind))
        h_mass_polm80_sig_hzqq=file_polm80_HZ_SignalHistos_.Get(dir_ind+"/h_jet1_mass")
        h_mass_polm80_sig_hzqq_AllEvents=file_polm80_HZ_SignalHistos_AllEvents_.Get(dir_ind+"/h_jet1_mass")
        h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get(dir_ind+"/h_jet1_mass")
        h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get(dir_ind+"/h_jet1_mass")
        h_mass_polm80_sig_ee_qqqqqq_BGHistos=file_polm80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_jet1_mass") 
        #print dir_ind,"integral of signal mass",h_mass_polm80_sig_hzqq.Integral(),h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_mass_polm80_sig_hzqq.Integral()>0:
            #print '#',dir_ind,"polm significance,purity ",h_mass_polm80_sig_hzqq.Integral()/sqrt(h_mass_polm80_sig_hzqq.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polm80_sig_hzqq.Integral()/(h_mass_polm80_sig_hzqq.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()),h_mass_polm80_sig_hzqq.Integral(),h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polm80_efficiency.append(100.*h_mass_polm80_sig_hzqq.Integral()/norm_polm)
            y_polm80_purity.append(100.*h_mass_polm80_sig_hzqq.Integral()/(h_mass_polm80_sig_hzqq.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()))
            y_polm80_significance.append(h_mass_polm80_sig_hzqq.Integral()/sqrt(h_mass_polm80_sig_hzqq.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()))
            print 'at polm ',dir_ind,'signif/pur/eff/events',h_mass_polm80_sig_hzqq.Integral()/sqrt(h_mass_polm80_sig_hzqq.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polm80_sig_hzqq.Integral()/(h_mass_polm80_sig_hzqq.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polm80_sig_hzqq.Integral()/norm_polm,h_mass_polm80_sig_hzqq.Integral()
        else:
            y_polm80_efficiency.append(0)
            y_polm80_purity.append(0)
            y_polm80_significance.append(0)
      #print dir_ind,"integral of signal mass",h_mass_polm80_sig_hzqq.Integral(),h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_mass_polm80_sig_hzqq_AllEvents.Integral()>0:
            #print '#',dir_ind,"polm significance,purity ",h_mass_polm80_sig_hzqq.Integral()/sqrt(h_mass_polm80_sig_hzqq.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polm80_sig_hzqq.Integral()/(h_mass_polm80_sig_hzqq.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()),h_mass_polm80_sig_hzqq.Integral(),h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polm80_AllEvents_efficiency.append(100.*h_mass_polm80_sig_hzqq_AllEvents.Integral()/norm_polm_AllEvents)
            y_polm80_AllEvents_purity.append(100.*h_mass_polm80_sig_hzqq_AllEvents.Integral()/(h_mass_polm80_sig_hzqq_AllEvents.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()))
            y_polm80_AllEvents_significance.append(h_mass_polm80_sig_hzqq_AllEvents.Integral()/sqrt(h_mass_polm80_sig_hzqq_AllEvents.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()))
            print 'at polm ',dir_ind,'signif/pur/eff/events all ',h_mass_polm80_sig_hzqq_AllEvents.Integral()/sqrt(h_mass_polm80_sig_hzqq_AllEvents.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polm80_sig_hzqq_AllEvents.Integral()/(h_mass_polm80_sig_hzqq_AllEvents.Integral()+h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polm80_sig_hzqq_AllEvents.Integral()/norm_polm,h_mass_polm80_sig_hzqq_AllEvents.Integral()
            print 'at polm ',dir_ind,'evt sig/sig all/qq/qqqq/qqqqqq',h_mass_polm80_sig_hzqq.Integral(),h_mass_polm80_sig_hzqq_AllEvents.Integral(),h_mass_polm80_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_polm80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_polm80_sig_ee_qqqqqq_BGHistos.Integral()
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
    canvas_polm80_BDT_significance.Print("hzqq_bbar_vs_totBG_polm80_sigificance_vs_BDT.eps")
    
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
    canvas_polm80_BDT_efficiency.Print("hzqq_bbar_vs_totBG_polm80_signal_efficiency_vs_BDT.eps")
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
    canvas_polm80_BDT_purity.Print("hzqq_bbar_vs_totBG_polm80_signal_purity_vs_BDT.eps")




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
    canvas_polm80_AllEvents_BDT_significance.Print("hzqq_bbar_vs_totBG_polm80_AllEvents_sigificance_vs_BDT.eps")
    
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
    canvas_polm80_AllEvents_BDT_efficiency.Print("hzqq_bbar_vs_totBG_polm80_AllEvents_signal_efficiency_vs_BDT.eps")
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
    canvas_polm80_AllEvents_BDT_purity.Print("hzqq_bbar_vs_totBG_polm80_AllEvents_signal_purity_vs_BDT.eps")


    #file_polp80_HZ_SignalHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVTrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj.root")  
    #file_polp80_HZ_SignalHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_all_jet1_and_jet2_subStruc_Vars_E_theta.root") 
    #file_polp80_HZ_SignalHistos_AllEvents_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_all_jet1_and_jet2_subStruc_Vars_E_theta_AllEvents.root") 
    #file_polp80_HZ_SignalHistos_AllEvents_noCorr_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_all_jet1_and_jet2_subStruc_Vars_E_theta_AllEvents.root")
    #file_polp80_ee_qq_mqq_1TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_all_jet1_and_jet2_subStruc_Vars_E_theta.root") 
    #file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_all_jet1_and_jet2_subStruc_Vars_E_theta.root") 
    #file_polp80_ee_qqqqqq_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_all_jet1_and_jet2_subStruc_Vars_E_theta.root") 


    #AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21_JESUnc0_99.root

    file_polp80_HZ_SignalHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root")  
    file_polp80_HZ_SignalHistos_AllEvents_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21_AllEvents.root")  
    file_polp80_HZ_SignalHistos_AllEventsParton_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/test_hzqq_noMass_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_withPartonHistos_AllEvents.root") 
    file_polp80_ee_qq_mqq_1TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root") 
    file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root") 
    file_polp80_ee_qqqqqq_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root") 
    #process_event(final_histo_name_,input_file_,files_weights_)

    h_2D_polp80_HZ_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polp80_HZ_SignalHistos_.Get("0.3/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polp80_HZ_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polp80_HZ_SignalHistos_.Get("0.3/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    print h_2D_polp80_HZ_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.GetBinContent(40,40)

    h_2D_polp80_HZ_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polp80_HZ_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")

    h_2D_polp80_HZ_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polp80_HZ_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
  
    h_2D_polp80_HZ_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polp80_HZ_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")

    h_2D_polp80_HZ_AllEvents_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_2D_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM")
    h_2D_polp80_HZ_AllEvents_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_2D_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM")
  

    h_2D_polp80_ee_qq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polp80_ee_qq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")

    h_2D_polp80_ee_qqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polp80_ee_qqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")

    h_2D_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_2D_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
    h_2D_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_2D_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM")
 
    A_pos_polp80_theta1_vs_theta2_signal=h_2D_polp80_HZ_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    B_neg_polp80_theta1_vs_theta2_signal=h_2D_polp80_HZ_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    alpha_polp80_theta1_vs_theta2_HZ=(A_pos_polp80_theta1_vs_theta2_signal+B_neg_polp80_theta1_vs_theta2_signal)/(abs(A_pos_polp80_theta1_vs_theta2_signal)+abs(B_neg_polp80_theta1_vs_theta2_signal))
    alpha_polp80_theta1_vs_theta2_HZ_error=sqrt(4.*abs(A_pos_polp80_theta1_vs_theta2_signal)*abs(B_neg_polp80_theta1_vs_theta2_signal)/pow((abs(A_pos_polp80_theta1_vs_theta2_signal)+abs(B_neg_polp80_theta1_vs_theta2_signal)),3))

    A_pos_polp80_theta1_vs_theta2_signal_all=h_2D_polp80_HZ_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    B_neg_polp80_theta1_vs_theta2_signal_all=h_2D_polp80_HZ_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    alpha_polp80_theta1_vs_theta2_HZ_all=(A_pos_polp80_theta1_vs_theta2_signal_all+B_neg_polp80_theta1_vs_theta2_signal_all)/(abs(A_pos_polp80_theta1_vs_theta2_signal_all)+abs(B_neg_polp80_theta1_vs_theta2_signal_all))
    alpha_polp80_theta1_vs_theta2_HZ_all_error=sqrt(4.*abs(A_pos_polp80_theta1_vs_theta2_signal_all)*abs(B_neg_polp80_theta1_vs_theta2_signal_all)/pow((abs(A_pos_polp80_theta1_vs_theta2_signal_all)+abs(B_neg_polp80_theta1_vs_theta2_signal_all)),3))
    
    A_pos_parton_polp80_theta1_vs_theta2_signal_all=h_2D_polp80_HZ_AllEvents_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Integral()
    B_neg_parton_polp80_theta1_vs_theta2_signal_all=h_2D_polp80_HZ_AllEvents_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Integral()
    alpha_parton_polp80_theta1_vs_theta2_HZ_all=(A_pos_parton_polp80_theta1_vs_theta2_signal_all+B_neg_parton_polp80_theta1_vs_theta2_signal_all)/(abs(A_pos_parton_polp80_theta1_vs_theta2_signal_all)+abs(B_neg_parton_polp80_theta1_vs_theta2_signal_all))
    alpha_parton_polp80_theta1_vs_theta2_HZ_all_error=sqrt(4.*abs(A_pos_parton_polp80_theta1_vs_theta2_signal_all)*abs(B_neg_parton_polp80_theta1_vs_theta2_signal_all)/pow((abs(A_pos_parton_polp80_theta1_vs_theta2_signal_all)+abs(B_neg_parton_polp80_theta1_vs_theta2_signal_all)),3))


    A_pos_polp80_theta1_vs_theta2_signal_BDTm020=h_2D_polp80_HZ_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020.Integral()
    B_neg_polp80_theta1_vs_theta2_signal_BDTm020=h_2D_polp80_HZ_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020.Integral()
    alpha_polp80_theta1_vs_theta2_HZ_BDTm020=(A_pos_polp80_theta1_vs_theta2_signal_BDTm020+B_neg_polp80_theta1_vs_theta2_signal_BDTm020)/(abs(A_pos_polp80_theta1_vs_theta2_signal_BDTm020)+abs(B_neg_polp80_theta1_vs_theta2_signal_BDTm020))
    alpha_polp80_theta1_vs_theta2_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polp80_theta1_vs_theta2_signal_BDTm020)*abs(B_neg_polp80_theta1_vs_theta2_signal_BDTm020)/pow((abs(A_pos_polp80_theta1_vs_theta2_signal_BDTm020)+abs(B_neg_polp80_theta1_vs_theta2_signal_BDTm020)),3))
    A_pos_polp80_theta1_vs_theta2_signal_all_BDTm020=h_2D_polp80_HZ_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020.Integral()
    B_neg_polp80_theta1_vs_theta2_signal_all_BDTm020=h_2D_polp80_HZ_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM_BDTm020.Integral()
    alpha_polp80_theta1_vs_theta2_HZ_all_BDTm020=(A_pos_polp80_theta1_vs_theta2_signal_all_BDTm020+B_neg_polp80_theta1_vs_theta2_signal_all_BDTm020)/(abs(A_pos_polp80_theta1_vs_theta2_signal_all_BDTm020)+abs(B_neg_polp80_theta1_vs_theta2_signal_all_BDTm020))
    alpha_polp80_theta1_vs_theta2_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polp80_theta1_vs_theta2_signal_all_BDTm020)*abs(B_neg_polp80_theta1_vs_theta2_signal_all_BDTm020)/pow((abs(A_pos_polp80_theta1_vs_theta2_signal_all_BDTm020)+abs(B_neg_polp80_theta1_vs_theta2_signal_all_BDTm020)),3))

    A_pos_polp80_theta1_vs_theta2_signal_all_withBG=h_2D_polp80_HZ_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polp80_ee_qq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polp80_ee_qqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    B_neg_polp80_theta1_vs_theta2_signal_all_withBG=h_2D_polp80_HZ_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polp80_ee_qq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polp80_ee_qqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    alpha_polp80_theta1_vs_theta2_HZ_all_withBG=(A_pos_polp80_theta1_vs_theta2_signal_all_withBG+B_neg_polp80_theta1_vs_theta2_signal_all_withBG)/(abs(A_pos_polp80_theta1_vs_theta2_signal_all_withBG)+abs(B_neg_polp80_theta1_vs_theta2_signal_all_withBG))
    alpha_polp80_theta1_vs_theta2_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polp80_theta1_vs_theta2_signal_all_withBG)*abs(B_neg_polp80_theta1_vs_theta2_signal_all_withBG)/pow((abs(A_pos_polp80_theta1_vs_theta2_signal_all_withBG)+abs(B_neg_polp80_theta1_vs_theta2_signal_all_withBG)),3))

    A_pos_polp80_theta1_vs_theta2_allBG=h_2D_polp80_ee_qq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polp80_ee_qqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    B_neg_polp80_theta1_vs_theta2_allBG=h_2D_polp80_ee_qq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polp80_ee_qqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()+h_2D_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgncostheta1_times_sgncostheta2_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Integral()
    alpha_polp80_theta1_vs_theta2_allBG=(A_pos_polp80_theta1_vs_theta2_allBG+B_neg_polp80_theta1_vs_theta2_allBG)/(abs(A_pos_polp80_theta1_vs_theta2_allBG)+abs(B_neg_polp80_theta1_vs_theta2_allBG))
    alpha_polp80_theta1_vs_theta2_allBG_error=sqrt(4.*abs(A_pos_polp80_theta1_vs_theta2_allBG)*abs(B_neg_polp80_theta1_vs_theta2_allBG)/pow((abs(A_pos_polp80_theta1_vs_theta2_allBG)+abs(B_neg_polp80_theta1_vs_theta2_allBG)),3))

    print 'alpha_polp80_theta1_vs_theta2_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_theta1_vs_theta2_HZ,alpha_polp80_theta1_vs_theta2_HZ_all,alpha_polp80_theta1_vs_theta2_HZ_BDTm020,alpha_polp80_theta1_vs_theta2_HZ_all_BDTm020,alpha_polp80_theta1_vs_theta2_HZ_all_withBG
    print 'errors of alpha_polp80_theta1_vs_theta2_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_theta1_vs_theta2_HZ_error,alpha_polp80_theta1_vs_theta2_HZ_all_error,alpha_polp80_theta1_vs_theta2_HZ_BDTm020_error,alpha_polp80_theta1_vs_theta2_HZ_all_BDTm020_error,alpha_polp80_theta1_vs_theta2_HZ_all_withBG_error


    print 'alpha_polp80_theta1_vs_theta2_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polp80_theta1_vs_theta2_HZ,alpha_polp80_theta1_vs_theta2_HZ_all,alpha_polp80_theta1_vs_theta2_allBG,alpha_polp80_theta1_vs_theta2_HZ_all_withBG,alpha_parton_polp80_theta1_vs_theta2_HZ_all
    print 'errors of alpha_polp80_theta1_vs_theta2_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polp80_theta1_vs_theta2_HZ_error,alpha_polp80_theta1_vs_theta2_HZ_all_error,alpha_polp80_theta1_vs_theta2_allBG_error,alpha_polp80_theta1_vs_theta2_HZ_all_withBG_error,alpha_parton_polp80_theta1_vs_theta2_HZ_all_error



    h_polp80_HZ_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_HZ_SignalHistos_.Get("0.3/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polp80_HZ_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_HZ_SignalHistos_.Get("0.3/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polp80_HZ_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polp80_HZ_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")

    h_polp80_HZ_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polp80_HZ_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polp80_HZ_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polp80_HZ_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")

    h_polp80_HZ_AllEvents_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom")
    h_polp80_HZ_AllEvents_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom")
  

    h_polp80_ee_qq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polp80_ee_qq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")

    h_polp80_ee_qqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polp80_ee_qqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")

    h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")

    A_pos_polp80_theta1_signal=h_polp80_HZ_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    B_neg_polp80_theta1_signal=h_polp80_HZ_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    alpha_polp80_theta1_HZ=(A_pos_polp80_theta1_signal+B_neg_polp80_theta1_signal)/(abs(A_pos_polp80_theta1_signal)+abs(B_neg_polp80_theta1_signal))
    alpha_polp80_theta1_HZ_error=sqrt(4.*abs(A_pos_polp80_theta1_signal)*abs(B_neg_polp80_theta1_signal)/pow((abs(A_pos_polp80_theta1_signal)+abs(B_neg_polp80_theta1_signal)),3))

    A_pos_polp80_theta1_signal_all=h_polp80_HZ_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    B_neg_polp80_theta1_signal_all=h_polp80_HZ_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    alpha_polp80_theta1_HZ_all=(A_pos_polp80_theta1_signal_all+B_neg_polp80_theta1_signal_all)/(abs(A_pos_polp80_theta1_signal_all)+abs(B_neg_polp80_theta1_signal_all))
    alpha_polp80_theta1_HZ_all_error=sqrt(4.*abs(A_pos_polp80_theta1_signal_all)*abs(B_neg_polp80_theta1_signal_all)/pow((abs(A_pos_polp80_theta1_signal_all)+abs(B_neg_polp80_theta1_signal_all)),3))
 
    A_pos_parton_polp80_theta1_signal_all=h_polp80_HZ_AllEvents_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral()
    B_neg_parton_polp80_theta1_signal_all=h_polp80_HZ_AllEvents_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral()
    alpha_parton_polp80_theta1_HZ_all=(A_pos_parton_polp80_theta1_signal_all+B_neg_parton_polp80_theta1_signal_all)/(abs(A_pos_parton_polp80_theta1_signal_all)+abs(B_neg_parton_polp80_theta1_signal_all))
    alpha_parton_polp80_theta1_HZ_all_error=sqrt(4.*abs(A_pos_parton_polp80_theta1_signal_all)*abs(B_neg_parton_polp80_theta1_signal_all)/pow((abs(A_pos_parton_polp80_theta1_signal_all)+abs(B_neg_parton_polp80_theta1_signal_all)),3))



    A_pos_polp80_theta1_signal_BDTm020=h_polp80_HZ_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020.Integral()
    B_neg_polp80_theta1_signal_BDTm020=h_polp80_HZ_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020.Integral()
    alpha_polp80_theta1_HZ_BDTm020=(A_pos_polp80_theta1_signal_BDTm020+B_neg_polp80_theta1_signal_BDTm020)/(abs(A_pos_polp80_theta1_signal_BDTm020)+abs(B_neg_polp80_theta1_signal_BDTm020))
    alpha_polp80_theta1_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polp80_theta1_signal_BDTm020)*abs(B_neg_polp80_theta1_signal_BDTm020)/pow((abs(A_pos_polp80_theta1_signal_BDTm020)+abs(B_neg_polp80_theta1_signal_BDTm020)),3))
    A_pos_polp80_theta1_signal_all_BDTm020=h_polp80_HZ_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020.Integral()
    B_neg_polp80_theta1_signal_all_BDTm020=h_polp80_HZ_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_BDTm020.Integral()
    alpha_polp80_theta1_HZ_all_BDTm020=(A_pos_polp80_theta1_signal_all_BDTm020+B_neg_polp80_theta1_signal_all_BDTm020)/(abs(A_pos_polp80_theta1_signal_all_BDTm020)+abs(B_neg_polp80_theta1_signal_all_BDTm020))
    alpha_polp80_theta1_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polp80_theta1_signal_all_BDTm020)*abs(B_neg_polp80_theta1_signal_all_BDTm020)/pow((abs(A_pos_polp80_theta1_signal_all_BDTm020)+abs(B_neg_polp80_theta1_signal_all_BDTm020)),3))

    A_pos_polp80_theta1_signal_all_withBG=h_polp80_HZ_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polp80_ee_qq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    B_neg_polp80_theta1_signal_all_withBG=h_polp80_HZ_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polp80_ee_qq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    alpha_polp80_theta1_HZ_all_withBG=(A_pos_polp80_theta1_signal_all_withBG+B_neg_polp80_theta1_signal_all_withBG)/(abs(A_pos_polp80_theta1_signal_all_withBG)+abs(B_neg_polp80_theta1_signal_all_withBG))
    alpha_polp80_theta1_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polp80_theta1_signal_all_withBG)*abs(B_neg_polp80_theta1_signal_all_withBG)/pow((abs(A_pos_polp80_theta1_signal_all_withBG)+abs(B_neg_polp80_theta1_signal_all_withBG)),3))

    A_pos_polp80_theta1_allBG=h_polp80_ee_qq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    B_neg_polp80_theta1_allBG=h_polp80_ee_qq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Integral()
    alpha_polp80_theta1_allBG=(A_pos_polp80_theta1_allBG+B_neg_polp80_theta1_allBG)/(abs(A_pos_polp80_theta1_allBG)+abs(B_neg_polp80_theta1_allBG))
    alpha_polp80_theta1_allBG_error=sqrt(4.*abs(A_pos_polp80_theta1_allBG)*abs(B_neg_polp80_theta1_allBG)/pow((abs(A_pos_polp80_theta1_allBG)+abs(B_neg_polp80_theta1_allBG)),3))
  
    print 'alpha_polp80_theta1_HZ/HZallevents/HZ BDTm020/HZ all BDTm020/ BG',alpha_polp80_theta1_HZ,alpha_polp80_theta1_HZ_all,alpha_polp80_theta1_HZ_BDTm020,alpha_polp80_theta1_HZ_all_BDTm020,alpha_polp80_theta1_HZ_all_withBG
    print 'errors of alpha_polp80_theta1_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_theta1_HZ_error,alpha_polp80_theta1_HZ_all_error,alpha_polp80_theta1_HZ_BDTm020_error,alpha_polp80_theta1_HZ_all_BDTm020_error,alpha_polp80_theta1_HZ_all_withBG_error

    print 'alpha_polp80_theta1_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polp80_theta1_HZ,alpha_polp80_theta1_HZ_all,alpha_polp80_theta1_allBG,alpha_polp80_theta1_HZ_all_withBG,alpha_parton_polp80_theta1_HZ_all
    print 'errors of alpha_polp80_theta1_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polp80_theta1_HZ_error,alpha_polp80_theta1_HZ_all_error,alpha_polp80_theta1_allBG_error,alpha_polp80_theta1_HZ_all_withBG_error,alpha_parton_polp80_theta1_HZ_all_error



    h_polp80_HZ_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_.Get("0.3/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_.Get("0.3/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_HZ_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_HZ_AllEvents_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep")
    h_polp80_HZ_AllEvents_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep")

    h_polp80_HZ_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")


    A_pos_polp80_phi1_signal=h_polp80_HZ_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi1_signal=h_polp80_HZ_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi1_HZ=(A_pos_polp80_phi1_signal+B_neg_polp80_phi1_signal)/(abs(A_pos_polp80_phi1_signal)+abs(B_neg_polp80_phi1_signal))
    alpha_polp80_phi1_HZ_error=sqrt(4.*abs(A_pos_polp80_phi1_signal)*abs(B_neg_polp80_phi1_signal)/pow((abs(A_pos_polp80_phi1_signal)+abs(B_neg_polp80_phi1_signal)),3))

    A_pos_polp80_phi1_signal_all=h_polp80_HZ_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi1_signal_all=h_polp80_HZ_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi1_HZ_all=(A_pos_polp80_phi1_signal_all+B_neg_polp80_phi1_signal_all)/(abs(A_pos_polp80_phi1_signal_all)+abs(B_neg_polp80_phi1_signal_all))
    alpha_polp80_phi1_HZ_all_error=sqrt(4.*abs(A_pos_polp80_phi1_signal_all)*abs(B_neg_polp80_phi1_signal_all)/pow((abs(A_pos_polp80_phi1_signal_all)+abs(B_neg_polp80_phi1_signal_all)),3))
    
    A_pos_parton_polp80_phi1_signal_all=h_polp80_HZ_AllEvents_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    B_neg_parton_polp80_phi1_signal_all=h_polp80_HZ_AllEvents_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    alpha_parton_polp80_phi1_HZ_all=(A_pos_parton_polp80_phi1_signal_all+B_neg_parton_polp80_phi1_signal_all)/(abs(A_pos_parton_polp80_phi1_signal_all)+abs(B_neg_parton_polp80_phi1_signal_all))
    alpha_parton_polp80_phi1_HZ_all_error=sqrt(4.*abs(A_pos_parton_polp80_phi1_signal_all)*abs(B_neg_parton_polp80_phi1_signal_all)/pow((abs(A_pos_parton_polp80_phi1_signal_all)+abs(B_neg_parton_polp80_phi1_signal_all)),3))



    A_pos_polp80_phi1_signal_BDTm020=h_polp80_HZ_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polp80_phi1_signal_BDTm020=h_polp80_HZ_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polp80_phi1_HZ_BDTm020=(A_pos_polp80_phi1_signal_BDTm020+B_neg_polp80_phi1_signal_BDTm020)/(abs(A_pos_polp80_phi1_signal_BDTm020)+abs(B_neg_polp80_phi1_signal_BDTm020))
    alpha_polp80_phi1_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polp80_phi1_signal_BDTm020)*abs(B_neg_polp80_phi1_signal_BDTm020)/pow((abs(A_pos_polp80_phi1_signal_BDTm020)+abs(B_neg_polp80_phi1_signal_BDTm020)),3))
    A_pos_polp80_phi1_signal_all_BDTm020=h_polp80_HZ_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polp80_phi1_signal_all_BDTm020=h_polp80_HZ_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polp80_phi1_HZ_all_BDTm020=(A_pos_polp80_phi1_signal_all_BDTm020+B_neg_polp80_phi1_signal_all_BDTm020)/(abs(A_pos_polp80_phi1_signal_all_BDTm020)+abs(B_neg_polp80_phi1_signal_all_BDTm020))
    alpha_polp80_phi1_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polp80_phi1_signal_all_BDTm020)*abs(B_neg_polp80_phi1_signal_all_BDTm020)/pow((abs(A_pos_polp80_phi1_signal_all_BDTm020)+abs(B_neg_polp80_phi1_signal_all_BDTm020)),3))

    A_pos_polp80_phi1_signal_all_withBG=h_polp80_HZ_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi1_signal_all_withBG=h_polp80_HZ_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi1_HZ_all_withBG=(A_pos_polp80_phi1_signal_all_withBG+B_neg_polp80_phi1_signal_all_withBG)/(abs(A_pos_polp80_phi1_signal_all_withBG)+abs(B_neg_polp80_phi1_signal_all_withBG))
    alpha_polp80_phi1_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polp80_phi1_signal_all_withBG)*abs(B_neg_polp80_phi1_signal_all_withBG)/pow((abs(A_pos_polp80_phi1_signal_all_withBG)+abs(B_neg_polp80_phi1_signal_all_withBG)),3))

  
    A_pos_polp80_phi1_allBG=h_polp80_ee_qq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi1_allBG=h_polp80_ee_qq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi1_allBG=(A_pos_polp80_phi1_allBG+B_neg_polp80_phi1_allBG)/(abs(A_pos_polp80_phi1_allBG)+abs(B_neg_polp80_phi1_allBG))
    alpha_polp80_phi1_allBG_error=sqrt(4.*abs(A_pos_polp80_phi1_allBG)*abs(B_neg_polp80_phi1_allBG)/pow((abs(A_pos_polp80_phi1_allBG)+abs(B_neg_polp80_phi1_allBG)),3))

    print 'alpha_polp80_phi1_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_phi1_HZ,alpha_polp80_phi1_HZ_all,alpha_polp80_phi1_HZ_BDTm020,alpha_polp80_phi1_HZ_all_BDTm020,alpha_polp80_phi1_HZ_all_withBG
    print 'errors of alpha_polp80_phi1_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_phi1_HZ_error,alpha_polp80_phi1_HZ_all_error,alpha_polp80_phi1_HZ_BDTm020_error,alpha_polp80_phi1_HZ_all_BDTm020_error,alpha_polp80_phi1_HZ_all_withBG_error

    print 'alpha_polp80_phi1_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polp80_phi1_HZ,alpha_polp80_phi1_HZ_all,alpha_polp80_phi1_allBG,alpha_polp80_phi1_HZ_all_withBG,alpha_parton_polp80_phi1_HZ_all
    print 'errors of alpha_polp80_phi1_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polp80_phi1_HZ_error,alpha_polp80_phi1_HZ_all_error,alpha_polp80_phi1_allBG_error,alpha_polp80_phi1_HZ_all_withBG_error,alpha_parton_polp80_phi1_HZ_all_error


    h_polp80_HZ_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_.Get("0.3/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_.Get("0.3/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_HZ_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_HZ_AllEvents_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep")
    h_polp80_HZ_AllEvents_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep")

    h_polp80_HZ_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qqqq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    A_pos_polp80_phi2_signal=h_polp80_HZ_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi2_signal=h_polp80_HZ_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi2_HZ=(A_pos_polp80_phi2_signal+B_neg_polp80_phi2_signal)/(abs(A_pos_polp80_phi2_signal)+abs(B_neg_polp80_phi2_signal))
    alpha_polp80_phi2_HZ_error=sqrt(4.*abs(A_pos_polp80_phi2_signal)*abs(B_neg_polp80_phi2_signal)/pow((abs(A_pos_polp80_phi2_signal)+abs(B_neg_polp80_phi2_signal)),3))

    A_pos_polp80_phi2_signal_all=h_polp80_HZ_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi2_signal_all=h_polp80_HZ_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi2_HZ_all=(A_pos_polp80_phi2_signal_all+B_neg_polp80_phi2_signal_all)/(abs(A_pos_polp80_phi2_signal_all)+abs(B_neg_polp80_phi2_signal_all))
    alpha_polp80_phi2_HZ_all_error=sqrt(4.*abs(A_pos_polp80_phi2_signal_all)*abs(B_neg_polp80_phi2_signal_all)/pow((abs(A_pos_polp80_phi2_signal_all)+abs(B_neg_polp80_phi2_signal_all)),3))
    
    A_pos_parton_polp80_phi2_signal_all=h_polp80_HZ_AllEvents_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    B_neg_parton_polp80_phi2_signal_all=h_polp80_HZ_AllEvents_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    alpha_parton_polp80_phi2_HZ_all=(A_pos_parton_polp80_phi2_signal_all+B_neg_parton_polp80_phi2_signal_all)/(abs(A_pos_parton_polp80_phi2_signal_all)+abs(B_neg_parton_polp80_phi2_signal_all))
    alpha_parton_polp80_phi2_HZ_all_error=sqrt(4.*abs(A_pos_parton_polp80_phi2_signal_all)*abs(B_neg_parton_polp80_phi2_signal_all)/pow((abs(A_pos_parton_polp80_phi2_signal_all)+abs(B_neg_parton_polp80_phi2_signal_all)),3))


    A_pos_polp80_phi2_signal_BDTm020=h_polp80_HZ_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polp80_phi2_signal_BDTm020=h_polp80_HZ_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polp80_phi2_HZ_BDTm020=(A_pos_polp80_phi2_signal_BDTm020+B_neg_polp80_phi2_signal_BDTm020)/(abs(A_pos_polp80_phi2_signal_BDTm020)+abs(B_neg_polp80_phi2_signal_BDTm020))
    alpha_polp80_phi2_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polp80_phi2_signal_BDTm020)*abs(B_neg_polp80_phi2_signal_BDTm020)/pow((abs(A_pos_polp80_phi2_signal_BDTm020)+abs(B_neg_polp80_phi2_signal_BDTm020)),3))
    A_pos_polp80_phi2_signal_all_BDTm020=h_polp80_HZ_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polp80_phi2_signal_all_BDTm020=h_polp80_HZ_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polp80_phi2_HZ_all_BDTm020=(A_pos_polp80_phi2_signal_all_BDTm020+B_neg_polp80_phi2_signal_all_BDTm020)/(abs(A_pos_polp80_phi2_signal_all_BDTm020)+abs(B_neg_polp80_phi2_signal_all_BDTm020))
    alpha_polp80_phi2_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polp80_phi2_signal_all_BDTm020)*abs(B_neg_polp80_phi2_signal_all_BDTm020)/pow((abs(A_pos_polp80_phi2_signal_all_BDTm020)+abs(B_neg_polp80_phi2_signal_all_BDTm020)),3))

    A_pos_polp80_phi2_signal_all_withBG=h_polp80_HZ_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi2_signal_all_withBG=h_polp80_HZ_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi2_HZ_all_withBG=(A_pos_polp80_phi2_signal_all_withBG+B_neg_polp80_phi2_signal_all_withBG)/(abs(A_pos_polp80_phi2_signal_all_withBG)+abs(B_neg_polp80_phi2_signal_all_withBG))
    alpha_polp80_phi2_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polp80_phi2_signal_all_withBG)*abs(B_neg_polp80_phi2_signal_all_withBG)/pow((abs(A_pos_polp80_phi2_signal_all_withBG)+abs(B_neg_polp80_phi2_signal_all_withBG)),3))

    A_pos_polp80_phi2_allBG=h_polp80_ee_qq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi2_allBG=h_polp80_ee_qq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_neg_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi2_allBG=(A_pos_polp80_phi2_allBG+B_neg_polp80_phi2_allBG)/(abs(A_pos_polp80_phi2_allBG)+abs(B_neg_polp80_phi2_allBG))
    alpha_polp80_phi2_allBG_error=sqrt(4.*abs(A_pos_polp80_phi2_allBG)*abs(B_neg_polp80_phi2_allBG)/pow((abs(A_pos_polp80_phi2_allBG)+abs(B_neg_polp80_phi2_allBG)),3))
  
    print 'alpha_polp80_phi2_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_phi2_HZ,alpha_polp80_phi2_HZ_all,alpha_polp80_phi2_HZ_BDTm020,alpha_polp80_phi2_HZ_all_BDTm020,alpha_polp80_phi2_HZ_all_withBG
    print 'errors of alpha_polp80_phi2_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_phi2_HZ_error,alpha_polp80_phi2_HZ_all_error,alpha_polp80_phi2_HZ_BDTm020_error,alpha_polp80_phi2_HZ_all_BDTm020_error,alpha_polp80_phi2_HZ_all_withBG_error

    print 'alpha_polp80_phi2_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polp80_phi2_HZ,alpha_polp80_phi2_HZ_all,alpha_polp80_phi2_allBG,alpha_polp80_phi2_HZ_all_withBG,alpha_parton_polp80_phi2_HZ_all
    print 'errors of alpha_polp80_phi2_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polp80_phi2_HZ_error,alpha_polp80_phi2_HZ_all_error,alpha_polp80_phi2_allBG_error,alpha_polp80_phi2_HZ_all_withBG_error,alpha_parton_polp80_phi2_HZ_all_error


    h_polp80_HZ_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_.Get("0.3/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_.Get("0.3/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_HZ_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_HZ_AllEvents_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep")
    h_polp80_HZ_AllEvents_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep")

    h_polp80_HZ_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qqqq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    A_pos_polp80_phi3_signal=h_polp80_HZ_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi3_signal=h_polp80_HZ_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi3_HZ=(A_pos_polp80_phi3_signal+B_neg_polp80_phi3_signal)/(abs(A_pos_polp80_phi3_signal)+abs(B_neg_polp80_phi3_signal))
    alpha_polp80_phi3_HZ_error=sqrt(4.*abs(A_pos_polp80_phi3_signal)*abs(B_neg_polp80_phi3_signal)/pow((abs(A_pos_polp80_phi3_signal)+abs(B_neg_polp80_phi3_signal)),3))

    A_pos_polp80_phi3_signal_all=h_polp80_HZ_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi3_signal_all=h_polp80_HZ_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi3_HZ_all=(A_pos_polp80_phi3_signal_all+B_neg_polp80_phi3_signal_all)/(abs(A_pos_polp80_phi3_signal_all)+abs(B_neg_polp80_phi3_signal_all))
    alpha_polp80_phi3_HZ_all_error=sqrt(4.*abs(A_pos_polp80_phi3_signal_all)*abs(B_neg_polp80_phi3_signal_all)/pow((abs(A_pos_polp80_phi3_signal_all)+abs(B_neg_polp80_phi3_signal_all)),3))
    
    A_pos_parton_polp80_phi3_signal_all=h_polp80_HZ_AllEvents_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    B_neg_parton_polp80_phi3_signal_all=h_polp80_HZ_AllEvents_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    alpha_parton_polp80_phi3_HZ_all=(A_pos_parton_polp80_phi3_signal_all+B_neg_parton_polp80_phi3_signal_all)/(abs(A_pos_parton_polp80_phi3_signal_all)+abs(B_neg_parton_polp80_phi3_signal_all))
    alpha_parton_polp80_phi3_HZ_all_error=sqrt(4.*abs(A_pos_parton_polp80_phi3_signal_all)*abs(B_neg_parton_polp80_phi3_signal_all)/pow((abs(A_pos_parton_polp80_phi3_signal_all)+abs(B_neg_parton_polp80_phi3_signal_all)),3))


    A_pos_polp80_phi3_signal_BDTm020=h_polp80_HZ_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polp80_phi3_signal_BDTm020=h_polp80_HZ_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polp80_phi3_HZ_BDTm020=(A_pos_polp80_phi3_signal_BDTm020+B_neg_polp80_phi3_signal_BDTm020)/(abs(A_pos_polp80_phi3_signal_BDTm020)+abs(B_neg_polp80_phi3_signal_BDTm020))
    alpha_polp80_phi3_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polp80_phi3_signal_BDTm020)*abs(B_neg_polp80_phi3_signal_BDTm020)/pow((abs(A_pos_polp80_phi3_signal_BDTm020)+abs(B_neg_polp80_phi3_signal_BDTm020)),3))
    A_pos_polp80_phi3_signal_all_BDTm020=h_polp80_HZ_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polp80_phi3_signal_all_BDTm020=h_polp80_HZ_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polp80_phi3_HZ_all_BDTm020=(A_pos_polp80_phi3_signal_all_BDTm020+B_neg_polp80_phi3_signal_all_BDTm020)/(abs(A_pos_polp80_phi3_signal_all_BDTm020)+abs(B_neg_polp80_phi3_signal_all_BDTm020))
    alpha_polp80_phi3_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polp80_phi3_signal_all_BDTm020)*abs(B_neg_polp80_phi3_signal_all_BDTm020)/pow((abs(A_pos_polp80_phi3_signal_all_BDTm020)+abs(B_neg_polp80_phi3_signal_all_BDTm020)),3))

    A_pos_polp80_phi3_signal_all_withBG=h_polp80_HZ_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi3_signal_all_withBG=h_polp80_HZ_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi3_HZ_all_withBG=(A_pos_polp80_phi3_signal_all_withBG+B_neg_polp80_phi3_signal_all_withBG)/(abs(A_pos_polp80_phi3_signal_all_withBG)+abs(B_neg_polp80_phi3_signal_all_withBG))
    alpha_polp80_phi3_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polp80_phi3_signal_all_withBG)*abs(B_neg_polp80_phi3_signal_all_withBG)/pow((abs(A_pos_polp80_phi3_signal_all_withBG)+abs(B_neg_polp80_phi3_signal_all_withBG)),3))

    A_pos_polp80_phi3_allBG=h_polp80_ee_qq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi3_allBG=h_polp80_ee_qq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_neg_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi3_allBG=(A_pos_polp80_phi3_allBG+B_neg_polp80_phi3_allBG)/(abs(A_pos_polp80_phi3_allBG)+abs(B_neg_polp80_phi3_allBG))
    alpha_polp80_phi3_allBG_error=sqrt(4.*abs(A_pos_polp80_phi3_allBG)*abs(B_neg_polp80_phi3_allBG)/pow((abs(A_pos_polp80_phi3_allBG)+abs(B_neg_polp80_phi3_allBG)),3))
  
    print 'alpha_polp80_phi3_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_phi3_HZ,alpha_polp80_phi3_HZ_all,alpha_polp80_phi3_HZ_BDTm020,alpha_polp80_phi3_HZ_all_BDTm020,alpha_polp80_phi3_HZ_all_withBG
    print 'errors of alpha_polp80_phi3_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_phi3_HZ_error,alpha_polp80_phi3_HZ_all_error,alpha_polp80_phi3_HZ_BDTm020_error,alpha_polp80_phi3_HZ_all_BDTm020_error,alpha_polp80_phi3_HZ_all_withBG_error


    print 'alpha_polp80_phi3_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polp80_phi3_HZ,alpha_polp80_phi3_HZ_all,alpha_polp80_phi3_allBG,alpha_polp80_phi3_HZ_all_withBG,alpha_parton_polp80_phi3_HZ_all
    print 'errors of alpha_polp80_phi3_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polp80_phi3_HZ_error,alpha_polp80_phi3_HZ_all_error,alpha_polp80_phi3_allBG_error,alpha_polp80_phi3_HZ_all_withBG_error,alpha_parton_polp80_phi3_HZ_all_error




    h_polp80_HZ_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_.Get("0.3/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_.Get("0.3/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_HZ_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_.Get("-0.2/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_HZ_AllEvents_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep")
    h_polp80_HZ_AllEvents_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep=file_polp80_HZ_SignalHistos_AllEventsParton_.Get("h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep")

    h_polp80_HZ_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_HZ_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qqqq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")

    A_pos_polp80_phi4_signal=h_polp80_HZ_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi4_signal=h_polp80_HZ_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi4_HZ=(A_pos_polp80_phi4_signal+B_neg_polp80_phi4_signal)/(abs(A_pos_polp80_phi4_signal)+abs(B_neg_polp80_phi4_signal))
    alpha_polp80_phi4_HZ_error=sqrt(4.*abs(A_pos_polp80_phi4_signal)*abs(B_neg_polp80_phi4_signal)/pow((abs(A_pos_polp80_phi4_signal)+abs(B_neg_polp80_phi4_signal)),3))

    A_pos_polp80_phi4_signal_all=h_polp80_HZ_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi4_signal_all=h_polp80_HZ_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi4_HZ_all=(A_pos_polp80_phi4_signal_all+B_neg_polp80_phi4_signal_all)/(abs(A_pos_polp80_phi4_signal_all)+abs(B_neg_polp80_phi4_signal_all))
    alpha_polp80_phi4_HZ_all_error=sqrt(4.*abs(A_pos_polp80_phi4_signal_all)*abs(B_neg_polp80_phi4_signal_all)/pow((abs(A_pos_polp80_phi4_signal_all)+abs(B_neg_polp80_phi4_signal_all)),3))
    
    A_pos_parton_polp80_phi4_signal_all=h_polp80_HZ_AllEvents_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    B_neg_parton_polp80_phi4_signal_all=h_polp80_HZ_AllEvents_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
    alpha_parton_polp80_phi4_HZ_all=(A_pos_parton_polp80_phi4_signal_all+B_neg_parton_polp80_phi4_signal_all)/(abs(A_pos_parton_polp80_phi4_signal_all)+abs(B_neg_parton_polp80_phi4_signal_all))
    alpha_parton_polp80_phi4_HZ_all_error=sqrt(4.*abs(A_pos_parton_polp80_phi4_signal_all)*abs(B_neg_parton_polp80_phi4_signal_all)/pow((abs(A_pos_parton_polp80_phi4_signal_all)+abs(B_neg_parton_polp80_phi4_signal_all)),3))


    A_pos_polp80_phi4_signal_BDTm020=h_polp80_HZ_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polp80_phi4_signal_BDTm020=h_polp80_HZ_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polp80_phi4_HZ_BDTm020=(A_pos_polp80_phi4_signal_BDTm020+B_neg_polp80_phi4_signal_BDTm020)/(abs(A_pos_polp80_phi4_signal_BDTm020)+abs(B_neg_polp80_phi4_signal_BDTm020))
    alpha_polp80_phi4_HZ_BDTm020_error=sqrt(4.*abs(A_pos_polp80_phi4_signal_BDTm020)*abs(B_neg_polp80_phi4_signal_BDTm020)/pow((abs(A_pos_polp80_phi4_signal_BDTm020)+abs(B_neg_polp80_phi4_signal_BDTm020)),3))
    A_pos_polp80_phi4_signal_all_BDTm020=h_polp80_HZ_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    B_neg_polp80_phi4_signal_all_BDTm020=h_polp80_HZ_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_BDTm020.Integral()
    alpha_polp80_phi4_HZ_all_BDTm020=(A_pos_polp80_phi4_signal_all_BDTm020+B_neg_polp80_phi4_signal_all_BDTm020)/(abs(A_pos_polp80_phi4_signal_all_BDTm020)+abs(B_neg_polp80_phi4_signal_all_BDTm020))
    alpha_polp80_phi4_HZ_all_BDTm020_error=sqrt(4.*abs(A_pos_polp80_phi4_signal_all_BDTm020)*abs(B_neg_polp80_phi4_signal_all_BDTm020)/pow((abs(A_pos_polp80_phi4_signal_all_BDTm020)+abs(B_neg_polp80_phi4_signal_all_BDTm020)),3))

    A_pos_polp80_phi4_signal_all_withBG=h_polp80_HZ_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi4_signal_all_withBG=h_polp80_HZ_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi4_HZ_all_withBG=(A_pos_polp80_phi4_signal_all_withBG+B_neg_polp80_phi4_signal_all_withBG)/(abs(A_pos_polp80_phi4_signal_all_withBG)+abs(B_neg_polp80_phi4_signal_all_withBG))
    alpha_polp80_phi4_HZ_all_withBG_error=sqrt(4.*abs(A_pos_polp80_phi4_signal_all_withBG)*abs(B_neg_polp80_phi4_signal_all_withBG)/pow((abs(A_pos_polp80_phi4_signal_all_withBG)+abs(B_neg_polp80_phi4_signal_all_withBG)),3))

    A_pos_polp80_phi4_allBG=h_polp80_ee_qq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_pos_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    B_neg_polp80_phi4_allBG=h_polp80_ee_qq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqq_AllEvents_recojet_neg_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()+h_polp80_ee_qqqqqq_AllEvents_recojet_neg_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Integral()
    alpha_polp80_phi4_allBG=(A_pos_polp80_phi4_allBG+B_neg_polp80_phi4_allBG)/(abs(A_pos_polp80_phi4_allBG)+abs(B_neg_polp80_phi4_allBG))
    alpha_polp80_phi4_allBG_error=sqrt(4.*abs(A_pos_polp80_phi4_allBG)*abs(B_neg_polp80_phi4_allBG)/pow((abs(A_pos_polp80_phi4_allBG)+abs(B_neg_polp80_phi4_allBG)),3))
  
    print 'alpha_polp80_phi4_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_phi4_HZ,alpha_polp80_phi4_HZ_all,alpha_polp80_phi4_HZ_BDTm020,alpha_polp80_phi4_HZ_all_BDTm020,alpha_polp80_phi4_HZ_all_withBG
    print 'errors of alpha_polp80_phi4_HZ/HZallevents/HZ BDT 020/HZ all BDTm020',alpha_polp80_phi4_HZ_error,alpha_polp80_phi4_HZ_all_error,alpha_polp80_phi4_HZ_BDTm020_error,alpha_polp80_phi4_HZ_all_BDTm020_error,alpha_polp80_phi4_HZ_all_withBG_error

    print 'alpha_polp80_phi4_HZ/HZallevents/BG/HZ+BG/HZ parton',alpha_polp80_phi4_HZ,alpha_polp80_phi4_HZ_all,alpha_polp80_phi4_allBG,alpha_polp80_phi4_HZ_all_withBG,alpha_parton_polp80_phi4_HZ_all
    print 'errors of alpha_polp80_phi4_HZ/HZallevents/BG/HZ+BG/HZ-parton',alpha_polp80_phi4_HZ_error,alpha_polp80_phi4_HZ_all_error,alpha_polp80_phi4_allBG_error,alpha_polp80_phi4_HZ_all_withBG_error,alpha_parton_polp80_phi4_HZ_all_error






    x_polp80_BDTScore = array( 'f' )
    y_polp80_significance = array( 'f' )
    y_polp80_purity = array( 'f' )
    y_polp80_efficiency = array( 'f' )
    norm_polp=file_polp80_HZ_SignalHistos_.Get("-0.2/h_BDT_output").Integral()

    y_polp80_AllEvents_significance = array( 'f' )
    y_polp80_AllEvents_purity = array( 'f' )
    y_polp80_AllEvents_efficiency = array( 'f' )
    norm_polp_AllEvents=file_polp80_HZ_SignalHistos_AllEvents_.Get("-0.2/h_BDT_output").Integral()
    #norm_polp=h_mass_sig_hzqq_norm_polp.Integral()
    #norm_polp = 10.
    print 'norm_polp',norm_polp

    for dir_ind in directory:
        x_polp80_BDTScore.append(float(dir_ind))
        h_mass_polp80_sig_hzqq=file_polp80_HZ_SignalHistos_.Get(dir_ind+"/h_jet1_mass")
        h_mass_polp80_sig_hzqq_AllEvents=file_polp80_HZ_SignalHistos_AllEvents_.Get(dir_ind+"/h_jet1_mass")
        h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get(dir_ind+"/h_jet1_mass")
        h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get(dir_ind+"/h_jet1_mass")
        h_mass_polp80_sig_ee_qqqqqq_BGHistos=file_polp80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_jet1_mass") 
        #print dir_ind,"integral of signal mass",h_mass_polp80_sig_hzqq.Integral(),h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_mass_polp80_sig_hzqq.Integral()>0:
            #print '#',dir_ind,"polp significance,purity ",h_mass_polp80_sig_hzqq.Integral()/sqrt(h_mass_polp80_sig_hzqq.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polp80_sig_hzqq.Integral()/(h_mass_polp80_sig_hzqq.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()),h_mass_polp80_sig_hzqq.Integral(),h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polp80_efficiency.append(100.*h_mass_polp80_sig_hzqq.Integral()/norm_polp)
            y_polp80_purity.append(100.*h_mass_polp80_sig_hzqq.Integral()/(h_mass_polp80_sig_hzqq.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()))
            y_polp80_significance.append(h_mass_polp80_sig_hzqq.Integral()/sqrt(h_mass_polp80_sig_hzqq.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()))
            print 'at polp ',dir_ind,'signif/pur/eff/events',h_mass_polp80_sig_hzqq.Integral()/sqrt(h_mass_polp80_sig_hzqq.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polp80_sig_hzqq.Integral()/(h_mass_polp80_sig_hzqq.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polp80_sig_hzqq.Integral()/norm_polp,h_mass_polp80_sig_hzqq.Integral()
        else:
            y_polp80_efficiency.append(0)
            y_polp80_purity.append(0)
            y_polp80_significance.append(0)
      #print dir_ind,"integral of signal mass",h_mass_polp80_sig_hzqq.Integral(),h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()
        if h_mass_polp80_sig_hzqq_AllEvents.Integral()>0:
            #print '#',dir_ind,"polp significance,purity ",h_mass_polp80_sig_hzqq.Integral()/sqrt(h_mass_polp80_sig_hzqq.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polp80_sig_hzqq.Integral()/(h_mass_polp80_sig_hzqq.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()),h_mass_polp80_sig_hzqq.Integral(),h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()
            y_polp80_AllEvents_efficiency.append(100.*h_mass_polp80_sig_hzqq_AllEvents.Integral()/norm_polp_AllEvents)
            y_polp80_AllEvents_purity.append(100.*h_mass_polp80_sig_hzqq_AllEvents.Integral()/(h_mass_polp80_sig_hzqq_AllEvents.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()))
            y_polp80_AllEvents_significance.append(h_mass_polp80_sig_hzqq_AllEvents.Integral()/sqrt(h_mass_polp80_sig_hzqq_AllEvents.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()))
            print 'at polp ',dir_ind,'signif/pur/eff/events all ',h_mass_polp80_sig_hzqq_AllEvents.Integral()/sqrt(h_mass_polp80_sig_hzqq_AllEvents.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polp80_sig_hzqq_AllEvents.Integral()/(h_mass_polp80_sig_hzqq_AllEvents.Integral()+h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_polp80_sig_hzqq_AllEvents.Integral()/norm_polp,h_mass_polp80_sig_hzqq_AllEvents.Integral()
            print 'at polp ',dir_ind,'evt sig/sig all/qq/qqqq/qqqqqq',h_mass_polp80_sig_hzqq.Integral(),h_mass_polp80_sig_hzqq_AllEvents.Integral(),h_mass_polp80_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_polp80_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_polp80_sig_ee_qqqqqq_BGHistos.Integral()
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
    l.DrawLatex(x2,y2,label3);
    canvas_polp80_BDT_significance.Print("hzqq_bbar_vs_totBG_polp80_sigificance_vs_BDT.eps")
    
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
    canvas_polp80_BDT_efficiency.Print("hzqq_bbar_vs_totBG_polp80_signal_efficiency_vs_BDT.eps")
    canvas_polp80_BDT_purity = setUpperCanvas("canvas_polp80_BDT_purity");
    canvas_polp80_BDT_purity.cd()
    graph_polp80_purity = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_purity )
    graph_polp80_purity.SetMarkerStyle(22)
    graph_polp80_purity.SetMarkerColor(2)
    graph_polp80_purity.GetXaxis().SetTitle('BDT Score')
    graph_polp80_purity.GetYaxis().SetTitle('purity [%]')
    graph_polp80_purity.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_polp80_BDT_purity.Print("hzqq_bbar_vs_totBG_polp80_signal_purity_vs_BDT.eps")




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
    l.DrawLatex(x2,y2,label3);
    canvas_polp80_AllEvents_BDT_significance.Print("hzqq_bbar_vs_totBG_polp80_AllEvents_sigificance_vs_BDT.eps")
    
    canvas_polp80_AllEvents_BDT_efficiency = setUpperCanvas("canvas_polp80_AllEvents_BDT_efficiency");
    canvas_polp80_AllEvents_BDT_efficiency.cd()
    graph_polp80_AllEvents_efficiency = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_AllEvents_efficiency )
    graph_polp80_AllEvents_efficiency.SetMarkerStyle(21)
    graph_polp80_AllEvents_efficiency.SetMarkerColor(1)
    graph_polp80_AllEvents_efficiency.GetXaxis().SetTitle('BDT Score')
    graph_polp80_AllEvents_efficiency.GetYaxis().SetTitle('signal efficiency [%]')
    graph_polp80_AllEvents_efficiency.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_polp80_AllEvents_BDT_efficiency.Print("hzqq_bbar_vs_totBG_polp80_AllEvents_signal_efficiency_vs_BDT.eps")
    canvas_polp80_AllEvents_BDT_purity = setUpperCanvas("canvas_polp80_AllEvents_BDT_purity");
    canvas_polp80_AllEvents_BDT_purity.cd()
    graph_polp80_AllEvents_purity = TGraph( n_graphs, x_polp80_BDTScore, y_polp80_AllEvents_purity )
    graph_polp80_AllEvents_purity.SetMarkerStyle(22)
    graph_polp80_AllEvents_purity.SetMarkerColor(2)
    graph_polp80_AllEvents_purity.GetXaxis().SetTitle('BDT Score')
    graph_polp80_AllEvents_purity.GetYaxis().SetTitle('purity [%]')
    graph_polp80_AllEvents_purity.Draw('AP')
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_polp80_AllEvents_BDT_purity.Print("hzqq_bbar_vs_totBG_polp80_AllEvents_signal_purity_vs_BDT.eps")



    #print 'A_theta1_parton',(h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral()+h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral())/(abs(h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral())+abs(h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral()))
    #print 'A_phi1_parton',(h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()+h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())/(abs(h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())+abs(h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()))
    #print 'A_phi2_parton',(h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()+h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())/(abs(h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())+abs(h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()))
    #print 'A_phi3_parton',(h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()+h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())/(abs(h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())+abs(h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()))
    #print 'A_phi4_parton',(h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()+h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())/(abs(h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())+abs(h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()))


    h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");
    h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetXaxis().SetTitle("#phi [#circ]");
    h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("Events");
    h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineWidth(3);
    h_ee_qq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");
    h_ee_qq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetXaxis().SetTitle("#phi [#circ]");
    h_ee_qq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("Events");
    h_ee_qq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetFillColor(kBlue);
    h_ee_qq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineColor(kBlue);
    h_ee_qqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");
    h_ee_qqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetXaxis().SetTitle("#phi [#circ]");
    h_ee_qqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("Events");
    h_ee_qqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetFillColor(kRed);
    h_ee_qqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineColor(kRed);
    h_ee_qqqqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");
    h_ee_qqqqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetXaxis().SetTitle("#phi [#circ]");
    h_ee_qqqqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("Events");
    h_ee_qqqqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetFillColor(kGreen-2);
    h_ee_qqqqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineColor(kGreen-2);
 
    h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Rebin(4)
    h_ee_qq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Rebin(4)
    h_ee_qqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Rebin(4)
    h_ee_qqqqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Rebin(4)

    hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep= THStack("hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep", "");
    hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_ee_qqqqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep);
    hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_ee_qqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep);
    hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_ee_qq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep);
    hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep);
    
    canvas_h_BDT_Signalpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack = setUpperCanvas("canvas_h_BDT_Signalpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack");
    canvas_h_BDT_Signalpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack.cd();
    hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Draw("hist");
    hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetXaxis().SetTitle("#phi [#circ]");
    hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("Events");
    hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetMaximum(120)
    canvas_h_BDT_Signalpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack.Modified();
    
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=TLegend(0.25,0.65,0.65,0.90);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetBorderSize(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetTextAlign(12);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetTextSize(0.050);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetTextFont(42);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetMargin(0.15);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineColor(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineStyle(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineWidth(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetFillColor(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetFillStyle(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.AddEntry(h_ee_qqqqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep,"ee#rightarrow qqqqqq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.AddEntry(h_ee_qqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep,"ee#rightarrow qqqq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.AddEntry(h_ee_qq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep,"ee#rightarrow qq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.AddEntry(h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep,"HZ");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);

    canvas_h_BDT_Signalpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack.Print("hzqq_and_totBG_polm80_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.eps")


    canvas_h_BG_over_Signal_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack = setUpperCanvas("canvas_h_BG_over_Signal_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack");
    canvas_h_BG_over_Signal_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack.cd();
    h_total_BG_polm80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=h_ee_qq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Clone("h_total_BG_polm80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_total_BG_polm80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_ee_qqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
    h_total_BG_polm80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_ee_qqqqqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
    hzqq_sig_BG_Stackpolm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Draw("hist");
    h_total_BG_polm80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("totBG/sig");
    h_total_BG_polm80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineColor(kBlack)
    h_total_BG_polm80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Divide(h_hzqq_BDT_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
    h_total_BG_polm80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Draw("PE")
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);

    canvas_h_BG_over_Signal_polm80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack.Print("totBG_vs_sig_hzqq_polm80_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.eps")


    h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");
    h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetXaxis().SetTitle("#phi [#circ]");
    h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("Events");
    h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineWidth(3);
    h_ee_qq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");
    h_ee_qq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetXaxis().SetTitle("#phi [#circ]");
    h_ee_qq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("Events");
    h_ee_qq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetFillColor(kBlue);
    h_ee_qq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineColor(kBlue);
    h_ee_qqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");
    h_ee_qqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetXaxis().SetTitle("#phi [#circ]");
    h_ee_qqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("Events");
    h_ee_qqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetFillColor(kRed);
    h_ee_qqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineColor(kRed);
    h_ee_qqqqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");
    h_ee_qqqqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetXaxis().SetTitle("#phi [#circ]");
    h_ee_qqqqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("Events");
    h_ee_qqqqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetFillColor(kGreen-2);
    h_ee_qqqqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineColor(kGreen-2);
 
    h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Rebin(4)
    h_ee_qq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Rebin(4)
    h_ee_qqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Rebin(4)
    h_ee_qqqqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Rebin(4)

    hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep= THStack("hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep", "");
    hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_ee_qqqqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep);
    hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_ee_qqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep);
    hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_ee_qq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep);
    hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep);
    
    canvas_h_BDT_Signalpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack = setUpperCanvas("canvas_h_BDT_Signalpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack");
    canvas_h_BDT_Signalpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack.cd();
    hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Draw("hist");
    hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetXaxis().SetTitle("#phi [#circ]");
    hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("Events");
    hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetMaximum(22)
    canvas_h_BDT_Signalpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack.Modified();
    
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=TLegend(0.25,0.65,0.65,0.90);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetBorderSize(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetTextAlign(12);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetTextSize(0.050);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetTextFont(42);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetMargin(0.15);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineColor(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineStyle(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineWidth(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetFillColor(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetFillStyle(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.AddEntry(h_ee_qqqqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep,"ee#rightarrow qqqqqq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.AddEntry(h_ee_qqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep,"ee#rightarrow qqqq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.AddEntry(h_ee_qq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep,"ee#rightarrow qq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.AddEntry(h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep,"HZ");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_BDT_Signalpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack.Print("hzqq_and_totBG_polp80_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.eps")

    
    canvas_h_BG_over_Signal_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack = setUpperCanvas("canvas_h_BG_over_Signal_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack");
    canvas_h_BG_over_Signal_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack.cd();
    h_total_BG_polp80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep=h_ee_qq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Clone("h_total_BG_polp80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep")
    h_total_BG_polp80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_ee_qqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
    h_total_BG_polp80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Add(h_ee_qqqqqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
    hzqq_sig_BG_Stackpolp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Draw("hist");
    h_total_BG_polp80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.GetYaxis().SetTitle("totBG/sig");
    h_total_BG_polp80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.SetLineColor(kBlack)
    h_total_BG_polp80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Divide(h_hzqq_BDT_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
    h_total_BG_polp80_recojetc_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Draw("PE")
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_BG_over_Signal_polp80_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_thstack.Print("totBG_vs_sig_hzqq_polp80_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.eps")


    h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("#theta_{1} [#circ]");
    h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineWidth(3);
    h_ee_qq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("#theta_{1} [#circ]");
    h_ee_qq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kBlue);
    h_ee_qq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kBlue);
    h_ee_qqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("#theta_{1} [#circ]");
    h_ee_qqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kRed);
    h_ee_qqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kRed);
    h_ee_qqqqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qqqqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("#theta_{1} [#circ]");
    h_ee_qqqqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qqqqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kGreen-2);
    h_ee_qqqqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kGreen-2);
 
    h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qqqqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)

    hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com= THStack("hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com", "");
    hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    
    canvas_h_BDT_Signalpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack = setUpperCanvas("canvas_h_BDT_Signalpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack");
    canvas_h_BDT_Signalpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.cd();
    hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("hist");
    hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("#theta_{1} [#circ]");
    hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMaximum(165)
    canvas_h_BDT_Signalpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Modified();
    
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=TLegend(0.20,0.65,0.60,0.90);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetBorderSize(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextAlign(12);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextSize(0.050);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextFont(42);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMargin(0.15);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineStyle(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineWidth(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillStyle(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qqqqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qqqqqq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qqqq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"HZ");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);

    canvas_h_BDT_Signalpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Print("hzqq_and_totBG_polm80_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.eps")


    canvas_h_BG_over_Signal_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack = setUpperCanvas("canvas_h_BG_over_Signal_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack");
    canvas_h_BG_over_Signal_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.cd();
    h_total_BG_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=h_ee_qq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Clone("h_total_BG_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_total_BG_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    h_total_BG_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqqqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    hzqq_sig_BG_Stackpolm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("hist");
    h_total_BG_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("totBG/sig");
    h_total_BG_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMaximum(1.5)
    h_total_BG_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kBlack)
    h_total_BG_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Divide(h_hzqq_BDT_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    h_total_BG_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("PE")
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);

    canvas_h_BG_over_Signal_polm80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Print("totBG_vs_sig_hzqq_polm80_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.eps")



    canvas_h_BDT_Signalpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack = setUpperCanvas("canvas_hzqq_and_totBG_polp80_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.eps")
    canvas_h_BDT_Signalpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.cd()
    h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("#theta_{1} [#circ]");
    h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineWidth(3);
    h_ee_qq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("#theta_{1} [#circ]");
    h_ee_qq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kBlue);
    h_ee_qq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kBlue);
    h_ee_qqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("#theta_{1} [#circ]");
    h_ee_qqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kRed);
    h_ee_qqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kRed);
    h_ee_qqqqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qqqqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("#theta_{1} [#circ]");
    h_ee_qqqqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qqqqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kGreen-2);
    h_ee_qqqqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kGreen-2);
 
    h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qqqqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)

    hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com= THStack("hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com", "");
    hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    
    canvas_h_BDT_Signalpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack = setUpperCanvas("canvas_h_BDT_Signalpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack");
    canvas_h_BDT_Signalpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.cd();
    hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("hist");
    hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("#theta_{1} [#circ]");
    hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMaximum(27)
    canvas_h_BDT_Signalpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Modified();
    
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=TLegend(0.20,0.65,0.60,0.90);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetBorderSize(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextAlign(12);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextSize(0.050);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextFont(42);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMargin(0.15);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineStyle(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineWidth(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillStyle(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qqqqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qqqqqq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qqqq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"HZ");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_BDT_Signalpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Print("hzqq_and_totBG_polp80_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.eps")

    canvas_h_BG_over_Signal_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack = setUpperCanvas("canvas_h_BG_over_Signal_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack");
    canvas_h_BG_over_Signal_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.cd();
    h_total_BG_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com=h_ee_qq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Clone("h_total_BG_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_total_BG_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    h_total_BG_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqqqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    hzqq_sig_BG_Stackpolp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("hist");
    h_total_BG_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("totBG/sig");
    h_total_BG_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMaximum(1.5)
    h_total_BG_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kBlack)
    h_total_BG_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Divide(h_hzqq_BDT_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    h_total_BG_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("PE")
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_BG_over_Signal_polp80_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Print("totBG_vs_sig_hzqq_polp80_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.eps")

    h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("cos#theta_{1}");
    h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineWidth(3);
    h_ee_qq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("cos#theta_{1}");
    h_ee_qq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kBlue);
    h_ee_qq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kBlue);
    h_ee_qqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("cos#theta_{1}");
    h_ee_qqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kRed);
    h_ee_qqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kRed);
    h_ee_qqqqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qqqqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("cos#theta_{1}");
    h_ee_qqqqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qqqqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kGreen-2);
    h_ee_qqqqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kGreen-2);
 
    h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qqqqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)

    hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com= THStack("hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com", "");
    hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    
    canvas_h_BDT_Signalpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack = setUpperCanvas("canvas_h_BDT_Signalpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack");
    canvas_h_BDT_Signalpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.cd();
    hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("hist");
    hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("cos#theta_{1}");
    hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMaximum(150)
    canvas_h_BDT_Signalpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Modified();
    
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=TLegend(0.20,0.65,0.60,0.90);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetBorderSize(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextAlign(12);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextSize(0.050);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextFont(42);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMargin(0.15);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineStyle(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineWidth(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillStyle(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qqqqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qqqqqq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qqqq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"HZ");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);

    canvas_h_BDT_Signalpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Print("hzqq_and_totBG_polm80_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.eps")

    canvas_h_BG_over_Signal_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack = setUpperCanvas("canvas_h_BG_over_Signal_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack");
    canvas_h_BG_over_Signal_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.cd();
    h_total_BG_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=h_ee_qq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Clone("h_total_BG_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_total_BG_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    h_total_BG_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqqqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    hzqq_sig_BG_Stackpolm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("hist");
    h_total_BG_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("totBG/sig");
    h_total_BG_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kBlack)
    h_total_BG_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMaximum(1.5)
    h_total_BG_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Divide(h_hzqq_BDT_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    h_total_BG_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("PE")
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);

    canvas_h_BG_over_Signal_polm80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Print("totBG_vs_sig_hzqq_polm80_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.eps")


    h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("cos#theta_{1}");
    h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineWidth(3);
    h_ee_qq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("cos#theta_{1}");
    h_ee_qq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kBlue);
    h_ee_qq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kBlue);
    h_ee_qqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("cos#theta_{1}");
    h_ee_qqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kRed);
    h_ee_qqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kRed);
    h_ee_qqqqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_ee_qqqqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("cos#theta_{1}");
    h_ee_qqqqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    h_ee_qqqqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(kGreen-2);
    h_ee_qqqqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kGreen-2);
 
    h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)
    h_ee_qqqqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Rebin(4)

    hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com= THStack("hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com", "");
    hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com);
    
    canvas_h_BDT_Signalpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack = setUpperCanvas("canvas_h_BDT_Signalpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack");
    canvas_h_BDT_Signalpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.cd();
    hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("hist");
    hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetXaxis().SetTitle("cos#theta_{1}");
    hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");
    hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMaximum(27)
    canvas_h_BDT_Signalpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Modified();
    
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=TLegend(0.20,0.65,0.60,0.90);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetBorderSize(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextAlign(12);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextSize(0.050);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetTextFont(42);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMargin(0.15);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineStyle(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineWidth(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillColor(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetFillStyle(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qqqqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qqqqqq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qqqq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_ee_qq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"ee#rightarrow qq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.AddEntry(h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com,"HZ");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_BDT_Signalpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Print("hzqq_and_totBG_polp80_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.eps")

    canvas_h_BG_over_Signal_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack = setUpperCanvas("canvas_h_BG_over_Signal_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack");
    canvas_h_BG_over_Signal_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.cd();
    h_total_BG_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com=h_ee_qq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Clone("h_total_BG_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com")
    h_total_BG_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    h_total_BG_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Add(h_ee_qqqqqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    hzqq_sig_BG_Stackpolp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("hist");
    h_total_BG_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("totBG/sig");
    h_total_BG_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetLineColor(kBlack)
    h_total_BG_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.SetMaximum(1.5)
    h_total_BG_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Divide(h_hzqq_BDT_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
    h_total_BG_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Draw("PE")
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_BG_over_Signal_polp80_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.Print("totBG_vs_sig_hzqq_polp80_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.eps")
























    h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_theta2_rj1_ep_E_totCOM");
    h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("#theta_{2} [#circ]");
    h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetLineWidth(3);
    h_ee_qq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_theta2_rj1_ep_E_totCOM");
    h_ee_qq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("#theta_{2} [#circ]");
    h_ee_qq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetFillColor(kBlue);
    h_ee_qq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetLineColor(kBlue);
    h_ee_qqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_theta2_rj1_ep_E_totCOM");
    h_ee_qqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("#theta_{2} [#circ]");
    h_ee_qqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetFillColor(kRed);
    h_ee_qqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetLineColor(kRed);
    h_ee_qqqqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_theta2_rj1_ep_E_totCOM");
    h_ee_qqqqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("#theta_{2} [#circ]");
    h_ee_qqqqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qqqqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetFillColor(kGreen-2);
    h_ee_qqqqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetLineColor(kGreen-2);
 
    h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qqqqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.Rebin(4)

    hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM= THStack("hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM", "");
    hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM.Add(h_ee_qqqqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM.Add(h_ee_qqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM.Add(h_ee_qq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM.Add(h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM);
    
    canvas_h_BDT_Signalpolm80_recojet_theta2_rj1_ep_E_totCOM_thstack = setUpperCanvas("canvas_h_BDT_Signalpolm80_recojet_theta2_rj1_ep_E_totCOM_thstack");
    canvas_h_BDT_Signalpolm80_recojet_theta2_rj1_ep_E_totCOM_thstack.cd();
    hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM.Draw("hist");
    hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("#theta_{2} [#circ]");
    hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM.SetMaximum(165)
    canvas_h_BDT_Signalpolm80_recojet_theta2_rj1_ep_E_totCOM_thstack.Modified();
    
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM=TLegend(0.20,0.65,0.60,0.90);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetBorderSize(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetTextAlign(12);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetTextSize(0.050);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetTextFont(42);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetMargin(0.15);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetLineColor(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetLineStyle(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetLineWidth(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetFillColor(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetFillStyle(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.AddEntry(h_ee_qqqqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM,"ee#rightarrow qqqqqq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.AddEntry(h_ee_qqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM,"ee#rightarrow qqqq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.AddEntry(h_ee_qq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM,"ee#rightarrow qq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.AddEntry(h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM,"HZ");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);

    canvas_h_BDT_Signalpolm80_recojet_theta2_rj1_ep_E_totCOM_thstack.Print("hzqq_and_totBG_polm80_theta2_rj1_ep_E_totCOM.eps")


    canvas_h_BG_over_Signal_polm80_recojet_theta2_rj1_ep_E_totCOM_thstack = setUpperCanvas("canvas_h_BG_over_Signal_polm80_recojet_theta2_rj1_ep_E_totCOM_thstack");
    canvas_h_BG_over_Signal_polm80_recojet_theta2_rj1_ep_E_totCOM_thstack.cd();
    h_total_BG_polm80_recojet_theta2_rj1_ep_E_totCOM=h_ee_qq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM.Clone("h_total_BG_polm80_recojet_theta2_rj1_ep_E_totCOM")
    h_total_BG_polm80_recojet_theta2_rj1_ep_E_totCOM.Add(h_ee_qqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM)
    h_total_BG_polm80_recojet_theta2_rj1_ep_E_totCOM.Add(h_ee_qqqqqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM)
    hzqq_sig_BG_Stackpolm80_recojet_theta2_rj1_ep_E_totCOM.Draw("hist");
    h_total_BG_polm80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("totBG/sig");
    h_total_BG_polm80_recojet_theta2_rj1_ep_E_totCOM.SetMaximum(1.5)
    h_total_BG_polm80_recojet_theta2_rj1_ep_E_totCOM.SetLineColor(kBlack)
    h_total_BG_polm80_recojet_theta2_rj1_ep_E_totCOM.Divide(h_hzqq_BDT_polm80_recojet_theta2_rj1_ep_E_totCOM)
    h_total_BG_polm80_recojet_theta2_rj1_ep_E_totCOM.Draw("PE")
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);

    canvas_h_BG_over_Signal_polm80_recojet_theta2_rj1_ep_E_totCOM_thstack.Print("totBG_vs_sig_hzqq_polm80_theta2_rj1_ep_E_totCOM.eps")



    canvas_h_BDT_Signalpolp80_recojet_theta2_rj1_ep_E_totCOM_thstack = setUpperCanvas("canvas_hzqq_and_totBG_polp80_theta2_rj1_ep_E_totCOM.eps")
    canvas_h_BDT_Signalpolp80_recojet_theta2_rj1_ep_E_totCOM_thstack.cd()
    h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_theta2_rj1_ep_E_totCOM");
    h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("#theta_{2} [#circ]");
    h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetLineWidth(3);
    h_ee_qq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_theta2_rj1_ep_E_totCOM");
    h_ee_qq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("#theta_{2} [#circ]");
    h_ee_qq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetFillColor(kBlue);
    h_ee_qq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetLineColor(kBlue);
    h_ee_qqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_theta2_rj1_ep_E_totCOM");
    h_ee_qqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("#theta_{2} [#circ]");
    h_ee_qqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetFillColor(kRed);
    h_ee_qqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetLineColor(kRed);
    h_ee_qqqqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_theta2_rj1_ep_E_totCOM");
    h_ee_qqqqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("#theta_{2} [#circ]");
    h_ee_qqqqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qqqqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetFillColor(kGreen-2);
    h_ee_qqqqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetLineColor(kGreen-2);
 
    h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qqqqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.Rebin(4)

    hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM= THStack("hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM", "");
    hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM.Add(h_ee_qqqqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM.Add(h_ee_qqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM.Add(h_ee_qq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM.Add(h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM);
    
    canvas_h_BDT_Signalpolp80_recojet_theta2_rj1_ep_E_totCOM_thstack = setUpperCanvas("canvas_h_BDT_Signalpolp80_recojet_theta2_rj1_ep_E_totCOM_thstack");
    canvas_h_BDT_Signalpolp80_recojet_theta2_rj1_ep_E_totCOM_thstack.cd();
    hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM.Draw("hist");
    hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("#theta_{2} [#circ]");
    hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM.SetMaximum(27)
    canvas_h_BDT_Signalpolp80_recojet_theta2_rj1_ep_E_totCOM_thstack.Modified();
    
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM=TLegend(0.20,0.65,0.60,0.90);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetBorderSize(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetTextAlign(12);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetTextSize(0.050);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetTextFont(42);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetMargin(0.15);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetLineColor(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetLineStyle(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetLineWidth(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetFillColor(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetFillStyle(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.AddEntry(h_ee_qqqqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM,"ee#rightarrow qqqqqq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.AddEntry(h_ee_qqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM,"ee#rightarrow qqqq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.AddEntry(h_ee_qq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM,"ee#rightarrow qq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.AddEntry(h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM,"HZ");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_BDT_Signalpolp80_recojet_theta2_rj1_ep_E_totCOM_thstack.Print("hzqq_and_totBG_polp80_theta2_rj1_ep_E_totCOM.eps")

    canvas_h_BG_over_Signal_polp80_recojet_theta2_rj1_ep_E_totCOM_thstack = setUpperCanvas("canvas_h_BG_over_Signal_polp80_recojet_theta2_rj1_ep_E_totCOM_thstack");
    canvas_h_BG_over_Signal_polp80_recojet_theta2_rj1_ep_E_totCOM_thstack.cd();
    h_total_BG_polp80_recojet_theta2_rj1_ep_E_totCOM=h_ee_qq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM.Clone("h_total_BG_polp80_recojet_theta2_rj1_ep_E_totCOM")
    h_total_BG_polp80_recojet_theta2_rj1_ep_E_totCOM.Add(h_ee_qqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM)
    h_total_BG_polp80_recojet_theta2_rj1_ep_E_totCOM.Add(h_ee_qqqqqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM)
    hzqq_sig_BG_Stackpolp80_recojet_theta2_rj1_ep_E_totCOM.Draw("hist");
    h_total_BG_polp80_recojet_theta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("totBG/sig");
    h_total_BG_polp80_recojet_theta2_rj1_ep_E_totCOM.SetMaximum(1.5)
    h_total_BG_polp80_recojet_theta2_rj1_ep_E_totCOM.SetLineColor(kBlack)
    h_total_BG_polp80_recojet_theta2_rj1_ep_E_totCOM.Divide(h_hzqq_BDT_polp80_recojet_theta2_rj1_ep_E_totCOM)
    h_total_BG_polp80_recojet_theta2_rj1_ep_E_totCOM.Draw("PE")
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_BG_over_Signal_polp80_recojet_theta2_rj1_ep_E_totCOM_thstack.Print("totBG_vs_sig_hzqq_polp80_theta2_rj1_ep_E_totCOM.eps")

    h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM=file_polm80_HZ_SignalHistos_AllEvents_.Get("0.375/h_recojet_costheta2_rj1_ep_E_totCOM");
    h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("cos#theta_{2}");
    h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetLineWidth(3);
    h_ee_qq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get("0.375/h_recojet_costheta2_rj1_ep_E_totCOM");
    h_ee_qq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("cos#theta_{2}");
    h_ee_qq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetFillColor(kBlue);
    h_ee_qq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetLineColor(kBlue);
    h_ee_qqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.375/h_recojet_costheta2_rj1_ep_E_totCOM");
    h_ee_qqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("cos#theta_{2}");
    h_ee_qqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetFillColor(kRed);
    h_ee_qqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetLineColor(kRed);
    h_ee_qqqqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM=file_polm80_ee_qqqqqq_BGHistos_.Get("0.375/h_recojet_costheta2_rj1_ep_E_totCOM");
    h_ee_qqqqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("cos#theta_{2}");
    h_ee_qqqqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qqqqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetFillColor(kGreen-2);
    h_ee_qqqqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetLineColor(kGreen-2);
 
    h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qqqqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.Rebin(4)

    hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM= THStack("hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM", "");
    hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_ee_qqqqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_ee_qqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_ee_qq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM);
    
    canvas_h_BDT_Signalpolm80_recojet_costheta2_rj1_ep_E_totCOM_thstack = setUpperCanvas("canvas_h_BDT_Signalpolm80_recojet_costheta2_rj1_ep_E_totCOM_thstack");
    canvas_h_BDT_Signalpolm80_recojet_costheta2_rj1_ep_E_totCOM_thstack.cd();
    hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM.Draw("hist");
    hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("cos#theta_{2}");
    hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM.SetMaximum(150)
    canvas_h_BDT_Signalpolm80_recojet_costheta2_rj1_ep_E_totCOM_thstack.Modified();
    
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM=TLegend(0.20,0.65,0.60,0.90);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetBorderSize(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetTextAlign(12);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetTextSize(0.050);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetTextFont(42);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetMargin(0.15);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetLineColor(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetLineStyle(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetLineWidth(1);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetFillColor(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetFillStyle(0);
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.AddEntry(h_ee_qqqqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM,"ee#rightarrow qqqqqq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.AddEntry(h_ee_qqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM,"ee#rightarrow qqqq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.AddEntry(h_ee_qq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM,"ee#rightarrow qq");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.AddEntry(h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM,"HZ");
    leg_h_BDT_h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);

    canvas_h_BDT_Signalpolm80_recojet_costheta2_rj1_ep_E_totCOM_thstack.Print("hzqq_and_totBG_polm80_costheta2_rj1_ep_E_totCOM.eps")

    canvas_h_BG_over_Signal_polm80_recojet_costheta2_rj1_ep_E_totCOM_thstack = setUpperCanvas("canvas_h_BG_over_Signal_polm80_recojet_costheta2_rj1_ep_E_totCOM_thstack");
    canvas_h_BG_over_Signal_polm80_recojet_costheta2_rj1_ep_E_totCOM_thstack.cd();
    h_total_BG_polm80_recojet_costheta2_rj1_ep_E_totCOM=h_ee_qq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM.Clone("h_total_BG_polm80_recojet_costheta2_rj1_ep_E_totCOM")
    h_total_BG_polm80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_ee_qqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM)
    h_total_BG_polm80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_ee_qqqqqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM)
    hzqq_sig_BG_Stackpolm80_recojet_costheta2_rj1_ep_E_totCOM.Draw("hist");
    h_total_BG_polm80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("totBG/sig");
    h_total_BG_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetLineColor(kBlack)
    h_total_BG_polm80_recojet_costheta2_rj1_ep_E_totCOM.SetMaximum(1.5)
    h_total_BG_polm80_recojet_costheta2_rj1_ep_E_totCOM.Divide(h_hzqq_BDT_polm80_recojet_costheta2_rj1_ep_E_totCOM)
    h_total_BG_polm80_recojet_costheta2_rj1_ep_E_totCOM.Draw("PE")
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);

    canvas_h_BG_over_Signal_polm80_recojet_costheta2_rj1_ep_E_totCOM_thstack.Print("totBG_vs_sig_hzqq_polm80_costheta2_rj1_ep_E_totCOM.eps")


    h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM=file_polp80_HZ_SignalHistos_AllEvents_.Get("0.3/h_recojet_costheta2_rj1_ep_E_totCOM");
    h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("cos#theta_{2}");
    h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetLineWidth(3);
    h_ee_qq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get("0.3/h_recojet_costheta2_rj1_ep_E_totCOM");
    h_ee_qq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("cos#theta_{2}");
    h_ee_qq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetFillColor(kBlue);
    h_ee_qq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetLineColor(kBlue);
    h_ee_qqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get("0.3/h_recojet_costheta2_rj1_ep_E_totCOM");
    h_ee_qqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("cos#theta_{2}");
    h_ee_qqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetFillColor(kRed);
    h_ee_qqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetLineColor(kRed);
    h_ee_qqqqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM=file_polp80_ee_qqqqqq_BGHistos_.Get("0.3/h_recojet_costheta2_rj1_ep_E_totCOM");
    h_ee_qqqqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("cos#theta_{2}");
    h_ee_qqqqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    h_ee_qqqqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetFillColor(kGreen-2);
    h_ee_qqqqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetLineColor(kGreen-2);
 
    h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.Rebin(4)
    h_ee_qqqqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.Rebin(4)

    hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM= THStack("hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM", "");
    hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_ee_qqqqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_ee_qqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_ee_qq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM);
    hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM);
    
    canvas_h_BDT_Signalpolp80_recojet_costheta2_rj1_ep_E_totCOM_thstack = setUpperCanvas("canvas_h_BDT_Signalpolp80_recojet_costheta2_rj1_ep_E_totCOM_thstack");
    canvas_h_BDT_Signalpolp80_recojet_costheta2_rj1_ep_E_totCOM_thstack.cd();
    hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM.Draw("hist");
    hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM.GetXaxis().SetTitle("cos#theta_{2}");
    hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("Events");
    hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM.SetMaximum(27)
    canvas_h_BDT_Signalpolp80_recojet_costheta2_rj1_ep_E_totCOM_thstack.Modified();
    
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM=TLegend(0.20,0.65,0.60,0.90);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetBorderSize(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetTextAlign(12);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetTextSize(0.050);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetTextFont(42);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetMargin(0.15);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetLineColor(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetLineStyle(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetLineWidth(1);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetFillColor(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetFillStyle(0);
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.AddEntry(h_ee_qqqqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM,"ee#rightarrow qqqqqq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.AddEntry(h_ee_qqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM,"ee#rightarrow qqqq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.AddEntry(h_ee_qq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM,"ee#rightarrow qq");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.AddEntry(h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM,"HZ");
    leg_h_BDT_h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_BDT_Signalpolp80_recojet_costheta2_rj1_ep_E_totCOM_thstack.Print("hzqq_and_totBG_polp80_costheta2_rj1_ep_E_totCOM.eps")

    canvas_h_BG_over_Signal_polp80_recojet_costheta2_rj1_ep_E_totCOM_thstack = setUpperCanvas("canvas_h_BG_over_Signal_polp80_recojet_costheta2_rj1_ep_E_totCOM_thstack");
    canvas_h_BG_over_Signal_polp80_recojet_costheta2_rj1_ep_E_totCOM_thstack.cd();
    h_total_BG_polp80_recojet_costheta2_rj1_ep_E_totCOM=h_ee_qq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM.Clone("h_total_BG_polp80_recojet_costheta2_rj1_ep_E_totCOM")
    h_total_BG_polp80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_ee_qqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM)
    h_total_BG_polp80_recojet_costheta2_rj1_ep_E_totCOM.Add(h_ee_qqqqqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM)
    hzqq_sig_BG_Stackpolp80_recojet_costheta2_rj1_ep_E_totCOM.Draw("hist");
    h_total_BG_polp80_recojet_costheta2_rj1_ep_E_totCOM.GetYaxis().SetTitle("totBG/sig");
    h_total_BG_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetLineColor(kBlack)
    h_total_BG_polp80_recojet_costheta2_rj1_ep_E_totCOM.SetMaximum(1.5)
    h_total_BG_polp80_recojet_costheta2_rj1_ep_E_totCOM.Divide(h_hzqq_BDT_polp80_recojet_costheta2_rj1_ep_E_totCOM)
    h_total_BG_polp80_recojet_costheta2_rj1_ep_E_totCOM.Draw("PE")
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_BG_over_Signal_polp80_recojet_costheta2_rj1_ep_E_totCOM_thstack.Print("totBG_vs_sig_hzqq_polp80_costheta2_rj1_ep_E_totCOM.eps")


    
    file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polm80=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingWeights_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees300NCuts_m1__hzqq__ee_qq_mqq_1TeV__ee_qqqq_mqqqq_2TeV__ee_qqqqqq_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root")
    file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polp80=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingWeights_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees250NCuts_m1__hzqq__ee_qq_mqq_1TeV__ee_qqqq_mqqqq_2TeV__ee_qqqqqq_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21.root")
    
    h_BDT_Signal_polm80_test=file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polm80.Get("dataset/Method_BDT/BDT/MVA_BDT_S");
    h_BDT_Signal_polm80_test.SetFillColor(kBlue-2);
    h_BDT_Signal_polm80_test.SetFillStyle(3001);
    h_BDT_Signal_polm80_test.GetXaxis().SetTitle("BDT Score");
    h_BDT_Signal_polm80_test.GetYaxis().SetTitle("Entries");
    
    h_BDT_Signal_polm80_train=file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polm80.Get("dataset/Method_BDT/BDT/MVA_BDT_Train_S");
    h_BDT_Signal_polm80_train.SetLineColor(kBlue);
    h_BDT_Signal_polm80_train.SetLineWidth(2);
    h_BDT_Signal_polm80_train.GetXaxis().SetTitle("BDT Score");
    h_BDT_Signal_polm80_train.GetYaxis().SetTitle("Entries");
    
    h_BDT_BG_polm80_test=file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polm80.Get("dataset/Method_BDT/BDT/MVA_BDT_B");
    h_BDT_BG_polm80_test.SetFillColor(kRed-2);
    h_BDT_BG_polm80_test.SetFillStyle(3002);
    h_BDT_BG_polm80_test.GetXaxis().SetTitle("BDT Score");
    h_BDT_BG_polm80_test.GetYaxis().SetTitle("Entries");
    
    h_BDT_BG_polm80_train=file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polm80.Get("dataset/Method_BDT/BDT/MVA_BDT_Train_B");
    h_BDT_BG_polm80_train.SetLineColor(kRed);
    h_BDT_BG_polm80_train.SetLineWidth(2);
    h_BDT_BG_polm80_train.GetXaxis().SetTitle("BDT Score");
    h_BDT_BG_polm80_train.GetYaxis().SetTitle("Entries");
    
    canvas_h_BDT_Signal_polm80_test = setUpperCanvas("canvas_h_BDT_Signal_polm80_test");
    canvas_h_BDT_Signal_polm80_test.cd();

    leg_h_BDT_Signal_polm80_test=TLegend(0.25,0.65,0.65,0.90);
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
    leg_h_BDT_Signal_polm80_test.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_Signal_polm80_test.AddEntry(h_BDT_Signal_polm80_test.DrawCopy("hist,e"),"Signal,test");
    leg_h_BDT_Signal_polm80_test.AddEntry(h_BDT_BG_polm80_test.DrawCopy("hist,e,same"),"BG,test");
    leg_h_BDT_Signal_polm80_test.AddEntry(h_BDT_Signal_polm80_train.DrawCopy("hist,e,same"),"Signal,Train");
    leg_h_BDT_Signal_polm80_test.AddEntry(h_BDT_BG_polm80_train.DrawCopy("hist,e,same"),"BG,Train");
    leg_h_BDT_Signal_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_BDT_Signal_polm80_test.Update()
    canvas_h_BDT_Signal_polm80_test.Modify()
    
    #canvas_h_BDT_Signal_polm80_test.Print("h_BDT_OverTrainingChecks_polm80.eps")
    #canvas_h_BDT_Signal_polm80_test.Print("h_BDT_OverTrainingChecks_polm80.pdf")
  
    h_BDT_Signal_polp80_test=file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polp80.Get("dataset/Method_BDT/BDT/MVA_BDT_S");
    h_BDT_Signal_polp80_test.SetFillColor(kBlue-2);
    h_BDT_Signal_polp80_test.SetFillStyle(3001);
    h_BDT_Signal_polp80_test.GetXaxis().SetTitle("BDT Score");
    h_BDT_Signal_polp80_test.GetYaxis().SetTitle("Entries");
    
    h_BDT_Signal_polp80_train=file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polp80.Get("dataset/Method_BDT/BDT/MVA_BDT_Train_S");
    h_BDT_Signal_polp80_train.SetLineColor(kBlue);
    h_BDT_Signal_polp80_train.SetLineWidth(2);
    h_BDT_Signal_polp80_train.GetXaxis().SetTitle("BDT Score");
    h_BDT_Signal_polp80_train.GetYaxis().SetTitle("Entries");
    
    h_BDT_BG_polp80_test=file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polp80.Get("dataset/Method_BDT/BDT/MVA_BDT_B");
    h_BDT_BG_polp80_test.SetFillColor(kRed-2);
    h_BDT_BG_polp80_test.SetFillStyle(3002);
    h_BDT_BG_polp80_test.GetXaxis().SetTitle("BDT Score");
    h_BDT_BG_polp80_test.GetYaxis().SetTitle("Entries");
    
    h_BDT_BG_polp80_train=file_BDT_HZqq_BG_qq_1TeV_qqqq_2TeV_qqqqqq_polp80.Get("dataset/Method_BDT/BDT/MVA_BDT_Train_B");
    h_BDT_BG_polp80_train.SetLineColor(kRed);
    h_BDT_BG_polp80_train.SetLineWidth(2);
    h_BDT_BG_polp80_train.GetXaxis().SetTitle("BDT Score");
    h_BDT_BG_polp80_train.GetYaxis().SetTitle("Entries");
    
    canvas_h_BDT_Signal_polp80_test = setUpperCanvas("canvas_h_BDT_Signal_polp80_test");
    canvas_h_BDT_Signal_polp80_test.cd();
    
    leg_h_BDT_Signal_polp80_test=TLegend(0.25,0.65,0.65,0.90);
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
    leg_h_BDT_Signal_polp80_test.SetHeader("#sqrt{s}>2500 GeV");
    leg_h_BDT_Signal_polp80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_BDT_Signal_polp80_test.Update()
    canvas_h_BDT_Signal_polp80_test.Modify()
    
    #canvas_h_BDT_Signal_polp80_test.Print("h_BDT_OverTrainingChecks_polp80.eps")
    #canvas_h_BDT_Signal_polp80_test.Print("h_BDT_OverTrainingChecks_polp80.pdf")
    
    file_polm80_HZ_noCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_noMassCut_Aug7/test_hzqq_AnalysisBaselineHistos_noCuts_EThetaVar_withSignalHistos.root")
    file_polm80_HZ_AllEvents_noCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_noMassCut_Aug7/test_hzqq_AnalysisBaselineHistos_noCuts_EThetaVar_withSignalHistos_AllEvents.root")
    file_polm80_ee_qq_noCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_noMassCut_Aug7/test_ee_qq_AnalysisBaselineHistos_noCuts_EThetaVar.root")
    file_polm80_ee_qqqq_noCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_noMassCut_Aug7/test_ee_qqqq_AnalysisBaselineHistos_noCuts_EThetaVar.root")
    file_polm80_ee_qqqqqq_noCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_noMassCut_Aug7/test_ee_qqqqqq_AnalysisBaselineHistos_noCuts_EThetaVar.root")

    h_mass_H_matched_orig_HZ_polm80_=file_polm80_HZ_noCuts_.Get("h_HZqq_signal_o_mass_H_matched")
    h_mass_H_matched_orig_HZ_polm80_.GetXaxis().SetTitle("m_{jet} [GeV]");
    h_mass_H_matched_orig_HZ_polm80_.SetLineWidth(3)
    h_mass_H_matched_orig_HZ_polm80_.SetLineColor(kRed)
    h_mass_H_matched_orig_HZ_polm80_.SetLineStyle(2)

    h_mass_H_matched_HZ_polm80_=file_polm80_HZ_noCuts_.Get("h_HZqq_signal_mass_H_matched")
    h_mass_H_matched_HZ_polm80_.GetXaxis().SetTitle("m_{jet} [GeV]");
    h_mass_H_matched_HZ_polm80_.SetLineWidth(3)
    h_mass_H_matched_HZ_polm80_.SetLineColor(kRed)

    h_mass_Z_matched_orig_HZ_polm80_=file_polm80_HZ_noCuts_.Get("h_HZqq_signal_o_mass_Z_matched")
    h_mass_Z_matched_orig_HZ_polm80_.GetXaxis().SetTitle("m_{jet} [GeV]");
    h_mass_Z_matched_orig_HZ_polm80_.SetLineWidth(3)
    h_mass_Z_matched_orig_HZ_polm80_.SetLineColor(kBlue)
    h_mass_Z_matched_orig_HZ_polm80_.SetLineStyle(2)

    h_mass_Z_matched_HZ_polm80_=file_polm80_HZ_noCuts_.Get("h_HZqq_signal_mass_Z_matched")
    h_mass_Z_matched_HZ_polm80_.GetXaxis().SetTitle("m_{jet} [GeV]");
    h_mass_Z_matched_HZ_polm80_.SetLineWidth(3)
    h_mass_Z_matched_HZ_polm80_.SetLineColor(kBlue)

    h_mass_H_matched_HZ_polm80_.SetMaximum(140)
    #h_mass_H_matched_HZ_polm80_.SetMaximum(180)

    canvas_HZqq_H_Z_mass_comp_polm80 = setUpperCanvas("canvas_HZqq_H_Z_mass_comp_polm80")
    canvas_HZqq_H_Z_mass_comp_polm80.cd()

    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_=TLegend(0.20,0.63,0.50,0.88);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetBorderSize(0);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetTextAlign(12);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetTextSize(0.050);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetTextFont(42);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetMargin(0.15);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetLineColor(1);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetLineStyle(1);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetLineWidth(1);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetFillColor(0);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetFillStyle(0);
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.AddEntry(h_mass_H_matched_HZ_polm80_.DrawCopy("hist"),"H jet, corr");
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.AddEntry(h_mass_H_matched_orig_HZ_polm80_.DrawCopy("hist,same"),"H jet, orig");
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.AddEntry(h_mass_Z_matched_HZ_polm80_.DrawCopy("hist,same"),"Z jet, corr");
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.AddEntry(h_mass_Z_matched_orig_HZ_polm80_.DrawCopy("hist,same"),"Z jet, orig");
    leg_hzqq_SIG_H_Z_mass_comp_polm80_noMassCuts_.Draw();


    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_HZqq_H_Z_mass_comp_polm80.Print("h_HZqq_H_Z_mass_comp_polm80.eps")



    h_2D_jet_mass_HZ_AllEvents=file_polm80_HZ_AllEvents_noCuts_.Get("h_2D_mass_jet1_vs_jet2");
    h_2D_jet_mass_HZ_AllEvents.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_2D_jet_mass_HZ_AllEvents.GetYaxis().SetTitle("jet2 mass [GeV]");


    canvas_h_2D_jet_mass_HZ_AllEvents_polm80_test = setUpperCanvas("canvas_h_2D_jet_mass_HZ_AllEvents_polm80_test");
    canvas_h_2D_jet_mass_HZ_AllEvents_polm80_test.cd();
    
    h_2D_jet_mass_HZ_AllEvents.Draw("col")

    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test=TLegend(0.25,0.80,0.65,0.85);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetBorderSize(0);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetTextAlign(12);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetTextSize(0.050);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetTextFont(42);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetMargin(0.15);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetLineColor(1);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetLineStyle(1);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetLineWidth(1);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetFillColor(0);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetFillStyle(0);
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.SetHeader("HZ, #sqrt{s}>2500 GeV");
    leg_h_2D_jet_mass_HZ_AllEvents_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_2D_jet_mass_HZ_AllEvents_polm80_test.Print("h_2D_HZ_AllEvents_jet1_vs_jet2_polm80.eps")


   
    h_2D_jet_mass_HZ=file_polm80_HZ_noCuts_.Get("h_2D_mass_jet1_vs_jet2");
    h_2D_jet_mass_HZ.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_2D_jet_mass_HZ.GetYaxis().SetTitle("jet2 mass [GeV]");


    canvas_h_2D_jet_mass_HZ_polm80_test = setUpperCanvas("canvas_h_2D_jet_mass_HZ_polm80_test");
    canvas_h_2D_jet_mass_HZ_polm80_test.cd();
    
    h_2D_jet_mass_HZ.Draw("col")

    leg_h_2D_jet_mass_HZ_polm80_test=TLegend(0.25,0.80,0.65,0.85);
    leg_h_2D_jet_mass_HZ_polm80_test.SetBorderSize(0);
    leg_h_2D_jet_mass_HZ_polm80_test.SetTextAlign(12);
    leg_h_2D_jet_mass_HZ_polm80_test.SetTextSize(0.050);
    leg_h_2D_jet_mass_HZ_polm80_test.SetTextFont(42);
    leg_h_2D_jet_mass_HZ_polm80_test.SetMargin(0.15);
    leg_h_2D_jet_mass_HZ_polm80_test.SetLineColor(1);
    leg_h_2D_jet_mass_HZ_polm80_test.SetLineStyle(1);
    leg_h_2D_jet_mass_HZ_polm80_test.SetLineWidth(1);
    leg_h_2D_jet_mass_HZ_polm80_test.SetFillColor(0);
    leg_h_2D_jet_mass_HZ_polm80_test.SetFillStyle(0);
    leg_h_2D_jet_mass_HZ_polm80_test.SetHeader("HZ with H#rightarrowb#bar{b}, #sqrt{s}>2500 GeV");
    leg_h_2D_jet_mass_HZ_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_2D_jet_mass_HZ_polm80_test.Print("h_2D_HZ_jet1_vs_jet2_polm80.eps")


    h_2D_jet_mass_ee_qq=file_polm80_ee_qq_noCuts_.Get("h_2D_mass_jet1_vs_jet2");
    h_2D_jet_mass_ee_qq.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_2D_jet_mass_ee_qq.GetYaxis().SetTitle("jet2 mass [GeV]");


    canvas_h_2D_jet_mass_ee_qq_polm80_test = setUpperCanvas("canvas_h_2D_jet_mass_ee_qq_polm80_test");
    canvas_h_2D_jet_mass_ee_qq_polm80_test.cd();
    
    h_2D_jet_mass_ee_qq.Draw("col")

    leg_h_2D_jet_mass_ee_qq_polm80_test=TLegend(0.25,0.80,0.65,0.85);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetBorderSize(0);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetTextAlign(12);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetTextSize(0.050);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetTextFont(42);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetMargin(0.15);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetLineColor(1);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetLineStyle(1);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetLineWidth(1);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetFillColor(0);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetFillStyle(0);
    leg_h_2D_jet_mass_ee_qq_polm80_test.SetHeader("ee#rightarrowqq, #sqrt{s}>2500 GeV");
    leg_h_2D_jet_mass_ee_qq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_2D_jet_mass_ee_qq_polm80_test.Print("h_2D_ee_qq_jet1_vs_jet2_polm80.eps")

    h_2D_jet_mass_ee_qqqq=file_polm80_ee_qqqq_noCuts_.Get("h_2D_mass_jet1_vs_jet2");
    h_2D_jet_mass_ee_qqqq.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_2D_jet_mass_ee_qqqq.GetYaxis().SetTitle("jet2 mass [GeV]");


    canvas_h_2D_jet_mass_ee_qqqq_polm80_test = setUpperCanvas("canvas_h_2D_jet_mass_ee_qqqq_polm80_test");
    canvas_h_2D_jet_mass_ee_qqqq_polm80_test.cd();
    
    h_2D_jet_mass_ee_qqqq.Draw("col")

    leg_h_2D_jet_mass_ee_qqqq_polm80_test=TLegend(0.25,0.80,0.65,0.85);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetBorderSize(0);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetTextAlign(12);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetTextSize(0.050);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetTextFont(42);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetMargin(0.15);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetLineColor(1);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetLineStyle(1);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetLineWidth(1);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetFillColor(0);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetFillStyle(0);
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.SetHeader("ee#rightarrowqqqq, #sqrt{s}>2500 GeV");
    leg_h_2D_jet_mass_ee_qqqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_2D_jet_mass_ee_qqqq_polm80_test.Print("h_2D_ee_qqqq_jet1_vs_jet2_polm80.eps")

    h_2D_jet_mass_ee_qqqqqq=file_polm80_ee_qqqqqq_noCuts_.Get("h_2D_mass_jet1_vs_jet2");
    h_2D_jet_mass_ee_qqqqqq.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_2D_jet_mass_ee_qqqqqq.GetYaxis().SetTitle("jet2 mass [GeV]");


    canvas_h_2D_jet_mass_ee_qqqqqq_polm80_test = setUpperCanvas("canvas_h_2D_jet_mass_ee_qqqqqq_polm80_test");
    canvas_h_2D_jet_mass_ee_qqqqqq_polm80_test.cd();
    
    h_2D_jet_mass_ee_qqqqqq.Draw("col")

    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test=TLegend(0.25,0.80,0.65,0.85);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetBorderSize(0);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetTextAlign(12);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetTextSize(0.050);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetTextFont(42);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetMargin(0.15);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetLineColor(1);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetLineStyle(1);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetLineWidth(1);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetFillColor(0);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetFillStyle(0);
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.SetHeader("ee#rightarrowqqqqqq, #sqrt{s}>2500 GeV");
    leg_h_2D_jet_mass_ee_qqqqqq_polm80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_2D_jet_mass_ee_qqqqqq_polm80_test.Print("h_2D_ee_qqqqqq_jet1_vs_jet2_polm80.eps")



    file_polp80_HZ_noCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_noMassCut_Aug7/test_hzqq_AnalysisBaselineHistos_noCuts_EThetaVar_withSignalHistos.root")
    file_polp80_HZ_AllEvents_noCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_noMassCut_Aug7/test_hzqq_AnalysisBaselineHistos_noCuts_EThetaVar_withSignalHistos_AllEvents.root")
    file_polp80_ee_qq_noCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_noMassCut_Aug7/test_ee_qq_AnalysisBaselineHistos_noCuts_EThetaVar.root")
    file_polp80_ee_qqqq_noCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_noMassCut_Aug7/test_ee_qqqq_AnalysisBaselineHistos_noCuts_EThetaVar.root")
    file_polp80_ee_qqqqqq_noCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_noMassCut_Aug7/test_ee_qqqqqq_AnalysisBaselineHistos_noCuts_EThetaVar.root")



    h_2D_jet_mass_HZ_AllEvents=file_polp80_HZ_AllEvents_noCuts_.Get("h_2D_mass_jet1_vs_jet2");
    h_2D_jet_mass_HZ_AllEvents.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_2D_jet_mass_HZ_AllEvents.GetYaxis().SetTitle("jet2 mass [GeV]");


    canvas_h_2D_jet_mass_HZ_AllEvents_polp80_test = setUpperCanvas("canvas_h_2D_jet_mass_HZ_AllEvents_polp80_test");
    canvas_h_2D_jet_mass_HZ_AllEvents_polp80_test.cd();
    
    h_2D_jet_mass_HZ_AllEvents.Draw("col")

    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test=TLegend(0.25,0.80,0.65,0.85);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetBorderSize(0);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetTextAlign(12);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetTextSize(0.050);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetTextFont(42);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetMargin(0.15);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetLineColor(1);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetLineStyle(1);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetLineWidth(1);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetFillColor(0);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetFillStyle(0);
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.SetHeader("HZ, #sqrt{s}>2500 GeV");
    leg_h_2D_jet_mass_HZ_AllEvents_polp80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_2D_jet_mass_HZ_AllEvents_polp80_test.Print("h_2D_HZ_AllEvents_jet1_vs_jet2_polp80.eps")


   
    h_2D_jet_mass_HZ_polp80=file_polp80_HZ_noCuts_.Get("h_2D_mass_jet1_vs_jet2");
    h_2D_jet_mass_HZ_polp80.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_2D_jet_mass_HZ_polp80.GetYaxis().SetTitle("jet2 mass [GeV]");


    canvas_h_2D_jet_mass_HZ_polp80_test = setUpperCanvas("canvas_h_2D_jet_mass_HZ_polp80_test");
    canvas_h_2D_jet_mass_HZ_polp80_test.cd();
    
    h_2D_jet_mass_HZ_polp80.Draw("col")

    leg_h_2D_jet_mass_HZ_polp80_test=TLegend(0.25,0.80,0.65,0.85);
    leg_h_2D_jet_mass_HZ_polp80_test.SetBorderSize(0);
    leg_h_2D_jet_mass_HZ_polp80_test.SetTextAlign(12);
    leg_h_2D_jet_mass_HZ_polp80_test.SetTextSize(0.050);
    leg_h_2D_jet_mass_HZ_polp80_test.SetTextFont(42);
    leg_h_2D_jet_mass_HZ_polp80_test.SetMargin(0.15);
    leg_h_2D_jet_mass_HZ_polp80_test.SetLineColor(1);
    leg_h_2D_jet_mass_HZ_polp80_test.SetLineStyle(1);
    leg_h_2D_jet_mass_HZ_polp80_test.SetLineWidth(1);
    leg_h_2D_jet_mass_HZ_polp80_test.SetFillColor(0);
    leg_h_2D_jet_mass_HZ_polp80_test.SetFillStyle(0);
    leg_h_2D_jet_mass_HZ_polp80_test.SetHeader("HZ with H#rightarrowb#bar{b}, #sqrt{s}>2500 GeV");
    leg_h_2D_jet_mass_HZ_polp80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_2D_jet_mass_HZ_polp80_test.Print("h_2D_HZ_jet1_vs_jet2_polp80.eps")


    h_2D_jet_mass_ee_qq_polp80=file_polp80_ee_qq_noCuts_.Get("h_2D_mass_jet1_vs_jet2");
    h_2D_jet_mass_ee_qq_polp80.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_2D_jet_mass_ee_qq_polp80.GetYaxis().SetTitle("jet2 mass [GeV]");


    canvas_h_2D_jet_mass_ee_qq_polp80_test = setUpperCanvas("canvas_h_2D_jet_mass_ee_qq_polp80_test");
    canvas_h_2D_jet_mass_ee_qq_polp80_test.cd();
    
    h_2D_jet_mass_ee_qq_polp80.Draw("col")

    leg_h_2D_jet_mass_ee_qq_polp80_test=TLegend(0.25,0.80,0.65,0.85);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetBorderSize(0);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetTextAlign(12);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetTextSize(0.050);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetTextFont(42);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetMargin(0.15);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetLineColor(1);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetLineStyle(1);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetLineWidth(1);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetFillColor(0);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetFillStyle(0);
    leg_h_2D_jet_mass_ee_qq_polp80_test.SetHeader("ee#rightarrowqq, #sqrt{s}>2500 GeV");
    leg_h_2D_jet_mass_ee_qq_polp80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_2D_jet_mass_ee_qq_polp80_test.Print("h_2D_ee_qq_jet1_vs_jet2_polp80.eps")

    h_2D_jet_mass_ee_qqqq_polp80=file_polp80_ee_qqqq_noCuts_.Get("h_2D_mass_jet1_vs_jet2");
    h_2D_jet_mass_ee_qqqq_polp80.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_2D_jet_mass_ee_qqqq_polp80.GetYaxis().SetTitle("jet2 mass [GeV]");


    canvas_h_2D_jet_mass_ee_qqqq_polp80_test = setUpperCanvas("canvas_h_2D_jet_mass_ee_qqqq_polp80_test");
    canvas_h_2D_jet_mass_ee_qqqq_polp80_test.cd();
    
    h_2D_jet_mass_ee_qqqq_polp80.Draw("col")

    leg_h_2D_jet_mass_ee_qqqq_polp80_test=TLegend(0.25,0.80,0.65,0.85);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetBorderSize(0);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetTextAlign(12);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetTextSize(0.050);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetTextFont(42);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetMargin(0.15);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetLineColor(1);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetLineStyle(1);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetLineWidth(1);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetFillColor(0);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetFillStyle(0);
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.SetHeader("ee#rightarrowqqqq, #sqrt{s}>2500 GeV");
    leg_h_2D_jet_mass_ee_qqqq_polp80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_2D_jet_mass_ee_qqqq_polp80_test.Print("h_2D_ee_qqqq_jet1_vs_jet2_polp80.eps")

    h_2D_jet_mass_ee_qqqqqq_polp80=file_polp80_ee_qqqqqq_noCuts_.Get("h_2D_mass_jet1_vs_jet2");
    h_2D_jet_mass_ee_qqqqqq_polp80.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_2D_jet_mass_ee_qqqqqq_polp80.GetYaxis().SetTitle("jet2 mass [GeV]");


    canvas_h_2D_jet_mass_ee_qqqqqq_polp80_test = setUpperCanvas("canvas_h_2D_jet_mass_ee_qqqqqq_polp80_test");
    canvas_h_2D_jet_mass_ee_qqqqqq_polp80_test.cd();
    
    h_2D_jet_mass_ee_qqqqqq_polp80.Draw("col")

    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test=TLegend(0.25,0.80,0.65,0.85);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetBorderSize(0);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetTextAlign(12);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetTextSize(0.050);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetTextFont(42);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetMargin(0.15);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetLineColor(1);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetLineStyle(1);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetLineWidth(1);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetFillColor(0);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetFillStyle(0);
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.SetHeader("ee#rightarrowqqqqqq, #sqrt{s}>2500 GeV");
    leg_h_2D_jet_mass_ee_qqqqqq_polp80_test.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_2D_jet_mass_ee_qqqqqq_polp80_test.Print("h_2D_ee_qqqqqq_jet1_vs_jet2_polp80.eps")






    file_polm80_HZ_massCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_hzqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar_withSignalHistos.root")
    file_polm80_HZ_AllEvents_massCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_hzqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar_withSignalHistos_AllEvents.root")
    file_polm80_ee_qq_massCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_ee_qq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root")
    file_polm80_ee_qqqq_massCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_ee_qqqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root")
    file_polm80_ee_qqqqqq_massCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_ee_qqqqqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root")

    h_jet1_mass_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_mass_jet1");
    h_jet1_mass_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_jet1_mass_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_mass_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_mass_jet1");
    h_jet1_mass_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_jet1_mass_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet1_mass_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet1_mass_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet1_mass_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet1_mass_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_mass_jet1");
    h_jet1_mass_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet1_mass_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet1_mass_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_mass_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_mass_jet1");
    h_jet1_mass_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet1_mass_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet1_mass_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_mass_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_mass_jet1");
    h_jet1_mass_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_mass_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_mass_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_tot_norm_jet1_mass_BG = h_jet1_mass_ee_qq_polm80_massCuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_jet1_mass_BG.Add(h_jet1_mass_ee_qqqq_polm80_massCuts_);
    h_tot_norm_jet1_mass_BG.Add(h_jet1_mass_ee_qqqqqq_polm80_massCuts_);
    norm_tot_BG_to_SIG=h_jet1_mass_HZ_polm80_massCuts_.Integral()/(h_jet1_mass_ee_qq_polm80_massCuts_.Integral()+h_jet1_mass_ee_qqqq_polm80_massCuts_.Integral()+h_jet1_mass_ee_qqqqqq_polm80_massCuts_.Integral())
    h_tot_norm_jet1_mass_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_jet1_mass_BG.SetLineColor(kBlack)
    h_tot_norm_jet1_mass_BG.SetFillColor(0)

    print 'scale or range ',h_tot_norm_jet1_mass_BG.Integral(),h_jet1_mass_HZ_polm80_massCuts_.Integral()

    h_jet1_mass_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_mass_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_mass_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)

    h_jet1_mass_ee_qq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    h_jet1_mass_ee_qqqq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    h_jet1_mass_ee_qqqqqq_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    h_jet1_mass_HZ_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)

    hzqq_BG_jet1_mass_polm80_massCuts_= THStack("hzqq_BG_jet1_mass_polm80_massCuts_", "");
    hzqq_BG_jet1_mass_polm80_massCuts_.Add(h_jet1_mass_ee_qq_polm80_massCuts_);
    hzqq_BG_jet1_mass_polm80_massCuts_.Add(h_jet1_mass_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet1_mass_polm80_massCuts_.Add(h_jet1_mass_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_mass_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_mass_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_mass_polm80_massCuts_thstack.cd();
    #h_jet1_mass_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_mass_BG.Draw("hist,e")
    #h_tot_norm_jet1_mass_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_mass_polm80_massCuts_.Draw("hist");
    hzqq_BG_jet1_mass_polm80_massCuts_.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_mass_polm80_massCuts_.GetXaxis().SetTitle("jet1 mass [GeV]");
    hzqq_BG_jet1_mass_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_mass_polm80_massCuts_.SetMaximum(155)
    h_jet1_mass_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_mass_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_mass_polm80_massCuts_=TLegend(0.35,0.63,0.75,0.88);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.AddEntry(h_jet1_mass_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.AddEntry(h_jet1_mass_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.AddEntry(h_jet1_mass_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.AddEntry(h_jet1_mass_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_mass_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet1_mass_polm80_massCuts_thstack.Print("h_jet1_mass_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet2_mass_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_mass_jet2");
    h_jet2_mass_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 mass [GeV]");
    h_jet2_mass_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_mass_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_mass_jet2");
    h_jet2_mass_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet2 mass [GeV]");
    h_jet2_mass_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet2_mass_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet2_mass_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet2_mass_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet2_mass_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_mass_jet2");
    h_jet2_mass_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet2_mass_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet2_mass_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_mass_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_mass_jet2");
    h_jet2_mass_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet2_mass_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet2_mass_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_mass_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_mass_jet2");
    h_jet2_mass_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_mass_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_mass_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet2_mass_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_mass_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_mass_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)

    hzqq_BG_jet2_mass_polm80_massCuts_= THStack("hzqq_BG_jet2_mass_polm80_massCuts_", "");
    hzqq_BG_jet2_mass_polm80_massCuts_.Add(h_jet2_mass_ee_qq_polm80_massCuts_);
    hzqq_BG_jet2_mass_polm80_massCuts_.Add(h_jet2_mass_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet2_mass_polm80_massCuts_.Add(h_jet2_mass_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_mass_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_mass_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_mass_polm80_massCuts_thstack.cd();
    #h_jet2_mass_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_mass_BG.Draw("hist,e")
    #h_tot_norm_jet2_mass_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_mass_polm80_massCuts_.Draw("hist");
    hzqq_BG_jet2_mass_polm80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet2_mass_polm80_massCuts_.GetXaxis().SetTitle("jet2 mass [GeV]");
    hzqq_BG_jet2_mass_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_mass_polm80_massCuts_.SetMaximum(160)
    h_jet2_mass_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_mass_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_mass_polm80_massCuts_=TLegend(0.20,0.63,0.55,0.88);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.AddEntry(h_jet2_mass_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.AddEntry(h_jet2_mass_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.AddEntry(h_jet2_mass_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.AddEntry(h_jet2_mass_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_mass_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet2_mass_polm80_massCuts_thstack.Print("h_jet2_mass_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")



    h_jet2_theta_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_theta_jet2");
    h_jet2_theta_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    h_jet2_theta_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_theta_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_theta_jet2");
    h_jet2_theta_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    h_jet2_theta_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet2_theta_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet2_theta_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet2_theta_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet2_theta_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_theta_jet2");
    h_jet2_theta_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet2_theta_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet2_theta_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_theta_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_theta_jet2");
    h_jet2_theta_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet2_theta_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet2_theta_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_theta_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_theta_jet2");
    h_jet2_theta_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_theta_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_theta_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet2_theta_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_theta_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_theta_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet2_theta_HZ_polm80_massCuts_.Rebin(2)
    h_jet2_theta_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet2_theta_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet2_theta_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet2_theta_polm80_massCuts_= THStack("hzqq_BG_jet2_theta_polm80_massCuts_", "");
    hzqq_BG_jet2_theta_polm80_massCuts_.Add(h_jet2_theta_ee_qq_polm80_massCuts_);
    hzqq_BG_jet2_theta_polm80_massCuts_.Add(h_jet2_theta_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet2_theta_polm80_massCuts_.Add(h_jet2_theta_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_theta_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_theta_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_theta_polm80_massCuts_thstack.cd();
    #h_jet2_theta_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_theta_BG.Draw("hist,e")
    #h_tot_norm_jet2_theta_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_theta_polm80_massCuts_.Draw("hist");
    #hzqq_BG_jet2_theta_polm80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet2_theta_polm80_massCuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    hzqq_BG_jet2_theta_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_theta_polm80_massCuts_.SetMaximum(95)
    h_jet2_theta_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_theta_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_theta_polm80_massCuts_=TLegend(0.20,0.63,0.55,0.88);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.AddEntry(h_jet2_theta_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.AddEntry(h_jet2_theta_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.AddEntry(h_jet2_theta_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.AddEntry(h_jet2_theta_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_theta_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet2_theta_polm80_massCuts_thstack.Print("h_jet2_theta_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet1_theta_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_theta_jet1");
    h_jet1_theta_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    h_jet1_theta_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_theta_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_theta_jet1");
    h_jet1_theta_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    h_jet1_theta_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet1_theta_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet1_theta_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet1_theta_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet1_theta_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_theta_jet1");
    h_jet1_theta_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet1_theta_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet1_theta_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_theta_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_theta_jet1");
    h_jet1_theta_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet1_theta_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet1_theta_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_theta_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_theta_jet1");
    h_jet1_theta_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_theta_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_theta_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet1_theta_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_theta_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_theta_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_theta_HZ_polm80_massCuts_.Rebin(2)
    h_jet1_theta_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet1_theta_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet1_theta_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet1_theta_polm80_massCuts_= THStack("hzqq_BG_jet1_theta_polm80_massCuts_", "");
    hzqq_BG_jet1_theta_polm80_massCuts_.Add(h_jet1_theta_ee_qq_polm80_massCuts_);
    hzqq_BG_jet1_theta_polm80_massCuts_.Add(h_jet1_theta_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet1_theta_polm80_massCuts_.Add(h_jet1_theta_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_theta_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_theta_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_theta_polm80_massCuts_thstack.cd();
    #h_jet1_theta_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_theta_BG.Draw("hist,e")
    #h_tot_norm_jet1_theta_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_theta_polm80_massCuts_.Draw("hist");
    #hzqq_BG_jet1_theta_polm80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet1_theta_polm80_massCuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    hzqq_BG_jet1_theta_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_theta_polm80_massCuts_.SetMaximum(95)
    h_jet1_theta_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_theta_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_theta_polm80_massCuts_=TLegend(0.20,0.63,0.55,0.88);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.AddEntry(h_jet1_theta_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.AddEntry(h_jet1_theta_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.AddEntry(h_jet1_theta_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.AddEntry(h_jet1_theta_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_theta_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet1_theta_polm80_massCuts_thstack.Print("h_jet1_theta_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    #defined as BTag of the larger BTagged subjet of jet1
    h_jet1_BTag_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_BTag_sj_BTagMax");
    h_jet1_BTag_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 BTag");
    h_jet1_BTag_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_BTag_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_BTag_sj_BTagMax");
    h_jet1_BTag_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet1 BTag");
    h_jet1_BTag_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet1_BTag_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet1_BTag_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet1_BTag_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet1_BTag_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_BTag_sj_BTagMax");
    h_jet1_BTag_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet1_BTag_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet1_BTag_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_BTag_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_BTag_sj_BTagMax");
    h_jet1_BTag_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet1_BTag_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet1_BTag_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_BTag_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_BTag_sj_BTagMax");
    h_jet1_BTag_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_BTag_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_BTag_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet1_BTag_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_BTag_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_BTag_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_BTag_HZ_polm80_massCuts_.Rebin(2)
    h_jet1_BTag_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet1_BTag_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet1_BTag_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet1_BTag_polm80_massCuts_= THStack("hzqq_BG_jet1_BTag_polm80_massCuts_", "");
    hzqq_BG_jet1_BTag_polm80_massCuts_.Add(h_jet1_BTag_ee_qq_polm80_massCuts_);
    hzqq_BG_jet1_BTag_polm80_massCuts_.Add(h_jet1_BTag_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet1_BTag_polm80_massCuts_.Add(h_jet1_BTag_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_BTag_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_BTag_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_BTag_polm80_massCuts_thstack.cd();
    #h_jet1_BTag_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_BTag_BG.Draw("hist,e")
    #h_tot_norm_jet1_BTag_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_BTag_polm80_massCuts_.Draw("hist");
    #hzqq_BG_jet1_BTag_polm80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet1_BTag_polm80_massCuts_.GetXaxis().SetTitle("jet1 BTag");
    hzqq_BG_jet1_BTag_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_BTag_polm80_massCuts_.SetMaximum(325)
    h_jet1_BTag_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_BTag_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_=TLegend(0.30,0.63,0.65,0.88);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.AddEntry(h_jet1_BTag_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.AddEntry(h_jet1_BTag_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.AddEntry(h_jet1_BTag_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.AddEntry(h_jet1_BTag_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_BTag_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet1_BTag_polm80_massCuts_thstack.Print("h_jet1_BTag_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")



    h_reco_y32_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_reco_y32");
    h_reco_y32_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("y_{23}");
    h_reco_y32_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_reco_y32_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_reco_y32");
    h_reco_y32_HZ_polm80_massCuts_.GetXaxis().SetTitle("y_{23}");
    h_reco_y32_HZ_polm80_massCuts_.SetLineWidth(3);
    h_reco_y32_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_reco_y32_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_reco_y32_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_reco_y32_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_reco_y32");
    h_reco_y32_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_reco_y32_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_reco_y32_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_reco_y32_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_reco_y32");
    h_reco_y32_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_reco_y32_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_reco_y32_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_reco_y32_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_reco_y32");
    h_reco_y32_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_reco_y32_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_reco_y32_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_reco_y32_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_reco_y32_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_reco_y32_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_reco_y32_HZ_polm80_massCuts_.Rebin(2)
    h_reco_y32_ee_qq_polm80_massCuts_.Rebin(2)
    h_reco_y32_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_reco_y32_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_reco_y32_polm80_massCuts_= THStack("hzqq_BG_reco_y32_polm80_massCuts_", "");
    hzqq_BG_reco_y32_polm80_massCuts_.Add(h_reco_y32_ee_qq_polm80_massCuts_);
    hzqq_BG_reco_y32_polm80_massCuts_.Add(h_reco_y32_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_reco_y32_polm80_massCuts_.Add(h_reco_y32_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_reco_y32_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_reco_y32_polm80_massCuts_thstack");
    canvas_h_SIG_BG_reco_y32_polm80_massCuts_thstack.cd();
    #h_reco_y32_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_reco_y32_BG.Draw("hist,e")
    #h_tot_norm_reco_y32_BG.GetXaxis().SetRangeUser(0,0.0045)
    hzqq_BG_reco_y32_polm80_massCuts_.Draw("hist");
    hzqq_BG_reco_y32_polm80_massCuts_.GetXaxis().SetTitle("y_{23}")
    hzqq_BG_reco_y32_polm80_massCuts_.GetXaxis().SetRangeUser(0,0.0045)
    hzqq_BG_reco_y32_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_reco_y32_polm80_massCuts_.SetMaximum(200)
    h_reco_y32_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_reco_y32_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_reco_y32_polm80_massCuts_=TLegend(0.50,0.63,0.85,0.88);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_reco_y32_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_reco_y32_polm80_massCuts_.AddEntry(h_reco_y32_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_reco_y32_polm80_massCuts_.AddEntry(h_reco_y32_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_reco_y32_polm80_massCuts_.AddEntry(h_reco_y32_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_reco_y32_polm80_massCuts_.AddEntry(h_reco_y32_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_reco_y32_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_reco_y32_polm80_massCuts_thstack.Print("h_reco_y32_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet2_beta1_D2_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_jet2_beta1_D2");
    h_jet2_beta1_D2_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 D_{2}^{(1)}");
    h_jet2_beta1_D2_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_beta1_D2_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_jet2_beta1_D2");
    h_jet2_beta1_D2_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet2 D_{2}^{(1)}");
    h_jet2_beta1_D2_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet2_beta1_D2_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet2_beta1_D2_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet2_beta1_D2_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet2_beta1_D2_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_jet2_beta1_D2");
    h_jet2_beta1_D2_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet2_beta1_D2_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet2_beta1_D2_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_D2_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_jet2_beta1_D2");
    h_jet2_beta1_D2_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet2_beta1_D2_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet2_beta1_D2_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_D2_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_jet2_beta1_D2");
    h_jet2_beta1_D2_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_beta1_D2_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_beta1_D2_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet2_beta1_D2_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_D2_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_D2_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet2_beta1_D2_HZ_polm80_massCuts_.Rebin(2)
    h_jet2_beta1_D2_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet2_beta1_D2_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet2_beta1_D2_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet2_beta1_D2_polm80_massCuts_= THStack("hzqq_BG_jet2_beta1_D2_polm80_massCuts_", "");
    hzqq_BG_jet2_beta1_D2_polm80_massCuts_.Add(h_jet2_beta1_D2_ee_qq_polm80_massCuts_);
    hzqq_BG_jet2_beta1_D2_polm80_massCuts_.Add(h_jet2_beta1_D2_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet2_beta1_D2_polm80_massCuts_.Add(h_jet2_beta1_D2_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_beta1_D2_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_beta1_D2_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_beta1_D2_polm80_massCuts_thstack.cd();
    #h_jet2_beta1_D2_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_beta1_D2_BG.Draw("hist,e")
    #h_tot_norm_jet2_beta1_D2_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_beta1_D2_polm80_massCuts_.Draw("hist");
    #hzqq_BG_jet2_beta1_D2_polm80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet2_beta1_D2_polm80_massCuts_.GetXaxis().SetTitle("jet2 D_{2}^{(1)}");
    hzqq_BG_jet2_beta1_D2_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetMaximum(200)
    h_jet2_beta1_D2_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_beta1_D2_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_=TLegend(0.40,0.63,0.75,0.88);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.AddEntry(h_jet2_beta1_D2_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.AddEntry(h_jet2_beta1_D2_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.AddEntry(h_jet2_beta1_D2_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.AddEntry(h_jet2_beta1_D2_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_beta1_D2_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet2_beta1_D2_polm80_massCuts_thstack.Print("h_jet2_beta1_D2_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet1_beta1_D2_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_jet1_beta1_D2");
    h_jet1_beta1_D2_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 D_{2}^{(1)}");
    h_jet1_beta1_D2_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_beta1_D2_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_jet1_beta1_D2");
    h_jet1_beta1_D2_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet1 D_{2}^{(1)}");
    h_jet1_beta1_D2_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet1_beta1_D2_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet1_beta1_D2_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet1_beta1_D2_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet1_beta1_D2_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_jet1_beta1_D2");
    h_jet1_beta1_D2_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet1_beta1_D2_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet1_beta1_D2_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_D2_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_jet1_beta1_D2");
    h_jet1_beta1_D2_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet1_beta1_D2_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet1_beta1_D2_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_D2_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_jet1_beta1_D2");
    h_jet1_beta1_D2_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_beta1_D2_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_beta1_D2_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet1_beta1_D2_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_D2_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_D2_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_beta1_D2_HZ_polm80_massCuts_.Rebin(2)
    h_jet1_beta1_D2_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet1_beta1_D2_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet1_beta1_D2_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet1_beta1_D2_polm80_massCuts_= THStack("hzqq_BG_jet1_beta1_D2_polm80_massCuts_", "");
    hzqq_BG_jet1_beta1_D2_polm80_massCuts_.Add(h_jet1_beta1_D2_ee_qq_polm80_massCuts_);
    hzqq_BG_jet1_beta1_D2_polm80_massCuts_.Add(h_jet1_beta1_D2_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet1_beta1_D2_polm80_massCuts_.Add(h_jet1_beta1_D2_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_beta1_D2_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_beta1_D2_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_beta1_D2_polm80_massCuts_thstack.cd();
    #h_jet1_beta1_D2_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_beta1_D2_BG.Draw("hist,e")
    #h_tot_norm_jet1_beta1_D2_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_beta1_D2_polm80_massCuts_.Draw("hist");
    #hzqq_BG_jet1_beta1_D2_polm80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet1_beta1_D2_polm80_massCuts_.GetXaxis().SetTitle("jet1 D_{2}^{(1)}");
    hzqq_BG_jet1_beta1_D2_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetMaximum(200)
    h_jet1_beta1_D2_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_beta1_D2_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_=TLegend(0.40,0.63,0.75,0.88);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.AddEntry(h_jet1_beta1_D2_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.AddEntry(h_jet1_beta1_D2_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.AddEntry(h_jet1_beta1_D2_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.AddEntry(h_jet1_beta1_D2_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_beta1_D2_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet1_beta1_D2_polm80_massCuts_thstack.Print("h_jet1_beta1_D2_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")







    h_jet2_beta1_C2_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_jet2_beta1_C2");
    h_jet2_beta1_C2_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 C_{2}^{(1)}");
    h_jet2_beta1_C2_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_beta1_C2_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_jet2_beta1_C2");
    h_jet2_beta1_C2_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet2 C_{2}^{(1)}");
    h_jet2_beta1_C2_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet2_beta1_C2_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet2_beta1_C2_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet2_beta1_C2_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet2_beta1_C2_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_jet2_beta1_C2");
    h_jet2_beta1_C2_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet2_beta1_C2_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet2_beta1_C2_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_C2_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_jet2_beta1_C2");
    h_jet2_beta1_C2_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet2_beta1_C2_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet2_beta1_C2_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_C2_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_jet2_beta1_C2");
    h_jet2_beta1_C2_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_beta1_C2_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_beta1_C2_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet2_beta1_C2_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_C2_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_C2_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet2_beta1_C2_HZ_polm80_massCuts_.Rebin(2)
    h_jet2_beta1_C2_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet2_beta1_C2_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet2_beta1_C2_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet2_beta1_C2_polm80_massCuts_= THStack("hzqq_BG_jet2_beta1_C2_polm80_massCuts_", "");
    hzqq_BG_jet2_beta1_C2_polm80_massCuts_.Add(h_jet2_beta1_C2_ee_qq_polm80_massCuts_);
    hzqq_BG_jet2_beta1_C2_polm80_massCuts_.Add(h_jet2_beta1_C2_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet2_beta1_C2_polm80_massCuts_.Add(h_jet2_beta1_C2_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_beta1_C2_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_beta1_C2_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_beta1_C2_polm80_massCuts_thstack.cd();
    #h_jet2_beta1_C2_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_beta1_C2_BG.Draw("hist,e")
    #h_tot_norm_jet2_beta1_C2_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_beta1_C2_polm80_massCuts_.Draw("hist");
    hzqq_BG_jet2_beta1_C2_polm80_massCuts_.GetXaxis().SetRangeUser(0,0.40)
    hzqq_BG_jet2_beta1_C2_polm80_massCuts_.GetXaxis().SetTitle("jet2 C_{2}^{(1)}");
    hzqq_BG_jet2_beta1_C2_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetMaximum(300)
    h_jet2_beta1_C2_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_beta1_C2_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_=TLegend(0.40,0.63,0.75,0.88);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.AddEntry(h_jet2_beta1_C2_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.AddEntry(h_jet2_beta1_C2_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.AddEntry(h_jet2_beta1_C2_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.AddEntry(h_jet2_beta1_C2_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_beta1_C2_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet2_beta1_C2_polm80_massCuts_thstack.Print("h_jet2_beta1_C2_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet1_beta1_C2_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_jet1_beta1_C2");
    h_jet1_beta1_C2_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 C_{2}^{(1)}");
    h_jet1_beta1_C2_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_beta1_C2_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_jet1_beta1_C2");
    h_jet1_beta1_C2_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet1 C_{2}^{(1)}");
    h_jet1_beta1_C2_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet1_beta1_C2_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet1_beta1_C2_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet1_beta1_C2_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet1_beta1_C2_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_jet1_beta1_C2");
    h_jet1_beta1_C2_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet1_beta1_C2_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet1_beta1_C2_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_C2_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_jet1_beta1_C2");
    h_jet1_beta1_C2_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet1_beta1_C2_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet1_beta1_C2_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_C2_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_jet1_beta1_C2");
    h_jet1_beta1_C2_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_beta1_C2_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_beta1_C2_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet1_beta1_C2_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_C2_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_C2_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_beta1_C2_HZ_polm80_massCuts_.Rebin(2)
    h_jet1_beta1_C2_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet1_beta1_C2_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet1_beta1_C2_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet1_beta1_C2_polm80_massCuts_= THStack("hzqq_BG_jet1_beta1_C2_polm80_massCuts_", "");
    hzqq_BG_jet1_beta1_C2_polm80_massCuts_.Add(h_jet1_beta1_C2_ee_qq_polm80_massCuts_);
    hzqq_BG_jet1_beta1_C2_polm80_massCuts_.Add(h_jet1_beta1_C2_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet1_beta1_C2_polm80_massCuts_.Add(h_jet1_beta1_C2_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_beta1_C2_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_beta1_C2_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_beta1_C2_polm80_massCuts_thstack.cd();
    #h_jet1_beta1_C2_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_beta1_C2_BG.Draw("hist,e")
    #h_tot_norm_jet1_beta1_C2_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_beta1_C2_polm80_massCuts_.Draw("hist");
    hzqq_BG_jet1_beta1_C2_polm80_massCuts_.GetXaxis().SetRangeUser(0,0.40)
    hzqq_BG_jet1_beta1_C2_polm80_massCuts_.GetXaxis().SetTitle("jet1 C_{2}^{(1)}");
    hzqq_BG_jet1_beta1_C2_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetMaximum(225)
    h_jet1_beta1_C2_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_beta1_C2_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_=TLegend(0.40,0.63,0.75,0.88);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.AddEntry(h_jet1_beta1_C2_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.AddEntry(h_jet1_beta1_C2_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.AddEntry(h_jet1_beta1_C2_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.AddEntry(h_jet1_beta1_C2_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_beta1_C2_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet1_beta1_C2_polm80_massCuts_thstack.Print("h_jet1_beta1_C2_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")







    h_jet2_beta1_C3_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_jet2_beta1_C3");
    h_jet2_beta1_C3_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 C_{3}^{(1)}");
    h_jet2_beta1_C3_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_beta1_C3_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_jet2_beta1_C3");
    h_jet2_beta1_C3_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet2 C_{3}^{(1)}");
    h_jet2_beta1_C3_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet2_beta1_C3_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet2_beta1_C3_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet2_beta1_C3_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet2_beta1_C3_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_jet2_beta1_C3");
    h_jet2_beta1_C3_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet2_beta1_C3_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet2_beta1_C3_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_C3_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_jet2_beta1_C3");
    h_jet2_beta1_C3_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet2_beta1_C3_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet2_beta1_C3_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_C3_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_jet2_beta1_C3");
    h_jet2_beta1_C3_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_beta1_C3_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_beta1_C3_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet2_beta1_C3_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_C3_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_C3_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet2_beta1_C3_HZ_polm80_massCuts_.Rebin(2)
    h_jet2_beta1_C3_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet2_beta1_C3_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet2_beta1_C3_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet2_beta1_C3_polm80_massCuts_= THStack("hzqq_BG_jet2_beta1_C3_polm80_massCuts_", "");
    hzqq_BG_jet2_beta1_C3_polm80_massCuts_.Add(h_jet2_beta1_C3_ee_qq_polm80_massCuts_);
    hzqq_BG_jet2_beta1_C3_polm80_massCuts_.Add(h_jet2_beta1_C3_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet2_beta1_C3_polm80_massCuts_.Add(h_jet2_beta1_C3_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_beta1_C3_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_beta1_C3_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_beta1_C3_polm80_massCuts_thstack.cd();
    #h_jet2_beta1_C3_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_beta1_C3_BG.Draw("hist,e")
    #h_tot_norm_jet2_beta1_C3_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_beta1_C3_polm80_massCuts_.Draw("hist");
    hzqq_BG_jet2_beta1_C3_polm80_massCuts_.GetXaxis().SetRangeUser(0.,0.625)
    hzqq_BG_jet2_beta1_C3_polm80_massCuts_.GetXaxis().SetTitle("jet2 C_{3}^{(1)}");
    hzqq_BG_jet2_beta1_C3_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetMaximum(150)
    h_jet2_beta1_C3_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_beta1_C3_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_=TLegend(0.40,0.63,0.75,0.88);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.AddEntry(h_jet2_beta1_C3_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.AddEntry(h_jet2_beta1_C3_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.AddEntry(h_jet2_beta1_C3_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.AddEntry(h_jet2_beta1_C3_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_beta1_C3_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet2_beta1_C3_polm80_massCuts_thstack.Print("h_jet2_beta1_C3_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet1_beta1_C3_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_jet1_beta1_C3");
    h_jet1_beta1_C3_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 C_{3}^{(1)}");
    h_jet1_beta1_C3_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_beta1_C3_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_jet1_beta1_C3");
    h_jet1_beta1_C3_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet1 C_{3}^{(1)}");
    h_jet1_beta1_C3_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet1_beta1_C3_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet1_beta1_C3_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet1_beta1_C3_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet1_beta1_C3_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_jet1_beta1_C3");
    h_jet1_beta1_C3_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet1_beta1_C3_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet1_beta1_C3_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_C3_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_jet1_beta1_C3");
    h_jet1_beta1_C3_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet1_beta1_C3_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet1_beta1_C3_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_C3_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_jet1_beta1_C3");
    h_jet1_beta1_C3_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_beta1_C3_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_beta1_C3_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet1_beta1_C3_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_C3_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_C3_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_beta1_C3_HZ_polm80_massCuts_.Rebin(2)
    h_jet1_beta1_C3_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet1_beta1_C3_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet1_beta1_C3_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet1_beta1_C3_polm80_massCuts_= THStack("hzqq_BG_jet1_beta1_C3_polm80_massCuts_", "");
    hzqq_BG_jet1_beta1_C3_polm80_massCuts_.Add(h_jet1_beta1_C3_ee_qq_polm80_massCuts_);
    hzqq_BG_jet1_beta1_C3_polm80_massCuts_.Add(h_jet1_beta1_C3_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet1_beta1_C3_polm80_massCuts_.Add(h_jet1_beta1_C3_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_beta1_C3_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_beta1_C3_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_beta1_C3_polm80_massCuts_thstack.cd();
    #h_jet1_beta1_C3_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_beta1_C3_BG.Draw("hist,e")
    #h_tot_norm_jet1_beta1_C3_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_beta1_C3_polm80_massCuts_.Draw("hist");
    hzqq_BG_jet1_beta1_C3_polm80_massCuts_.GetXaxis().SetRangeUser(0,0.625)
    hzqq_BG_jet1_beta1_C3_polm80_massCuts_.GetXaxis().SetTitle("jet1 C_{3}^{(1)}");
    hzqq_BG_jet1_beta1_C3_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetMaximum(150)
    h_jet1_beta1_C3_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_beta1_C3_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_=TLegend(0.50,0.63,0.85,0.88);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.AddEntry(h_jet1_beta1_C3_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.AddEntry(h_jet1_beta1_C3_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.AddEntry(h_jet1_beta1_C3_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.AddEntry(h_jet1_beta1_C3_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_beta1_C3_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet1_beta1_C3_polm80_massCuts_thstack.Print("h_jet1_beta1_C3_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")



    h_jet2_tau21_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_jet2_tau21");
    h_jet2_tau21_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 C_{3}^{(1)}");
    h_jet2_tau21_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_tau21_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_jet2_tau21");
    h_jet2_tau21_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet2 C_{3}^{(1)}");
    h_jet2_tau21_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet2_tau21_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet2_tau21_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet2_tau21_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet2_tau21_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_jet2_tau21");
    h_jet2_tau21_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet2_tau21_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet2_tau21_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_tau21_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_jet2_tau21");
    h_jet2_tau21_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet2_tau21_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet2_tau21_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet2_tau21_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_jet2_tau21");
    h_jet2_tau21_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_tau21_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_tau21_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet2_tau21_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_tau21_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_tau21_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet2_tau21_HZ_polm80_massCuts_.Rebin(2)
    h_jet2_tau21_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet2_tau21_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet2_tau21_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet2_tau21_polm80_massCuts_= THStack("hzqq_BG_jet2_tau21_polm80_massCuts_", "");
    hzqq_BG_jet2_tau21_polm80_massCuts_.Add(h_jet2_tau21_ee_qq_polm80_massCuts_);
    hzqq_BG_jet2_tau21_polm80_massCuts_.Add(h_jet2_tau21_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet2_tau21_polm80_massCuts_.Add(h_jet2_tau21_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_tau21_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_tau21_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_tau21_polm80_massCuts_thstack.cd();
    #h_jet2_tau21_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_tau21_BG.Draw("hist,e")
    #h_tot_norm_jet2_tau21_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_tau21_polm80_massCuts_.Draw("hist");
    #hzqq_BG_jet2_tau21_polm80_massCuts_.GetXaxis().SetRangeUser(0.,0.625)
    hzqq_BG_jet2_tau21_polm80_massCuts_.GetXaxis().SetTitle("jet2 #tau_{21}");
    hzqq_BG_jet2_tau21_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_tau21_polm80_massCuts_.SetMaximum(100)
    h_jet2_tau21_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_tau21_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_=TLegend(0.60,0.63,0.95,0.88);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.AddEntry(h_jet2_tau21_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.AddEntry(h_jet2_tau21_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.AddEntry(h_jet2_tau21_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.AddEntry(h_jet2_tau21_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_tau21_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet2_tau21_polm80_massCuts_thstack.Print("h_jet2_tau21_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet1_tau21_HZ_polm80_AllEvents_massCuts_=file_polm80_HZ_AllEvents_massCuts_.Get("h_jet1_tau21");
    h_jet1_tau21_HZ_polm80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 #tau_{21}");
    h_jet1_tau21_HZ_polm80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_tau21_HZ_polm80_massCuts_=file_polm80_HZ_massCuts_.Get("h_jet1_tau21");
    h_jet1_tau21_HZ_polm80_massCuts_.GetXaxis().SetTitle("jet1 #tau_{21}");
    h_jet1_tau21_HZ_polm80_massCuts_.SetLineWidth(3);
    h_jet1_tau21_HZ_polm80_massCuts_.SetFillColor(kWhite)
    h_jet1_tau21_HZ_polm80_massCuts_.SetLineColor(kBlack)
    #h_jet1_tau21_HZ_polm80_massCuts_.SetFillStyle(3001)

    h_jet1_tau21_ee_qq_polm80_massCuts_=file_polm80_ee_qq_massCuts_.Get("h_jet1_tau21");
    h_jet1_tau21_ee_qq_polm80_massCuts_.SetFillColor(kBlue);
    h_jet1_tau21_ee_qq_polm80_massCuts_.SetLineColor(kBlue);
    h_jet1_tau21_ee_qq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_tau21_ee_qqqq_polm80_massCuts_=file_polm80_ee_qqqq_massCuts_.Get("h_jet1_tau21");
    h_jet1_tau21_ee_qqqq_polm80_massCuts_.SetFillColor(kRed);
    h_jet1_tau21_ee_qqqq_polm80_massCuts_.SetLineColor(kRed);
    h_jet1_tau21_ee_qqqq_polm80_massCuts_.SetFillStyle(3002);
    h_jet1_tau21_ee_qqqqqq_polm80_massCuts_=file_polm80_ee_qqqqqq_massCuts_.Get("h_jet1_tau21");
    h_jet1_tau21_ee_qqqqqq_polm80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_tau21_ee_qqqqqq_polm80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_tau21_ee_qqqqqq_polm80_massCuts_.SetFillStyle(3002);

    h_jet1_tau21_ee_qq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_tau21_ee_qqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_tau21_ee_qqqqqq_polm80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_tau21_HZ_polm80_massCuts_.Rebin(2)
    h_jet1_tau21_ee_qq_polm80_massCuts_.Rebin(2)
    h_jet1_tau21_ee_qqqq_polm80_massCuts_.Rebin(2)
    h_jet1_tau21_ee_qqqqqq_polm80_massCuts_.Rebin(2)

    hzqq_BG_jet1_tau21_polm80_massCuts_= THStack("hzqq_BG_jet1_tau21_polm80_massCuts_", "");
    hzqq_BG_jet1_tau21_polm80_massCuts_.Add(h_jet1_tau21_ee_qq_polm80_massCuts_);
    hzqq_BG_jet1_tau21_polm80_massCuts_.Add(h_jet1_tau21_ee_qqqqqq_polm80_massCuts_);
    hzqq_BG_jet1_tau21_polm80_massCuts_.Add(h_jet1_tau21_ee_qqqq_polm80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_tau21_polm80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_tau21_polm80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_tau21_polm80_massCuts_thstack.cd();
    #h_jet1_tau21_HZ_polm80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_tau21_BG.Draw("hist,e")
    #h_tot_norm_jet1_tau21_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_tau21_polm80_massCuts_.Draw("hist");
    #hzqq_BG_jet1_tau21_polm80_massCuts_.GetXaxis().SetRangeUser(0,0.625)
    hzqq_BG_jet1_tau21_polm80_massCuts_.GetXaxis().SetTitle("jet1 #tau_{21}");
    hzqq_BG_jet1_tau21_polm80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_tau21_polm80_massCuts_.SetMaximum(100)
    h_jet1_tau21_HZ_polm80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_tau21_polm80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_=TLegend(0.60,0.63,0.90,0.88);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.AddEntry(h_jet1_tau21_ee_qqqqqq_polm80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.AddEntry(h_jet1_tau21_ee_qqqq_polm80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.AddEntry(h_jet1_tau21_ee_qq_polm80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.AddEntry(h_jet1_tau21_HZ_polm80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_tau21_polm80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_h_SIG_BG_jet1_tau21_polm80_massCuts_thstack.Print("h_jet1_tau21_polm80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")

    #now the BDT polp plots
    file_polp80_HZ_massCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/test_hzqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar_withSignalHistos.root")
    file_polp80_HZ_AllEvents_massCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/test_hzqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar_withSignalHistos_AllEvents.root")
    file_polp80_ee_qq_massCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/test_ee_qq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root")
    file_polp80_ee_qqqq_massCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/test_ee_qqqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root")
    file_polp80_ee_qqqqqq_massCuts_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/test_ee_qqqqqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root")

    h_jet1_mass_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_mass_jet1");
    h_jet1_mass_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_jet1_mass_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_mass_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_mass_jet1");
    h_jet1_mass_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet1 mass [GeV]");
    h_jet1_mass_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet1_mass_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet1_mass_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet1_mass_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet1_mass_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_mass_jet1");
    h_jet1_mass_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet1_mass_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet1_mass_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_mass_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_mass_jet1");
    h_jet1_mass_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet1_mass_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet1_mass_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_mass_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_mass_jet1");
    h_jet1_mass_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_mass_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_mass_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_tot_norm_jet1_mass_BG = h_jet1_mass_ee_qq_polp80_massCuts_.Clone("h_tot_BG_normalisation")
    h_tot_norm_jet1_mass_BG.Add(h_jet1_mass_ee_qqqq_polp80_massCuts_);
    h_tot_norm_jet1_mass_BG.Add(h_jet1_mass_ee_qqqqqq_polp80_massCuts_);
    norm_tot_BG_to_SIG=h_jet1_mass_HZ_polp80_massCuts_.Integral()/(h_jet1_mass_ee_qq_polp80_massCuts_.Integral()+h_jet1_mass_ee_qqqq_polp80_massCuts_.Integral()+h_jet1_mass_ee_qqqqqq_polp80_massCuts_.Integral())
    h_tot_norm_jet1_mass_BG.Scale(norm_tot_BG_to_SIG)
    h_tot_norm_jet1_mass_BG.SetLineColor(kBlack)
    h_tot_norm_jet1_mass_BG.SetFillColor(0)

    print 'scale or range ',h_tot_norm_jet1_mass_BG.Integral(),h_jet1_mass_HZ_polp80_massCuts_.Integral()

    h_jet1_mass_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_mass_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_mass_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)

    h_jet1_mass_ee_qq_polp80_massCuts_.GetXaxis().SetRangeUser(91,161)
    h_jet1_mass_ee_qqqq_polp80_massCuts_.GetXaxis().SetRangeUser(91,161)
    h_jet1_mass_ee_qqqqqq_polp80_massCuts_.GetXaxis().SetRangeUser(91,161)
    h_jet1_mass_HZ_polp80_massCuts_.GetXaxis().SetRangeUser(91,161)

    hzqq_BG_jet1_mass_polp80_massCuts_= THStack("hzqq_BG_jet1_mass_polp80_massCuts_", "");
    hzqq_BG_jet1_mass_polp80_massCuts_.Add(h_jet1_mass_ee_qq_polp80_massCuts_);
    hzqq_BG_jet1_mass_polp80_massCuts_.Add(h_jet1_mass_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet1_mass_polp80_massCuts_.Add(h_jet1_mass_ee_qqqq_polp80_massCuts_);
 
    h_jet1_mass_HZ_polp80_massCuts_.Rebin(2)
    h_jet1_mass_ee_qq_polp80_massCuts_.Rebin(2)
    h_jet1_mass_ee_qqqq_polp80_massCuts_.Rebin(2)
    h_jet1_mass_ee_qqqqqq_polp80_massCuts_.Rebin(2)


    canvas_h_SIG_BG_jet1_mass_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_mass_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_mass_polp80_massCuts_thstack.cd();
    #h_jet1_mass_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_mass_BG.Draw("hist,e")
    #h_tot_norm_jet1_mass_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_mass_polp80_massCuts_.Draw("hist");
    hzqq_BG_jet1_mass_polp80_massCuts_.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_mass_polp80_massCuts_.GetXaxis().SetTitle("jet1 mass [GeV]");
    hzqq_BG_jet1_mass_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_mass_polp80_massCuts_.SetMaximum(50)
    h_jet1_mass_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_mass_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_mass_polp80_massCuts_=TLegend(0.30,0.63,0.70,0.88);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.AddEntry(h_jet1_mass_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.AddEntry(h_jet1_mass_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.AddEntry(h_jet1_mass_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.AddEntry(h_jet1_mass_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_mass_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet1_mass_polp80_massCuts_thstack.Print("h_jet1_mass_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet2_mass_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_mass_jet2");
    h_jet2_mass_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 mass [GeV]");
    h_jet2_mass_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_mass_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_mass_jet2");
    h_jet2_mass_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet2 mass [GeV]");
    h_jet2_mass_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet2_mass_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet2_mass_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet2_mass_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet2_mass_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_mass_jet2");
    h_jet2_mass_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet2_mass_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet2_mass_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_mass_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_mass_jet2");
    h_jet2_mass_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet2_mass_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet2_mass_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_mass_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_mass_jet2");
    h_jet2_mass_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_mass_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_mass_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet2_mass_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_mass_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_mass_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)

    hzqq_BG_jet2_mass_polp80_massCuts_= THStack("hzqq_BG_jet2_mass_polp80_massCuts_", "");
    hzqq_BG_jet2_mass_polp80_massCuts_.Add(h_jet2_mass_ee_qq_polp80_massCuts_);
    hzqq_BG_jet2_mass_polp80_massCuts_.Add(h_jet2_mass_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet2_mass_polp80_massCuts_.Add(h_jet2_mass_ee_qqqq_polp80_massCuts_);
 
    h_jet2_mass_HZ_polp80_massCuts_.Rebin(2)
    h_jet2_mass_ee_qq_polp80_massCuts_.Rebin(2)
    h_jet2_mass_ee_qqqq_polp80_massCuts_.Rebin(2)
    h_jet2_mass_ee_qqqqqq_polp80_massCuts_.Rebin(2)


    canvas_h_SIG_BG_jet2_mass_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_mass_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_mass_polp80_massCuts_thstack.cd();
    #h_jet2_mass_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_mass_BG.Draw("hist,e")
    #h_tot_norm_jet2_mass_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_mass_polp80_massCuts_.Draw("hist");
    hzqq_BG_jet2_mass_polp80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet2_mass_polp80_massCuts_.GetXaxis().SetTitle("jet2 mass [GeV]");
    hzqq_BG_jet2_mass_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_mass_polp80_massCuts_.SetMaximum(50)
    h_jet2_mass_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_mass_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_mass_polp80_massCuts_=TLegend(0.20,0.63,0.55,0.88);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.AddEntry(h_jet2_mass_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.AddEntry(h_jet2_mass_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.AddEntry(h_jet2_mass_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.AddEntry(h_jet2_mass_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_mass_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet2_mass_polp80_massCuts_thstack.Print("h_jet2_mass_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")



    h_jet2_theta_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_theta_jet2");
    h_jet2_theta_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    h_jet2_theta_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_theta_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_theta_jet2");
    h_jet2_theta_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    h_jet2_theta_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet2_theta_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet2_theta_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet2_theta_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet2_theta_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_theta_jet2");
    h_jet2_theta_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet2_theta_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet2_theta_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_theta_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_theta_jet2");
    h_jet2_theta_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet2_theta_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet2_theta_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_theta_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_theta_jet2");
    h_jet2_theta_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_theta_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_theta_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet2_theta_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_theta_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_theta_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet2_theta_HZ_polp80_massCuts_.Rebin(4)
    h_jet2_theta_ee_qq_polp80_massCuts_.Rebin(4)
    h_jet2_theta_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_jet2_theta_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_jet2_theta_polp80_massCuts_= THStack("hzqq_BG_jet2_theta_polp80_massCuts_", "");
    hzqq_BG_jet2_theta_polp80_massCuts_.Add(h_jet2_theta_ee_qq_polp80_massCuts_);
    hzqq_BG_jet2_theta_polp80_massCuts_.Add(h_jet2_theta_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet2_theta_polp80_massCuts_.Add(h_jet2_theta_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_theta_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_theta_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_theta_polp80_massCuts_thstack.cd();
    #h_jet2_theta_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_theta_BG.Draw("hist,e")
    #h_tot_norm_jet2_theta_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_theta_polp80_massCuts_.Draw("hist");
    #hzqq_BG_jet2_theta_polp80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet2_theta_polp80_massCuts_.GetXaxis().SetTitle("jet2 #theta [#circ]");
    hzqq_BG_jet2_theta_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_theta_polp80_massCuts_.SetMaximum(40)
    h_jet2_theta_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_theta_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_theta_polp80_massCuts_=TLegend(0.20,0.63,0.55,0.88);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.AddEntry(h_jet2_theta_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.AddEntry(h_jet2_theta_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.AddEntry(h_jet2_theta_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.AddEntry(h_jet2_theta_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_theta_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet2_theta_polp80_massCuts_thstack.Print("h_jet2_theta_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet1_theta_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_theta_jet1");
    h_jet1_theta_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    h_jet1_theta_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_theta_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_theta_jet1");
    h_jet1_theta_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    h_jet1_theta_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet1_theta_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet1_theta_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet1_theta_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet1_theta_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_theta_jet1");
    h_jet1_theta_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet1_theta_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet1_theta_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_theta_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_theta_jet1");
    h_jet1_theta_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet1_theta_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet1_theta_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_theta_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_theta_jet1");
    h_jet1_theta_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_theta_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_theta_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet1_theta_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_theta_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_theta_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_theta_HZ_polp80_massCuts_.Rebin(4)
    h_jet1_theta_ee_qq_polp80_massCuts_.Rebin(4)
    h_jet1_theta_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_jet1_theta_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_jet1_theta_polp80_massCuts_= THStack("hzqq_BG_jet1_theta_polp80_massCuts_", "");
    hzqq_BG_jet1_theta_polp80_massCuts_.Add(h_jet1_theta_ee_qq_polp80_massCuts_);
    hzqq_BG_jet1_theta_polp80_massCuts_.Add(h_jet1_theta_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet1_theta_polp80_massCuts_.Add(h_jet1_theta_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_theta_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_theta_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_theta_polp80_massCuts_thstack.cd();
    #h_jet1_theta_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_theta_BG.Draw("hist,e")
    #h_tot_norm_jet1_theta_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_theta_polp80_massCuts_.Draw("hist");
    #hzqq_BG_jet1_theta_polp80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet1_theta_polp80_massCuts_.GetXaxis().SetTitle("jet1 #theta [#circ]");
    hzqq_BG_jet1_theta_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_theta_polp80_massCuts_.SetMaximum(40)
    h_jet1_theta_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_theta_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_theta_polp80_massCuts_=TLegend(0.20,0.63,0.55,0.88);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.AddEntry(h_jet1_theta_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.AddEntry(h_jet1_theta_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.AddEntry(h_jet1_theta_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.AddEntry(h_jet1_theta_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_theta_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet1_theta_polp80_massCuts_thstack.Print("h_jet1_theta_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    #defined as BTag of the larger BTagged subjet of jet1
    h_jet1_BTag_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_BTag_sj_BTagMax");
    h_jet1_BTag_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 BTag");
    h_jet1_BTag_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_BTag_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_BTag_sj_BTagMax");
    h_jet1_BTag_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet1 BTag");
    h_jet1_BTag_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet1_BTag_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet1_BTag_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet1_BTag_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet1_BTag_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_BTag_sj_BTagMax");
    h_jet1_BTag_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet1_BTag_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet1_BTag_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_BTag_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_BTag_sj_BTagMax");
    h_jet1_BTag_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet1_BTag_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet1_BTag_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_BTag_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_BTag_sj_BTagMax");
    h_jet1_BTag_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_BTag_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_BTag_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet1_BTag_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_BTag_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_BTag_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_BTag_HZ_polp80_massCuts_.Rebin(2)
    h_jet1_BTag_ee_qq_polp80_massCuts_.Rebin(2)
    h_jet1_BTag_ee_qqqq_polp80_massCuts_.Rebin(2)
    h_jet1_BTag_ee_qqqqqq_polp80_massCuts_.Rebin(2)

    hzqq_BG_jet1_BTag_polp80_massCuts_= THStack("hzqq_BG_jet1_BTag_polp80_massCuts_", "");
    hzqq_BG_jet1_BTag_polp80_massCuts_.Add(h_jet1_BTag_ee_qq_polp80_massCuts_);
    hzqq_BG_jet1_BTag_polp80_massCuts_.Add(h_jet1_BTag_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet1_BTag_polp80_massCuts_.Add(h_jet1_BTag_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_BTag_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_BTag_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_BTag_polp80_massCuts_thstack.cd();
    #h_jet1_BTag_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_BTag_BG.Draw("hist,e")
    #h_tot_norm_jet1_BTag_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_BTag_polp80_massCuts_.Draw("hist");
    #hzqq_BG_jet1_BTag_polp80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet1_BTag_polp80_massCuts_.GetXaxis().SetTitle("jet1 BTag");
    hzqq_BG_jet1_BTag_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_BTag_polp80_massCuts_.SetMaximum(60)
    h_jet1_BTag_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_BTag_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_=TLegend(0.30,0.63,0.65,0.88);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.AddEntry(h_jet1_BTag_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.AddEntry(h_jet1_BTag_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.AddEntry(h_jet1_BTag_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.AddEntry(h_jet1_BTag_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_BTag_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet1_BTag_polp80_massCuts_thstack.Print("h_jet1_BTag_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")



    h_reco_y32_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_reco_y32");
    h_reco_y32_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("y_{23}");
    h_reco_y32_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_reco_y32_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_reco_y32");
    h_reco_y32_HZ_polp80_massCuts_.GetXaxis().SetTitle("y_{23}");
    h_reco_y32_HZ_polp80_massCuts_.SetLineWidth(3);
    h_reco_y32_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_reco_y32_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_reco_y32_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_reco_y32_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_reco_y32");
    h_reco_y32_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_reco_y32_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_reco_y32_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_reco_y32_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_reco_y32");
    h_reco_y32_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_reco_y32_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_reco_y32_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_reco_y32_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_reco_y32");
    h_reco_y32_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_reco_y32_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_reco_y32_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_reco_y32_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_reco_y32_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_reco_y32_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_reco_y32_HZ_polp80_massCuts_.Rebin(4)
    h_reco_y32_ee_qq_polp80_massCuts_.Rebin(4)
    h_reco_y32_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_reco_y32_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_reco_y32_polp80_massCuts_= THStack("hzqq_BG_reco_y32_polp80_massCuts_", "");
    hzqq_BG_reco_y32_polp80_massCuts_.Add(h_reco_y32_ee_qq_polp80_massCuts_);
    hzqq_BG_reco_y32_polp80_massCuts_.Add(h_reco_y32_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_reco_y32_polp80_massCuts_.Add(h_reco_y32_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_reco_y32_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_reco_y32_polp80_massCuts_thstack");
    canvas_h_SIG_BG_reco_y32_polp80_massCuts_thstack.cd();
    #h_reco_y32_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_reco_y32_BG.Draw("hist,e")
    #h_tot_norm_reco_y32_BG.GetXaxis().SetRangeUser(0,0.0045)
    hzqq_BG_reco_y32_polp80_massCuts_.Draw("hist");
    hzqq_BG_reco_y32_polp80_massCuts_.GetXaxis().SetTitle("y_{23}")
    hzqq_BG_reco_y32_polp80_massCuts_.GetXaxis().SetRangeUser(0,0.0045)
    hzqq_BG_reco_y32_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_reco_y32_polp80_massCuts_.SetMaximum(60)
    h_reco_y32_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_reco_y32_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_reco_y32_polp80_massCuts_=TLegend(0.50,0.63,0.85,0.88);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_reco_y32_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_reco_y32_polp80_massCuts_.AddEntry(h_reco_y32_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_reco_y32_polp80_massCuts_.AddEntry(h_reco_y32_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_reco_y32_polp80_massCuts_.AddEntry(h_reco_y32_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_reco_y32_polp80_massCuts_.AddEntry(h_reco_y32_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_reco_y32_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_reco_y32_polp80_massCuts_thstack.Print("h_reco_y32_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet2_beta1_D2_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_jet2_beta1_D2");
    h_jet2_beta1_D2_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 D_{2}^{(1)}");
    h_jet2_beta1_D2_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_beta1_D2_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_jet2_beta1_D2");
    h_jet2_beta1_D2_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet2 D_{2}^{(1)}");
    h_jet2_beta1_D2_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet2_beta1_D2_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet2_beta1_D2_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet2_beta1_D2_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet2_beta1_D2_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_jet2_beta1_D2");
    h_jet2_beta1_D2_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet2_beta1_D2_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet2_beta1_D2_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_D2_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_jet2_beta1_D2");
    h_jet2_beta1_D2_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet2_beta1_D2_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet2_beta1_D2_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_D2_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_jet2_beta1_D2");
    h_jet2_beta1_D2_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_beta1_D2_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_beta1_D2_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet2_beta1_D2_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_D2_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_D2_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet2_beta1_D2_HZ_polp80_massCuts_.Rebin(4)
    h_jet2_beta1_D2_ee_qq_polp80_massCuts_.Rebin(4)
    h_jet2_beta1_D2_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_jet2_beta1_D2_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_jet2_beta1_D2_polp80_massCuts_= THStack("hzqq_BG_jet2_beta1_D2_polp80_massCuts_", "");
    hzqq_BG_jet2_beta1_D2_polp80_massCuts_.Add(h_jet2_beta1_D2_ee_qq_polp80_massCuts_);
    hzqq_BG_jet2_beta1_D2_polp80_massCuts_.Add(h_jet2_beta1_D2_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet2_beta1_D2_polp80_massCuts_.Add(h_jet2_beta1_D2_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_beta1_D2_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_beta1_D2_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_beta1_D2_polp80_massCuts_thstack.cd();
    #h_jet2_beta1_D2_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_beta1_D2_BG.Draw("hist,e")
    #h_tot_norm_jet2_beta1_D2_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_beta1_D2_polp80_massCuts_.Draw("hist");
    #hzqq_BG_jet2_beta1_D2_polp80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet2_beta1_D2_polp80_massCuts_.GetXaxis().SetTitle("jet2 D_{2}^{(1)}");
    hzqq_BG_jet2_beta1_D2_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetMaximum(80)
    h_jet2_beta1_D2_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_beta1_D2_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_=TLegend(0.40,0.63,0.75,0.88);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.AddEntry(h_jet2_beta1_D2_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.AddEntry(h_jet2_beta1_D2_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.AddEntry(h_jet2_beta1_D2_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.AddEntry(h_jet2_beta1_D2_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_beta1_D2_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet2_beta1_D2_polp80_massCuts_thstack.Print("h_jet2_beta1_D2_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet1_beta1_D2_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_jet1_beta1_D2");
    h_jet1_beta1_D2_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 D_{2}^{(1)}");
    h_jet1_beta1_D2_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_beta1_D2_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_jet1_beta1_D2");
    h_jet1_beta1_D2_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet1 D_{2}^{(1)}");
    h_jet1_beta1_D2_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet1_beta1_D2_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet1_beta1_D2_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet1_beta1_D2_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet1_beta1_D2_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_jet1_beta1_D2");
    h_jet1_beta1_D2_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet1_beta1_D2_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet1_beta1_D2_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_D2_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_jet1_beta1_D2");
    h_jet1_beta1_D2_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet1_beta1_D2_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet1_beta1_D2_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_D2_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_jet1_beta1_D2");
    h_jet1_beta1_D2_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_beta1_D2_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_beta1_D2_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet1_beta1_D2_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_D2_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_D2_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_beta1_D2_HZ_polp80_massCuts_.Rebin(4)
    h_jet1_beta1_D2_ee_qq_polp80_massCuts_.Rebin(4)
    h_jet1_beta1_D2_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_jet1_beta1_D2_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_jet1_beta1_D2_polp80_massCuts_= THStack("hzqq_BG_jet1_beta1_D2_polp80_massCuts_", "");
    hzqq_BG_jet1_beta1_D2_polp80_massCuts_.Add(h_jet1_beta1_D2_ee_qq_polp80_massCuts_);
    hzqq_BG_jet1_beta1_D2_polp80_massCuts_.Add(h_jet1_beta1_D2_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet1_beta1_D2_polp80_massCuts_.Add(h_jet1_beta1_D2_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_beta1_D2_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_beta1_D2_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_beta1_D2_polp80_massCuts_thstack.cd();
    #h_jet1_beta1_D2_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_beta1_D2_BG.Draw("hist,e")
    #h_tot_norm_jet1_beta1_D2_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_beta1_D2_polp80_massCuts_.Draw("hist");
    #hzqq_BG_jet1_beta1_D2_polp80_massCuts_.GetXaxis().SetRangeUser(57,127)
    hzqq_BG_jet1_beta1_D2_polp80_massCuts_.GetXaxis().SetTitle("jet1 D_{2}^{(1)}");
    hzqq_BG_jet1_beta1_D2_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetMaximum(80)
    h_jet1_beta1_D2_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_beta1_D2_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_=TLegend(0.40,0.63,0.75,0.88);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.AddEntry(h_jet1_beta1_D2_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.AddEntry(h_jet1_beta1_D2_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.AddEntry(h_jet1_beta1_D2_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.AddEntry(h_jet1_beta1_D2_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_beta1_D2_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet1_beta1_D2_polp80_massCuts_thstack.Print("h_jet1_beta1_D2_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")







    h_jet2_beta1_C2_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_jet2_beta1_C2");
    h_jet2_beta1_C2_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 C_{2}^{(1)}");
    h_jet2_beta1_C2_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_beta1_C2_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_jet2_beta1_C2");
    h_jet2_beta1_C2_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet2 C_{2}^{(1)}");
    h_jet2_beta1_C2_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet2_beta1_C2_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet2_beta1_C2_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet2_beta1_C2_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet2_beta1_C2_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_jet2_beta1_C2");
    h_jet2_beta1_C2_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet2_beta1_C2_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet2_beta1_C2_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_C2_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_jet2_beta1_C2");
    h_jet2_beta1_C2_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet2_beta1_C2_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet2_beta1_C2_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_C2_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_jet2_beta1_C2");
    h_jet2_beta1_C2_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_beta1_C2_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_beta1_C2_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet2_beta1_C2_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_C2_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_C2_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet2_beta1_C2_HZ_polp80_massCuts_.Rebin(4)
    h_jet2_beta1_C2_ee_qq_polp80_massCuts_.Rebin(4)
    h_jet2_beta1_C2_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_jet2_beta1_C2_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_jet2_beta1_C2_polp80_massCuts_= THStack("hzqq_BG_jet2_beta1_C2_polp80_massCuts_", "");
    hzqq_BG_jet2_beta1_C2_polp80_massCuts_.Add(h_jet2_beta1_C2_ee_qq_polp80_massCuts_);
    hzqq_BG_jet2_beta1_C2_polp80_massCuts_.Add(h_jet2_beta1_C2_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet2_beta1_C2_polp80_massCuts_.Add(h_jet2_beta1_C2_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_beta1_C2_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_beta1_C2_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_beta1_C2_polp80_massCuts_thstack.cd();
    #h_jet2_beta1_C2_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_beta1_C2_BG.Draw("hist,e")
    #h_tot_norm_jet2_beta1_C2_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_beta1_C2_polp80_massCuts_.Draw("hist");
    hzqq_BG_jet2_beta1_C2_polp80_massCuts_.GetXaxis().SetRangeUser(0,0.40)
    hzqq_BG_jet2_beta1_C2_polp80_massCuts_.GetXaxis().SetTitle("jet2 C_{2}^{(1)}");
    hzqq_BG_jet2_beta1_C2_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetMaximum(100)
    h_jet2_beta1_C2_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_beta1_C2_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_=TLegend(0.40,0.63,0.75,0.88);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.AddEntry(h_jet2_beta1_C2_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.AddEntry(h_jet2_beta1_C2_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.AddEntry(h_jet2_beta1_C2_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.AddEntry(h_jet2_beta1_C2_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_beta1_C2_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet2_beta1_C2_polp80_massCuts_thstack.Print("h_jet2_beta1_C2_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet1_beta1_C2_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_jet1_beta1_C2");
    h_jet1_beta1_C2_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 C_{2}^{(1)}");
    h_jet1_beta1_C2_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_beta1_C2_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_jet1_beta1_C2");
    h_jet1_beta1_C2_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet1 C_{2}^{(1)}");
    h_jet1_beta1_C2_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet1_beta1_C2_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet1_beta1_C2_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet1_beta1_C2_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet1_beta1_C2_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_jet1_beta1_C2");
    h_jet1_beta1_C2_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet1_beta1_C2_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet1_beta1_C2_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_C2_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_jet1_beta1_C2");
    h_jet1_beta1_C2_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet1_beta1_C2_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet1_beta1_C2_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_C2_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_jet1_beta1_C2");
    h_jet1_beta1_C2_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_beta1_C2_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_beta1_C2_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet1_beta1_C2_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_C2_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_C2_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_beta1_C2_HZ_polp80_massCuts_.Rebin(4)
    h_jet1_beta1_C2_ee_qq_polp80_massCuts_.Rebin(4)
    h_jet1_beta1_C2_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_jet1_beta1_C2_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_jet1_beta1_C2_polp80_massCuts_= THStack("hzqq_BG_jet1_beta1_C2_polp80_massCuts_", "");
    hzqq_BG_jet1_beta1_C2_polp80_massCuts_.Add(h_jet1_beta1_C2_ee_qq_polp80_massCuts_);
    hzqq_BG_jet1_beta1_C2_polp80_massCuts_.Add(h_jet1_beta1_C2_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet1_beta1_C2_polp80_massCuts_.Add(h_jet1_beta1_C2_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_beta1_C2_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_beta1_C2_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_beta1_C2_polp80_massCuts_thstack.cd();
    #h_jet1_beta1_C2_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_beta1_C2_BG.Draw("hist,e")
    #h_tot_norm_jet1_beta1_C2_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_beta1_C2_polp80_massCuts_.Draw("hist");
    hzqq_BG_jet1_beta1_C2_polp80_massCuts_.GetXaxis().SetRangeUser(0,0.40)
    hzqq_BG_jet1_beta1_C2_polp80_massCuts_.GetXaxis().SetTitle("jet1 C_{2}^{(1)}");
    hzqq_BG_jet1_beta1_C2_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetMaximum(80)
    h_jet1_beta1_C2_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_beta1_C2_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_=TLegend(0.40,0.63,0.75,0.88);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.AddEntry(h_jet1_beta1_C2_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.AddEntry(h_jet1_beta1_C2_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.AddEntry(h_jet1_beta1_C2_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.AddEntry(h_jet1_beta1_C2_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_beta1_C2_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet1_beta1_C2_polp80_massCuts_thstack.Print("h_jet1_beta1_C2_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")







    h_jet2_beta1_C3_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_jet2_beta1_C3");
    h_jet2_beta1_C3_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 C_{3}^{(1)}");
    h_jet2_beta1_C3_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_beta1_C3_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_jet2_beta1_C3");
    h_jet2_beta1_C3_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet2 C_{3}^{(1)}");
    h_jet2_beta1_C3_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet2_beta1_C3_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet2_beta1_C3_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet2_beta1_C3_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet2_beta1_C3_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_jet2_beta1_C3");
    h_jet2_beta1_C3_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet2_beta1_C3_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet2_beta1_C3_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_C3_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_jet2_beta1_C3");
    h_jet2_beta1_C3_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet2_beta1_C3_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet2_beta1_C3_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_beta1_C3_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_jet2_beta1_C3");
    h_jet2_beta1_C3_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_beta1_C3_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_beta1_C3_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet2_beta1_C3_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_C3_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_beta1_C3_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet2_beta1_C3_HZ_polp80_massCuts_.Rebin(4)
    h_jet2_beta1_C3_ee_qq_polp80_massCuts_.Rebin(4)
    h_jet2_beta1_C3_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_jet2_beta1_C3_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_jet2_beta1_C3_polp80_massCuts_= THStack("hzqq_BG_jet2_beta1_C3_polp80_massCuts_", "");
    hzqq_BG_jet2_beta1_C3_polp80_massCuts_.Add(h_jet2_beta1_C3_ee_qq_polp80_massCuts_);
    hzqq_BG_jet2_beta1_C3_polp80_massCuts_.Add(h_jet2_beta1_C3_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet2_beta1_C3_polp80_massCuts_.Add(h_jet2_beta1_C3_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_beta1_C3_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_beta1_C3_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_beta1_C3_polp80_massCuts_thstack.cd();
    #h_jet2_beta1_C3_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_beta1_C3_BG.Draw("hist,e")
    #h_tot_norm_jet2_beta1_C3_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_beta1_C3_polp80_massCuts_.Draw("hist");
    hzqq_BG_jet2_beta1_C3_polp80_massCuts_.GetXaxis().SetRangeUser(0.,0.625)
    hzqq_BG_jet2_beta1_C3_polp80_massCuts_.GetXaxis().SetTitle("jet2 C_{3}^{(1)}");
    hzqq_BG_jet2_beta1_C3_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetMaximum(60)
    h_jet2_beta1_C3_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_beta1_C3_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_=TLegend(0.40,0.63,0.75,0.88);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.AddEntry(h_jet2_beta1_C3_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.AddEntry(h_jet2_beta1_C3_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.AddEntry(h_jet2_beta1_C3_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.AddEntry(h_jet2_beta1_C3_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_beta1_C3_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet2_beta1_C3_polp80_massCuts_thstack.Print("h_jet2_beta1_C3_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet1_beta1_C3_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_jet1_beta1_C3");
    h_jet1_beta1_C3_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 C_{3}^{(1)}");
    h_jet1_beta1_C3_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_beta1_C3_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_jet1_beta1_C3");
    h_jet1_beta1_C3_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet1 C_{3}^{(1)}");
    h_jet1_beta1_C3_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet1_beta1_C3_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet1_beta1_C3_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet1_beta1_C3_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet1_beta1_C3_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_jet1_beta1_C3");
    h_jet1_beta1_C3_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet1_beta1_C3_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet1_beta1_C3_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_C3_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_jet1_beta1_C3");
    h_jet1_beta1_C3_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet1_beta1_C3_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet1_beta1_C3_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_beta1_C3_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_jet1_beta1_C3");
    h_jet1_beta1_C3_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_beta1_C3_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_beta1_C3_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet1_beta1_C3_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_C3_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_beta1_C3_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_beta1_C3_HZ_polp80_massCuts_.Rebin(4)
    h_jet1_beta1_C3_ee_qq_polp80_massCuts_.Rebin(4)
    h_jet1_beta1_C3_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_jet1_beta1_C3_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_jet1_beta1_C3_polp80_massCuts_= THStack("hzqq_BG_jet1_beta1_C3_polp80_massCuts_", "");
    hzqq_BG_jet1_beta1_C3_polp80_massCuts_.Add(h_jet1_beta1_C3_ee_qq_polp80_massCuts_);
    hzqq_BG_jet1_beta1_C3_polp80_massCuts_.Add(h_jet1_beta1_C3_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet1_beta1_C3_polp80_massCuts_.Add(h_jet1_beta1_C3_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_beta1_C3_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_beta1_C3_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_beta1_C3_polp80_massCuts_thstack.cd();
    #h_jet1_beta1_C3_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_beta1_C3_BG.Draw("hist,e")
    #h_tot_norm_jet1_beta1_C3_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_beta1_C3_polp80_massCuts_.Draw("hist");
    hzqq_BG_jet1_beta1_C3_polp80_massCuts_.GetXaxis().SetRangeUser(0,0.625)
    hzqq_BG_jet1_beta1_C3_polp80_massCuts_.GetXaxis().SetTitle("jet1 C_{3}^{(1)}");
    hzqq_BG_jet1_beta1_C3_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetMaximum(60)
    h_jet1_beta1_C3_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_beta1_C3_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_=TLegend(0.50,0.63,0.85,0.88);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.AddEntry(h_jet1_beta1_C3_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.AddEntry(h_jet1_beta1_C3_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.AddEntry(h_jet1_beta1_C3_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.AddEntry(h_jet1_beta1_C3_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_beta1_C3_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet1_beta1_C3_polp80_massCuts_thstack.Print("h_jet1_beta1_C3_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")



    h_jet2_tau21_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_jet2_tau21");
    h_jet2_tau21_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet2 C_{3}^{(1)}");
    h_jet2_tau21_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet2_tau21_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_jet2_tau21");
    h_jet2_tau21_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet2 C_{3}^{(1)}");
    h_jet2_tau21_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet2_tau21_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet2_tau21_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet2_tau21_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet2_tau21_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_jet2_tau21");
    h_jet2_tau21_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet2_tau21_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet2_tau21_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_tau21_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_jet2_tau21");
    h_jet2_tau21_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet2_tau21_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet2_tau21_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet2_tau21_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_jet2_tau21");
    h_jet2_tau21_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet2_tau21_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet2_tau21_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet2_tau21_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_tau21_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet2_tau21_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet2_tau21_HZ_polp80_massCuts_.Rebin(4)
    h_jet2_tau21_ee_qq_polp80_massCuts_.Rebin(4)
    h_jet2_tau21_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_jet2_tau21_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_jet2_tau21_polp80_massCuts_= THStack("hzqq_BG_jet2_tau21_polp80_massCuts_", "");
    hzqq_BG_jet2_tau21_polp80_massCuts_.Add(h_jet2_tau21_ee_qq_polp80_massCuts_);
    hzqq_BG_jet2_tau21_polp80_massCuts_.Add(h_jet2_tau21_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet2_tau21_polp80_massCuts_.Add(h_jet2_tau21_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet2_tau21_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet2_tau21_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet2_tau21_polp80_massCuts_thstack.cd();
    #h_jet2_tau21_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet2_tau21_BG.Draw("hist,e")
    #h_tot_norm_jet2_tau21_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet2_tau21_polp80_massCuts_.Draw("hist");
    #hzqq_BG_jet2_tau21_polp80_massCuts_.GetXaxis().SetRangeUser(0.,0.625)
    hzqq_BG_jet2_tau21_polp80_massCuts_.GetXaxis().SetTitle("jet2 #tau_{21}");
    hzqq_BG_jet2_tau21_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet2_tau21_polp80_massCuts_.SetMaximum(40)
    h_jet2_tau21_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet2_tau21_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_=TLegend(0.60,0.63,0.95,0.88);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.AddEntry(h_jet2_tau21_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.AddEntry(h_jet2_tau21_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.AddEntry(h_jet2_tau21_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.AddEntry(h_jet2_tau21_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet2_tau21_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);
    canvas_h_SIG_BG_jet2_tau21_polp80_massCuts_thstack.Print("h_jet2_tau21_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")


    h_jet1_tau21_HZ_polp80_AllEvents_massCuts_=file_polp80_HZ_AllEvents_massCuts_.Get("h_jet1_tau21");
    h_jet1_tau21_HZ_polp80_AllEvents_massCuts_.GetXaxis().SetTitle("jet1 #tau_{21}");
    h_jet1_tau21_HZ_polp80_AllEvents_massCuts_.SetLineWidth(3);
    h_jet1_tau21_HZ_polp80_massCuts_=file_polp80_HZ_massCuts_.Get("h_jet1_tau21");
    h_jet1_tau21_HZ_polp80_massCuts_.GetXaxis().SetTitle("jet1 #tau_{21}");
    h_jet1_tau21_HZ_polp80_massCuts_.SetLineWidth(3);
    h_jet1_tau21_HZ_polp80_massCuts_.SetFillColor(kWhite)
    h_jet1_tau21_HZ_polp80_massCuts_.SetLineColor(kBlack)
    #h_jet1_tau21_HZ_polp80_massCuts_.SetFillStyle(3001)

    h_jet1_tau21_ee_qq_polp80_massCuts_=file_polp80_ee_qq_massCuts_.Get("h_jet1_tau21");
    h_jet1_tau21_ee_qq_polp80_massCuts_.SetFillColor(kBlue);
    h_jet1_tau21_ee_qq_polp80_massCuts_.SetLineColor(kBlue);
    h_jet1_tau21_ee_qq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_tau21_ee_qqqq_polp80_massCuts_=file_polp80_ee_qqqq_massCuts_.Get("h_jet1_tau21");
    h_jet1_tau21_ee_qqqq_polp80_massCuts_.SetFillColor(kRed);
    h_jet1_tau21_ee_qqqq_polp80_massCuts_.SetLineColor(kRed);
    h_jet1_tau21_ee_qqqq_polp80_massCuts_.SetFillStyle(3002);
    h_jet1_tau21_ee_qqqqqq_polp80_massCuts_=file_polp80_ee_qqqqqq_massCuts_.Get("h_jet1_tau21");
    h_jet1_tau21_ee_qqqqqq_polp80_massCuts_.SetFillColor(kGreen-2);
    h_jet1_tau21_ee_qqqqqq_polp80_massCuts_.SetLineColor(kGreen-2);
    h_jet1_tau21_ee_qqqqqq_polp80_massCuts_.SetFillStyle(3002);

    h_jet1_tau21_ee_qq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_tau21_ee_qqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
    h_jet1_tau21_ee_qqqqqq_polp80_massCuts_.Scale(norm_tot_BG_to_SIG)
 
    h_jet1_tau21_HZ_polp80_massCuts_.Rebin(4)
    h_jet1_tau21_ee_qq_polp80_massCuts_.Rebin(4)
    h_jet1_tau21_ee_qqqq_polp80_massCuts_.Rebin(4)
    h_jet1_tau21_ee_qqqqqq_polp80_massCuts_.Rebin(4)

    hzqq_BG_jet1_tau21_polp80_massCuts_= THStack("hzqq_BG_jet1_tau21_polp80_massCuts_", "");
    hzqq_BG_jet1_tau21_polp80_massCuts_.Add(h_jet1_tau21_ee_qq_polp80_massCuts_);
    hzqq_BG_jet1_tau21_polp80_massCuts_.Add(h_jet1_tau21_ee_qqqqqq_polp80_massCuts_);
    hzqq_BG_jet1_tau21_polp80_massCuts_.Add(h_jet1_tau21_ee_qqqq_polp80_massCuts_);
 
    
    canvas_h_SIG_BG_jet1_tau21_polp80_massCuts_thstack = setUpperCanvas("canvas_h_SIG_BG_jet1_tau21_polp80_massCuts_thstack");
    canvas_h_SIG_BG_jet1_tau21_polp80_massCuts_thstack.cd();
    #h_jet1_tau21_HZ_polp80_massCuts_.Draw("hist,e")
    #h_tot_norm_jet1_tau21_BG.Draw("hist,e")
    #h_tot_norm_jet1_tau21_BG.GetXaxis().SetRangeUser(91,161)
    hzqq_BG_jet1_tau21_polp80_massCuts_.Draw("hist");
    #hzqq_BG_jet1_tau21_polp80_massCuts_.GetXaxis().SetRangeUser(0,0.625)
    hzqq_BG_jet1_tau21_polp80_massCuts_.GetXaxis().SetTitle("jet1 #tau_{21}");
    hzqq_BG_jet1_tau21_polp80_massCuts_.GetYaxis().SetTitle("A.U.");
    hzqq_BG_jet1_tau21_polp80_massCuts_.SetMaximum(40)
    h_jet1_tau21_HZ_polp80_massCuts_.Draw("hist,e,same")
    canvas_h_SIG_BG_jet1_tau21_polp80_massCuts_thstack.Modified();
    
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_=TLegend(0.60,0.63,0.90,0.88);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetBorderSize(0);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetTextAlign(12);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetTextSize(0.050);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetTextFont(42);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetMargin(0.15);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetLineColor(1);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetLineStyle(1);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetLineWidth(1);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetFillColor(0);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetFillStyle(0);
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.SetHeader("#sqrt{s}>2500 GeV");
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.AddEntry(h_jet1_tau21_ee_qqqqqq_polp80_massCuts_,"ee#rightarrow qqqqqq");
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.AddEntry(h_jet1_tau21_ee_qqqq_polp80_massCuts_,"ee#rightarrow qqqq");
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.AddEntry(h_jet1_tau21_ee_qq_polp80_massCuts_,"ee#rightarrow qq");
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.AddEntry(h_jet1_tau21_HZ_polp80_massCuts_,"HZ");
    leg_hzqq_BG_jet1_tau21_polp80_massCuts_.Draw();
    
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label3);

    canvas_h_SIG_BG_jet1_tau21_polp80_massCuts_thstack.Print("h_jet1_tau21_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")
    
    file_polm80_HZ_BDT_Unfolding_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21_Unfolding.root")
    h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_=file_polm80_HZ_BDT_Unfolding_.Get("0.375/h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM");
    h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetXaxis().SetTitle("cos #theta_{2}^{parton}")
    h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetYaxis().SetTitle("cos #theta_{2}^{reco}")


    canvas_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80 = setUpperCanvas("parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80");
    canvas_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80.cd();
    #gPad.SetRightMargin(0.18)
    h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.Draw("col")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi,y2_noLumi,label2_noLumi);
    canvas_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80.Print("h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80.eps")

    #normed transpose matrix
    h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_=file_polm80_HZ_BDT_Unfolding_.Get("0.375/h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM");

    norm_polm80_costheta2_H_epgen=[]
    norm_polm80_costheta2_H_eptest=[]
    for i in range(0,h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetNbinsX()+1):
        norm_polm80_costheta2_H_epgen.append(0.)
        norm_polm80_costheta2_H_eptest.append(0)
        for j in range(0,h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetNbinsX()+1):
            norm_polm80_costheta2_H_epgen[i]+=h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetBinContent(i,j)
            h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.SetBinContent(j,i,h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetBinContent(i,j))
        print("norm costheta2",i,norm_polm80_costheta2_H_epgen[i])


    for i in range(0,h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetNbinsX()+1):  
        for j in range(0,h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetNbinsX()+1):
            if norm_polm80_costheta2_H_epgen[j]!=0:
                h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.SetBinContent(i,j,h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetBinContent(i,j)/norm_polm80_costheta2_H_epgen[j])
            norm_polm80_costheta2_H_eptest[j]+=h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetBinContent(i,j)

    print("norm test costheta2",norm_polm80_costheta2_H_eptest)

    #h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_

    h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetXaxis().SetTitle("cos #theta_{2}^{reco}")
    h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.GetYaxis().SetTitle("cos #theta_{2}^{parton}")

    canvas_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80 = setUpperCanvasWColz("normed_parton_recojet_costheta2_rj1_ep_E_totCOM_vs_costheta2_H_ep_approx_HZ_COM__polm80");
    canvas_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80.cd();
    #gPad.SetRightMargin(0.18)
    h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80_.Draw("colz")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi_normed,y2_noLumi,label2_noLumi);
    canvas_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polm80.Print("h_recojet_costheta2_rj1_ep_E_totCOM_vs_parton_costheta2_H_ep_approx_HZ_COM__polm80_normed.eps")



    h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_=file_polm80_HZ_BDT_Unfolding_.Get("0.375/h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetXaxis().SetTitle("cos #theta_{1}^{parton}")
    h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetYaxis().SetTitle("cos #theta_{1}^{reco}")
    canvas_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80 = setUpperCanvas("parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80");
    canvas_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80.cd();
    #gPad.SetRightMargin(0.18)
    h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.Draw("col")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi,y2_noLumi,label2_noLumi);
    canvas_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80.Print("h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80.eps")

    #now the normed transposed transfer matrix
    h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_=file_polm80_HZ_BDT_Unfolding_.Get("0.375/h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");


    norm_polm80_costheta1_Z_qpos_Zcomgen=[]
    norm_polm80_costheta1_Z_qpos_Zcomtest=[]
    for i in range(0,h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetNbinsX()+1):
        norm_polm80_costheta1_Z_qpos_Zcomgen.append(0.)
        norm_polm80_costheta1_Z_qpos_Zcomtest.append(0)
        for j in range(0,h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetNbinsX()+1):
            norm_polm80_costheta1_Z_qpos_Zcomgen[i]+=h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetBinContent(i,j)
            h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.SetBinContent(j,i,h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetBinContent(i,j))
        print("norm costheta2",i,norm_polm80_costheta1_Z_qpos_Zcomgen[i])


    for i in range(0,h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetNbinsX()+1):  
        for j in range(0,h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetNbinsX()+1):
            if norm_polm80_costheta1_Z_qpos_Zcomgen[j]!=0:
                h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.SetBinContent(i,j,h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetBinContent(i,j)/norm_polm80_costheta1_Z_qpos_Zcomgen[j])
            norm_polm80_costheta1_Z_qpos_Zcomtest[j]+=h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetBinContent(i,j)

    print("norm test costheta2",norm_polm80_costheta1_Z_qpos_Zcomtest)

    #h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_

    h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetXaxis().SetTitle("cos #theta_{1}^{reco}")
    h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.GetYaxis().SetTitle("cos #theta_{1}^{parton}")

    canvas_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80 = setUpperCanvasWColz("normed_parton_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2comvs_costheta1_Z_qpos_Zcom__polm80");
    canvas_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80.cd();
    #gPad.SetRightMargin(0.18)
    h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.SetMaximum(0.40);
    h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80_.Draw("colz")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi_normed,y2_noLumi,label2_noLumi);
    canvas_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polm80.Print("h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2comvs_parton_costheta1_Z_qpos_Zcom__polm80_normed.eps")




    h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_=file_polm80_HZ_BDT_Unfolding_.Get("0.375/h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");
    h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetXaxis().SetTitle("#phi^{parton} [#circ]")
    h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetYaxis().SetTitle("#phi^{reco} [#circ]")
    canvas_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80 = setUpperCanvas("parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80");
    canvas_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80.cd();
    #gPad.SetRightMargin(0.18)
    h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.Draw("col")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi,y2_noLumi,label2_noLumi);
    canvas_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80.Print("h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80.eps")


    #now the normed transposed transfer matrix
    h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_=file_polm80_HZ_BDT_Unfolding_.Get("0.375/h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");


    norm_polm80_costheta1_Z_qpos_Zcomgen=[]
    norm_polm80_costheta1_Z_qpos_Zcomtest=[]
    for i in range(0,h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetNbinsX()+1):
        norm_polm80_costheta1_Z_qpos_Zcomgen.append(0.)
        norm_polm80_costheta1_Z_qpos_Zcomtest.append(0)
        for j in range(0,h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetNbinsX()+1):
            norm_polm80_costheta1_Z_qpos_Zcomgen[i]+=h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetBinContent(i,j)
            h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.SetBinContent(j,i,h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetBinContent(i,j))
        print("norm phi",i,norm_polm80_costheta1_Z_qpos_Zcomgen[i])


    for i in range(0,h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetNbinsX()+1):  
        for j in range(0,h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetNbinsX()+1):
            if norm_polm80_costheta1_Z_qpos_Zcomgen[j]!=0:
                h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.SetBinContent(i,j,h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetBinContent(i,j)/norm_polm80_costheta1_Z_qpos_Zcomgen[j])
            norm_polm80_costheta1_Z_qpos_Zcomtest[j]+=h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetBinContent(i,j)

    print("norm test phi",norm_polm80_costheta1_Z_qpos_Zcomtest)

    #h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_

    h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetXaxis().SetTitle("#phi^{reco} [#circ]")
    h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.GetYaxis().SetTitle("#phi^{parton} [#circ]")

    canvas_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80 = setUpperCanvasWColz("normed_parton_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_epvs_costheta1_Z_qpos_Zcom__polm80");
    canvas_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80.cd();
    #gPad.SetRightMargin(0.18)
    h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80_.Draw("colz")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi_normed,y2_noLumi,label2_noLumi);
    canvas_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polm80.Print("h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_vs_parton_phi_plane_Z_qpos_vs_plane_H_ep__polm80_normed.eps")





    file_polp80_HZ_BDT_Unfolding_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts_m1_AnglesMETProj_eprealcalc_METOverSqrtSOrig_1_00_sqrtSOrig_0_BDT_y32_jet1_D2_C2_C3_tau21_jet2_D2_C2_C3_tau21_Unfolding.root")
    h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_=file_polp80_HZ_BDT_Unfolding_.Get("0.3/h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM");
    h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetXaxis().SetTitle("cos #theta_{2}^{parton}");
    h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetYaxis().SetTitle("cos #theta_{2}^{reco}");


    canvas_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80 = setUpperCanvas("parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80");
    canvas_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80.cd();
    #gPad.SetRightMargin(0.05)
    h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.Draw("col")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi,y2_noLumi,label3_noLumi);
    canvas_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80.Print("h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80.eps")

    #normed transpose matrix
    h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_=file_polp80_HZ_BDT_Unfolding_.Get("0.375/h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM");


    norm_polp80_costheta2_H_epgen=[]
    norm_polp80_costheta2_H_eptest=[]
    for i in range(0,h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetNbinsX()+1):
        norm_polp80_costheta2_H_epgen.append(0.)
        norm_polp80_costheta2_H_eptest.append(0)
        for j in range(0,h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetNbinsX()+1):
            norm_polp80_costheta2_H_epgen[i]+=h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetBinContent(i,j)
            h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.SetBinContent(j,i,h_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetBinContent(i,j))
        print("norm costheta2",i,norm_polp80_costheta2_H_epgen[i])


    for i in range(0,h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetNbinsX()+1):  
        for j in range(0,h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetNbinsX()+1):
            if norm_polp80_costheta2_H_epgen[j]!=0:
                h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.SetBinContent(i,j,h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetBinContent(i,j)/norm_polp80_costheta2_H_epgen[j])
            norm_polp80_costheta2_H_eptest[j]+=h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetBinContent(i,j)

    print("norm test costheta2",norm_polp80_costheta2_H_eptest)

    #h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_

    h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetXaxis().SetTitle("cos #theta_{2}^{reco}")
    h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.GetYaxis().SetTitle("cos #theta_{2}^{parton}")

    canvas_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80 = setUpperCanvasWColz("normed_parton_recojet_costheta2_rj1_ep_E_totCOM_vs_costheta2_H_ep_approx_HZ_COM__polp80");
    canvas_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80.cd();
    #gPad.SetRightMargin(0.18)
    h_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80_.Draw("colz")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi_normed,y2_noLumi,label3_noLumi);
    canvas_normed_parton_costheta2_H_ep_approx_HZ_COM_vs_recojet_costheta2_rj1_ep_E_totCOM_polp80.Print("h_recojet_costheta2_rj1_ep_E_totCOM_vs_parton_costheta2_H_ep_approx_HZ_COM__polp80_normed.eps")






    h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_=file_polp80_HZ_BDT_Unfolding_.Get("0.3/h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");
    h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetXaxis().SetTitle("cos #theta_{1}^{parton}");
    h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetYaxis().SetTitle("cos #theta_{1}^{reco}");
    canvas_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80 = setUpperCanvas("parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80");
    canvas_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80.cd();
    #gPad.SetRightMargin(0.05)
    h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.Draw("col")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi,y2_noLumi,label3_noLumi);
    canvas_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80.Print("h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80.eps")


    #now the normed transposed transfer matrix
    h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_=file_polp80_HZ_BDT_Unfolding_.Get("0.375/h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com");


    norm_polp80_costheta1_Z_qpos_Zcomgen=[]
    norm_polp80_costheta1_Z_qpos_Zcomtest=[]
    for i in range(0,h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetNbinsX()+1):
        norm_polp80_costheta1_Z_qpos_Zcomgen.append(0.)
        norm_polp80_costheta1_Z_qpos_Zcomtest.append(0)
        for j in range(0,h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetNbinsX()+1):
            norm_polp80_costheta1_Z_qpos_Zcomgen[i]+=h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetBinContent(i,j)
            h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.SetBinContent(j,i,h_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetBinContent(i,j))
        print("norm costheta2",i,norm_polp80_costheta1_Z_qpos_Zcomgen[i])


    for i in range(0,h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetNbinsX()+1):  
        for j in range(0,h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetNbinsX()+1):
            if norm_polp80_costheta1_Z_qpos_Zcomgen[j]!=0:
                h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.SetBinContent(i,j,h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetBinContent(i,j)/norm_polp80_costheta1_Z_qpos_Zcomgen[j])
            norm_polp80_costheta1_Z_qpos_Zcomtest[j]+=h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetBinContent(i,j)

    print("norm test costheta2",norm_polp80_costheta1_Z_qpos_Zcomtest)

    #h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_

    h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetXaxis().SetTitle("cos #theta_{1}^{reco}")
    h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.GetYaxis().SetTitle("cos #theta_{1}^{parton}")

    canvas_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80 = setUpperCanvasWColz("normed_parton_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2comvs_costheta1_Z_qpos_Zcom__polp80");
    canvas_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80.cd();
    #gPad.SetRightMargin(0.18)
    h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.SetMaximum(0.40);
    h_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80_.Draw("colz")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi_normed,y2_noLumi,label3_noLumi);
    canvas_normed_parton_costheta1_Z_qpos_Zcom_vs_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_polp80.Print("h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2comvs_parton_costheta1_Z_qpos_Zcom__polp80_normed.eps")


    h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_=file_polp80_HZ_BDT_Unfolding_.Get("0.3/h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");
    h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetXaxis().SetTitle("#phi^{parton} [#circ]");
    h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetYaxis().SetTitle("#phi^{reco} [#circ]");
    canvas_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80 = setUpperCanvas("parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80");
    canvas_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80.cd();
    #gPad.SetRightMargin(0.05)
    h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.Draw("col")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi,y2_noLumi,label3_noLumi);
    canvas_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80.Print("h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80.eps")



    #now the normed transposed transfer matrix
    h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_=file_polp80_HZ_BDT_Unfolding_.Get("0.375/h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep");


    norm_polp80_costheta1_Z_qpos_Zcomgen=[]
    norm_polp80_costheta1_Z_qpos_Zcomtest=[]
    for i in range(0,h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetNbinsX()+1):
        norm_polp80_costheta1_Z_qpos_Zcomgen.append(0.)
        norm_polp80_costheta1_Z_qpos_Zcomtest.append(0)
        for j in range(0,h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetNbinsX()+1):
            norm_polp80_costheta1_Z_qpos_Zcomgen[i]+=h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetBinContent(i,j)
            h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.SetBinContent(j,i,h_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetBinContent(i,j))
        print("norm phi",i,norm_polp80_costheta1_Z_qpos_Zcomgen[i])


    for i in range(0,h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetNbinsX()+1):  
        for j in range(0,h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetNbinsX()+1):
            if norm_polp80_costheta1_Z_qpos_Zcomgen[j]!=0:
                h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.SetBinContent(i,j,h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetBinContent(i,j)/norm_polp80_costheta1_Z_qpos_Zcomgen[j])
            norm_polp80_costheta1_Z_qpos_Zcomtest[j]+=h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetBinContent(i,j)

    print("norm test phi",norm_polp80_costheta1_Z_qpos_Zcomtest)

    #h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_

    h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetXaxis().SetTitle("#phi^{reco} [#circ]")
    h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.GetYaxis().SetTitle("#phi^{parton} [#circ]")

    canvas_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80 = setUpperCanvasWColz("normed_parton_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_epvs_costheta1_Z_qpos_Zcom__polp80");
    canvas_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80.cd();
    #gPad.SetRightMargin(0.18)
    h_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80_.Draw("colz")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2_noLumi_normed,y2_noLumi,label3_noLumi);
    canvas_normed_parton_phi_plane_Z_qpos_vs_plane_H_ep_vs_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_polp80.Print("h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep_vs_parton_phi_plane_Z_qpos_vs_plane_H_ep__polp80_normed.eps")




    file_polm80_HZ_full_withParton_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_hzqq_AnalysisBaselineHistos_noCuts_EThetaVar_withSignalHistos.root")
    h_HZqq_parton_sqrtS_polm80_=file_polm80_HZ_full_withParton_.Get("h_HZqq_signal_1D_sqrtS_part");
    h_HZqq_parton_sqrtS_polm80_.GetXaxis().SetTitle("#sqrt{s}_{parton} [GeV]");
    h_HZqq_parton_sqrtS_polm80_.GetYaxis().SetTitle("Events");
    h_HZqq_parton_sqrtS_polm80_.SetLineColor(kBlack)
    h_HZqq_parton_sqrtS_polm80_.SetLineWidth(3)


    canvas_HZqq_parton_sqrtS_polm80 = setUpperCanvas("canvas_HZqq_parton_sqrtS_polm80")
    canvas_HZqq_parton_sqrtS_polm80.cd()
    h_HZqq_parton_sqrtS_polm80_.Rebin(2)
    h_HZqq_parton_sqrtS_polm80_.Draw("hist,e")
    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_HZqq_parton_sqrtS_polm80.Print("h_HZqq_AllEvents_parton_sqrtS_polm80.eps")


    h_HZqq_reco_j1_j2_EMiss_high_sqrtS_polm80_=file_polm80_HZ_full_withParton_.Get("h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss_high_sqrtS");
    h_HZqq_reco_j1_j2_EMiss_high_sqrtS_polm80_.GetXaxis().SetTitle("#sqrt{s}_{reco} [GeV]");
    h_HZqq_reco_j1_j2_EMiss_high_sqrtS_polm80_.GetYaxis().SetTitle("Events");
    h_HZqq_reco_j1_j2_EMiss_high_sqrtS_polm80_.SetLineColor(kRed)
    h_HZqq_reco_j1_j2_EMiss_high_sqrtS_polm80_.SetLineWidth(3)
    h_HZqq_reco_j1_j2_high_sqrtS_polm80_=file_polm80_HZ_full_withParton_.Get("h_HZqq_signal_1D_sqrtS_reco_j1_j2_high_sqrtS");
    h_HZqq_reco_j1_j2_high_sqrtS_polm80_.GetXaxis().SetTitle("#sqrt{s}_{reco} [GeV]");
    h_HZqq_reco_j1_j2_high_sqrtS_polm80_.SetLineColor(kRed)
    h_HZqq_reco_j1_j2_high_sqrtS_polm80_.SetLineWidth(3)
    h_HZqq_reco_j1_j2_high_sqrtS_polm80_.SetLineStyle(2)

    h_HZqq_recoTot_isoPh_EMiss_high_sqrtS_polm80_=file_polm80_HZ_full_withParton_.Get("h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss_high_sqrtS");
    h_HZqq_recoTot_isoPh_EMiss_high_sqrtS_polm80_.GetXaxis().SetTitle("#sqrt{s}_{reco} [GeV]");
    h_HZqq_recoTot_isoPh_EMiss_high_sqrtS_polm80_.SetLineColor(kBlue)
    h_HZqq_recoTot_isoPh_EMiss_high_sqrtS_polm80_.SetLineWidth(3)

    h_HZqq_recoTot_isoPh_high_sqrtS_polm80_=file_polm80_HZ_full_withParton_.Get("h_HZqq_signal_1D_sqrtS_recoTot_isoPh_high_sqrtS");
    h_HZqq_recoTot_isoPh_high_sqrtS_polm80_.GetXaxis().SetTitle("#sqrt{s}_{reco} [GeV]");
    h_HZqq_recoTot_isoPh_high_sqrtS_polm80_ .GetYaxis().SetTitle("Events");
    h_HZqq_recoTot_isoPh_high_sqrtS_polm80_.SetLineColor(kBlue)
    h_HZqq_recoTot_isoPh_high_sqrtS_polm80_.SetLineWidth(3)
    h_HZqq_recoTot_isoPh_high_sqrtS_polm80_.SetLineStyle(2)

    h_HZqq_reco_j1_j2_EMiss_high_sqrtS_polm80_.Rebin(2)
    h_HZqq_reco_j1_j2_high_sqrtS_polm80_.Rebin(2)
    
    h_HZqq_recoTot_isoPh_high_sqrtS_polm80_.Rebin(2)
    h_HZqq_recoTot_isoPh_EMiss_high_sqrtS_polm80_.Rebin(2)
    #175 for bbar only
    #h_HZqq_reco_j1_j2_EMiss_high_sqrtS_polm80_.SetMaximum(250)

    canvas_HZqq_reco_sqrtS_polm80 = setUpperCanvas("canvas_HZqq_reco_sqrtS_polm80")
    canvas_HZqq_reco_sqrtS_polm80.cd()

    h_HZqq_recoTot_isoPh_high_sqrtS_polm80_.SetMaximum(175)

    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_=TLegend(0.60,0.63,0.90,0.88);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetBorderSize(0);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetTextAlign(12);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetTextSize(0.050);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetTextFont(42);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetMargin(0.15);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetLineColor(1);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetLineStyle(1);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetLineWidth(1);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetFillColor(0);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetFillStyle(0);
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.SetHeader("#sqrt{s}_{part}>2500 GeV");
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.AddEntry(h_HZqq_recoTot_isoPh_high_sqrtS_polm80_.DrawCopy("hist"),"all PFOs, orig");
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.AddEntry(h_HZqq_recoTot_isoPh_EMiss_high_sqrtS_polm80_.DrawCopy("hist,same"),"all PFOs, corr");
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.AddEntry(h_HZqq_reco_j1_j2_high_sqrtS_polm80_.DrawCopy("hist,same"),"jets, orig");
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.AddEntry(h_HZqq_reco_j1_j2_EMiss_high_sqrtS_polm80_.DrawCopy("hist,same"),"jets, corr");
    leg_hzqq_SIG_sqrtS_reco_polm80_noMassCuts_.Draw();


    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_HZqq_reco_sqrtS_polm80.Print("h_HZqq_reco_sqrtS_polm80.eps")


    h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500_polm80_=file_polm80_HZ_full_withParton_.Get("h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500");
    h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500_polm80_.GetXaxis().SetTitle("(#sqrt{s}_{part}-#sqrt{s}_{reco})/#sqrt{s}_{part}");
    h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500_polm80_.GetYaxis().SetTitle("Events");
    h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500_polm80_.SetLineColor(kRed)
    h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500_polm80_.SetLineWidth(3)
    h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500_polm80_=file_polm80_HZ_full_withParton_.Get("h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500");
    h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500_polm80_.GetXaxis().SetTitle("(#sqrt{s}_{part}-#sqrt{s}_{reco})/#sqrt{s}_{part}");
    h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500_polm80_.SetLineColor(kRed)
    h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500_polm80_.SetLineWidth(3)
    h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500_polm80_.SetLineStyle(2)

    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_EMiss_high_sqrtS_polm80_=file_polm80_HZ_full_withParton_.Get("h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_EMiss_part2500");
    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_EMiss_high_sqrtS_polm80_.GetXaxis().SetTitle("(#sqrt{s}_{part}-#sqrt{s}_{reco})/#sqrt{s}_{part}");
    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_EMiss_high_sqrtS_polm80_.SetLineColor(kBlue)
    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_EMiss_high_sqrtS_polm80_.SetLineWidth(3)

    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_high_sqrtS_polm80_=file_polm80_HZ_full_withParton_.Get("h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_part2500");
    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_high_sqrtS_polm80_.GetXaxis().SetTitle("(#sqrt{s}_{part}-#sqrt{s}_{reco})/#sqrt{s}_{part}");
    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_high_sqrtS_polm80_.GetYaxis().SetTitle("Events");
    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_high_sqrtS_polm80_.SetLineColor(kBlue)
    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_high_sqrtS_polm80_.SetLineWidth(3)
    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_high_sqrtS_polm80_.SetLineStyle(2)

    #h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500_polm80_.Rebin(2)
    #h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500_polm80_.Rebin(2)
    
    #h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_high_sqrtS_polm80_.Rebin(2)
    #h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_EMiss_high_sqrtS_polm80_.Rebin(2)
    
    canvas_HZqq_delta_rel_sqrtS_polm80 = setUpperCanvas("canvas_HZqq_delta_rel_sqrtS_polm80")
    canvas_HZqq_delta_rel_sqrtS_polm80.cd()

    h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_high_sqrtS_polm80_.SetMaximum(200)
    #h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_high_sqrtS_polm80_.SetMaximum(h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500_polm80_.GetMaximum())

    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_=TLegend(0.60,0.63,0.90,0.88);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetBorderSize(0);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetTextAlign(12);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetTextSize(0.050);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetTextFont(42);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetMargin(0.15);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetLineColor(1);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetLineStyle(1);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetLineWidth(1);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetFillColor(0);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetFillStyle(0);
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.SetHeader("#sqrt{s}_{part}>2500 GeV");
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.AddEntry(h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_high_sqrtS_polm80_.DrawCopy("hist"),"all PFOs, orig");
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.AddEntry(h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_EMiss_high_sqrtS_polm80_.DrawCopy("hist,same"),"all PFOs, corr");
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.AddEntry(h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500_polm80_.DrawCopy("hist,same"),"jets, orig");
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.AddEntry(h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500_polm80_.DrawCopy("hist,same"),"jets, corr");
    leg_hzqq_SIG_rel_sqrtS_reco_polm80_noMassCuts_.Draw();


    l.DrawLatex(x,y,label);
    l.DrawLatex(x2,y2,label2);
    canvas_HZqq_delta_rel_sqrtS_polm80.Print("h_HZqq_delta_rel_sqrtS_polm80.eps")

    """   canvas_h_SIG_BG_jet1_tau21_polp80_massCuts_thstack.Print("h_jet1_tau21_polp80_ee_qqqqqq_qqqq_qq_and_hzqq.eps")    file_polp80_ee_qqqqqqqqqqqq_SignalHistos_partonInfo_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/ee_qqqqqqqqqqqqAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_hzqqqqqqqqqqqq_AllEvents_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_AnglesMETProj_AllEvents.root")
    h_hzqqqqqqqqqqqq_AllEvents_BDT_polp80_recojet_BTag1_rj2_pos_sj_jetChargeE_0_30_rj2com_fake=file_polp80_ee_qqqqqqqqqqqq_SignalHistos_AllEvents_.Get("0.25/h_parton_cosBTag1_Z_qpos_Zcom_miss_HbbSelection");
    h_hzqqqqqqqqqqqq_AllEvents_BDT_polp80_recojet_BTag1_rj2_pos_sj_jetChargeE_0_30_rj2com_miss=file_polp80_ee_qqqqqqqqqqqq_SignalHistos_AllEvents_.Get("0.25/h_recojet_cosBTag1_rj2_pos_sj_jetChargeE_0_30_rj2com_fake");
    h_hzqqqqqqqqqqqq_AllEvents_BDT_polp80_recojet_BTag1_rj2_pos_sj_jetChargeE_0_30_rj2com_miss.GetXaxis().SetTitle("#BTag_{1} [#circ]");
    h_ee_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq_BDT_polp80_recojet_BTag1_rj2_pos_sj_jetChargeE_0_30_rj2com.GetYaxis().SetTitle("Events");

    
    canvas_h_2D_jet_mass_ee_qqqqqqqqqqqq_AllEventspolp80_recojet_BTag1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack = setUpperCanvas("canvas_h_2D_jet_mass_ee_qqqqqqqqqqqq_AllEventspolp80_recojet_BTag1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack");
    canvas_h_2D_jet_mass_ee_qqqqqqqqqqqq_AllEventspolp80_recojet_BTag1_rj2_pos_sj_jetChargeE_0_30_rj2com_thstack.cd();
    """
    return None


process_files()
root.gApplication.Run()

#
