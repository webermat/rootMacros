from ROOT import gROOT, TCanvas, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor,TMVA,TCut,TString,TDirectory
import ROOT as root
from math import cos, sin, pi, degrees, radians, pow, sqrt,acos
from array import array

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
  #c1.SetTopMargin(0.055);
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


CLICdpStyle()

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

#file_polm80_HZ_SignalHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVTrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20.root")  
file_polm80_HZ_SignalHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20.root") 
file_polm80_HZ_SignalHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AllEvents.root") 
file_polm80_ee_qq_mqq_1TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20.root") 
file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20.root") 
file_polm80_ee_qqqqqq_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20.root") 

x_polm80_BDTScore = array( 'f' )
y_polm80_significance = array( 'f' )
y_polm80_purity = array( 'f' )
y_polm80_efficiency = array( 'f' )
norm=file_polm80_HZ_SignalHistos_.Get("-0.2/h_jet1_mass").Integral()
#norm=h_mass_sig_hzqq_norm.Integral()
#norm = 10.
print 'norm',norm

for dir_ind in directory:
    x_polm80_BDTScore.append(float(dir_ind))
    h_mass_sig_hzqq=file_polm80_HZ_SignalHistos_.Get(dir_ind+"/h_jet1_mass")
    h_mass_sig_ee_qq_mqq_1TeV_BGHistos=file_polm80_ee_qq_mqq_1TeV_BGHistos_.Get(dir_ind+"/h_jet1_mass")
    h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos=file_polm80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get(dir_ind+"/h_jet1_mass")
    h_mass_sig_ee_qqqqqq_BGHistos=file_polm80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_jet1_mass") 
    #print dir_ind,"integral of signal mass",h_mass_sig_hzqq.Integral(),h_mass_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_sig_ee_qqqqqq_BGHistos.Integral()
    if h_mass_sig_hzqq.Integral()>0:
        print '#',dir_ind,"polm significance,purity ",h_mass_sig_hzqq.Integral()/sqrt(h_mass_sig_hzqq.Integral()+h_mass_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_sig_hzqq.Integral()/(h_mass_sig_hzqq.Integral()+h_mass_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_sig_ee_qqqqqq_BGHistos.Integral()),h_mass_sig_hzqq.Integral(),h_mass_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_sig_ee_qqqqqq_BGHistos.Integral()
        y_polm80_efficiency.append(100.*h_mass_sig_hzqq.Integral()/norm)
        y_polm80_purity.append(100.*h_mass_sig_hzqq.Integral()/(h_mass_sig_hzqq.Integral()+h_mass_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_sig_ee_qqqqqq_BGHistos.Integral()))
        y_polm80_significance.append(h_mass_sig_hzqq.Integral()/sqrt(h_mass_sig_hzqq.Integral()+h_mass_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_sig_ee_qqqqqq_BGHistos.Integral()))
    else:
        y_polm80_efficiency.append(0)
        y_polm80_purity.append(0)
        y_polm80_significance.append(0)

canvas_polm80_BDT_significance = setUpperCanvas("canvas_polm80_BDT_significance");
canvas_polm80_BDT_significance.cd()
graph_polm80_significance = TGraph( n_graphs, x_polm80_BDTScore, y_polm80_significance )
graph_polm80_significance.SetMarkerStyle(21)
graph_polm80_significance.SetMarkerColor(0)
graph_polm80_significance.GetXaxis().SetTitle('BDT Score')
graph_polm80_significance.GetYaxis().SetTitle('significance [#sigma]')
graph_polm80_significance.Draw('ACP')
canvas_polm80_BDT_significance.Update()
canvas_polm80_BDT_significance.GetFrame().SetFillColor( 21 )
canvas_polm80_BDT_significance.GetFrame().SetBorderSize( 12 )
canvas_polm80_BDT_significance.Modified()
canvas_polm80_BDT_significance.Update()

canvas_polm80_BDT_efficiency = setUpperCanvas("canvas_polm80_BDT_efficiency");
canvas_polm80_BDT_efficiency.cd()
graph_polm80_efficiency = TGraph( n_graphs, x_polm80_BDTScore, y_polm80_efficiency )
graph_polm80_efficiency.SetMarkerStyle(21)
graph_polm80_efficiency.SetMarkerColor(0)
graph_polm80_efficiency.GetXaxis().SetTitle('BDT Score')
graph_polm80_efficiency.GetYaxis().SetTitle('signal efficiency [%]')
graph_polm80_purity = TGraph( n_graphs, x_polm80_BDTScore, y_polm80_purity )
graph_polm80_purity.SetMarkerStyle(22)
graph_polm80_purity.SetMarkerColor(1)
graph_polm80_purity.GetXaxis().SetTitle('BDT Score')
graph_polm80_purity.GetYaxis().SetTitle('purity [%]')

#file_polp80_HZ_SignalHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVTrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20.root")  
file_polp80_HZ_SignalHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20.root")  
file_polp80_HZ_SignalHistosAllEvents_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AllEvents.root")  
file_polp80_ee_qq_mqq_1TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20.root") 
file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20.root") 
file_polp80_ee_qqqqqq_BGHistos_=root.TFile.Open("/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20.root") 
#process_event(final_histo_name_,input_file_,files_weights_)

for dir_ind in directory:
    h_mass_sig_hzqq=file_polp80_HZ_SignalHistos_.Get(dir_ind+"/h_jet1_mass")
    h_mass_sig_ee_qq_mqq_1TeV_BGHistos=file_polp80_ee_qq_mqq_1TeV_BGHistos_.Get(dir_ind+"/h_jet1_mass")
    h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos=file_polp80_ee_qqqq_mqqqq_2TeV_BGHistos_.Get(dir_ind+"/h_jet1_mass")
    h_mass_sig_ee_qqqqqq_BGHistos=file_polp80_ee_qqqqqq_BGHistos_.Get(dir_ind+"/h_jet1_mass")
    if h_mass_sig_hzqq.Integral()>0:
        #print dir_ind,"integral of signal mass",h_mass_sig_hzqq.Integral(),h_mass_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_sig_ee_qqqqqq_BGHistos.Integral()
        print '#',dir_ind,"polp significance,purity ",h_mass_sig_hzqq.Integral()/sqrt(h_mass_sig_hzqq.Integral()+h_mass_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_sig_ee_qqqqqq_BGHistos.Integral()),100.*h_mass_sig_hzqq.Integral()/(h_mass_sig_hzqq.Integral()+h_mass_sig_ee_qq_mqq_1TeV_BGHistos.Integral()+h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral()+h_mass_sig_ee_qqqqqq_BGHistos.Integral()),h_mass_sig_hzqq.Integral(),h_mass_sig_ee_qq_mqq_1TeV_BGHistos.Integral(),h_mass_sig_ee_qqqq_mqqqq_2TeV_BGHistos.Integral(),h_mass_sig_ee_qqqqqq_BGHistos.Integral()


