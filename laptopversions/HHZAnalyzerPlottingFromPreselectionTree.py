from ROOT import gROOT, TCanvas, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor,TMVA,TCut,TString,TDirectory,TDatabasePDG
import ROOT as root
from math import cos, sin, pi, degrees, radians, pow, sqrt,acos,log,tan,atan2
from array import array

#done using angles
def DeltaPhi( Phi1, Phi2 ):
    deltaphi=abs(Phi1-Phi2)
    if (deltaphi>pi):
        deltaphi=2*pi-deltaphi        
    return deltaphi

def DeltaPhiDir( Phi1, Phi2 ):
    deltaphi=Phi1-Phi2
    if(deltaphi>pi):
        deltaphi=deltaphi-2*pi   
    if(deltaphi<(-pi)):
        deltaphi=2*pi+deltaphi
    return deltaphi

def DeltaPhiAngles( Phi1, Phi2 ):
    deltaphi=abs(Phi1-Phi2)
    pi_angle = 180.
    if (deltaphi>pi_angle):
        deltaphi=2*pi_angle-deltaphi        
    return deltaphi

def DeltaPhiDirAngles( Phi1, Phi2 ):
    deltaphi=Phi1-Phi2    
    pi_angle=180.
    if(deltaphi>pi_angle):
        deltaphi=deltaphi-2*pi_angle   
    if(deltaphi<(-pi_angle)):
        deltaphi=2*pi_angle+deltaphi
    return deltaphi

class RecoJet(TLorentzVector):
    def __init__(self,x,y,z,t,tau21,C2_beta1,D2_beta1,BTag_rfj_BTagMax,subjetE_ratio,subjet_DeltaPhi):
        TLorentzVector.__init__(self,x,y,z,t)
        self.tau21 = tau21
        self.C2_beta1 = C2_beta1
        self.D2_beta1 = D2_beta1
        self.BTag_rfj_BTagMax = BTag_rfj_BTagMax
        self.subjetE_ratio = subjetE_ratio
        self.subjet_DeltaPhi = subjet_DeltaPhi

class RecoSubJet(TLorentzVector):
    def __init__(self,jetChargePt_kappa_0_25,jetChargePt_kappa_0_30,jetChargeE_kappa_0_25,jetChargeE_kappa_0_30,nTracks,chFrac):
        self.jetChargePt_kappa_0_25 = jetChargePt_kappa_0_25
        self.jetChargePt_kappa_0_30 = jetChargePt_kappa_0_30
        self.jetChargeE_kappa_0_25 = jetChargeE_kappa_0_25
        self.jetChargeE_kappa_0_30 = jetChargeE_kappa_0_30


def setUpperCanvas(canvas_name) :
    c1= TCanvas(canvas_name,canvas_name,10,50,600,500)
    c1.cd()
    gPad.SetTopMargin(0.06)
    return c1
 
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


    
def process_event(i_final_histo_directory_,i_input_file_name_,fileout):

    n_bins_mass_high=250
    lim_mass_low=0
    lim_mass_high=250

    lim_E_low=0
    lim_E1_high=1650
    n_bins_E1 = 55
    
    lim_E2_high=1400
    n_bins_E2 = 56
    
    lim_E3_high=750
    n_bins_E3 = 75
    
    lim_E4_high=600
    n_bins_E4 = 60

    n_bins_theta = 90
    
    lim_theta_low=-0.0
    lim_theta_high=180.0
    
    n_bins_BTag_maxsum_low=150
    lim_BTag_maxsum_low=0.000
    lim_BTag_maxsum_high=3.000



    #fileout = root.TFile(i_final_histo_name_,"RECREATE")
    tdirectory=fileout.GetDirectory("")
    fileout.mkdir(str(i_final_histo_directory_))
    fileout.cd(str(i_final_histo_directory_))

    h_comb_jet1_mass = TH1F( "h_comb_jet1_mass", "", n_bins_mass_high, lim_mass_low, lim_mass_high)
    h_comb_jet2_mass = TH1F( "h_comb_jet2_mass", "", n_bins_mass_high, lim_mass_low, lim_mass_high)
    h_comb_jet3_mass = TH1F( "h_comb_jet3_mass", "", n_bins_mass_high, lim_mass_low, lim_mass_high)
    h_jet1_theta = TH1F( "h_jet1_theta", "", n_bins_theta, lim_theta_low, lim_theta_high)
    h_jet2_theta = TH1F( "h_jet2_theta", "", n_bins_theta, lim_theta_low, lim_theta_high)

    h_jet1_E = TH1F( "h_jet1_E", "", n_bins_E1,lim_E_low,lim_E1_high)
    h_jet2_E = TH1F( "h_jet2_E", "", n_bins_E2,lim_E_low,lim_E2_high)
    h_jet3_E = TH1F( "h_jet3_E", "", n_bins_E3,lim_E_low,lim_E3_high)
    h_jet4_E = TH1F( "h_jet4_E", "", n_bins_E4,lim_E_low,lim_E4_high)
    h_BTag_sum_max3 = TH1F( "h_BTag_sum_max3", "", n_bins_BTag_maxsum_low, lim_BTag_maxsum_low, lim_BTag_maxsum_high)

    input_file_=root.TFile.Open(i_input_file_name_)
    tree = input_file_.Get("MVATrainingVariables")

    print 'start processing',i_input_file_name_
    for ientry in tree:
        weight = ientry.weight
        h_comb_jet1_mass.Fill(ientry.comb_jet1_mass,weight)
        h_comb_jet2_mass.Fill(ientry.comb_jet2_mass,weight)
        h_comb_jet3_mass.Fill(ientry.comb_jet3_mass,weight)
        h_jet1_theta.Fill(degrees(ientry.jet1_theta),weight)
        h_jet2_theta.Fill(degrees(ientry.jet1_theta),weight)

        h_jet1_E.Fill(ientry.jet1_E,weight)
        h_jet2_E.Fill(ientry.jet2_E,weight)
        h_jet3_E.Fill(ientry.jet3_E,weight)
        h_jet4_E.Fill(ientry.jet4_E,weight)
   
        h_BTag_sum_max3.Fill(ientry.BTag_sum_max3,weight)
                 

    fileout.Write()
    return None

def process_files():

    #for negative polarisation, NTrees 300 gives best results, for positive polarisation change to NTrees 250

    #GradShrink200
    #AdaBoostBeta020
    #noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTag_noLTagMins --> has also no noj3CMax
    #noEratios_noCosThetaHel
    i_final_histo_name_="testfile_preselectionhistos_polm80.root"
    fileout=root.TFile(i_final_histo_name_,"RECREATE")
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis/ntuple_HHZ_hhqq_14364_polm80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="hhqq"

    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis/ntuple_HHZ_hhqq_14364_polm80_3TeV_wO_CLIC_o3_v14_AllEvents.root"
    i_final_histo_directory_="hhqq_AllEvents"

    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis/ntuple_HHZ_hzqq_13391_polm80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="hzqq"

    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis/ntuple_HHZ_ee_qq_13399_to_13402_polm80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="ee_qq"
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis/ntuple_HHZ_ee_qqqq_13394_to_13397_polm80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="ee_qqqq"

    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis/ntuple_HHZ_ee_qqqqqq_polm80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="ee_qqqqqq"
    
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)

    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis/ntuple_HHZ_WWH_qqqqH_14734_polm80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="WWH_qqqqH"
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)

    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis/ntuple_HHZ_ZZH_qqqqH_14726_polm80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="ZZH_qqqqH"
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    fileout.Close()
    
    i_final_histo_name_="testfile_preselectionhistos_polp80.root"
    fileout=root.TFile(i_final_histo_name_,"RECREATE")
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis/ntuple_HHZ_hhqq_14365_polp80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="hhqq"

    process_event(i_final_histo_directory_,i_input_file_name_,fileout)

    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis/ntuple_HHZ_hhqq_14365_polp80_3TeV_wO_CLIC_o3_v14_AllEvents.root"
    i_final_histo_directory_="hhqq_AllEvents"
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis/ntuple_HHZ_hzqq_13392_polp80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="hzqq"
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis/ntuple_HHZ_ee_qq_13398_polp80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="ee_qq"
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis/ntuple_HHZ_ee_qqqq_13393_polp80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="ee_qqqq"
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis/ntuple_HHZ_ee_qqqqqq_polp80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="ee_qqqqqq"
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)

    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis/ntuple_HHZ_WWH_qqqqH_14735_polp80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="WWH_qqqqH"
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    
    i_input_file_name_ = "/Users/matthiasweber/rootfilesHHZ/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis/ntuple_HHZ_ZZH_qqqqH_14727_polp80_3TeV_wO_CLIC_o3_v14.root"
    i_final_histo_directory_="ZZH_qqqqH"
    
    process_event(i_final_histo_directory_,i_input_file_name_,fileout)
    fileout.Close()


    return None

process_files()



