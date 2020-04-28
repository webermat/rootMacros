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


    
def process_event(i_final_histo_name_,i_input_file_name_,i_weight_file,i_isSignalData):

    reader = root.TMVA.Reader("V:Color:!Silent")


    #noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTag_noLTagMins --> also no CTagMax 3
    #noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins__

    v_comb_jet1_mass = array('f',[0])
    v_comb_jet2_mass = array('f',[0])
    v_comb_jet3_mass = array('f',[0])
    v_comb_jet1_dalpha = array('f',[0])
    v_comb_jet2_dalpha = array('f',[0])
    v_comb_jet3_dalpha = array('f',[0])
    v_comb_jet1_E1_over_Etot = array('f',[0])
    v_comb_jet2_E1_over_Etot = array('f',[0])
    v_comb_jet3_E1_over_Etot = array('f',[0])
    v_comb_jet1_BTagMax = array('f',[0])
    v_comb_jet2_BTagMax = array('f',[0])
    v_comb_jet3_BTagMax = array('f',[0])
    v_comb_jet1_cosThetaHel_absmin = array('f',[0])
    v_comb_jet2_cosThetaHel_absmin = array('f',[0])
    v_comb_jet3_cosThetaHel_absmin = array('f',[0])
    v_comb_jet1_E = array('f',[0])
    v_comb_jet2_E = array('f',[0])
    v_comb_jet3_E = array('f',[0])
    v_comb_jet1_theta = array('f',[0])
    v_comb_jet2_theta = array('f',[0])
    v_comb_jet3_theta = array('f',[0])
    v_BTag_sum_max2 = array('f',[0])
    v_BTag_sum_max3 = array('f',[0])
    v_BTag_sum_max4 = array('f',[0])
    v_BTag_sum_all = array('f',[0])
    #new variables for tight mass cut samples
    v_CTag_sum_max2 = array('f',[0])    
    v_CTag_sum_max3 = array('f',[0])
    v_LTag_sum_all = array('f',[0])
    v_comb_jet1_CTagMax = array('f',[0])
    v_comb_jet2_CTagMax = array('f',[0])
    v_comb_jet3_CTagMax = array('f',[0])
    v_comb_jet1_LTagMax = array('f',[0])
    v_comb_jet2_LTagMax = array('f',[0])
    v_comb_jet3_LTagMax = array('f',[0])
    v_comb_jet1_LTagMin = array('f',[0])
    v_comb_jet2_LTagMin = array('f',[0])
    v_comb_jet3_LTagMin = array('f',[0])
    v_jet1_theta = array('f',[0])
    v_jet2_theta = array('f',[0])
    v_jet3_theta = array('f',[0])
    v_jet4_theta = array('f',[0])
    v_jet5_theta = array('f',[0])
    v_jet6_theta = array('f',[0])
    v_jet1_E = array('f',[0])
    v_jet2_E = array('f',[0])
    v_jet3_E = array('f',[0])
    v_jet4_E = array('f',[0])
    v_jet5_E = array('f',[0])
    v_jet6_E = array('f',[0])
    v_jet1_BTag = array('f',[0])
    v_jet2_BTag = array('f',[0])
    v_jet3_BTag = array('f',[0])
    v_jet4_BTag = array('f',[0])
    v_jet5_BTag = array('f',[0])
    v_jet6_BTag = array('f',[0])
    v_y12 = array('f',[0])
    v_y23 = array('f',[0])
    v_y34 = array('f',[0])
    v_y45 = array('f',[0])
    v_y56 = array('f',[0])
    

    reader.AddVariable("comb_jet1_mass", v_comb_jet1_mass)
    reader.AddVariable("comb_jet2_mass", v_comb_jet2_mass)
    reader.AddVariable("comb_jet3_mass", v_comb_jet3_mass)
    reader.AddVariable("comb_jet1_dalpha", v_comb_jet1_dalpha)
    reader.AddVariable("comb_jet2_dalpha", v_comb_jet2_dalpha)
    reader.AddSpectator("comb_jet3_dalpha", v_comb_jet3_dalpha)
    reader.AddSpectator("comb_jet1_E1_over_Etot", v_comb_jet1_E1_over_Etot)
    reader.AddSpectator("comb_jet2_E1_over_Etot", v_comb_jet2_E1_over_Etot)
    reader.AddSpectator("comb_jet3_E1_over_Etot", v_comb_jet3_E1_over_Etot)
    reader.AddVariable("comb_jet1_BTagMax", v_comb_jet1_BTagMax)
    reader.AddVariable("comb_jet2_BTagMax", v_comb_jet2_BTagMax)
    reader.AddVariable("comb_jet3_BTagMax", v_comb_jet3_BTagMax)
    reader.AddSpectator("comb_jet1_cosThetaHel_absmin", v_comb_jet1_cosThetaHel_absmin)
    reader.AddSpectator("comb_jet2_cosThetaHel_absmin", v_comb_jet2_cosThetaHel_absmin)
    reader.AddSpectator("comb_jet3_cosThetaHel_absmin", v_comb_jet3_cosThetaHel_absmin)
    reader.AddSpectator("comb_jet1_E", v_comb_jet1_E)
    reader.AddSpectator("comb_jet2_E", v_comb_jet2_E)
    reader.AddSpectator("comb_jet3_E", v_comb_jet3_E)
    reader.AddSpectator("comb_jet1_theta", v_comb_jet1_theta)
    reader.AddSpectator("comb_jet2_theta", v_comb_jet2_theta)
    reader.AddSpectator("comb_jet3_theta", v_comb_jet3_theta)
    reader.AddVariable("BTag_sum_max2", v_BTag_sum_max2)
    reader.AddVariable("BTag_sum_max3", v_BTag_sum_max3)
    reader.AddSpectator("BTag_sum_max4", v_BTag_sum_max4)
    reader.AddSpectator("BTag_sum_all", v_BTag_sum_all)
    reader.AddSpectator("CTag_sum_max2", v_CTag_sum_max2)
    reader.AddSpectator("CTag_sum_max3", v_CTag_sum_max3)
    reader.AddVariable("LTag_sum_all", v_LTag_sum_all)
    reader.AddVariable("comb_jet1_CTagMax", v_comb_jet1_CTagMax)
    reader.AddVariable("comb_jet2_CTagMax", v_comb_jet2_CTagMax)
    reader.AddSpectator("comb_jet3_CTagMax", v_comb_jet3_CTagMax)
    reader.AddVariable("comb_jet1_LTagMax", v_comb_jet1_LTagMax)
    reader.AddVariable("comb_jet2_LTagMax", v_comb_jet2_LTagMax)
    reader.AddVariable("comb_jet3_LTagMax", v_comb_jet3_LTagMax)
    reader.AddSpectator("comb_jet1_LTagMin", v_comb_jet1_LTagMin)
    reader.AddSpectator("comb_jet2_LTagMin", v_comb_jet2_LTagMin)
    reader.AddSpectator("comb_jet3_LTagMin", v_comb_jet3_LTagMin)
    reader.AddVariable("jet1_theta", v_jet1_theta)
    reader.AddVariable("jet2_theta", v_jet2_theta)
    reader.AddVariable("jet3_theta", v_jet3_theta)
    reader.AddVariable("jet4_theta", v_jet4_theta)
    reader.AddVariable("jet5_theta", v_jet5_theta)
    reader.AddVariable("jet6_theta", v_jet6_theta)
    reader.AddVariable("jet1_E", v_jet1_E)
    reader.AddVariable("jet2_E", v_jet2_E)
    reader.AddVariable("jet3_E", v_jet3_E)
    reader.AddVariable("jet4_E", v_jet4_E)
    reader.AddVariable("jet5_E", v_jet5_E)
    reader.AddVariable("jet6_E", v_jet6_E)
    reader.AddSpectator("jet1_BTag", v_jet1_BTag)
    reader.AddSpectator("jet2_BTag", v_jet2_BTag)
    reader.AddSpectator("jet3_BTag", v_jet3_BTag)
    reader.AddSpectator("jet4_BTag", v_jet4_BTag)
    reader.AddSpectator("jet5_BTag", v_jet5_BTag)
    reader.AddSpectator("jet6_BTag", v_jet6_BTag)
    reader.AddVariable("y12", v_y12)
    reader.AddVariable("y23", v_y23)
    reader.AddVariable("y34", v_y34)
    reader.AddVariable("y45", v_y45)
    reader.AddVariable("y56", v_y56)



    v_CTag_sum_max4 = array('f',[0])
    v_CTag_sum_all = array('f',[0])
    input_file_=root.TFile.Open(i_input_file_name_)
    tree = input_file_.Get("MVATrainingVariables")
    tree.SetBranchAddress('comb_jet1_mass',v_comb_jet1_mass)
    tree.SetBranchAddress('comb_jet2_mass',v_comb_jet2_mass)
    tree.SetBranchAddress('comb_jet3_mass',v_comb_jet3_mass)
    tree.SetBranchAddress('comb_jet1_dalpha',v_comb_jet1_dalpha)
    tree.SetBranchAddress('comb_jet2_dalpha',v_comb_jet2_dalpha)
    tree.SetBranchAddress('comb_jet3_dalpha',v_comb_jet3_dalpha)
    tree.SetBranchAddress('comb_jet1_E1_over_Etot',v_comb_jet1_E1_over_Etot)
    tree.SetBranchAddress('comb_jet2_E1_over_Etot',v_comb_jet2_E1_over_Etot)
    tree.SetBranchAddress('comb_jet3_E1_over_Etot',v_comb_jet3_E1_over_Etot)
    tree.SetBranchAddress('comb_jet1_BTagMax',v_comb_jet1_BTagMax)
    tree.SetBranchAddress('comb_jet2_BTagMax',v_comb_jet2_BTagMax)
    tree.SetBranchAddress('comb_jet3_BTagMax',v_comb_jet3_BTagMax)
    tree.SetBranchAddress('comb_jet1_cosThetaHel_absmin',v_comb_jet1_cosThetaHel_absmin)
    tree.SetBranchAddress('comb_jet2_cosThetaHel_absmin',v_comb_jet2_cosThetaHel_absmin)
    tree.SetBranchAddress('comb_jet3_cosThetaHel_absmin',v_comb_jet3_cosThetaHel_absmin)
    tree.SetBranchAddress('comb_jet1_E',v_comb_jet1_E)
    tree.SetBranchAddress('comb_jet2_E',v_comb_jet2_E)
    tree.SetBranchAddress('comb_jet3_E',v_comb_jet3_E)
    tree.SetBranchAddress('comb_jet1_theta',v_comb_jet1_theta)
    tree.SetBranchAddress('comb_jet2_theta',v_comb_jet2_theta)
    tree.SetBranchAddress('comb_jet3_theta',v_comb_jet3_theta)
    tree.SetBranchAddress('BTag_sum_max2',v_BTag_sum_max2)
    tree.SetBranchAddress('BTag_sum_max3',v_BTag_sum_max3)
    tree.SetBranchAddress('BTag_sum_max4',v_BTag_sum_max4)
    tree.SetBranchAddress('BTag_sum_all',v_BTag_sum_all)
    #extra variables only there for tight mass selections
    tree.SetBranchAddress('comb_jet1_CTagMax',v_comb_jet1_CTagMax)
    tree.SetBranchAddress('comb_jet2_CTagMax',v_comb_jet2_CTagMax)
    tree.SetBranchAddress('comb_jet3_CTagMax',v_comb_jet3_CTagMax)
    tree.SetBranchAddress('comb_jet1_LTagMax',v_comb_jet1_LTagMax)
    tree.SetBranchAddress('comb_jet2_LTagMax',v_comb_jet2_LTagMax)
    tree.SetBranchAddress('comb_jet3_LTagMax',v_comb_jet3_LTagMax)
    tree.SetBranchAddress('comb_jet1_LTagMin',v_comb_jet1_LTagMin)
    tree.SetBranchAddress('comb_jet2_LTagMin',v_comb_jet2_LTagMin)
    tree.SetBranchAddress('comb_jet3_LTagMin',v_comb_jet3_LTagMin)
    tree.SetBranchAddress('CTag_sum_max2',v_CTag_sum_max2)
    tree.SetBranchAddress('CTag_sum_max3',v_CTag_sum_max3)
    tree.SetBranchAddress('LTag_sum_all',v_LTag_sum_all)
    #extra variables done
    tree.SetBranchAddress('jet1_theta',v_jet1_theta)
    tree.SetBranchAddress('jet2_theta',v_jet2_theta)
    tree.SetBranchAddress('jet3_theta',v_jet3_theta)
    tree.SetBranchAddress('jet4_theta',v_jet4_theta)
    tree.SetBranchAddress('jet5_theta',v_jet5_theta)
    tree.SetBranchAddress('jet6_theta',v_jet6_theta)
    tree.SetBranchAddress('jet1_E',v_jet1_E)
    tree.SetBranchAddress('jet2_E',v_jet2_E)
    tree.SetBranchAddress('jet3_E',v_jet3_E)
    tree.SetBranchAddress('jet4_E',v_jet4_E)
    tree.SetBranchAddress('jet5_E',v_jet5_E)
    tree.SetBranchAddress('jet6_E',v_jet6_E)
    tree.SetBranchAddress('jet1_BTag',v_jet1_BTag)
    tree.SetBranchAddress('jet2_BTag',v_jet2_BTag)
    tree.SetBranchAddress('jet3_BTag',v_jet3_BTag)
    tree.SetBranchAddress('jet4_BTag',v_jet4_BTag)
    tree.SetBranchAddress('jet5_BTag',v_jet5_BTag)
    tree.SetBranchAddress('jet6_BTag',v_jet6_BTag)
    tree.SetBranchAddress('y12',v_y12)
    tree.SetBranchAddress('y23',v_y23)
    tree.SetBranchAddress('y34',v_y34)
    tree.SetBranchAddress('y45',v_y45)
    tree.SetBranchAddress('y56',v_y56)


    v_weight = array('f',[0])
    tree.SetBranchAddress('weight',v_weight)
    v_sqrtS_parton = array('f',[0])
    if i_isSignalData:
        tree.SetBranchAddress('sqrtS_parton',v_sqrtS_parton)
    v_sqrtS_gen = array('f',[0])
    tree.SetBranchAddress('sqrtS_gen',v_sqrtS_gen)
    v_sqrtS_genjets = array('f',[0])
    tree.SetBranchAddress('sqrtS_genjets',v_sqrtS_genjets)
    v_sqrtS = array('f',[0])
    tree.SetBranchAddress('sqrtS',v_sqrtS)
    v_sqrtS_jets = array('f',[0])
    tree.SetBranchAddress('sqrtS_jets',v_sqrtS_jets)
    v_MET = array('f',[0])
    tree.SetBranchAddress('MET',v_MET)
    v_MET_gen = array('f',[0])
    tree.SetBranchAddress('MET_gen',v_MET_gen)
    v_MHT = array('f',[0])
    tree.SetBranchAddress('MHT',v_MHT)
    v_MHT_gen = array('f',[0])
    tree.SetBranchAddress('MHT_gen',v_MHT_gen)


    reader.BookMVA('BDT',TString(i_weight_file))
    
        
    n_bins_high=100
    lim_BDT_low=-1.0
    lim_BDT_high=1.0
    lim_mass_low=0
    lim_mass_high=250
    lim_theta_low=-0.5
    lim_theta_high=180.5
    n_bins_very_high=500
    lim_E_low=0
    lim_E_high=2000
    lim_MET_low=0
    lim_MET_high=600

    lim_theta_low=-0.0
    lim_theta_high=180.0

    lim_dtheta_low=-180.5
    lim_dtheta_high=180.5

    lim_phi_low=-0.0
    lim_phi_high=360.0

    lim_BTag_low=-0.005
    lim_BTag_high=1.005

    lim_dAlpha_sj12_low=0
    lim_dAlpha_sj12_high=25

    lim_Efrac_low=0.0
    lim_Efrac_high=0.5

    lim_sqrtS_low=0.0
    lim_sqrtS_high=3500

    lim_jetE_ratio_low=0.5
    lim_jetE_ratio_high=1.0

    lim_y_ij_low=0.0
    lim_y12_high=0.25
    lim_y23_high=0.1
    lim_y34_high=0.04
    lim_y45_high=0.02
    lim_y56_high=0.015

    BDT_cuts=[]
    #print 'do i get here maybe 1'
    BDT_cuts.append(-0.200)
    #print 'do i get here maybe 2'
    #BDT_cuts.append(-0.150)
    BDT_cuts.append(-0.100)
    #BDT_cuts.append(-0.050)
    BDT_cuts.append(0.000)
    #BDT_cuts.append(0.050)
    BDT_cuts.append(0.100)
    BDT_cuts.append(0.150)
    #BDT_cuts.append(0.175)
    BDT_cuts.append(0.200)
    #BDT_cuts.append(0.225)
    BDT_cuts.append(0.250)
    #BDT_cuts.append(0.275)
    BDT_cuts.append(0.300)
    BDT_cuts.append(0.325)
    BDT_cuts.append(0.350)
    BDT_cuts.append(0.375)
    BDT_cuts.append(0.385)
    BDT_cuts.append(0.400)
    BDT_cuts.append(0.410)
    BDT_cuts.append(0.425)
    BDT_cuts.append(0.435)
    BDT_cuts.append(0.440)
    BDT_cuts.append(0.450)
    BDT_cuts.append(0.460)
    BDT_cuts.append(0.470)
    BDT_cuts.append(0.480)
    BDT_cuts.append(0.490)
    BDT_cuts.append(0.500)
    BDT_cuts.append(0.510)
    BDT_cuts.append(0.525)
    BDT_cuts.append(0.550)
    BDT_cuts.append(0.575)
    BDT_cuts.append(0.600)
    BDT_cuts.append(0.625)
    BDT_cuts.append(0.630)
    BDT_cuts.append(0.640)
    BDT_cuts.append(0.650)
    BDT_cuts.append(0.675)
    BDT_cuts.append(0.700)
    BDT_cuts.append(0.725)
    BDT_cuts.append(0.735)
    BDT_cuts.append(0.740)
    BDT_cuts.append(0.745)
    BDT_cuts.append(0.750)
    BDT_cuts.append(0.760)
    BDT_cuts.append(0.770)
    BDT_cuts.append(0.775)
    BDT_cuts.append(0.800)
    BDT_cuts.append(0.810)
    BDT_cuts.append(0.825)
    BDT_cuts.append(0.830)
    BDT_cuts.append(0.850)
    BDT_cuts.append(0.875)
    BDT_cuts.append(0.900)
    BDT_cuts.append(0.925)
    BDT_cuts.append(0.940)
    BDT_cuts.append(0.950)
    BDT_cuts.append(0.960)
    BDT_cuts.append(0.975)
    BDT_cuts.append(1.000)
    fileout = root.TFile(i_final_histo_name_,"RECREATE")
    for bdt_value in BDT_cuts:
        tdirectory=fileout.GetDirectory(i_final_histo_name_)
        fileout.mkdir(str(bdt_value))
        fileout.cd(str(bdt_value))

        h_BDT_output = TH1F( "h_BDT_output", "", n_bins_high, lim_BDT_low, lim_BDT_high)
        h_BDT_output_Eff = TH1F( "h_BDT_output_Eff", "", n_bins_high, lim_BDT_low, lim_BDT_high)
        h_comb_jet1_mass = TH1F( "h_comb_jet1_mass", "", n_bins_high, lim_mass_low, lim_mass_high)
        h_comb_jet2_mass = TH1F( "h_comb_jet2_mass", "", n_bins_high, lim_mass_low, lim_mass_high)
        h_comb_jet3_mass = TH1F( "h_comb_jet3_mass", "", n_bins_high, lim_mass_low, lim_mass_high)
        h_comb_jet1_theta = TH1F( "h_comb_jet1_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_comb_jet2_theta = TH1F( "h_comb_jet2_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_comb_jet3_theta = TH1F( "h_comb_jet3_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_comb_jet1_E = TH1F( "h_comb_jet1_E", "", n_bins_very_high, lim_E_low, lim_E_high)
        h_comb_jet2_E = TH1F( "h_comb_jet2_E", "", n_bins_very_high, lim_E_low, lim_E_high)
        h_comb_jet3_E = TH1F( "h_comb_jet3_E", "", n_bins_very_high, lim_E_low, lim_E_high)
        h_comb_jet1_BTagMax = TH1F( "h_comb_jet1_BTagMax", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet2_BTagMax = TH1F( "h_comb_jet2_BTagMax", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet3_BTagMax = TH1F( "h_comb_jet3_BTagMax", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet1_cosThetaHel_absmin = TH1F( "h_comb_jet1_cosThetaHel_absmin", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet2_cosThetaHel_absmin = TH1F( "h_comb_jet2_cosThetaHel_absmin", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet3_cosThetaHel_absmin = TH1F( "h_comb_jet3_cosThetaHel_absmin", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet1_dalpha = TH1F( "h_comb_jet1_dalpha", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_comb_jet2_dalpha = TH1F( "h_comb_jet2_dalpha", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_comb_jet3_dalpha = TH1F( "h_comb_jet3_dalpha", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_comb_jet1_E1_over_Etot = TH1F( "h_comb_jet1_E1_over_Etot", "", n_bins_high, lim_jetE_ratio_low, lim_jetE_ratio_high)
        h_comb_jet2_E1_over_Etot = TH1F( "h_comb_jet2_E1_over_Etot", "", n_bins_high, lim_jetE_ratio_low, lim_jetE_ratio_high)
        h_comb_jet3_E1_over_Etot = TH1F( "h_comb_jet3_E1_over_Etot", "", n_bins_high, lim_jetE_ratio_low, lim_jetE_ratio_high)
        h_BTag_sum_max2 = TH1F( "h_BTag_sum_max2", "", n_bins_high, lim_BTag_low, 2.*lim_BTag_high)
        h_BTag_sum_max3 = TH1F( "h_BTag_sum_max3", "", n_bins_high, lim_BTag_low, 3.*lim_BTag_high)
        h_BTag_sum_max4 = TH1F( "h_BTag_sum_max4", "", n_bins_high, lim_BTag_low, 4.*lim_BTag_high)
        h_BTag_sum_all = TH1F( "h_BTag_sum_all", "", n_bins_high, lim_BTag_low, 6.*lim_BTag_high)
        #for tight mass cuts available
        h_comb_jet1_CTagMax = TH1F( "h_comb_jet1_CTagMax", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet2_CTagMax = TH1F( "h_comb_jet2_CTagMax", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet3_CTagMax = TH1F( "h_comb_jet3_CTagMax", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet1_LTagMax = TH1F( "h_comb_jet1_LTagMax", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet2_LTagMax = TH1F( "h_comb_jet2_LTagMax", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet3_LTagMax = TH1F( "h_comb_jet3_LTagMax", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet1_LTagMin = TH1F( "h_comb_jet1_LTagMin", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet2_LTagMin = TH1F( "h_comb_jet2_LTagMin", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_comb_jet3_LTagMin = TH1F( "h_comb_jet3_LTagMin", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_CTag_sum_max2 = TH1F( "h_CTag_sum_max2", "", n_bins_high, lim_BTag_low, 2.*lim_BTag_high)
        h_CTag_sum_max3 = TH1F( "h_CTag_sum_max3", "", n_bins_high, lim_BTag_low, 3.*lim_BTag_high)
        h_CTag_sum_max4 = TH1F( "h_CTag_sum_max4", "", n_bins_high, lim_BTag_low, 4.*lim_BTag_high)
        h_CTag_sum_all = TH1F( "h_CTag_sum_all", "", n_bins_high, lim_BTag_low, 6.*lim_BTag_high)
        h_LTag_sum_all = TH1F( "h_LTag_sum_all", "", n_bins_high, lim_BTag_low, 6.*lim_BTag_high)
        #extra histos are done here
        h_jet1_theta = TH1F( "h_jet1_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_jet2_theta = TH1F( "h_jet2_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_jet3_theta = TH1F( "h_jet3_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_jet4_theta = TH1F( "h_jet4_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_jet5_theta = TH1F( "h_jet5_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_jet6_theta = TH1F( "h_jet6_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_jet1_E = TH1F( "h_jet1_E", "", n_bins_very_high, lim_E_low, lim_E_high)
        h_jet2_E = TH1F( "h_jet2_E", "", n_bins_very_high, lim_E_low, lim_E_high)
        h_jet3_E = TH1F( "h_jet3_E", "", n_bins_very_high, lim_E_low, lim_E_high)
        h_jet4_E = TH1F( "h_jet4_E", "", n_bins_very_high, lim_E_low, lim_E_high)
        h_jet5_E = TH1F( "h_jet5_E", "", n_bins_very_high, lim_E_low, lim_E_high)
        h_jet6_E = TH1F( "h_jet6_E", "", n_bins_very_high, lim_E_low, lim_E_high)
        h_jet1_BTag = TH1F( "h_jet1_BTag", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_jet2_BTag = TH1F( "h_jet2_BTag", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_jet3_BTag = TH1F( "h_jet3_BTag", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_jet4_BTag = TH1F( "h_jet4_BTag", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_jet5_BTag = TH1F( "h_jet5_BTag", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_jet6_BTag = TH1F( "h_jet6_BTag", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_y12 = TH1F( "h_y12", "", n_bins_high, lim_y_ij_low, lim_y12_high)
        h_y23 = TH1F( "h_y23", "", n_bins_high, lim_y_ij_low, lim_y23_high)
        h_y34 = TH1F( "h_y34", "", n_bins_high, lim_y_ij_low, lim_y34_high)
        h_y45 = TH1F( "h_y45", "", n_bins_high, lim_y_ij_low, lim_y45_high)
        h_y56 = TH1F( "h_y56", "", n_bins_high, lim_y_ij_low, lim_y56_high)
        h_sqrtS = TH1F( "h_sqrtS", "", n_bins_high, lim_sqrtS_low, lim_sqrtS_high)
        h_sqrtS_jets = TH1F( "h_sqrtS_jets", "", n_bins_high, lim_sqrtS_low, lim_sqrtS_high)
        h_sqrtS_gen = TH1F( "h_sqrtS_gen", "", n_bins_high, lim_sqrtS_low, lim_sqrtS_high)
        h_sqrtS_genjets = TH1F( "h_sqrtS_genjet", "", n_bins_high, lim_sqrtS_low, lim_sqrtS_high)
        h_sqrtS_parton = TH1F( "h_sqrtS_parton", "", n_bins_high, lim_sqrtS_low, lim_sqrtS_high)
        h_MET = TH1F( "h_MET", "", n_bins_high, lim_MET_low, lim_MET_high)
        h_MET_gen = TH1F( "h_MET_gen", "", n_bins_high, lim_MET_low, lim_MET_high)
        h_MHT = TH1F( "h_MHT", "", n_bins_high, lim_MET_low, lim_MET_high)
        h_MHT_gen = TH1F( "h_MHT_gen", "", n_bins_high, lim_MET_low, lim_MET_high)

        num_entry=-1
        #BDTCut=0.50
        total_weight=0
        for ientry in range(tree.GetEntries()):
            num_entry+=1;
            #print "sig BG in entry ",num_entry,bdt_value             
            tree.GetEntry(ientry)
            mvaValue=reader.EvaluateMVA("BDT")
            if num_entry%(int(tree.GetEntries()/5.)) == 0:
                print "sig BG in entry/bdt-folder/mva/tot_weight ",num_entry,bdt_value,mvaValue,total_weight
            #we only get here, if we pass reco style cuts
            h_BDT_output.Fill(mvaValue,v_weight[0])
            if(mvaValue>bdt_value):
                total_weight+=v_weight[0]
                #print 'parton stuff',v_sqrtS_parton[0]
                h_BDT_output_Eff.Fill(mvaValue,v_weight[0])
                if i_isSignalData:
                    h_sqrtS_parton.Fill(v_sqrtS_parton[0],v_weight[0])
                h_sqrtS_gen.Fill(v_sqrtS_gen[0],v_weight[0])
                h_sqrtS_genjets.Fill(v_sqrtS_genjets[0],v_weight[0])
                h_sqrtS.Fill(v_sqrtS[0],v_weight[0])
                h_sqrtS_jets.Fill(v_sqrtS_jets[0],v_weight[0])
                h_MET.Fill(v_MET[0],v_weight[0])
                h_MET_gen.Fill(v_MET_gen[0],v_weight[0])
                h_MHT.Fill(v_MHT[0],v_weight[0])
                h_MHT_gen.Fill(v_MHT_gen[0],v_weight[0])
                h_comb_jet1_mass.Fill(v_comb_jet1_mass[0],v_weight[0])
                h_comb_jet2_mass.Fill(v_comb_jet2_mass[0],v_weight[0])
                h_comb_jet3_mass.Fill(v_comb_jet3_mass[0],v_weight[0])
                h_comb_jet1_theta.Fill(degrees(v_comb_jet1_theta[0]),v_weight[0])
                h_comb_jet2_theta.Fill(degrees(v_comb_jet2_theta[0]),v_weight[0])
                h_comb_jet3_theta.Fill(degrees(v_comb_jet3_theta[0]),v_weight[0])
                h_comb_jet1_E.Fill(v_comb_jet1_E[0],v_weight[0])
                h_comb_jet2_E.Fill(v_comb_jet2_E[0],v_weight[0])
                h_comb_jet3_E.Fill(v_comb_jet3_E[0],v_weight[0])
                h_comb_jet1_BTagMax.Fill(v_comb_jet1_BTagMax[0],v_weight[0])
                h_comb_jet2_BTagMax.Fill(v_comb_jet2_BTagMax[0],v_weight[0])
                h_comb_jet3_BTagMax.Fill(v_comb_jet3_BTagMax[0],v_weight[0])
                h_comb_jet1_cosThetaHel_absmin.Fill(v_comb_jet1_cosThetaHel_absmin[0],v_weight[0])
                h_comb_jet2_cosThetaHel_absmin.Fill(v_comb_jet1_cosThetaHel_absmin[0],v_weight[0])
                h_comb_jet3_cosThetaHel_absmin.Fill(v_comb_jet1_cosThetaHel_absmin[0],v_weight[0])
                h_comb_jet1_dalpha.Fill(degrees(v_comb_jet1_dalpha[0]),v_weight[0])
                h_comb_jet2_dalpha.Fill(degrees(v_comb_jet2_dalpha[0]),v_weight[0])
                h_comb_jet3_dalpha.Fill(degrees(v_comb_jet3_dalpha[0]),v_weight[0])
                h_comb_jet1_E1_over_Etot.Fill(v_comb_jet1_E1_over_Etot[0],v_weight[0])
                h_comb_jet2_E1_over_Etot.Fill(v_comb_jet2_E1_over_Etot[0],v_weight[0])
                h_comb_jet3_E1_over_Etot.Fill(v_comb_jet3_E1_over_Etot[0],v_weight[0])
                h_BTag_sum_max2.Fill(v_BTag_sum_max2[0],v_weight[0])
                h_BTag_sum_max3.Fill(v_BTag_sum_max3[0],v_weight[0])
                h_BTag_sum_max4.Fill(v_BTag_sum_max4[0],v_weight[0])
                h_BTag_sum_all.Fill(v_BTag_sum_all[0],v_weight[0])
                #extra plots for tight mass selection
                h_comb_jet1_CTagMax.Fill(v_comb_jet1_CTagMax[0],v_weight[0])
                h_comb_jet2_CTagMax.Fill(v_comb_jet2_CTagMax[0],v_weight[0])
                h_comb_jet3_CTagMax.Fill(v_comb_jet3_CTagMax[0],v_weight[0])
                h_comb_jet1_LTagMax.Fill(v_comb_jet1_LTagMax[0],v_weight[0])
                h_comb_jet2_LTagMax.Fill(v_comb_jet2_LTagMax[0],v_weight[0])
                h_comb_jet3_LTagMax.Fill(v_comb_jet3_LTagMax[0],v_weight[0])
                h_comb_jet1_LTagMin.Fill(v_comb_jet1_LTagMin[0],v_weight[0])
                h_comb_jet2_LTagMin.Fill(v_comb_jet2_LTagMin[0],v_weight[0])
                h_comb_jet3_LTagMin.Fill(v_comb_jet3_LTagMin[0],v_weight[0])
                h_CTag_sum_max2.Fill(v_CTag_sum_max2[0],v_weight[0])
                h_CTag_sum_max3.Fill(v_CTag_sum_max3[0],v_weight[0])
                h_CTag_sum_max4.Fill(v_CTag_sum_max4[0],v_weight[0])
                h_CTag_sum_all.Fill(v_CTag_sum_all[0],v_weight[0])
                h_LTag_sum_all.Fill(v_LTag_sum_all[0],v_weight[0])
                #finish extra plots for tight mass selection
                h_jet1_theta.Fill(degrees(v_jet1_theta[0]),v_weight[0])
                h_jet2_theta.Fill(degrees(v_jet2_theta[0]),v_weight[0])
                h_jet3_theta.Fill(degrees(v_jet3_theta[0]),v_weight[0])
                h_jet4_theta.Fill(degrees(v_jet4_theta[0]),v_weight[0])
                h_jet5_theta.Fill(degrees(v_jet5_theta[0]),v_weight[0])
                h_jet6_theta.Fill(degrees(v_jet6_theta[0]),v_weight[0])
                h_jet1_E.Fill(v_jet1_E[0],v_weight[0])
                h_jet2_E.Fill(v_jet2_E[0],v_weight[0])
                h_jet3_E.Fill(v_jet3_E[0],v_weight[0])
                h_jet4_E.Fill(v_jet4_E[0],v_weight[0])
                h_jet5_E.Fill(v_jet5_E[0],v_weight[0])
                h_jet6_E.Fill(v_jet6_E[0],v_weight[0])
                h_jet1_BTag.Fill(v_jet1_BTag[0],v_weight[0])
                h_jet2_BTag.Fill(v_jet2_BTag[0],v_weight[0])
                h_jet3_BTag.Fill(v_jet3_BTag[0],v_weight[0])
                h_jet4_BTag.Fill(v_jet4_BTag[0],v_weight[0])
                h_jet5_BTag.Fill(v_jet5_BTag[0],v_weight[0])
                h_jet6_BTag.Fill(v_jet6_BTag[0],v_weight[0])
                h_y12.Fill(v_y12[0],v_weight[0])
                h_y23.Fill(v_y23[0],v_weight[0])
                h_y34.Fill(v_y34[0],v_weight[0])
                h_y45.Fill(v_y45[0],v_weight[0])
                h_y56.Fill(v_y56[0],v_weight[0])

        print 'at end of file processing of ',i_final_histo_name_,bdt_value,h_jet1_theta.Integral()
                 
        fileout.Write()
    fileout.Close()

 
    return None

def process_files():

    #for negative polarisation, NTrees 300 gives best results, for positive polarisation change to NTrees 250

    #GradShrink200
    #AdaBoostBeta020
    #noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTag_noLTagMins --> has also no noj3CMax
    #noEratios_noCosThetaHel


    print 'polm80 VLC7 files'
    files_weights_polm80_VLC7_='/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq/dataset/weights/TMVAClassification_BDT.weights.xml'

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14364.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC7_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14364_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC7_,isSignalData_)
   

    #in all events weight files, actually run with -1 in NCuts
    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14343.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC7_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14343_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC7_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hzqq_13391.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hzqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC7_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qq_13399_to_13402.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC7_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqq_13394_to_13397.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC7_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqqqq.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC7_,isSignalData_)


    print 'start VLC7 polp80 now'
    files_weights_polp80_VLC7_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq/dataset/weights/TMVAClassification_BDT.weights.xml"
                                 

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14365.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC7_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14365_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC7_,isSignalData_)
   

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14344.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC7_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14344_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC7_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hzqq_13392.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hzqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC7_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qq_13398.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC7_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqq_13393.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC7_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqqqq.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC7_,isSignalData_)





    print 'polm80 VLC11 files'
    files_weights_polm80_VLC11_='/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq/dataset/weights/TMVAClassification_BDT.weights.xml'

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14364.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC11_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14364_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC11_,isSignalData_)
   

    #in all events weight files, actually run with -1 in NCuts
    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14343.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC11_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14343_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC11_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hzqq_13391.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hzqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC11_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qq_13399_to_13402.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC11_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqq_13394_to_13397.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC11_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqqqq.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC11_,isSignalData_)


    print 'start VLC11 polp80 now'
    files_weights_polp80_VLC11_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq/dataset/weights/TMVAClassification_BDT.weights.xml"
                                 

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14365.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC11_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14365_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC11_,isSignalData_)
   

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14344.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC11_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14344_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC11_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hzqq_13392.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hzqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC11_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qq_13398.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC11_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqq_13393.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC11_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqqqq.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC11_,isSignalData_)





    print 'polm80 VLC14 files'
    files_weights_polm80_VLC14_='/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq/dataset/weights/TMVAClassification_BDT.weights.xml'

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14364.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC14_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14364_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC14_,isSignalData_)
   

    #in all events weight files, actually run with -1 in NCuts
    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14343.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC14_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14343_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC14_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hzqq_13391.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hzqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC14_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qq_13399_to_13402.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC14_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqq_13394_to_13397.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC14_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqqqq.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polm80_VLC14_,isSignalData_)


    print 'start VLC14 polp80 now'
    files_weights_polp80_VLC14_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq/dataset/weights/TMVAClassification_BDT.weights.xml"
                                 

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14365.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC14_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhqq_14365_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC14_,isSignalData_)
   

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14344.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC14_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hhz_14344_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hhz__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_AllEvents.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC14_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_hzqq_13392.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__hzqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC14_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qq_13398.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC14_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqq_13393.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC14_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_M2_75_M3_50_150_ee_qqqqqq.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag_tight_Mass_Cuts/MVATrainingReader_BDT__ee_qqqqqq__HHZ__ee_qqqqqq__histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins.root"
    process_event(final_histo_name_,input_file_,files_weights_polp80_VLC14_,isSignalData_)


    return None

process_files()



