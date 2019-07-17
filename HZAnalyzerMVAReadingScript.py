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

    v_jet1_mass = array('f',[0])
    reader.AddVariable('jet1_mass',v_jet1_mass)
    v_jet2_mass = array('f',[0])
    reader.AddVariable('jet2_mass',v_jet2_mass)
    v_jet1_theta = array('f',[0])
    reader.AddVariable('jet1_theta',v_jet1_theta)
    v_jet2_theta = array('f',[0])
    reader.AddVariable('jet2_theta',v_jet2_theta)
    #v_deltatheta = array('f',[0])
    #reader.AddVariable('deltatheta:=jet1_theta-jet2_theta',v_deltatheta)
    v_jet1_D2_beta1 = array('f',[0])
    reader.AddVariable('jet1_D2_beta1',v_jet1_D2_beta1)
    v_jet2_D2_beta1 = array('f',[0])
    reader.AddVariable('jet2_D2_beta1',v_jet2_D2_beta1)
    v_jet1_BTag_rfj_BTagMax = array('f',[0])
    reader.AddVariable('jet1_BTag_rfj_BTagMax', v_jet1_BTag_rfj_BTagMax)

    v_jet1_C2_beta1 = array('f',[0])
    reader.AddVariable('jet1_C2_beta1',v_jet1_C2_beta1)
    v_jet2_C2_beta1 = array('f',[0])
    reader.AddVariable('jet2_C2_beta1',v_jet2_C2_beta1)
    v_jet1_tau21 = array('f',[0])
    reader.AddVariable('jet1_tau21',v_jet1_tau21)
    v_jet2_tau21 = array('f',[0])
    reader.AddVariable('jet2_tau21',v_jet2_tau21)


    #s_jet1_C2_beta1 = array('f',[0])
    #reader.AddSpectator('jet1_C2_beta1',s_jet1_C2_beta1)
    #s_jet2_C2_beta1 = array('f',[0])
    #reader.AddSpectator('jet2_C2_beta1',s_jet2_C2_beta1)
    #s_jet1_tau21 = array('f',[0])
    #reader.AddSpectator('jet1_tau21',s_jet1_tau21)
    #s_jet2_tau21 = array('f',[0])
    #reader.AddSpectator('jet2_tau21',s_jet2_tau21)

    input_file_=root.TFile.Open(i_input_file_name_)
    tree = input_file_.Get("MVATrainingVariables")
    tree.SetBranchAddress('jet1_mass',v_jet1_mass)
    tree.SetBranchAddress('jet2_mass',v_jet2_mass)
    tree.SetBranchAddress('jet1_theta',v_jet1_theta)
    tree.SetBranchAddress('jet2_theta',v_jet2_theta)
    tree.SetBranchAddress('jet1_D2_beta1',v_jet1_D2_beta1)
    tree.SetBranchAddress('jet2_D2_beta1',v_jet2_D2_beta1)
    tree.SetBranchAddress('jet1_BTag_rfj_BTagMax', v_jet1_BTag_rfj_BTagMax)
    tree.SetBranchAddress('jet1_C2_beta1',v_jet1_C2_beta1)
    tree.SetBranchAddress('jet2_C2_beta1',v_jet2_C2_beta1)
    tree.SetBranchAddress('jet1_tau21',v_jet1_tau21)
    tree.SetBranchAddress('jet2_tau21',v_jet2_tau21)
    v_eventWeight = array('f',[0])
    tree.SetBranchAddress('eventWeight',v_eventWeight)

    reader.BookMVA('BDT',TString(i_weight_file))

    v_jet1_E = array('f',[0])
    v_jet2_E = array('f',[0])
    tree.SetBranchAddress('jet1_E',v_jet1_E)
    tree.SetBranchAddress('jet2_E',v_jet2_E)

    v_jet1_phi = array('f',[0])
    v_jet2_phi = array('f',[0])
    tree.SetBranchAddress('jet1_phi',v_jet1_phi)
    tree.SetBranchAddress('jet2_phi',v_jet2_phi)

    v_jet1_Pt = array('f',[0])
    v_jet2_Pt = array('f',[0])
    tree.SetBranchAddress('jet1_Pt',v_jet1_Pt)
    tree.SetBranchAddress('jet2_Pt',v_jet2_Pt)

    v_jet1_sj1_E = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_E',v_jet1_sj1_E)
    v_jet1_sj2_E = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_E',v_jet1_sj2_E)
    v_jet1_sj1_Px = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_Px',v_jet1_sj1_Px)
    v_jet1_sj2_Px = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_Px',v_jet1_sj2_Px)
    v_jet1_sj1_Py = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_Py',v_jet1_sj1_Py)
    v_jet1_sj2_Py = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_Py',v_jet1_sj2_Py)
    v_jet1_sj1_Pz = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_Pz',v_jet1_sj1_Pz)
    v_jet1_sj2_Pz = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_Pz',v_jet1_sj2_Pz)

    v_jet1_sj1_nTracks = array('i',[0])
    tree.SetBranchAddress('jet1_sj1_nTracks',v_jet1_sj1_nTracks)
    v_jet1_sj2_nTracks = array('i',[0])
    tree.SetBranchAddress('jet1_sj2_nTracks',v_jet1_sj2_nTracks)
    v_jet2_sj1_nTracks = array('i',[0])
    tree.SetBranchAddress('jet2_sj1_nTracks',v_jet2_sj1_nTracks)
    v_jet2_sj2_nTracks = array('i',[0])
    tree.SetBranchAddress('jet2_sj2_nTracks',v_jet2_sj2_nTracks)

    v_jet1_sj1_chFrac = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_chFrac',v_jet1_sj1_chFrac)
    v_jet1_sj2_chFrac = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_chFrac',v_jet1_sj2_chFrac)
    v_jet2_sj1_chFrac = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_chFrac',v_jet2_sj1_chFrac)
    v_jet2_sj2_chFrac = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_chFrac',v_jet2_sj2_chFrac)
    v_jet1_dAlpha_sj1sj2 = array('f',[0])
    tree.SetBranchAddress('jet1_dAlpha_sj1sj2',v_jet1_dAlpha_sj1sj2)
    v_jet2_dAlpha_sj1sj2 = array('f',[0])
    tree.SetBranchAddress('jet2_dAlpha_sj1sj2',v_jet2_dAlpha_sj1sj2)

    v_jet2_sj1_E = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_E',v_jet2_sj1_E)
    v_jet2_sj2_E = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_E',v_jet2_sj2_E)
    v_jet2_sj1_Px = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_Px',v_jet2_sj1_Px)
    v_jet2_sj2_Px = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_Px',v_jet2_sj2_Px)
    v_jet2_sj1_Py = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_Py',v_jet2_sj1_Py)
    v_jet2_sj2_Py = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_Py',v_jet2_sj2_Py)
    v_jet2_sj1_Pz = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_Pz',v_jet2_sj1_Pz)
    v_jet2_sj2_Pz = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_Pz',v_jet2_sj2_Pz)


    v_jet1_sj1_jetChargePt_kappa_0_25 = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_jetChargePt_kappa_0_25',v_jet1_sj1_jetChargePt_kappa_0_25)
    v_jet1_sj2_jetChargePt_kappa_0_25 = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_jetChargePt_kappa_0_25',v_jet1_sj2_jetChargePt_kappa_0_25)
    v_jet2_sj1_jetChargePt_kappa_0_25 = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_jetChargePt_kappa_0_25',v_jet2_sj1_jetChargePt_kappa_0_25)
    v_jet2_sj2_jetChargePt_kappa_0_25 = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_jetChargePt_kappa_0_25',v_jet2_sj2_jetChargePt_kappa_0_25)
    v_jet1_sj1_jetChargeE_kappa_0_25 = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_jetChargeE_kappa_0_25',v_jet1_sj1_jetChargeE_kappa_0_25)
    v_jet1_sj2_jetChargeE_kappa_0_25 = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_jetChargeE_kappa_0_25',v_jet1_sj2_jetChargeE_kappa_0_25)
    v_jet2_sj1_jetChargeE_kappa_0_25 = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_jetChargeE_kappa_0_25',v_jet2_sj1_jetChargeE_kappa_0_25)
    v_jet2_sj2_jetChargeE_kappa_0_25 = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_jetChargeE_kappa_0_25',v_jet2_sj2_jetChargeE_kappa_0_25)

    v_jet1_sj1_jetChargePt_kappa_0_30 = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_jetChargePt_kappa_0_30',v_jet1_sj1_jetChargePt_kappa_0_30)
    v_jet1_sj2_jetChargePt_kappa_0_30 = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_jetChargePt_kappa_0_30',v_jet1_sj2_jetChargePt_kappa_0_30)
    v_jet2_sj1_jetChargePt_kappa_0_30 = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_jetChargePt_kappa_0_30',v_jet2_sj1_jetChargePt_kappa_0_30)
    v_jet2_sj2_jetChargePt_kappa_0_30 = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_jetChargePt_kappa_0_30',v_jet2_sj2_jetChargePt_kappa_0_30)

    v_jet1_sj1_jetChargeE_kappa_0_30 = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_jetChargeE_kappa_0_30',v_jet1_sj1_jetChargeE_kappa_0_30)
    v_jet1_sj2_jetChargeE_kappa_0_30 = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_jetChargeE_kappa_0_30',v_jet1_sj2_jetChargeE_kappa_0_30)
    v_jet2_sj1_jetChargeE_kappa_0_30 = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_jetChargeE_kappa_0_30',v_jet2_sj1_jetChargeE_kappa_0_30)
    v_jet2_sj2_jetChargeE_kappa_0_30 = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_jetChargeE_kappa_0_30',v_jet2_sj2_jetChargeE_kappa_0_30)

    v_jet1_sj1_Angle_closestMatch = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_Angle_closestMatch',v_jet1_sj1_Angle_closestMatch)
    v_jet1_sj2_Angle_closestMatch = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_Angle_closestMatch',v_jet1_sj2_Angle_closestMatch)
    v_jet1_sj1_Angle_decMatch = array('f',[0])
    tree.SetBranchAddress('jet1_sj1_Angle_decMatch',v_jet1_sj1_Angle_decMatch)
    v_jet1_sj2_Angle_decMatch = array('f',[0])
    tree.SetBranchAddress('jet1_sj2_Angle_decMatch',v_jet1_sj2_Angle_decMatch)
    v_jet2_sj1_Angle_closestMatch = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_Angle_closestMatch',v_jet2_sj1_Angle_closestMatch)
    v_jet2_sj2_Angle_closestMatch = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_Angle_closestMatch',v_jet2_sj2_Angle_closestMatch)
    v_jet2_sj1_Angle_decMatch = array('f',[0])
    tree.SetBranchAddress('jet2_sj1_Angle_decMatch',v_jet2_sj1_Angle_decMatch)
    v_jet2_sj2_Angle_decMatch = array('f',[0])
    tree.SetBranchAddress('jet2_sj2_Angle_decMatch',v_jet2_sj2_Angle_decMatch)

    v_jet1_sj1_closestMatch = array('i',[0])
    tree.SetBranchAddress('jet1_sj1_closestMatch',v_jet1_sj1_closestMatch)
    v_jet1_sj2_closestMatch = array('i',[0])
    tree.SetBranchAddress('jet1_sj2_closestMatch',v_jet1_sj2_closestMatch)
    v_jet1_sj1_decMatch = array('i',[0])
    tree.SetBranchAddress('jet1_sj1_decMatch',v_jet1_sj1_decMatch)
    v_jet1_sj2_decMatch = array('i',[0])
    tree.SetBranchAddress('jet1_sj2_decMatch',v_jet1_sj2_decMatch)
    v_jet2_sj1_closestMatch = array('i',[0])
    tree.SetBranchAddress('jet2_sj1_closestMatch',v_jet2_sj1_closestMatch)
    v_jet2_sj2_closestMatch = array('i',[0])
    tree.SetBranchAddress('jet2_sj2_closestMatch',v_jet2_sj2_closestMatch)
    v_jet2_sj1_decMatch = array('i',[0])
    tree.SetBranchAddress('jet2_sj1_decMatch',v_jet2_sj1_decMatch)
    v_jet2_sj2_decMatch = array('i',[0])
    tree.SetBranchAddress('jet2_sj2_decMatch',v_jet2_sj2_decMatch)

    v_parton_H_E = array('f',[0])
    v_parton_H_Px = array('f',[0])  
    v_parton_H_Py = array('f',[0])
    v_parton_H_Pz = array('f',[0])

    v_parton_Z_qneg_E = array('f',[0])
    v_parton_Z_qneg_Px = array('f',[0])  
    v_parton_Z_qneg_Py = array('f',[0])
    v_parton_Z_qneg_Pz = array('f',[0])
    v_parton_Z_qneg_PDGID = array('i',[0])
    
    v_parton_Z_qpos_E = array('f',[0])
    v_parton_Z_qpos_Px = array('f',[0])  
    v_parton_Z_qpos_Py = array('f',[0])
    v_parton_Z_qpos_Pz = array('f',[0])
    v_parton_Z_qpos_PDGID = array('i',[0])

    v_parton_em_E = array('f',[0])
    v_parton_em_Px = array('f',[0])  
    v_parton_em_Py = array('f',[0])
    v_parton_em_Pz = array('f',[0])
    
    v_parton_ep_E = array('f',[0])
    v_parton_ep_Px = array('f',[0])  
    v_parton_ep_Py = array('f',[0])
    v_parton_ep_Pz = array('f',[0])

    v_genjet1_sj1_E = array('f',[0])
    v_genjet1_sj1_Px = array('f',[0])  
    v_genjet1_sj1_Py = array('f',[0])
    v_genjet1_sj1_Pz = array('f',[0])
    v_genjet1_sj1_jetChargeE_kappa_0_30 = array('f',[0])
    v_genjet1_sj1_jetChargePt_kappa_0_30 = array('f',[0])
    
    v_genjet1_sj2_E = array('f',[0])
    v_genjet1_sj2_Px = array('f',[0])  
    v_genjet1_sj2_Py = array('f',[0])
    v_genjet1_sj2_Pz = array('f',[0])
    v_genjet1_sj2_jetChargeE_kappa_0_30 = array('f',[0])
    v_genjet1_sj2_jetChargePt_kappa_0_30 = array('f',[0])

    v_genjet2_sj1_E = array('f',[0])
    v_genjet2_sj1_Px = array('f',[0])  
    v_genjet2_sj1_Py = array('f',[0])
    v_genjet2_sj1_Pz = array('f',[0])
    v_genjet2_sj1_jetChargeE_kappa_0_30 = array('f',[0])
    v_genjet2_sj1_jetChargePt_kappa_0_30 = array('f',[0])
    
    v_genjet2_sj2_E = array('f',[0])
    v_genjet2_sj2_Px = array('f',[0])  
    v_genjet2_sj2_Py = array('f',[0])
    v_genjet2_sj2_Pz = array('f',[0])
    v_genjet2_sj2_jetChargeE_kappa_0_30 = array('f',[0])
    v_genjet2_sj2_jetChargePt_kappa_0_30 = array('f',[0])


    if(i_isSignalData):
        tree.SetBranchAddress('parton_H_E',v_parton_H_E)  
        tree.SetBranchAddress('parton_H_Px',v_parton_H_Px)
        tree.SetBranchAddress('parton_H_Py',v_parton_H_Py)
        tree.SetBranchAddress('parton_H_Pz',v_parton_H_Pz)

        tree.SetBranchAddress('parton_Z_qpos_E',v_parton_Z_qpos_E)  
        tree.SetBranchAddress('parton_Z_qpos_Px',v_parton_Z_qpos_Px)
        tree.SetBranchAddress('parton_Z_qpos_Py',v_parton_Z_qpos_Py)
        tree.SetBranchAddress('parton_Z_qpos_Pz',v_parton_Z_qpos_Pz)
        tree.SetBranchAddress('parton_Z_qpos_PDGID',v_parton_Z_qpos_PDGID)

        tree.SetBranchAddress('parton_Z_qneg_E',v_parton_Z_qneg_E)  
        tree.SetBranchAddress('parton_Z_qneg_Px',v_parton_Z_qneg_Px)
        tree.SetBranchAddress('parton_Z_qneg_Py',v_parton_Z_qneg_Py)
        tree.SetBranchAddress('parton_Z_qneg_Pz',v_parton_Z_qneg_Pz)
        tree.SetBranchAddress('parton_Z_qneg_PDGID',v_parton_Z_qneg_PDGID)

        tree.SetBranchAddress('parton_ep_E',v_parton_ep_E)  
        tree.SetBranchAddress('parton_ep_Px',v_parton_ep_Px)
        tree.SetBranchAddress('parton_ep_Py',v_parton_ep_Py)
        tree.SetBranchAddress('parton_ep_Pz',v_parton_ep_Pz)

        tree.SetBranchAddress('parton_em_E',v_parton_em_E)  
        tree.SetBranchAddress('parton_em_Px',v_parton_em_Px)
        tree.SetBranchAddress('parton_em_Py',v_parton_em_Py)
        tree.SetBranchAddress('parton_em_Pz',v_parton_em_Pz)

        tree.SetBranchAddress('genjet1_sj1_E',v_genjet1_sj1_E)  
        tree.SetBranchAddress('genjet1_sj1_Px',v_genjet1_sj1_Px)
        tree.SetBranchAddress('genjet1_sj1_Py',v_genjet1_sj1_Py)
        tree.SetBranchAddress('genjet1_sj1_Pz',v_genjet1_sj1_Pz)
        tree.SetBranchAddress('genjet1_sj1_jetChargeE_kappa_0_30',v_genjet1_sj1_jetChargeE_kappa_0_30)
        tree.SetBranchAddress('genjet1_sj1_jetChargeE_kappa_0_30',v_genjet1_sj1_jetChargeE_kappa_0_30)

        tree.SetBranchAddress('genjet1_sj2_E',v_genjet1_sj2_E)  
        tree.SetBranchAddress('genjet1_sj2_Px',v_genjet1_sj2_Px)
        tree.SetBranchAddress('genjet1_sj2_Py',v_genjet1_sj2_Py)
        tree.SetBranchAddress('genjet1_sj2_Pz',v_genjet1_sj2_Pz)
        tree.SetBranchAddress('genjet1_sj2_jetChargeE_kappa_0_30',v_genjet1_sj2_jetChargeE_kappa_0_30)
        tree.SetBranchAddress('genjet1_sj2_jetChargeE_kappa_0_30',v_genjet1_sj2_jetChargeE_kappa_0_30)

        tree.SetBranchAddress('genjet2_sj1_E',v_genjet2_sj1_E)  
        tree.SetBranchAddress('genjet2_sj1_Px',v_genjet2_sj1_Px)
        tree.SetBranchAddress('genjet2_sj1_Py',v_genjet2_sj1_Py)
        tree.SetBranchAddress('genjet2_sj1_Pz',v_genjet2_sj1_Pz)
        tree.SetBranchAddress('genjet2_sj1_jetChargeE_kappa_0_30',v_genjet2_sj1_jetChargeE_kappa_0_30)
        tree.SetBranchAddress('genjet2_sj1_jetChargeE_kappa_0_30',v_genjet2_sj1_jetChargeE_kappa_0_30)
        tree.SetBranchAddress('genjet2_sj1_jetChargePt_kappa_0_30',v_genjet2_sj1_jetChargePt_kappa_0_30)
        tree.SetBranchAddress('genjet2_sj1_jetChargePt_kappa_0_30',v_genjet2_sj1_jetChargePt_kappa_0_30)

        tree.SetBranchAddress('genjet2_sj2_E',v_genjet2_sj2_E)  
        tree.SetBranchAddress('genjet2_sj2_Px',v_genjet2_sj2_Px)
        tree.SetBranchAddress('genjet2_sj2_Py',v_genjet2_sj2_Py)
        tree.SetBranchAddress('genjet2_sj2_Pz',v_genjet2_sj2_Pz)
        tree.SetBranchAddress('genjet2_sj2_jetChargeE_kappa_0_30',v_genjet2_sj2_jetChargeE_kappa_0_30)
        tree.SetBranchAddress('genjet2_sj2_jetChargeE_kappa_0_30',v_genjet2_sj2_jetChargeE_kappa_0_30)
        tree.SetBranchAddress('genjet2_sj2_jetChargePt_kappa_0_30',v_genjet2_sj2_jetChargePt_kappa_0_30)
        tree.SetBranchAddress('genjet2_sj2_jetChargePt_kappa_0_30',v_genjet2_sj2_jetChargePt_kappa_0_30)
        
    n_bins_high=100
    lim_BDT_low=-1.0
    lim_BDT_high=1.0
    lim_mass_low=0
    lim_mass_high=200
    lim_theta_low=-0.5
    lim_theta_high=180.5

    lim_cosTheta_low=-1.0
    lim_cosTheta_high=1.0

    lim_theta_low=-0.0
    lim_theta_high=180.0

    lim_dtheta_low=-180.5
    lim_dtheta_high=180.5

    lim_phi_low=-0.0
    lim_phi_high=360.0

    
    lim_BTag_low=0
    lim_BTag_high=0

    lim_D2_low=0
    lim_D2_high=10


    lim_C2_low=0
    lim_C2_high=0.55

    lim_tau21_low=0
    lim_tau21_high=1.3

    n_bins_tracks=31
    lim_ntracks_low=-0.5
    lim_ntracks_high=30.5

    lim_chfrac_low=0.0
    lim_chfrac_high=1.0

    lim_jetCharge_low=-2.5
    lim_jetCharge_high=2.5

    lim_dAlpha_sj12_low=0
    lim_dAlpha_sj12_high=25

    lim_Efrac_low=0.0
    lim_Efrac_high=0.5

    lim_jetE_ratio_low=0.5
    lim_jetE_ratio_high=1.5

    lim_jetP_ratio_low=0.5
    lim_jetP_ratio_high=1.5

    nbins_charge=20
    lim_charge_low=-1.50
    lim_charge_high=1.50

    lim_cosProd_low=-1.05
    lim_cosProd_high=1.05

    lim_cosProdNorm_low=-1.55
    lim_cosProdNorm_high=1.55

    n_bins_high_2D=40


    BDT_cuts=[]
    #print 'do i get here maybe 1'
    BDT_cuts.append(-0.200)
    #print 'do i get here maybe 2'
    BDT_cuts.append(-0.150)
    BDT_cuts.append(-0.100)
    BDT_cuts.append(-0.050)
    BDT_cuts.append(0.000)
    BDT_cuts.append(0.050)
    BDT_cuts.append(0.100)
    BDT_cuts.append(0.150)
    BDT_cuts.append(0.175)
    BDT_cuts.append(0.200)
    BDT_cuts.append(0.225)
    BDT_cuts.append(0.250)
    BDT_cuts.append(0.275)
    BDT_cuts.append(0.300)
    BDT_cuts.append(0.325)
    BDT_cuts.append(0.350)
    BDT_cuts.append(0.375)
    BDT_cuts.append(0.400)
    BDT_cuts.append(0.425)
    BDT_cuts.append(0.450)
    BDT_cuts.append(0.500)
    BDT_cuts.append(0.550)
    BDT_cuts.append(0.600)
    BDT_cuts.append(0.650)
    BDT_cuts.append(0.700)
    BDT_cuts.append(0.750)
    BDT_cuts.append(0.800)
    BDT_cuts.append(0.850)
    BDT_cuts.append(0.875)
    BDT_cuts.append(0.900)
    BDT_cuts.append(0.925)
    BDT_cuts.append(0.950)
    BDT_cuts.append(0.975)
    BDT_cuts.append(1.000)
    fileout = root.TFile(i_final_histo_name_,"RECREATE")
    for bdt_value in BDT_cuts:
        tdirectory=fileout.GetDirectory(i_final_histo_name_)
        fileout.mkdir(str(bdt_value))
        fileout.cd(str(bdt_value))

        h_BDT_output = TH1F( "h_BDT_output", "", n_bins_high, lim_BDT_low, lim_BDT_high)
        h_jet1_mass = TH1F( "h_jet1_mass", "", n_bins_high, lim_mass_low, lim_mass_high)
        h_jet2_mass = TH1F( "h_jet2_mass", "", n_bins_high, lim_mass_low, lim_mass_high)
        h_jet1_theta = TH1F( "h_jet1_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_jet2_theta = TH1F( "h_jet2_theta", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_jet1_min_jet2_theta = TH1F( "h_jet1_min_jet2_theta", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
        h_jet1_BTag = TH1F( "h_jet1_Btag", "", n_bins_high, lim_BTag_low, lim_BTag_high)
        h_jet1_D2 = TH1F( "h_jet1_D2_beta1", "", n_bins_high, lim_D2_low, lim_D2_high)
        h_jet2_D2 = TH1F( "h_jet2_D2_beta1", "", n_bins_high, lim_D2_low, lim_D2_high)
        h_jet1_C2 = TH1F( "h_jet1_C2_beta1", "", n_bins_high, lim_C2_low, lim_C2_high)
        h_jet2_C2 = TH1F( "h_jet2_C2_beta1", "", n_bins_high, lim_C2_low, lim_C2_high)
        h_jet1_tau21 = TH1F( "h_jet1_tau21", "", n_bins_high, lim_tau21_low, lim_tau21_high)
        h_jet2_tau21 = TH1F( "h_jet2_tau21", "", n_bins_high, lim_tau21_low, lim_tau21_high)
        h_BDT_output_Eff = TH1F( "h_BDT_output_wCut", "", n_bins_high, lim_BDT_low, lim_BDT_high)
 
        
        h_signal_background_1D_hist_list=[]
        h_signal_background_1D_hist_list.append(h_BDT_output)
        h_signal_background_1D_hist_list.append(h_jet1_mass)
        h_signal_background_1D_hist_list.append(h_jet2_mass)
        h_signal_background_1D_hist_list.append(h_jet1_theta)
        h_signal_background_1D_hist_list.append(h_jet2_theta)
        h_signal_background_1D_hist_list.append(h_jet1_min_jet2_theta)
        h_signal_background_1D_hist_list.append(h_jet1_BTag)
        h_signal_background_1D_hist_list.append(h_jet1_D2)
        h_signal_background_1D_hist_list.append(h_jet2_D2)
        h_signal_background_1D_hist_list.append(h_jet1_C2)
        h_signal_background_1D_hist_list.append(h_jet2_C2)
        h_signal_background_1D_hist_list.append(h_jet1_tau21)
        h_signal_background_1D_hist_list.append(h_jet2_tau21)
        h_signal_background_1D_hist_list.append(h_BDT_output_Eff)


        h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com = TH1F( "h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_rj2com = TH1F( "h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_rj2com", "", n_bins_high, lim_theta_low, lim_theta_high)

        h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com = TH1F( "h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)
        h_recojet_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com = TH1F( "h_recojet_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)

        h_recojet_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com = TH1F( "h_recojet_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)
        h_recojet_sgncos2theta1_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com = TH1F( "h_recojet_sgncos2theta1_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)

        h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_dm_0_20_rj2com = TH1F( "h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_dm_0_20_rj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_dm_0_20_rj2com = TH1F( "h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_dm_0_20_rj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_dm_0_30_rj2com = TH1F( "h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_dm_0_30_rj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_dm_0_30_rj2com = TH1F( "h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_dm_0_30_rj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
        
        h_recojet_theta_rj1_rj2_E_totCOM = TH1F("h_recojet_theta_rj1_rj2_E_totCOM", "", n_bins_high, lim_theta_low, lim_theta_high)
        h_recojet_theta2_rj1_ep_E_totCOM = TH1F( "h_recojet_theta2_rj1_ep_E_totCOM", "", n_bins_high, lim_theta_low, lim_theta_high)

        h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
        h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)

        h_recojet_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
        h_recojet_sgnsinphi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_sgnsinphi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
        h_recojet_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
        h_recojet_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
        h_recojet_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
        h_recojet_sgncosphi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_sgncosphi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
        h_recojet_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
        h_recojet_sgncos2phi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_sgncos2phi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)

        h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_dm_0_20_vs_plane_rj1_ep = TH1F( "h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_dm_0_20_vs_plane_rj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
        h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_dm_0_20_vs_plane_rj1_ep = TH1F( "h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_dm_0_20_vs_plane_rj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
        h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_dm_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_dm_0_30_vs_plane_rj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
        h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_dm_0_30_vs_plane_rj1_ep = TH1F( "h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_dm_0_30_vs_plane_rj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
        h_recojet_theta_sj1_sj2_rj2__rj2COM = TH1F("h_recojet_theta_sj1_sj2_rj2__rj2COM", "", n_bins_high, lim_theta_low, lim_theta_high)

        h_signal_background_1D_hist_list.append(h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
        h_signal_background_1D_hist_list.append(h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_rj2com)
        h_signal_background_1D_hist_list.append(h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_dm_0_20_rj2com)
        h_signal_background_1D_hist_list.append(h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_dm_0_20_rj2com)
        h_signal_background_1D_hist_list.append(h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_dm_0_30_rj2com)
        h_signal_background_1D_hist_list.append(h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_dm_0_30_rj2com)
        
        h_signal_background_1D_hist_list.append(h_recojet_theta_rj1_rj2_E_totCOM)
        h_signal_background_1D_hist_list.append(h_recojet_theta2_rj1_ep_E_totCOM)

        h_signal_background_1D_hist_list.append(h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_dm_0_20_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_dm_0_20_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_dm_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_dm_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_theta_sj1_sj2_rj2__rj2COM)

        h_signal_background_1D_hist_list.append(h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
        h_signal_background_1D_hist_list.append(h_recojet_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com)
        h_signal_background_1D_hist_list.append(h_recojet_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com)
        h_signal_background_1D_hist_list.append(h_recojet_sgncos2theta1_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com)

        h_signal_background_1D_hist_list.append(h_recojet_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_sgnsinphi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_sgncosphi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep)
        h_signal_background_1D_hist_list.append(h_recojet_sgncos2phi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep)

        for hist in h_signal_background_1D_hist_list:
            hist.Sumw2()


        h_2D_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM = TH2F( "h_2D_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
        h_2D_recojet_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM = TH2F( "h_2D_recojet_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
        h_2D_recojet_sgncostheta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_sgncostheta2_costheta2_rj1_ep_E_totCOM = TH2F( "h_2D_recojet_sgncostheta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_sgncostheta2_costheta2_rj1_ep_E_totCOM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
        h_2D_recojet_sgncostheta1_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com_vs_sgncostheta2_costheta2_rj1_ep_E_totCOM = TH2F( "h_2D_recojet_sgncostheta1_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com_vs_sgncostheta2_costheta2_rj1_ep_E_totCOM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
        h_2D_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Sumw2()
        h_2D_recojet_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Sumw2()
        h_2D_recojet_sgncostheta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_sgncostheta2_costheta2_rj1_ep_E_totCOM.Sumw2()
        h_2D_recojet_sgncostheta1_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com_vs_sgncostheta2_costheta2_rj1_ep_E_totCOM.Sumw2()

        h_signal_jetChargeHistos_1D_hist_list=[]
        if(i_isSignalData):    
            print 'histos should be defined'
            h_jet2_mass_H_matched = TH1F( "h_jet2_mass_H_matched", "",n_bins_high, lim_mass_low, lim_mass_high)
            h_jet2_mass_Z_matched = TH1F( "h_jet2_mass_Z_matched", "",n_bins_high, lim_mass_low, lim_mass_high)

            h_jetE_reco_over_parton_H_matched_orig = TH1F( "h_jetE_reco_over_parton_H_matched_orig", "",n_bins_high,lim_jetE_ratio_low, lim_jetE_ratio_high)
            h_jetE_reco_over_parton_H_matched_corr = TH1F( "h_jetE_reco_over_parton_H_matched_corr", "",n_bins_high,lim_jetE_ratio_low, lim_jetE_ratio_high)

            h_jetP_reco_over_parton_H_matched_orig = TH1F( "h_jetP_reco_over_parton_H_matched_orig", "",n_bins_high,lim_jetP_ratio_low, lim_jetP_ratio_high)
            h_jetP_reco_over_parton_H_matched_corr = TH1F( "h_jetP_reco_over_parton_H_matched_corr", "",n_bins_high,lim_jetP_ratio_low, lim_jetP_ratio_high)

            h_jetE_reco_over_parton_Z_matched_orig = TH1F( "h_jetE_reco_over_parton_Z_matched_orig", "",n_bins_high,lim_jetE_ratio_low, lim_jetE_ratio_high)
            h_jetE_reco_over_parton_Z_matched_corr = TH1F( "h_jetE_reco_over_parton_Z_matched_corr", "",n_bins_high,lim_jetE_ratio_low, lim_jetE_ratio_high)

            h_jetP_reco_over_parton_Z_matched_orig = TH1F( "h_jetP_reco_over_parton_Z_matched_orig", "",n_bins_high,lim_jetP_ratio_low, lim_jetP_ratio_high)
            h_jetP_reco_over_parton_Z_matched_corr = TH1F( "h_jetP_reco_over_parton_Z_matched_corr", "",n_bins_high,lim_jetP_ratio_low, lim_jetP_ratio_high)

            h_jet2_sj_b_cM_jetChargeE_kappa_0_30= TH1F( "h_jet2_sj_b_cM_jetChargeE_kappa_0_30", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_bbar_cM_jetChargeE_kappa_0_30= TH1F( "h_jet2_sj_bbar_cM_jetChargeE_kappa_0_30", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_b_dM_jetChargeE_kappa_0_30= TH1F( "h_jet2_sj_b_dM_jetChargeE_kappa_0_30", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_bbar_dM_jetChargeE_kappa_0_30= TH1F( "h_jet2_sj_bbar_dM_jetChargeE_kappa_0_30", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            
            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30= TH1F( "h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30= TH1F( "h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30= TH1F( "h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30= TH1F( "h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            
            #decided and closest the same for BOTH subjets
            h_jet2_sj_b_cM_jetChargeE_kappa_0_30_dec_also_close= TH1F( "h_jet2_sj_b_cM_jetChargeE_kappa_0_30_dec_also_close", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_bbar_cM_jetChargeE_kappa_0_30_dec_also_close= TH1F( "h_jet2_sj_bbar_cM_jetChargeE_kappa_0_30_dec_also_close", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close= TH1F( "h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30_dec_also_close", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close= TH1F( "h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30_dec_also_close", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            
            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_oppCharge= TH1F( "h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30_oppCharge", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_oppCharge= TH1F( "h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30_oppCharge", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            
            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_chFrac= TH1F( "h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_chFrac", "",n_bins_high, lim_chfrac_low, lim_chfrac_high)
            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_chFrac= TH1F( "h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_chFrac", "",n_bins_high, lim_chfrac_low, lim_chfrac_high)
            
            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_nTrack= TH1F( "h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_nTrack", "",n_bins_high, lim_ntracks_low, lim_ntracks_high)
            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_nTrack= TH1F( "h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_nTrack", "",n_bins_high, lim_ntracks_low, lim_ntracks_high)
            
            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_chFrac= TH1F( "h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_chFrac", "",n_bins_high, lim_chfrac_low, lim_chfrac_high)
            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_chFrac= TH1F( "h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_chFrac", "",n_bins_high, lim_chfrac_low, lim_chfrac_high)
            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_nTrack= TH1F( "h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_nTrack", "",n_bins_high, lim_ntracks_low, lim_ntracks_high)
            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_nTrack= TH1F( "h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_nTrack", "",n_bins_high, lim_ntracks_low, lim_ntracks_high)

            h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qneg= TH1F( "h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qneg", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg= TH1F( "h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qpos= TH1F( "h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qpos", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qpos= TH1F( "h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qpos", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos= TH1F( "h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_qpos", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos= TH1F( "h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_qpos", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)


            h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qneg= TH1F( "h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qneg", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg= TH1F( "h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qpos= TH1F( "h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qpos", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qpos= TH1F( "h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qpos", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)


            h_jet2_E_sj2_over_E_sum_dec_also_close= TH1F( "h_jet2_E_sj2_over_E_sum_dec_also_close", "",n_bins_high, lim_Efrac_low, lim_Efrac_high)            
            h_jet2_E_sj2_over_E_sum_jetChargeE_kappa_0_30_dec_also_close_rightCharges= TH1F( "h_jet2_E_sj2_over_E_sum_jetChargeE_kappa_0_30_dec_also_close_rightCharges", "",n_bins_high, lim_Efrac_low, lim_Efrac_high)
            h_jet2_E_sj2_over_E_sum_jetChargeE_kappa_0_30_dec_also_close_wrongCharges= TH1F( "h_jet2_E_sj2_over_E_sum_jetChargeE_kappa_0_30_dec_also_close_wrongCharges", "",n_bins_high, lim_Efrac_low, lim_Efrac_high)  
            h_jet2_dAlpha_sj1sj2_dec_also_close= TH1F( "h_jet2_dAlpha_sj1sj2_sum_dec_also_close", "",n_bins_high, lim_dAlpha_sj12_low, lim_dAlpha_sj12_high)    
            h_jet2_dAlpha_sj1sj2_jetChargeE_kappa_0_30_dec_also_close_rightCharges= TH1F( "h_jet2_dAlpha_sj1sj2_jetChargeE_kappa_0_30_dec_also_close_rightCharges", "",n_bins_high, lim_dAlpha_sj12_low, lim_dAlpha_sj12_high)
            h_jet2_dAlpha_sj1sj2_jetChargeE_kappa_0_30_dec_also_close_wrongCharges= TH1F( "h_jet2_dAlpha_sj1sj2_jetChargeE_kappa_0_30_dec_also_close_wrongCharges", "",n_bins_high, lim_dAlpha_sj12_low, lim_dAlpha_sj12_high)

            h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20= TH1F( "h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_jetCharge_0_20", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20= TH1F( "h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_jetCharge_0_20", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
 
            h_jet2_sj_qneg_dM_ChargeHemisphere_kappa_0_30_dm_jetCharge_0_20= TH1F("h_jet2_sj_qneg_dM_ChargeHemisphere_kappa_0_30_dm_jetCharge_0_20", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high)

            h_parton_theta_Z_qneg_Zcom = TH1F( "h_parton_theta_Z_qneg_Zcom", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_parton_theta1_Z_qpos_Zcom = TH1F( "h_parton_theta1_Z_qpos_Zcom", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_parton_costheta1_Z_qpos_Zcom = TH1F( "h_parton_costheta1_Z_qpos_Zcom", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)
            h_parton_sgncos2theta1_costheta1_Z_qpos_Zcom = TH1F( "h_parton_sgncos2theta1_costheta1_Z_qpos_Zcom", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)

            h_parton_theta_Z_qpos_Z_qneg_Zcom = TH1F( "h_parton_theta_Z_qpos_Z_neg_Zcom", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_parton_theta_H_Z_E_totCOM = TH1F("h_parton_theta_H_Z_E_totCOM", "", n_bins_high, lim_theta_low, lim_theta_high)


            h_parton_cosTheta_Z_qneg_Zcom = TH1F( "h_parton_cosTheta_Z_qneg_Zcom", "", n_bins_high, lim_cosTheta_low, lim_cosTheta_high)
            h_parton_cosTheta_Z_qpos_Zcom = TH1F( "h_parton_cosTheta_Z_qpos_Zcom", "", n_bins_high, lim_cosTheta_low, lim_cosTheta_high)

            h_parton_theta2_H_ep_HZ_COM = TH1F( "h_parton_theta2_H_ep_HZ_COM", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_parton_theta2_H_ep_approx_HZ_COM = TH1F( "h_parton_theta2_H_ep_approx_HZ_COM", "", n_bins_high, lim_theta_low, lim_theta_high)

            h_parton_polarAngle_Z_qneg_Zcom = TH1F( "h_parton_polarAngle_Z_qneg_Zcom", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_parton_polarAngle_Z_qpos_Zcom = TH1F( "h_parton_polarAngle_Z_qpos_Zcom", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_parton_cosPolarAngle_Z_qneg_Zcom = TH1F( "h_parton_cosPolarAngle_Z_qneg_Zcom", "", n_bins_high, lim_cosTheta_low, lim_cosTheta_high)
            h_parton_cosPolarAngle_Z_qpos_Zcom = TH1F( "h_parton_cosPolarAngle_Z_qpos_Zcom", "", n_bins_high, lim_cosTheta_low, lim_cosTheta_high)
            h_parton_phi_plane_Z_qpos_vs_plane_H_ep = TH1F( "h_parton_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
         
            h_parton_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            h_parton_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            h_parton_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            h_parton_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            
            h_jet2_sj_qneg_dM_vs_qpos_jetChargeE_kappa_0_30_oppCharge= TH2F( "h_jet2_sj_qneg_cM_vs_qpos_dM_jetChargeE_kappa_0_30_oppCharge", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high,n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qneg_dM_vs_qpos_jetChargeE_kappa_0_30= TH2F( "h_jet2_sj_qneg_cM_vs_qpos_dM_jetChargeE_kappa_0_30", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high,n_bins_high, lim_jetCharge_low, lim_jetCharge_high) 
            h_jet2_sj_qneg_dM_vs_qpos_jetChargeE_kappa_0_30_oppCharge.Sumw2()
            h_jet2_sj_qneg_dM_vs_qpos_jetChargeE_kappa_0_30.Sumw2()

 
 


            h_jet2_sj_qneg_dM_jetChargeE_vs_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20=TH2F( "h_jet2_sj_qneg_dM_jetChargeE_vs_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high,n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qneg_dM_jetChargeE_vs_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20.Sumw2()
            h_jet2_sj_qneg_dM_jetChargeE_vs_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_25=TH2F( "h_jet2_sj_qneg_dM_jetChargeE_vs_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_25", "",n_bins_high, lim_jetCharge_low, lim_jetCharge_high,n_bins_high, lim_jetCharge_low, lim_jetCharge_high)
            h_jet2_sj_qneg_dM_jetChargeE_vs_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_25.Sumw2()


            h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_gj2com = TH1F( "h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_gj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_gj2com = TH1F( "h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_gj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_genjet_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com = TH1F( "h_genjet_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)
            h_genjet_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com = TH1F( "h_genjet_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)
            h_genjet_sgncos2theta1_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com = TH1F( "h_genjet_sgncos2theta1_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)
            h_genjet_sgncos2theta1_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com = TH1F( "h_genjet_sgncos2theta1_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)

            h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_dm_0_20_gj2com = TH1F( "h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_dm_0_20_gj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_dm_0_20_gj2com = TH1F( "h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_dm_0_20_gj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_dm_0_30_gj2com = TH1F( "h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_dm_0_30_gj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_dm_0_30_gj2com = TH1F( "h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_dm_0_30_gj2com", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_deltaAngle_part_gj_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_gj_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_gj_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)

            h_cosAngle_Product_gj_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_gj_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)

            h_cosAngle_Product_gj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_gj_E1_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_E1_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_gj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_gj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_gj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_gj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_gj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)

            h_cosAngle_normCosPart_gj_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_gj_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)

            h_cosAngle_normCosPart_gj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_gj_E1_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_E1_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_gj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_gj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_gj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_gj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_gj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)

            h_partCharge_gj_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_jcPt_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_gj_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_jcE_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)

            h_partCharge_gj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_gj_E1_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_E1_jcE_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_gj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_gj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_gj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_gj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_gj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)


            h_genjet_theta_gj1_gj2_E_totCOM = TH1F("h_genjet_theta_gj1_gj2_E_totCOM", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_genjet_theta2_gj1_ep_E_totCOM = TH1F( "h_genjet_theta2_gj1_ep_E_totCOM", "", n_bins_high, lim_theta_low, lim_theta_high)
            h_deltaAngle_part_gj_theta2_H_ep_E_totCOM = TH1F( "h_deltaAngle_part_gj_theta2_H_ep_E_totCOM", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)

            h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
            h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)

            h_genjet_sgnsinphi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_sgnsinphi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            h_genjet_sgnsinphi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_sgnsinphi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            h_genjet_sgnsin2phi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_sgnsin2phi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            h_genjet_sgnsin2phi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_sgnsin2phi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            h_genjet_sgncosphi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_sgncosphi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            h_genjet_sgncosphi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_sgncosphi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            h_genjet_sgncos2phi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_sgncos2phi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)
            h_genjet_sgncos2phi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_sgncos2phi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep", "", n_bins_high, -lim_phi_high, lim_phi_high)


            h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_dm_0_20_vs_plane_gj1_ep = TH1F( "h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_dm_0_20_vs_plane_gj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
            h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_dm_0_20_vs_plane_gj1_ep = TH1F( "h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_dm_0_20_vs_plane_gj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
            h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_dm_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_dm_0_30_vs_plane_gj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
            h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_dm_0_30_vs_plane_gj1_ep = TH1F( "h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_dm_0_30_vs_plane_gj1_ep", "", n_bins_high, lim_phi_low, lim_phi_high)

            h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep = TH1F( "h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep = TH1F( "h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_dm_0_20_vs_plane_H_ep = TH1F( "h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_dm_0_20_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_dm_0_20_vs_plane_H_ep = TH1F( "h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_dm_0_20_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_dm_0_30_vs_plane_H_ep = TH1F( "h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_dm_0_30_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_dm_0_30_vs_plane_H_ep = TH1F( "h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_dm_0_30_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)

            h_deltaAngle_part_gj_jcE_0_30_theta1_Z_q_pos_Zcom_oppJC  = TH1F( "h_deltaAngle_part_gj_jcE_0_30_theta1_Z_q_pos_Zcom_oppJC ", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_jcPt_0_30_theta1_Z_q_pos_Zcom_oppJC  = TH1F( "h_deltaAngle_part_gj_jcPt_0_30_theta1_Z_q_pos_Zcom_oppJC ", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep_oppJC = TH1F( "h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep_oppJC", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep_oppJC = TH1F( "h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep_oppJC", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)

            h_genjet_theta_sj1_sj2_gj2__gj2COM = TH1F("h_genjet_theta_sj1_sj2_gj2__gj2COM", "", n_bins_high, lim_theta_low, lim_theta_high)

            #compare recojets and parton level, only relevant for signal
            h_cosAngle_Product_rj_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)

            h_cosAngle_Product_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)

            h_cosAngle_Product_rj_nTrack_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_nTrack_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_chFrac_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_chFrac_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_chEnergy_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_chEnergy_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            h_cosAngle_Product_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_Product_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProd_low,lim_cosProd_high)
            
            h_cosAngle_normCosPart_rj_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)

            h_cosAngle_normCosPart_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)

            h_cosAngle_normCosPart_rj_nTrack_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_nTrack_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_chFrac_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_chFrac_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_chEnergy_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_chEnergy_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
            h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high,lim_cosProdNorm_low,lim_cosProdNorm_high)
    
            h_partCharge_rj_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_jcPt_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_rj_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_jcE_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)

            h_partCharge_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)
            h_partCharge_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_partCharge_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", nbins_charge,lim_charge_low,lim_charge_high)

            h_deltaAngle_part_rj_jcE_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_rj_jcE_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_jcPt_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_rj_jcPt_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom = TH1F( "h_deltaAngle_part_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)       
            h_deltaAngle_part_rj_theta2_H_ep_E_totCOM = TH1F( "h_deltaAngle_part_rj_theta2_H_ep_E_totCOM", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            
            h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep = TH1F( "h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep = TH1F( "h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_dm_0_20_vs_plane_H_ep = TH1F( "h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_dm_0_20_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_dm_0_20_vs_plane_H_ep = TH1F( "h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_dm_0_20_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_dm_0_30_vs_plane_H_ep = TH1F( "h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_dm_0_30_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_dm_0_30_vs_plane_H_ep = TH1F( "h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_dm_0_30_vs_plane_H_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            
            h_deltaAngle_part_rj_jcE_0_30_theta1_Z_q_pos_Zcom_oppJC  = TH1F( "h_deltaAngle_part_rj_jcE_0_30_theta1_Z_q_pos_Zcom_oppJC ", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_jcPt_0_30_theta1_Z_q_pos_Zcom_oppJC  = TH1F( "h_deltaAngle_part_rj_jcPt_0_30_theta1_Z_q_pos_Zcom_oppJC ", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep_oppJC = TH1F( "h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep_oppJC", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep_oppJC = TH1F( "h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep_oppJC", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)

            #compare recojets and parton level, only relevant for signal
            h_deltaAngle_gj_rj_jcE_0_30_theta_gj_jcE_0_30 = TH1F( "h_deltaAngle_gj_rj_jcE_0_30_theta_gj_jcE_0_30", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_jcPt_0_30_theta_gj_jcPt_0_30 = TH1F( "h_deltaAngle_gj_rj_jcPt_0_30_theta_gj_jcPt_0_30", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_jcE_0_30_dm_0_20_theta_gj_jcE_0_30 = TH1F( "h_deltaAngle_gj_rj_jcE_0_30_dm_0_20_theta_gj_jcE_0_30", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_jcPt_0_30_dm_0_20_theta_gj_jcPt_0_30 = TH1F( "h_deltaAngle_gj_rj_jcPt_0_30_dm_0_20_theta_gj_jcPt_0_30", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_jcE_0_30_dm_0_30_theta_gj_jcE_0_30 = TH1F( "h_deltaAngle_gj_rj_jcE_0_30_dm_0_30_theta_gj_jcE_0_30", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_jcPt_0_30_dm_0_30_theta_gj_jcPt_0_30 = TH1F( "h_deltaAngle_gj_rj_jcPt_0_30_dm_0_30_theta_gj_jcPt_0_30", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)       
            h_deltaAngle_gj_rj_phi_gj1_ep_E_totCOM = TH1F( "h_deltaAngle_gj_rj_phi_gj1_ep_E_totCOM", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            
            h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_vs_plane_gj1_ep = TH1F( "h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_vs_plane_gj1_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_vs_plane_gj1_ep = TH1F( "h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_vs_plane_gj1_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_dm_0_20_vs_plane_gj1_ep = TH1F( "h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_dm_0_20_vs_plane_gj1_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_dm_0_20_vs_plane_gj1_ep = TH1F( "h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_dm_0_20_vs_plane_gj1_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_dm_0_30_vs_plane_gj1_ep = TH1F( "h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_dm_0_30_vs_plane_gj1_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_dm_0_30_vs_plane_gj1_ep = TH1F( "h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_dm_0_30_vs_plane_gj1_ep", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            
            h_deltaAngle_gj_rj_jcE_0_30_theta_gj_jcE_0_30_oppJC  = TH1F( "h_deltaAngle_gj_rj_jcE_0_30_theta_gj_jcE_0_30_oppJC ", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_jcPt_0_30_theta_gj_jcPt_0_30_oppJC  = TH1F( "h_deltaAngle_gj_rj_jcPt_0_30_theta_gj_jcPt_0_30_oppJC ", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_vs_plane_gj1_ep_oppJC = TH1F( "h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_vs_plane_gj1_ep_oppJC", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)
            h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_vs_plane_gj1_ep_oppJC = TH1F( "h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_vs_plane_gj1_ep_oppJC", "", n_bins_high, lim_dtheta_low, lim_dtheta_high)

            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_mass_H_matched)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_mass_Z_matched)

            h_signal_jetChargeHistos_1D_hist_list.append(h_jetE_reco_over_parton_H_matched_orig)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jetE_reco_over_parton_H_matched_corr)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jetP_reco_over_parton_H_matched_orig)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jetP_reco_over_parton_H_matched_corr)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jetE_reco_over_parton_Z_matched_orig)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jetE_reco_over_parton_Z_matched_corr)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jetP_reco_over_parton_Z_matched_orig)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jetP_reco_over_parton_Z_matched_corr)

            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_b_cM_jetChargeE_kappa_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_bbar_cM_jetChargeE_kappa_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_b_dM_jetChargeE_kappa_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_bbar_dM_jetChargeE_kappa_0_30)
            
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30)
            
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_b_cM_jetChargeE_kappa_0_30_dec_also_close)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_bbar_cM_jetChargeE_kappa_0_30_dec_also_close)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close)
            
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_oppCharge)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_oppCharge)
            
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_chFrac)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_chFrac)     
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_nTrack)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_nTrack)

            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_chFrac)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_chFrac)            
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_nTrack)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_nTrack)

            
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_E_sj2_over_E_sum_dec_also_close)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_E_sj2_over_E_sum_jetChargeE_kappa_0_30_dec_also_close_rightCharges)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_E_sj2_over_E_sum_jetChargeE_kappa_0_30_dec_also_close_wrongCharges)     
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_dAlpha_sj1sj2_dec_also_close)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_dAlpha_sj1sj2_jetChargeE_kappa_0_30_dec_also_close_rightCharges)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_dAlpha_sj1sj2_jetChargeE_kappa_0_30_dec_also_close_wrongCharges)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qneg)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg)

            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qpos)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qpos)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qneg)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg)

            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qpos)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qpos)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20)
            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20)


            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_theta_Z_qneg_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_theta1_Z_qpos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_costheta1_Z_qpos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_sgncos2theta1_costheta1_Z_qpos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_theta_Z_qpos_Z_qneg_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_cosTheta_Z_qneg_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_cosTheta_Z_qpos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_polarAngle_Z_qneg_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_polarAngle_Z_qpos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_cosPolarAngle_Z_qneg_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_cosPolarAngle_Z_qpos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_theta2_H_ep_HZ_COM)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_phi_plane_Z_qpos_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_theta_H_Z_E_totCOM)
            h_signal_jetChargeHistos_1D_hist_list.append(h_parton_theta2_H_ep_approx_HZ_COM)


            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_gj2com)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_gj2com)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_sgncos2theta1_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_sgncos2theta1_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_dm_0_20_gj2com)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_dm_0_20_gj2com)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_dm_0_30_gj2com)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_dm_0_30_gj2com)

            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_E1_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_gj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)


            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_E1_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_gj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_E1_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_gj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_jet2_sj_qneg_dM_ChargeHemisphere_kappa_0_30_dm_jetCharge_0_20)

            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_theta_gj1_gj2_E_totCOM)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_theta2_gj1_ep_E_totCOM)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_theta2_H_ep_E_totCOM)

            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep)

            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_sgnsinphi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_sgnsinphi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_sgnsin2phi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_sgnsin2phi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_sgncosphi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_sgncosphi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_sgncos2phi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_sgncos2phi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep)

            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_dm_0_20_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_dm_0_20_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_dm_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_dm_0_30_vs_plane_gj1_ep)

            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_dm_0_20_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_dm_0_20_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_dm_0_30_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_dm_0_30_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_genjet_theta_sj1_sj2_gj2__gj2COM)

            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_jcE_0_30_theta1_Z_q_pos_Zcom_oppJC)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_jcPt_0_30_theta1_Z_q_pos_Zcom_oppJC)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep_oppJC)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep_oppJC)

            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_partCharge_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_theta2_H_ep_E_totCOM)
            
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_dm_0_20_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_dm_0_20_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_dm_0_30_vs_plane_H_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_dm_0_30_vs_plane_H_ep)
            
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_jcE_0_30_theta1_Z_q_pos_Zcom_oppJC)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_jcPt_0_30_theta1_Z_q_pos_Zcom_oppJC)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep_oppJC)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep_oppJC)

            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_jcE_0_30_theta_gj_jcE_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_jcPt_0_30_theta_gj_jcPt_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_jcE_0_30_dm_0_20_theta_gj_jcE_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_jcPt_0_30_dm_0_20_theta_gj_jcPt_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_jcE_0_30_dm_0_30_theta_gj_jcE_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_jcPt_0_30_dm_0_30_theta_gj_jcPt_0_30)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_phi_gj1_ep_E_totCOM)
            
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_dm_0_20_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_dm_0_20_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_dm_0_30_vs_plane_gj1_ep)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_dm_0_30_vs_plane_gj1_ep)
            
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_jcE_0_30_theta_gj_jcE_0_30_oppJC)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_jcPt_0_30_theta_gj_jcPt_0_30_oppJC)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_vs_plane_gj1_ep_oppJC)
            h_signal_jetChargeHistos_1D_hist_list.append(h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_vs_plane_gj1_ep_oppJC)

            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom)

            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_nTrack_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_nTrack_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom)
            
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_chFrac_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_chFrac_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_chEnergy_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_Product_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_chEnergy_jcPt_0_30_theta1_Z_q_pos_Zcom)
            h_signal_jetChargeHistos_1D_hist_list.append(h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom)


            for hist in h_signal_jetChargeHistos_1D_hist_list:
                hist.Sumw2()
 
            h_2D_genjet_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com_vs_costheta2_gj1_ep_E_totCOM = TH2F( "h_2D_genjet_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com_vs_costheta2_gj1_ep_E_totCOM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
            h_2D_genjet_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com_vs_costheta2_gj1_ep_E_totCOM = TH2F( "h_2D_genjet_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com_vs_costheta2_gj1_ep_E_totCOM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
            h_2D_genjet_sgncostheta1_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com_vs_sgncostheta2_costheta2_gj1_ep_E_totCOM = TH2F( "h_2D_genjet_sgncostheta1_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com_vs_sgncostheta2_costheta2_gj1_ep_E_totCOM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
            h_2D_genjet_sgncostheta1_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com_vs_sgncostheta2_costheta2_gj1_ep_E_totCOM = TH2F( "h_2D_genjet_sgncostheta1_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com_vs_sgncostheta2_costheta2_gj1_ep_E_totCOM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
            h_2D_genjet_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com_vs_costheta2_gj1_ep_E_totCOM.Sumw2()
            h_2D_genjet_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com_vs_costheta2_gj1_ep_E_totCOM.Sumw2()
            h_2D_genjet_sgncostheta1_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com_vs_sgncostheta2_costheta2_gj1_ep_E_totCOM.Sumw2()
            h_2D_genjet_sgncostheta1_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com_vs_sgncostheta2_costheta2_gj1_ep_E_totCOM.Sumw2()

            h_2D_parton_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM = TH2F( "h_2D_parton_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
            h_2D_parton_sgncostheta1_costheta1_costheta1_Z_qpos_Zcom_vs_sgncostheta2_costheta2_H_ep_HZ_COM = TH2F( "h_2D_parton_sgncostheta1_costheta1_Z_qpos_Zcom_vs_sgncostheta2_costheta2_H_ep_HZ_COM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
            h_2D_parton_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Sumw2()
            h_2D_parton_sgncostheta1_costheta1_costheta1_Z_qpos_Zcom_vs_sgncostheta2_costheta2_H_ep_HZ_COM.Sumw2()


        num_entry=-1
        #BDTCut=0.50
        for ientry in range(tree.GetEntries()):
            num_entry+=1;
            if num_entry%(int(tree.GetEntries()/5.)) == 0:
                print "sig BG in entry ",num_entry,bdt_value
            #print "sig BG in entry ",num_entry,bdt_value             
            tree.GetEntry(ientry)
            mvaValue=reader.EvaluateMVA("BDT")
            #print 'input variables ',v_jet1_mass,v_jet2_mass,v_jet1_theta,v_jet2_theta,v_deltatheta,v_jet1_D2_beta1,v_jet2_D2_beta1,v_jet1_BTag_rfj_BTagMax,v_eventWeight,mvaValue
            #print 'input variables ',mvaValue
            h_signal_background_1D_hist_list[0].Fill(mvaValue,v_eventWeight[0]);
            if(mvaValue>bdt_value):
                h_BDT_output_Eff.Fill(mvaValue,v_eventWeight[0])
                h_signal_background_1D_hist_list[1].Fill(v_jet1_mass[0],v_eventWeight[0])
                h_signal_background_1D_hist_list[2].Fill(v_jet2_mass[0],v_eventWeight[0])
                h_jet1_theta.Fill(v_jet1_theta[0],v_eventWeight[0])
                h_jet2_theta.Fill(v_jet2_theta[0],v_eventWeight[0])
                h_jet1_min_jet2_theta.Fill(v_jet1_theta[0]-v_jet2_theta[0],v_eventWeight[0])
                h_jet1_BTag.Fill(v_jet1_BTag_rfj_BTagMax[0],v_eventWeight[0])
                h_jet1_D2.Fill(v_jet1_D2_beta1[0],v_eventWeight[0])
                h_jet2_D2.Fill(v_jet2_D2_beta1[0],v_eventWeight[0])
                h_jet1_C2.Fill(v_jet1_C2_beta1[0],v_eventWeight[0])
                h_jet2_C2.Fill(v_jet2_C2_beta1[0],v_eventWeight[0])
                h_jet1_tau21.Fill(v_jet1_tau21[0],v_eventWeight[0])
                h_jet2_tau21.Fill(v_jet2_tau21[0],v_eventWeight[0])

                obj_jet1= RecoJet(0,0,0,0,v_jet1_tau21[0],v_jet1_C2_beta1[0],v_jet1_D2_beta1[0],v_jet1_BTag_rfj_BTagMax[0],0,0)

                jet_p=sqrt(v_jet1_E[0]*v_jet1_E[0]-v_jet1_mass[0]*v_jet1_mass[0])
                obj_jet1.SetPtEtaPhiM(v_jet1_Pt[0],-log(tan(radians(v_jet1_theta[0]/2.))),radians(v_jet1_phi[0]),v_jet1_mass[0])
                obj_jet2= RecoJet(0,0,0,0,v_jet2_tau21[0],v_jet2_C2_beta1[0],v_jet2_D2_beta1[0],0,0,0)
                jet_p=sqrt(v_jet2_E[0]*v_jet2_E[0]-v_jet2_mass[0]*v_jet2_mass[0])
                obj_jet2.SetPtEtaPhiM(v_jet2_Pt[0],-log(tan(radians(v_jet2_theta[0]/2.))),radians(v_jet2_phi[0]),v_jet2_mass[0])




                if(abs(degrees(obj_jet1.Phi())-v_jet1_phi[0])>0.001 or abs(v_jet1_theta[0]-degrees(obj_jet1.Theta()))>0.001 or (abs(v_jet1_E[0]-obj_jet1.E())/obj_jet1.E())>0.001 or (abs(v_jet1_mass[0]-obj_jet1.M())/obj_jet1.M())>0.001):
                    print 'rather large differences',abs(degrees(obj_jet1.Phi())-v_jet1_phi[0]),abs(v_jet1_theta[0]-degrees(obj_jet1.Theta())),(abs(v_jet1_E[0]-obj_jet1.E())/obj_jet1.E()),(abs(v_jet1_mass[0]-obj_jet1.M())/obj_jet1.M())
                    print 'jettry phi/theta/mass/E ',degrees(obj_jet1.Phi()),degrees(obj_jet1.Theta()),obj_jet1.M(),obj_jet1.E()
                    print 'jet single phi/theta/mass/E ',v_jet1_phi[0],v_jet1_theta[0],v_jet1_mass[0],v_jet1_E[0]
                #elif num_entry%(int(tree.GetEntries()/5.)) == 0:
                #    print 'things are alright',abs(degrees(obj_jet1.Phi())-v_jet1_phi[0]),abs(v_jet1_theta[0]-degrees(obj_jet1.Theta())),(abs(v_jet1_E[0]-obj_jet1.E())/obj_jet1.E()),(abs(v_jet1_mass[0]-obj_jet1.M())/obj_jet1.M())
                #    print 'jettry phi/theta/mass/E ',degrees(obj_jet1.Phi()),degrees(obj_jet1.Theta()),obj_jet1.M(),obj_jet1.E()
                #    print 'jet single phi/theta/mass/E ',v_jet1_phi[0],v_jet1_theta[0],v_jet1_mass[0],v_jet1_E[0]
                #    print 'other variables',obj_jet1.tau21, obj_jet1.C2_beta1,obj_jet1.D2_beta1,obj_jet1.BTag_rfj_BTagMax,obj_jet1.subjetE_ratio,obj_jet1.subjet_DeltaPhi
                #    print 'compared to',v_jet1_tau21[0],v_jet1_C2_beta1[0],v_jet1_D2_beta1[0],v_jet1_BTag_rfj_BTagMax[0],0,0

                temp_rj1=TLorentzVector(0,0,0,0)
                temp_rj1.SetPtEtaPhiM(v_jet1_Pt[0],-log(tan(radians(v_jet1_theta[0]/2.))),radians(v_jet1_phi[0]),v_jet1_mass[0])
                temp_rj2=TLorentzVector(0,0,0,0)
                temp_rj2.SetPtEtaPhiM(v_jet2_Pt[0],-log(tan(radians(v_jet2_theta[0]/2.))),radians(v_jet2_phi[0]),v_jet2_mass[0])
                temp_rj2_orig_corr=TLorentzVector(0,0,0,0)
                temp_rj2_orig_corr.SetPtEtaPhiM(v_jet2_Pt[0],-log(tan(radians(v_jet2_theta[0]/2.))),radians(v_jet2_phi[0]),v_jet2_mass[0])
                temp_rj2_orig=TLorentzVector(0,0,0,0)
                temp_rj2_orig.SetPxPyPzE(v_jet2_sj1_Px[0]+v_jet2_sj2_Px[0],v_jet2_sj1_Py[0]+v_jet2_sj2_Py[0],v_jet2_sj1_Pz[0]+v_jet2_sj2_Pz[0],v_jet2_sj1_E[0]+v_jet2_sj2_E[0])
                temp_rj1_orig=TLorentzVector(0,0,0,0)
                temp_rj1_orig.SetPxPyPzE(v_jet1_sj1_Px[0]+v_jet1_sj2_Px[0],v_jet1_sj1_Py[0]+v_jet1_sj2_Py[0],v_jet1_sj1_Pz[0]+v_jet1_sj2_Pz[0],v_jet1_sj1_E[0]+v_jet1_sj2_E[0])
                if(i_isSignalData):    
                   tempVecH=TLorentzVector(0,0,0,0)
                   tempVecH.SetPxPyPzE(v_parton_H_Px[0],v_parton_H_Py[0],v_parton_H_Pz[0],v_parton_H_E[0])
                   tempVecZ_qneg=TLorentzVector(0,0,0,0)
                   tempVecZ_qneg.SetPxPyPzE(v_parton_Z_qneg_Px[0],v_parton_Z_qneg_Py[0],v_parton_Z_qneg_Pz[0],v_parton_Z_qneg_E[0])
                   tempVecZ_qpos=TLorentzVector(0,0,0,0)
                   tempVecZ_qpos.SetPxPyPzE(v_parton_Z_qpos_Px[0],v_parton_Z_qpos_Py[0],v_parton_Z_qpos_Pz[0],v_parton_Z_qpos_E[0])
                   tempVecZ=TLorentzVector(0,0,0,0)
                   tempVecZ=tempVecZ_qneg+tempVecZ_qpos
                   if (temp_rj1.Angle(tempVecZ.Vect())-temp_rj1_orig.Angle(tempVecZ.Vect())>1.e-5):
                       print 'angle should be so similarr j1',temp_rj1.Angle(tempVecZ.Vect()),temp_rj1_orig.Angle(tempVecZ.Vect())
                   if (temp_rj2.Angle(tempVecZ.Vect())-temp_rj2_orig.Angle(tempVecZ.Vect())>1.e-5):
                       print 'angle should be so similar rj2',temp_rj2.Angle(tempVecZ.Vect()),temp_rj2_orig.Angle(tempVecZ.Vect())
                   if temp_rj1.Angle(tempVecZ.Vect())<temp_rj1.Angle(tempVecH.Vect()):
                       h_jetE_reco_over_parton_Z_matched_orig.Fill(temp_rj1_orig.E()/tempVecZ.E(),v_eventWeight[0])
                       h_jetE_reco_over_parton_Z_matched_corr.Fill(temp_rj1.E()/tempVecZ.E(),v_eventWeight[0])
                       h_jetE_reco_over_parton_H_matched_orig.Fill(temp_rj2_orig.E()/tempVecH.E(),v_eventWeight[0])
                       h_jetE_reco_over_parton_H_matched_corr.Fill(temp_rj2.E()/tempVecH.E(),v_eventWeight[0])
                       h_jetP_reco_over_parton_Z_matched_orig.Fill(temp_rj1_orig.P()/tempVecZ.P(),v_eventWeight[0])
                       h_jetP_reco_over_parton_Z_matched_corr.Fill(temp_rj1.P()/tempVecZ.P(),v_eventWeight[0])
                       h_jetP_reco_over_parton_H_matched_orig.Fill(temp_rj2_orig.P()/tempVecH.P(),v_eventWeight[0])
                       h_jetP_reco_over_parton_H_matched_corr.Fill(temp_rj2.P()/tempVecH.P(),v_eventWeight[0])
                   else:
                       h_jetE_reco_over_parton_Z_matched_orig.Fill(temp_rj2_orig.E()/tempVecZ.E(),v_eventWeight[0])
                       h_jetE_reco_over_parton_Z_matched_corr.Fill(temp_rj2.E()/tempVecZ.E(),v_eventWeight[0])
                       h_jetE_reco_over_parton_H_matched_orig.Fill(temp_rj1_orig.E()/tempVecH.E(),v_eventWeight[0])
                       h_jetE_reco_over_parton_H_matched_corr.Fill(temp_rj1.E()/tempVecH.E(),v_eventWeight[0])
                       h_jetP_reco_over_parton_Z_matched_orig.Fill(temp_rj2_orig.P()/tempVecZ.P(),v_eventWeight[0])
                       h_jetP_reco_over_parton_Z_matched_corr.Fill(temp_rj2.P()/tempVecZ.P(),v_eventWeight[0])
                       h_jetP_reco_over_parton_H_matched_orig.Fill(temp_rj1_orig.P()/tempVecH.P(),v_eventWeight[0])
                       h_jetP_reco_over_parton_H_matched_corr.Fill(temp_rj1.P()/tempVecH.P(),v_eventWeight[0])
                #reco level
                temp_tot_Event_rj_orig=TLorentzVector(0,0,0,0)
                temp_tot_Event_rj_orig.SetPxPyPzE(v_jet1_sj1_Px[0]+v_jet1_sj2_Px[0]+v_jet2_sj1_Px[0]+v_jet2_sj2_Px[0],v_jet1_sj1_Py[0]+v_jet1_sj2_Py[0]+v_jet2_sj1_Py[0]+v_jet2_sj2_Py[0],v_jet1_sj1_Pz[0]+v_jet1_sj2_Pz[0]+v_jet2_sj1_Pz[0]+v_jet2_sj2_Pz[0],v_jet1_sj1_E[0]+v_jet1_sj2_E[0]+v_jet2_sj1_E[0]+v_jet2_sj2_E[0])
                temp_tot_Event_rj=TLorentzVector(0,0,0,0)
                temp_tot_Event_rj=temp_rj1+temp_rj2
                if(v_jet1_sj2_E[0]>v_jet1_sj1_E[0]):
                    print 'sj1 E<sj2 E', v_jet1_sj1_E[0],v_jet1_sj2_E[0]
                sf_rj1=obj_jet1.Px()/(v_jet1_sj1_Px[0]+v_jet1_sj2_Px[0])
                if abs(sf_rj1-1.)<1.e-6 :
                    sf_rj1=1.
                    #print 'wtf scale factor s1 ',sf_rj1,obj_jet1.Px()/(v_jet1_sj1_Px[0]+v_jet1_sj2_Px[0])
                temp_rj1_sj1_orig=TLorentzVector(0,0,0,0)
                temp_rj1_sj1_orig.SetPxPyPzE(v_jet1_sj1_Px[0],v_jet1_sj1_Py[0],v_jet1_sj1_Pz[0],v_jet1_sj1_E[0])
                temp_rj1_sj2_orig=TLorentzVector(0,0,0,0)
                temp_rj1_sj2_orig.SetPxPyPzE(v_jet1_sj2_Px[0],v_jet1_sj2_Py[0],v_jet1_sj2_Pz[0],v_jet1_sj2_E[0])
                temp_rj1_sj1=TLorentzVector(0,0,0,0)
                if sf_rj1==1.:
                    temp_rj1_sj1=temp_rj1_sj1_orig
                else:
                    temp_rj1_sj1.SetPxPyPzE(sf_rj1*v_jet1_sj1_Px[0],sf_rj1*v_jet1_sj1_Py[0],sf_rj1*v_jet1_sj1_Pz[0],v_jet1_sj1_E[0]+v_jet1_sj1_E[0]/(v_jet1_sj1_E[0]+v_jet1_sj2_E[0])*(obj_jet1.E()-(v_jet1_sj1_E[0]+v_jet1_sj2_E[0])))
                temp_rj1_sj2=TLorentzVector(0,0,0,0)
                if sf_rj1==1.:
                    temp_rj1_sj2=temp_rj1_sj2_orig
                else:
                   temp_rj1_sj2.SetPxPyPzE(sf_rj1*v_jet1_sj2_Px[0],sf_rj1*v_jet1_sj2_Py[0],sf_rj1*v_jet1_sj2_Pz[0],v_jet1_sj2_E[0]+v_jet1_sj2_E[0]/(v_jet1_sj1_E[0]+v_jet1_sj2_E[0])*(obj_jet1.E()-(v_jet1_sj1_E[0]+v_jet1_sj2_E[0])))
                sf_rj2=obj_jet2.Px()/(v_jet2_sj1_Px[0]+v_jet2_sj2_Px[0])
                if abs(sf_rj2-1.)<1.e-6 :
                    sf_rj2=1.
                temp_rj2_sj1_orig=TLorentzVector(0,0,0,0)
                temp_rj2_sj1_orig.SetPxPyPzE(v_jet2_sj1_Px[0],v_jet2_sj1_Py[0],v_jet2_sj1_Pz[0],v_jet2_sj1_E[0])
                temp_rj2_sj2_orig=TLorentzVector(0,0,0,0)
                temp_rj2_sj2_orig.SetPxPyPzE(v_jet2_sj2_Px[0],v_jet2_sj2_Py[0],v_jet2_sj2_Pz[0],v_jet2_sj2_E[0])
                temp_rj2_sj1=TLorentzVector(0,0,0,0)
                if(sf_rj2==1):
                    temp_rj2_sj1=temp_rj2_sj1_orig
                else:
                    temp_rj2_sj1.SetPxPyPzE(sf_rj2*v_jet2_sj1_Px[0],sf_rj2*v_jet2_sj1_Py[0],sf_rj2*v_jet2_sj1_Pz[0],v_jet2_sj1_E[0]+v_jet2_sj1_E[0]/(v_jet2_sj1_E[0]+v_jet2_sj2_E[0])*(obj_jet2.E()-(v_jet2_sj1_E[0]+v_jet2_sj2_E[0])))
                temp_rj2_sj2=TLorentzVector(0,0,0,0)
                if(sf_rj2==1):
                    temp_rj2_sj2=temp_rj2_sj2_orig
                else:
                    temp_rj2_sj2.SetPxPyPzE(sf_rj2*v_jet2_sj2_Px[0],sf_rj2*v_jet2_sj2_Py[0],sf_rj2*v_jet2_sj2_Pz[0],v_jet2_sj2_E[0]+v_jet2_sj2_E[0]/(v_jet2_sj1_E[0]+v_jet2_sj2_E[0])*(obj_jet2.E()-(v_jet2_sj1_E[0]+v_jet2_sj2_E[0])))
                temp_ep_approx_rj=TLorentzVector(0,0,0,0)
                temp_ep_approx_E_rj=sqrt(pow(0.5*temp_tot_Event_rj.E()-0.5*temp_tot_Event_rj.Pz(),2)+pow(root.TDatabasePDG.Instance().GetParticle(11).Mass(),2))
                temp_ep_approx_rj.SetPxPyPzE(0.,0.,-0.5*temp_tot_Event_rj.E()+0.5*temp_tot_Event_rj.Pz(),temp_ep_approx_E_rj)
                boostE_tot_COM_rj=TVector3(0,0,0) 
                boostE_tot_COM_rj=-temp_tot_Event_rj.BoostVector()
                temp_ep_approx_rj.Boost(boostE_tot_COM_rj)
                if temp_ep_approx_rj.Pz()>0 :
                    print 'pz should have been negative per principle ',0.5*temp_tot_Event_rj.E()-0.5*temp_tot_Event_rj.Pz(),temp_ep_approx_rj.Pz(),temp_tot_Event_rj.Pz(),temp_tot_Event_rj.Px(),temp_tot_Event_rj.Py()
                if((temp_rj1_sj1+temp_rj1_sj2).Angle(obj_jet1.Vect())<0.01 and abs((temp_rj1_sj1+temp_rj1_sj2).E()-obj_jet1.E())>0.01) :
                    print 'energy of new calc/rescaled rj1 big',  degrees((temp_rj1_sj1+temp_rj1_sj2).Angle(obj_jet1.Vect())),obj_jet1.E(),(temp_rj1_sj1+temp_rj1_sj2).E()

                if((temp_rj1_sj1+temp_rj1_sj2).Angle(obj_jet1.Vect()))>0.01 :
                    print 'delta angle new calc/rescaled rj1 big',  degrees((temp_rj1_sj1+temp_rj1_sj2).Angle(obj_jet1.Vect())),obj_jet1.E(),(temp_rj1_sj1+temp_rj1_sj2).E()
                #if temp_rj1_sj1.M()<=0 or temp_rj1_sj2.M()<=0 :
                #    print 'wtf with masses rj1 ',temp_rj1_sj1.M(),temp_rj1_sj1.P(),temp_rj1_sj1.E(),temp_rj1_sj2.M(),temp_rj1_sj2.P(),temp_rj1_sj2.E(),sf_rj1,temp_tot_Event_rj.M()/temp_tot_Event_rj_orig.M()
                if((temp_rj2_sj1+temp_rj2_sj2).Angle(obj_jet2.Vect())<0.01 and abs((temp_rj2_sj1+temp_rj2_sj2).E()-obj_jet2.E())>0.01) :
                    print 'energy of new calc/rescaled rj2 big',  degrees((temp_rj2_sj1+temp_rj2_sj2).Angle(obj_jet2.Vect())),obj_jet2.E(),(temp_rj2_sj1+temp_rj2_sj2).E()
                if((temp_rj2_sj1+temp_rj2_sj2).Angle(obj_jet2.Vect()))>0.01 :
                    print 'delta angle new calc/rescaled rj2 big',  degrees((temp_rj2_sj1+temp_rj2_sj2).Angle(obj_jet2.Vect())),obj_jet2.E(),(temp_rj2_sj1+temp_rj2_sj2).E()
                #defined after rescaling with MET projection above
                temp_rj1.Boost(boostE_tot_COM_rj)
                plane_rj1_ep_approx_rj_boost=TVector3(0,0,0)
                plane_rj1_ep_approx_rj_boost=temp_rj1.Vect().Cross(temp_ep_approx_rj.Vect()).Unit()
                #if temp_rj2_sj1.M()<=0 or temp_rj2_sj2.M()<=0 :
                #    print 'wtf with masses rj2 ',temp_rj2_sj1.M(),temp_rj2_sj1.P(),temp_rj2_sj1.E(),temp_rj2_sj2.M(),temp_rj2_sj2.P(),temp_rj2_sj2.E(),sf_rj2,temp_tot_Event_rj.M()/temp_tot_Event_rj_orig.M(),temp_tot_Event_rj.M()
                temp_rj2.Boost(boostE_tot_COM_rj)
                boostrj2_COM=TVector3(0,0,0) 
                boostrj2_COM=-temp_rj2_orig_corr.BoostVector()
                #boost subjets into rj2 COM
                temp_rj2_sj1.Boost(boostrj2_COM)
                temp_rj2_sj2.Boost(boostrj2_COM)

                h_recojet_theta_rj1_rj2_E_totCOM.Fill(degrees(temp_rj1.Angle(temp_rj2.Vect())),v_eventWeight[0])
                h_recojet_theta_sj1_sj2_rj2__rj2COM.Fill(degrees(temp_rj2_sj1.Angle(temp_rj2_sj2.Vect())),v_eventWeight[0])
                h_recojet_theta2_rj1_ep_E_totCOM.Fill(degrees(temp_rj1.Angle(temp_ep_approx_rj.Vect())),v_eventWeight[0])
                
                #print 'ep approx part',temp_ep_approx.Px(),temp_ep_approx.Py(),temp_ep_approx.Pz(),temp_ep_approx.E()
                #print 'ep approx rj',temp_ep_approx_rj.Px(),temp_ep_approx_rj.Py(),temp_ep_approx_rj.Pz(),temp_ep_approx_rj.E()
                ind_rj2_sj_pos_jcE=-1
                if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0])> abs(v_jet2_sj2_jetChargeE_kappa_0_30[0]):
                    if v_jet2_sj1_jetChargeE_kappa_0_30[0]>0:
                        ind_rj2_sj_pos_jcE=1
                    else:
                        ind_rj2_sj_pos_jcE=2
                else:
                    if v_jet2_sj2_jetChargeE_kappa_0_30[0]>0:
                        ind_rj2_sj_pos_jcE=2
                    else:
                        ind_rj2_sj_pos_jcE=1
                if ind_rj2_sj_pos_jcE==-1:
                    print 'sth wrong in jcE rj subjet index'
                #subjets are already boosted at that point
                temp_rj2_sj_jcE_pos=TLorentzVector(0,0,0,0)
                if ind_rj2_sj_pos_jcE==1:
                    temp_rj2_sj_jcE_pos=temp_rj2_sj1
                else:
                    temp_rj2_sj_jcE_pos=temp_rj2_sj2
          

                #print 'ep approx part',temp_ep_approx.Px(),temp_ep_approx.Py(),temp_ep_approx.Pz(),temp_ep_approx.E()
                #print 'ep approx rj',temp_ep_approx_rj.Px(),temp_ep_approx_rj.Py(),temp_ep_approx_rj.Pz(),temp_ep_approx_rj.E()
                ind_rj2_sj_pos_jcPt=-1
                if abs(v_jet2_sj1_jetChargePt_kappa_0_30[0])> abs(v_jet2_sj2_jetChargePt_kappa_0_30[0]):
                    if v_jet2_sj1_jetChargePt_kappa_0_30[0]>0:
                        ind_rj2_sj_pos_jcPt=1
                    else:
                        ind_rj2_sj_pos_jcPt=2
                else:
                    if v_jet2_sj2_jetChargePt_kappa_0_30[0]>0:
                        ind_rj2_sj_pos_jcPt=2
                    else:
                        ind_rj2_sj_pos_jcPt=1
                if ind_rj2_sj_pos_jcPt==-1:
                    print 'sth wrong in jcPt rj subjet index'
                    
                temp_rj2_sj_jcPt_pos=TLorentzVector(0,0,0,0)
                if ind_rj2_sj_pos_jcPt==1:
                    temp_rj2_sj_jcPt_pos=temp_rj2_sj1
                    #print "jcPt 1 is pos 1/2 ",v_jet2_sj1_jetChargePt_kappa_0_30[0],v_jet2_sj2_jetChargePt_kappa_0_30[0]
                else:
                    temp_rj2_sj_jcPt_pos=temp_rj2_sj2
                    #print "jcPt 2 is pos 2/1 ",v_jet2_sj2_jetChargePt_kappa_0_30[0],v_jet2_sj1_jetChargePt_kappa_0_30[0]

                #subjet indices ordered by energy
                ind_rj2_sj_pos_jcE_sj_E1=-1  
                temp_rj2_sj_jcE_pos_sj_E1=TLorentzVector(0,0,0,0)
                if v_jet2_sj1_jetChargeE_kappa_0_30[0]>0 :
                    ind_rj2_sj_pos_jcE_sj_E1=1
                    temp_rj2_sj_jcE_pos_sj_E1=temp_rj2_sj1
                else:
                    ind_rj2_sj_pos_jcE_sj_E1=2
                    temp_rj2_sj_jcE_pos_sj_E1=temp_rj2_sj2
                if ind_rj2_sj_pos_jcE_sj_E1==-1:
                    print 'sth wrong in jcE rj subjet index, case E1'

                #subjet indices ordered by energy
                ind_rj2_sj_pos_jcPt_sj_E1=-1  
                temp_rj2_sj_jcPt_pos_sj_E1=TLorentzVector(0,0,0,0)
                if v_jet2_sj1_jetChargeE_kappa_0_30[0]>0 :
                    ind_rj2_sj_pos_jcPt_sj_E1=1
                    temp_rj2_sj_jcPt_pos_sj_E1=temp_rj2_sj1
                else:
                    ind_rj2_sj_pos_jcPt_sj_E1=2
                    temp_rj2_sj_jcPt_pos_sj_E1=temp_rj2_sj2
                if ind_rj2_sj_pos_jcPt_sj_E1==-1:
                    print 'sth wrong in jcPt rj subjet index, case E1'

                #subjet indices ordered by length of tracks
                ind_rj2_sj_pos_jcE_sj_nTracksMax=-1  
                if v_jet2_sj1_nTracks[0]> v_jet2_sj2_nTracks[0] :
                    if v_jet2_sj1_jetChargeE_kappa_0_30[0]>0:
                        ind_rj2_sj_pos_jcE_sj_nTracksMax=1
                    else:
                        ind_rj2_sj_pos_jcE_sj_nTracksMax=2
                elif v_jet2_sj2_nTracks[0]> v_jet2_sj1_nTracks[0] :
                    if v_jet2_sj2_jetChargeE_kappa_0_30[0]>0:
                        ind_rj2_sj_pos_jcE_sj_nTracksMax=2
                    else:
                        ind_rj2_sj_pos_jcE_sj_nTracksMax=1
                else:
                    if v_jet2_sj1_jetChargeE_kappa_0_30[0]>0 :
                        ind_rj2_sj_pos_jcE_sj_nTracksMax=1
                    else:
                        ind_rj2_sj_pos_jcE_sj_nTracksMax=2
                if ind_rj2_sj_pos_jcE_sj_nTracksMax==-1:
                    print 'sth wrong in jcE rj subjet index, case nTracksMax'
                temp_rj2_sj_jcE_pos_sj_nTracksMax=TLorentzVector(0,0,0,0)
                if ind_rj2_sj_pos_jcE_sj_nTracksMax==1:
                    temp_rj2_sj_jcE_pos_sj_nTracksMax=temp_rj2_sj1
                else:
                    temp_rj2_sj_jcE_pos_sj_nTracksMax=temp_rj2_sj2

                #subjet indices ordered by length of tracks
                ind_rj2_sj_pos_jcPt_sj_nTracksMax=-1  
                if v_jet2_sj1_nTracks[0]> v_jet2_sj2_nTracks[0] :
                    if v_jet2_sj1_jetChargePt_kappa_0_30[0]>0:
                        ind_rj2_sj_pos_jcPt_sj_nTracksMax=1
                    else:
                        ind_rj2_sj_pos_jcPt_sj_nTracksMax=2
                elif v_jet2_sj2_nTracks[0]> v_jet2_sj1_nTracks[0] :
                    if v_jet2_sj2_jetChargeE_kappa_0_30[0]>0:
                        ind_rj2_sj_pos_jcPt_sj_nTracksMax=2
                    else:
                        ind_rj2_sj_pos_jcPt_sj_nTracksMax=1
                else:
                    if v_jet2_sj1_jetChargePt_kappa_0_30[0]>0 :
                        ind_rj2_sj_pos_jcPt_sj_nTracksMax=1
                    else:
                        ind_rj2_sj_pos_jcPt_sj_nTracksMax=2
                if ind_rj2_sj_pos_jcPt_sj_nTracksMax==-1:
                    print 'sth wrong in jcPt rj subjet index, case nTracksMax'
                temp_rj2_sj_jcPt_pos_sj_nTracksMax=TLorentzVector(0,0,0,0)
                if ind_rj2_sj_pos_jcPt_sj_nTracksMax==1:
                    temp_rj2_sj_jcPt_pos_sj_nTracksMax=temp_rj2_sj1
                else:
                    temp_rj2_sj_jcPt_pos_sj_nTracksMax=temp_rj2_sj2

                #subjet indices ordered by chargedFraction
                ind_rj2_sj_chFrac=-1  
                if v_jet2_sj1_chFrac[0]>v_jet2_sj2_chFrac[0]:
                    ind_rj2_sj_chFrac=1
                else:
                    ind_rj2_sj_chFrac=2
                if ind_rj2_sj_chFrac==-1:
                    print 'sth wrong in rj subjet index, case chFrac'
                temp_rj2_sj_jcE_pos_sj_chFrac=TLorentzVector(0,0,0,0)
                temp_rj2_sj_jcPt_pos_sj_chFrac=TLorentzVector(0,0,0,0)
                #larger charged fraction in subjet 1
                if ind_rj2_sj_chFrac==1:
                    if v_jet2_sj1_jetChargeE_kappa_0_30[0]>0 :
                       temp_rj2_sj_jcE_pos_sj_chFrac =temp_rj2_sj1
                    else:
                        temp_rj2_sj_jcE_pos_sj_chFrac=temp_rj2_sj2
                    if v_jet2_sj1_jetChargePt_kappa_0_30[0]>0 :
                       temp_rj2_sj_jcPt_pos_sj_chFrac =temp_rj2_sj1
                    else:
                        temp_rj2_sj_jcPt_pos_sj_chFrac=temp_rj2_sj2
                #larger charged fraction in subjet 2
                else:
                    if v_jet2_sj2_jetChargeE_kappa_0_30[0]>0 :
                       temp_rj2_sj_jcE_pos_sj_chFrac =temp_rj2_sj2
                    else:
                        temp_rj2_sj_jcE_pos_sj_chFrac=temp_rj2_sj1
                    if v_jet2_sj2_jetChargePt_kappa_0_30[0]>0 :
                       temp_rj2_sj_jcPt_pos_sj_chFrac =temp_rj2_sj2
                    else:
                        temp_rj2_sj_jcPt_pos_sj_chFrac=temp_rj2_sj1
 

                #subjet indices ordered by chargedEnergy = chFrac*subjetEnergy (even after Neutrino correction, charged energy REMAINS the same)
                ind_rj2_sj_chEnergy=-1  
                if (v_jet2_sj1_chFrac[0]*v_jet2_sj1_E[0])>(v_jet2_sj2_chFrac[0]*v_jet2_sj2_E[0]):
                    ind_rj2_sj_chEnergy=1
                else:
                    ind_rj2_sj_chEnergy=2
                if ind_rj2_sj_chEnergy==-1:
                    print 'sth wrong in rj subjet index, case chEnergy'
                temp_rj2_sj_jcE_pos_sj_chEnergy=TLorentzVector(0,0,0,0)
                temp_rj2_sj_jcPt_pos_sj_chEnergy=TLorentzVector(0,0,0,0)
                #larger charged fraction in subjet 1
                if ind_rj2_sj_chEnergy==1:
                    if v_jet2_sj1_jetChargeE_kappa_0_30[0]>0 :
                       temp_rj2_sj_jcE_pos_sj_chEnergy =temp_rj2_sj1
                    else:
                        temp_rj2_sj_jcE_pos_sj_chEnergy=temp_rj2_sj2
                    if v_jet2_sj1_jetChargePt_kappa_0_30[0]>0 :
                       temp_rj2_sj_jcPt_pos_sj_chEnergy =temp_rj2_sj1
                    else:
                        temp_rj2_sj_jcPt_pos_sj_chEnergy=temp_rj2_sj2
                #larger charged fraction in subjet 2
                else:
                    if v_jet2_sj2_jetChargeE_kappa_0_30[0]>0 :
                       temp_rj2_sj_jcE_pos_sj_chEnergy =temp_rj2_sj2
                    else:
                        temp_rj2_sj_jcE_pos_sj_chEnergy=temp_rj2_sj1
                    if v_jet2_sj2_jetChargePt_kappa_0_30[0]>0 :
                       temp_rj2_sj_jcPt_pos_sj_chEnergy =temp_rj2_sj2
                    else:
                        temp_rj2_sj_jcPt_pos_sj_chEnergy=temp_rj2_sj1


                plane_rj2_rj2_sj_jcPt_pos_boost=TVector3(0,0,0)
                plane_rj2_rj2_sj_jcPt_pos=temp_rj2.Vect().Cross(temp_rj2_sj_jcPt_pos.Vect()).Unit()
              

                plane_rj2_rj2_sj_jcE_pos_boost=TVector3(0,0,0)
                plane_rj2_rj2_sj_jcE_pos=temp_rj2.Vect().Cross(temp_rj2_sj_jcE_pos.Vect()).Unit()

                h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect())),v_eventWeight[0])
                phi_plane_rj1_ep_rj2_sj_jcE_pos=degrees(acos(plane_rj1_ep_approx_rj_boost.Dot(plane_rj2_rj2_sj_jcE_pos)))
                #plane is already a TVector3
                if(degrees(temp_rj2_sj_jcE_pos.Angle(plane_rj1_ep_approx_rj_boost))>90.):
                    phi_plane_rj1_ep_rj2_sj_jcE_pos=degrees(acos(plane_rj1_ep_approx_rj_boost.Dot(plane_rj2_rj2_sj_jcE_pos)))+180.
                h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Fill(phi_plane_rj1_ep_rj2_sj_jcE_pos,v_eventWeight[0])
                if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.20:
                    h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_dm_0_20_rj2com.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect())),v_eventWeight[0])
                    h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_dm_0_20_vs_plane_rj1_ep.Fill(phi_plane_rj1_ep_rj2_sj_jcE_pos,v_eventWeight[0])
                    if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.30:
                        h_recojet_theta1_rj2_pos_sj_jetChargeE_0_30_dm_0_30_rj2com.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect())),v_eventWeight[0])
                        h_recojet_phi_plane_rj2_pos_sj_jetChargeE_0_30_dm_0_30_vs_plane_rj1_ep.Fill(phi_plane_rj1_ep_rj2_sj_jcE_pos,v_eventWeight[0])


                h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_rj2com.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect())),v_eventWeight[0])

                phi_plane_rj1_ep_rj2_sj_jcPt_pos=degrees(acos(plane_rj1_ep_approx_rj_boost.Dot(plane_rj2_rj2_sj_jcPt_pos)))
                #plane is already a TVector3
                if(degrees(temp_rj2_sj_jcPt_pos.Angle(plane_rj1_ep_approx_rj_boost))>90.):
                    phi_plane_rj1_ep_rj2_sj_jcPt_pos=degrees(acos(plane_rj1_ep_approx_rj_boost.Dot(plane_rj2_rj2_sj_jcPt_pos)))+180.
                h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep.Fill(phi_plane_rj1_ep_rj2_sj_jcPt_pos,v_eventWeight[0])
                if abs(v_jet2_sj1_jetChargePt_kappa_0_30[0]-v_jet2_sj2_jetChargePt_kappa_0_30[0])>0.20:
                    h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_dm_0_20_rj2com.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect())),v_eventWeight[0])
                    h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_dm_0_20_vs_plane_rj1_ep.Fill(phi_plane_rj1_ep_rj2_sj_jcPt_pos,v_eventWeight[0])
                    if abs(v_jet2_sj1_jetChargePt_kappa_0_30[0]-v_jet2_sj2_jetChargePt_kappa_0_30[0])>0.30:
                        h_recojet_theta1_rj2_pos_sj_jetChargePt_0_30_dm_0_30_rj2com.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect())),v_eventWeight[0])
                        h_recojet_phi_plane_rj2_pos_sj_jetChargePt_0_30_dm_0_30_vs_plane_rj1_ep.Fill(phi_plane_rj1_ep_rj2_sj_jcPt_pos,v_eventWeight[0])

                h_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Fill(cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect())),v_eventWeight[0])
                h_recojet_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com.Fill(cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect())),v_eventWeight[0])
                if cos(2.*temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))!=0:
                    h_recojet_sgncos2theta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com.Fill(cos(2.*temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))/abs(cos(2.*temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect())))*cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect())),v_eventWeight[0])
                if cos(2.*temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))!=0:
                    h_recojet_sgncos2theta1_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com.Fill(cos(2.*temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))/abs(cos(2.*temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect())))*cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect())),v_eventWeight[0])


                h_2D_recojet_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Fill(cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect())),cos(temp_rj1.Angle(temp_ep_approx_rj.Vect())),v_eventWeight[0])
                h_2D_recojet_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com_vs_costheta2_rj1_ep_E_totCOM.Fill(cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect())),cos(temp_rj1.Angle(temp_ep_approx_rj.Vect())),v_eventWeight[0])
                if cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))!=0 and cos(temp_rj1.Angle(temp_ep_approx_rj.Vect()))!=0:
                    h_2D_recojet_sgncostheta1_costheta1_rj2_pos_sj_jetChargeE_0_30_rj2com_vs_sgncostheta2_costheta2_rj1_ep_E_totCOM.Fill(cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))/abs(cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect())))*cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect())),cos(temp_rj1.Angle(temp_ep_approx_rj.Vect()))/abs(cos(temp_rj1.Angle(temp_ep_approx_rj.Vect())))*cos(temp_rj1.Angle(temp_ep_approx_rj.Vect())),v_eventWeight[0])
                if cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))!=0 and cos(temp_rj1.Angle(temp_ep_approx_rj.Vect()))!=0:
                    h_2D_recojet_sgncostheta1_costheta1_rj2_pos_sj_jetChargePt_0_30_rj2com_vs_sgncostheta2_costheta2_rj1_ep_E_totCOM.Fill(cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))/abs(cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect())))*cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect())),cos(temp_rj1.Angle(temp_ep_approx_rj.Vect()))/abs(cos(temp_rj1.Angle(temp_ep_approx_rj.Vect())))*cos(temp_rj1.Angle(temp_ep_approx_rj.Vect())),v_eventWeight[0])

                #print 'degrees phi plane ',degrees(phi_plane_rj1_ep_rj2_sj_jcE_pos)
                phi_plane_rj1_ep_rj2_sj_jcE_pos_rad = radians(phi_plane_rj1_ep_rj2_sj_jcE_pos)
                if sin(phi_plane_rj1_ep_rj2_sj_jcE_pos_rad)!=0:
                    h_recojet_sgnsinphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Fill(sin(phi_plane_rj1_ep_rj2_sj_jcE_pos_rad)/abs(sin(phi_plane_rj1_ep_rj2_sj_jcE_pos_rad))*phi_plane_rj1_ep_rj2_sj_jcE_pos,v_eventWeight[0])
                #print 'degrees phi plane,sin phi, sign sin phi/fill ',phi_plane_rj1_ep_rj2_sj_jcE_pos,sin(phi_plane_rj1_ep_rj2_sj_jcE_pos),sin(phi_plane_rj1_ep_rj2_sj_jcE_pos)/abs(sin(phi_plane_rj1_ep_rj2_sj_jcE_pos)),sin(phi_plane_rj1_ep_rj2_sj_jcE_pos)/abs(sin(phi_plane_rj1_ep_rj2_sj_jcE_pos))*phi_plane_rj1_ep_rj2_sj_jcE_pos
                if sin(2.*phi_plane_rj1_ep_rj2_sj_jcE_pos_rad)!=0:
                    h_recojet_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Fill(sin(2.*phi_plane_rj1_ep_rj2_sj_jcE_pos_rad)/abs(sin(2.*phi_plane_rj1_ep_rj2_sj_jcE_pos_rad))*phi_plane_rj1_ep_rj2_sj_jcE_pos,v_eventWeight[0])
                if cos(phi_plane_rj1_ep_rj2_sj_jcE_pos_rad)!=0:
                    h_recojet_sgncosphi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Fill(cos(phi_plane_rj1_ep_rj2_sj_jcE_pos_rad)/abs(cos(phi_plane_rj1_ep_rj2_sj_jcE_pos_rad))*phi_plane_rj1_ep_rj2_sj_jcE_pos,v_eventWeight[0])
                if cos(2.*phi_plane_rj1_ep_rj2_sj_jcE_pos_rad)!=0:
                    h_recojet_sgncos2phi_phi_plane_rj2_pos_sj_jetChargeE_0_30_vs_plane_rj1_ep.Fill(cos(2.*phi_plane_rj1_ep_rj2_sj_jcE_pos_rad)/abs(cos(2.*phi_plane_rj1_ep_rj2_sj_jcE_pos_rad))*phi_plane_rj1_ep_rj2_sj_jcE_pos,v_eventWeight[0])
                #print 'degrees phi plane ',degrees(phi_plane_rj1_ep_rj2_sj_jcPt_pos)
                phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad = radians(phi_plane_rj1_ep_rj2_sj_jcPt_pos)
                if sin(phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad)!=0:
                    h_recojet_sgnsinphi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep.Fill(sin(phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad)/abs(sin(phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad))*phi_plane_rj1_ep_rj2_sj_jcPt_pos,v_eventWeight[0])
                #print 'degrees phi plane,sin phi, sign sin phi/fill ',phi_plane_rj1_ep_rj2_sj_jcPt_pos,sin(phi_plane_rj1_ep_rj2_sj_jcPt_pos),sin(phi_plane_rj1_ep_rj2_sj_jcPt_pos)/abs(sin(phi_plane_rj1_ep_rj2_sj_jcPt_pos)),sin(phi_plane_rj1_ep_rj2_sj_jcPt_pos)/abs(sin(phi_plane_rj1_ep_rj2_sj_jcPt_pos))*phi_plane_rj1_ep_rj2_sj_jcPt_pos
                if sin(2.*phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad)!=0:
                    h_recojet_sgnsin2phi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep.Fill(sin(2.*phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad)/abs(sin(2.*phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad))*phi_plane_rj1_ep_rj2_sj_jcPt_pos,v_eventWeight[0])
                if cos(phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad)!=0:
                    h_recojet_sgncosphi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep.Fill(cos(phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad)/abs(cos(phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad))*phi_plane_rj1_ep_rj2_sj_jcPt_pos,v_eventWeight[0])
                if cos(2.*phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad)!=0:
                    h_recojet_sgncos2phi_phi_plane_rj2_pos_sj_jetChargePt_0_30_vs_plane_rj1_ep.Fill(cos(2.*phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad)/abs(cos(2.*phi_plane_rj1_ep_rj2_sj_jcPt_pos_rad))*phi_plane_rj1_ep_rj2_sj_jcPt_pos,v_eventWeight[0])
 
                #jet charge histos for matching etc only relevant for signal if len(h_signal_jetChargeHistos_1D_hist_list)>0:
                if len(h_signal_jetChargeHistos_1D_hist_list)>0:
                    temp_tot_Event=TLorentzVector(0,0,0,0)
                    temp_tot_Event.SetPxPyPzE(v_parton_em_Px[0]+v_parton_ep_Px[0],v_parton_em_Py[0]+v_parton_ep_Py[0],v_parton_em_Pz[0]+v_parton_ep_Pz[0],v_parton_em_E[0]+v_parton_ep_E[0])
                    boostE_tot_COM=TVector3(0,0,0) 
                    boostE_tot_COM=-temp_tot_Event.BoostVector()
                    temp_ep_orig=TLorentzVector(0,0,0,0)
                    temp_ep_orig.SetPxPyPzE(v_parton_ep_Px[0],v_parton_ep_Py[0],v_parton_ep_Pz[0],v_parton_ep_E[0])
                    temp_ep=TLorentzVector(0,0,0,0)
                    temp_ep.SetPxPyPzE(v_parton_ep_Px[0],v_parton_ep_Py[0],v_parton_ep_Pz[0],v_parton_ep_E[0])
                    temp_ep.Boost(boostE_tot_COM)
                    temp_em=TLorentzVector(0,0,0,0)
                    temp_em.SetPxPyPzE(v_parton_em_Px[0],v_parton_em_Py[0],v_parton_em_Pz[0],v_parton_em_E[0])
                    temp_em.Boost(boostE_tot_COM)
                    temp_H=TLorentzVector(0,0,0,0)
                    temp_H.SetPxPyPzE(v_parton_H_Px[0],v_parton_H_Py[0],v_parton_H_Pz[0],v_parton_H_E[0])
                    temp_H.Boost(boostE_tot_COM)
                    temp_Z_qneg=TLorentzVector(0,0,0,0)
                    temp_Z_qneg.SetPxPyPzE(v_parton_Z_qneg_Px[0],v_parton_Z_qneg_Py[0],v_parton_Z_qneg_Pz[0],v_parton_Z_qneg_E[0])
                    temp_Z_qpos=TLorentzVector(0,0,0,0)
                    temp_Z_qpos.SetPxPyPzE(v_parton_Z_qpos_Px[0],v_parton_Z_qpos_Py[0],v_parton_Z_qpos_Pz[0],v_parton_Z_qpos_E[0])
                    temp_Z=TLorentzVector(0,0,0,0)
                    temp_Z=temp_Z_qneg+temp_Z_qpos
                    #print 'Z before COM boost ',temp_Z.M(),temp_Z.Px(),temp_Z.Py(),temp_Z.Pz(),temp_Z.E()
                    temp_Z.Boost(boostE_tot_COM)
                    #print 'Z after COM boost ',temp_Z.M(),temp_Z.Px(),temp_Z.Py(),temp_Z.Pz(),temp_Z.E()
                    #print 'event boost',temp_tot_Event.M(),temp_tot_Event.Px(),temp_tot_Event.Py(),temp_tot_Event.Pz(),temp_tot_Event.E()
                    boostZ_COM=TVector3(0,0,0) 
                    boostZ_COM=-temp_Z.BoostVector()
                    temp_Z_qpos_boostZ_COM=TLorentzVector(0,0,0,0)
                    temp_Z_qpos_boostZ_COM.SetPxPyPzE(v_parton_Z_qpos_Px[0],v_parton_Z_qpos_Py[0],v_parton_Z_qpos_Pz[0],v_parton_Z_qpos_E[0])
                    temp_Z_qpos_boostZ_COM.Boost(boostE_tot_COM)
                    temp_Z_qpos_boostZ_COM.Boost(boostZ_COM)
                    temp_Z_qneg_boostZ_COM=TLorentzVector(0,0,0,0)
                    temp_Z_qneg_boostZ_COM.SetPxPyPzE(v_parton_Z_qneg_Px[0],v_parton_Z_qneg_Py[0],v_parton_Z_qneg_Pz[0],v_parton_Z_qneg_E[0])
                    temp_Z_qneg_boostZ_COM.Boost(boostE_tot_COM)
                    temp_Z_qneg_boostZ_COM.Boost(boostZ_COM)
                    temp_Z_sum=TLorentzVector(0,0,0,0)
                    temp_Z_sum=temp_Z_qpos_boostZ_COM+temp_Z_qneg_boostZ_COM
                    temp_Z_COM=TLorentzVector(0,0,0,0)
                    temp_Z_COM.SetPxPyPzE(temp_Z.Px(),temp_Z.Py(),temp_Z.Pz(),temp_Z.E())
                    temp_Z_COM.Boost(boostZ_COM)
                    #print 'Z COM', temp_Z_COM.M(),temp_Z_COM.Px(),temp_Z_COM.Py(),temp_Z_COM.Pz(),temp_Z_COM.E()
                    #print 'Z sum', temp_Z_sum.M(),temp_Z_sum.Px(),temp_Z_sum.Py(),temp_Z_sum.Pz(),temp_Z_sum.E()
                    h_parton_theta2_H_ep_HZ_COM.Fill(degrees(temp_H.Angle(temp_ep.Vect())),v_eventWeight[0])
                    h_parton_theta_Z_qneg_Zcom.Fill(degrees(temp_Z_qneg_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_parton_theta1_Z_qpos_Zcom.Fill(degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_parton_costheta1_Z_qpos_Zcom.Fill(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))*degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    if cos(2.*temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0 :
                        h_parton_sgncos2theta1_costheta1_Z_qpos_Zcom.Fill(cos(2.*temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(2.*temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))*degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_2D_parton_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Fill(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),cos(temp_H.Angle(temp_ep.Vect())),v_eventWeight[0])
                    if cos(temp_H.Angle(temp_ep.Vect()))!=0 and cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0:
                        h_2D_parton_sgncostheta1_costheta1_costheta1_Z_qpos_Zcom_vs_sgncostheta2_costheta2_H_ep_HZ_COM.Fill(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),cos(temp_H.Angle(temp_ep.Vect()))/abs(cos(temp_H.Angle(temp_ep.Vect())))*cos(temp_H.Angle(temp_ep.Vect())),v_eventWeight[0])
 
                    h_parton_theta_Z_qpos_Z_qneg_Zcom.Fill(degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z_qneg_boostZ_COM.Vect())),v_eventWeight[0])
                    h_parton_cosTheta_Z_qneg_Zcom.Fill(cos(temp_Z_qneg_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_parton_cosTheta_Z_qpos_Zcom.Fill(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_parton_polarAngle_Z_qneg_Zcom.Fill(degrees(temp_Z_qneg_boostZ_COM.Theta()),v_eventWeight[0])
                    h_parton_polarAngle_Z_qpos_Zcom.Fill(degrees(temp_Z_qpos_boostZ_COM.Theta()),v_eventWeight[0])
                    h_parton_cosPolarAngle_Z_qneg_Zcom.Fill(cos(temp_Z_qneg_boostZ_COM.Theta()),v_eventWeight[0])
                    h_parton_cosPolarAngle_Z_qpos_Zcom.Fill(cos(temp_Z_qpos_boostZ_COM.Theta()),v_eventWeight[0])
                    plane_Z_Z_qpos_boost=TVector3(0,0,0)
                    plane_Z_Z_qpos_boost=temp_Z.Vect().Cross(temp_Z_qpos_boostZ_COM.Vect()).Unit()
                    #print 'plane vector Z, pos q ',plane_Z_Z_qpos_boost.Px(),plane_Z_Z_qpos_boost.Py(),plane_Z_Z_qpos_boost.Pz()
                    tempZ=TLorentzVector(0,0,0,0)
                    tempZ=temp_Z_qneg+temp_Z_qpos
                    plane_Z_Z_qneg_orig=TVector3(0,0,0)
                    plane_Z_Z_qneg_orig=tempZ.Vect().Cross(temp_Z_qneg.Vect()).Unit()
                    #print 'plane vector Z, neg q orig ',plane_Z_Z_qneg_orig.Px(),plane_Z_Z_qneg_orig.Py(),plane_Z_Z_qneg_orig.Pz()
                    tempH=TLorentzVector(0,0,0,0)
                    tempH.SetPxPyPzE(v_parton_H_Px[0],v_parton_H_Py[0],v_parton_H_Pz[0],v_parton_H_E[0])
                    temp_tot_Event_V2=TLorentzVector(0,0,0,0)
                    temp_tot_Event_V2=tempH+tempZ
                    #print 'tot event A',temp_tot_Event.Px(),temp_tot_Event.Py(),temp_tot_Event.Pz(),temp_tot_Event.E()
                    #print 'tot event B',temp_tot_Event_V2.Px(),temp_tot_Event_V2.Py(),temp_tot_Event_V2.Pz(),temp_tot_Event_V2.E()
                    #print 'tot event C',v_parton_H_Px[0]+v_parton_Z_qneg_Px[0]+v_parton_Z_qpos_Px[0],v_parton_H_Py[0]+v_parton_Z_qneg_Py[0]+v_parton_Z_qpos_Py[0],v_parton_H_Pz[0]+v_parton_Z_qneg_Pz[0]+v_parton_Z_qpos_Pz[0],v_parton_H_E[0]+v_parton_Z_qneg_E[0]+v_parton_Z_qpos_E[0]
                    #print 'ep real ',v_parton_ep_Px[0],v_parton_ep_Py[0],v_parton_ep_Pz[0]
                    #print 'ep calc ',-0.01*(0.5*temp_tot_Event_V2.Pz()-50.*temp_tot_Event_V2.Px()),0,0.5*temp_tot_Event_V2.Pz()-50.*temp_tot_Event_V2.Px()
                    #print 'em real ',v_parton_em_Px[0],v_parton_em_Py[0],v_parton_em_Pz[0],v_parton_em_E[0]
                    #print 'em calc ',0.01*(0.5*temp_tot_Event_V2.Pz()+50.*temp_tot_Event_V2.Px()),0,0.5*temp_tot_Event_V2.Pz()+50.*temp_tot_Event_V2.Px()
                    plane_H_ep_boost=TVector3(0,0,0)
                    plane_H_ep_boost=temp_H.Vect().Cross(temp_ep.Vect()).Unit()
                    #print 'plane vector H,ep boost',plane_H_ep_boost.Px(),plane_H_ep_boost.Py(),plane_H_ep_boost.Pz()
                    plane_H_ep_orig=TVector3(0,0,0)
                    plane_H_ep_orig=tempH.Vect().Cross(temp_ep_orig.Vect()).Unit()
                    #print 'plane vector H,ep orig',plane_H_ep_orig.Px(),plane_H_ep_orig.Py(),plane_H_ep_orig.Pz()
                    #print 'dot product of both planes ',degrees(acos(plane_H_ep_boost.Dot (plane_Z_Z_qpos_boost)))
                    #print 'dot product of both planes orig ',degrees(acos(plane_H_ep_orig.Dot (plane_Z_Z_qneg_orig)))

                    phi_plane_parton_H_ep_Z_Z_qpos = degrees(acos(plane_H_ep_boost.Dot(plane_Z_Z_qpos_boost)))
                    #check if positive quark direction points towards normal vector of H-ep plane, or in opposite direction
                    #plane is already a TVector3
                    if degrees(temp_Z_qpos_boostZ_COM.Angle(plane_H_ep_boost))>90.:
                        phi_plane_parton_H_ep_Z_Z_qpos=degrees(acos(plane_H_ep_boost.Dot(plane_Z_Z_qpos_boost)))+180.
                    h_parton_phi_plane_Z_qpos_vs_plane_H_ep.Fill(phi_plane_rj1_ep_rj2_sj_jcE_pos,v_eventWeight[0])


                    phi_plane_parton_H_ep_Z_Z_qpos_rad = radians(phi_plane_parton_H_ep_Z_Z_qpos)
                    if sin(phi_plane_parton_H_ep_Z_Z_qpos_rad)!=0:
                        h_parton_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(sin(phi_plane_parton_H_ep_Z_Z_qpos_rad)/abs(sin(phi_plane_parton_H_ep_Z_Z_qpos_rad))*phi_plane_parton_H_ep_Z_Z_qpos,v_eventWeight[0])
                    if sin(2.*phi_plane_parton_H_ep_Z_Z_qpos_rad)!=0:
                        h_parton_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(sin(2.*phi_plane_parton_H_ep_Z_Z_qpos_rad)/abs(sin(2.*phi_plane_parton_H_ep_Z_Z_qpos_rad))*phi_plane_parton_H_ep_Z_Z_qpos,v_eventWeight[0])
                    if cos(phi_plane_parton_H_ep_Z_Z_qpos_rad)!=0:
                        h_parton_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(cos(phi_plane_parton_H_ep_Z_Z_qpos_rad)/abs(cos(phi_plane_parton_H_ep_Z_Z_qpos_rad))*phi_plane_parton_H_ep_Z_Z_qpos,v_eventWeight[0])
                    if cos(2.*phi_plane_parton_H_ep_Z_Z_qpos_rad)!=0:
                        h_parton_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(cos(2.*phi_plane_parton_H_ep_Z_Z_qpos_rad)/abs(cos(2.*phi_plane_parton_H_ep_Z_Z_qpos_rad))*phi_plane_parton_H_ep_Z_Z_qpos,v_eventWeight[0])
                    h_parton_theta_H_Z_E_totCOM.Fill(degrees(temp_H.Angle(temp_Z.Vect())),v_eventWeight[0])

                    temp_ep_approx=TLorentzVector(0,0,0,0)
                    #below is the actual calculation, but that won't work on reconstructed events, due to large smearing and clustering effects
                    #temp_ep_approx_E=sqrt(pow(-0.01*(0.5*temp_tot_Event_V2.Pz()-50.*temp_tot_Event_V2.Px()),2)+pow(0.5*temp_tot_Event_V2.Pz()-50.*temp_tot_Event_V2.Px(),2)+pow(root.TDatabasePDG.Instance().GetParticle(11).Mass(),2))
                    #temp_ep_approx.SetPxPyPzE(-0.01*(0.5*temp_tot_Event_V2.Pz()-50.*temp_tot_Event_V2.Px()),0,0.5*temp_tot_Event_V2.Pz()-50.*temp_tot_Event_V2.Px(),temp_ep_approx_E)
                    temp_ep_approx_E=sqrt(pow(0.5*temp_tot_Event_V2.E()-0.5*temp_tot_Event_V2.Pz(),2)+pow(root.TDatabasePDG.Instance().GetParticle(11).Mass(),2))
                    temp_ep_approx.SetPxPyPzE(0.,0.,-0.5*temp_tot_Event_V2.E()+0.5*temp_tot_Event_V2.Pz(),temp_ep_approx_E)

                    #temp_ep_approx.Boost(boostE_tot_COM)
 
                    h_parton_theta2_H_ep_approx_HZ_COM.Fill(degrees(temp_H.Angle(temp_ep_approx.Vect())),v_eventWeight[0])
                    temp_ep_approx.Boost(boostE_tot_COM)
                    #tiny change for approx after/before
                    #print 'ep approx after',temp_ep_approx.Px(),temp_ep_approx.Py(),temp_ep_approx.Pz(),temp_ep_approx.E(),degrees(temp_H.Angle(temp_ep_approx.Vect())),root.TDatabasePDG.Instance().GetParticle(11).Mass()

                    #gen level
                    temp_tot_Event_gj=TLorentzVector(0,0,0,0)
                    temp_tot_Event_gj.SetPxPyPzE(v_genjet1_sj1_Px[0]+v_genjet1_sj2_Px[0]+v_genjet2_sj1_Px[0]+v_genjet2_sj2_Px[0],v_genjet1_sj1_Py[0]+v_genjet1_sj2_Py[0]+v_genjet2_sj1_Py[0]+v_genjet2_sj2_Py[0],v_genjet1_sj1_Pz[0]+v_genjet1_sj2_Pz[0]+v_genjet2_sj1_Pz[0]+v_genjet2_sj2_Pz[0],v_genjet1_sj1_E[0]+v_genjet1_sj2_E[0]+v_genjet2_sj1_E[0]+v_genjet2_sj2_E[0])
                    temp_ep_approx_gj=TLorentzVector(0,0,0,0)
                    temp_ep_approx_E_gj=sqrt(pow(0.5*temp_tot_Event_gj.E()-0.5*temp_tot_Event_gj.Pz(),2)+pow(root.TDatabasePDG.Instance().GetParticle(11).Mass(),2))
                    temp_ep_approx_gj.SetPxPyPzE(0.,0.,-0.5*temp_tot_Event_gj.E()+0.5*temp_tot_Event_gj.Pz(),temp_ep_approx_E_gj)
                    boostE_tot_COM_gj=TVector3(0,0,0) 
                    boostE_tot_COM_gj=-temp_tot_Event_gj.BoostVector()
                    temp_ep_approx_gj.Boost(boostE_tot_COM_gj)
                    if temp_ep_approx_gj.Pz()>0 :
                        print 'pz should have been negative per principle ',0.5*temp_tot_Event_gj.E()-0.5*temp_tot_Event_gj.Pz(),temp_ep_approx_gj.Pz(),temp_tot_Event_gj.Pz(),temp_tot_Event_gj.Px(),temp_tot_Event_gj.Py()
                    #if (temp_ep_approx_gj.E()/temp_ep_approx.E())<0.10 or (temp_ep_approx_gj.E()/temp_ep_approx.E())>1.10 :
                    #    print 'huge difference in energy ',temp_tot_Event_gj.E(),temp_tot_Event.E()
                    #if degrees(temp_ep_approx_gj.Angle(temp_ep_approx.Vect()))>25.:
                    #    print 'huge difference in angle ',degrees(temp_ep_approx_gj.Angle(temp_ep_approx.Vect())),temp_ep_approx_gj.Vect().Unit().Px(),temp_ep_approx_gj.Vect().Unit().Py(),temp_ep_approx_gj.Vect().Unit().Pz(),temp_ep_approx.Vect().Unit().Px(),temp_ep_approx.Vect().Unit().Py(),temp_ep_approx.Vect().Unit().Pz()
                    #if degrees(temp_ep_approx_gj.Angle(temp_ep.Vect()))>20.:
                    #    print 'huge difference in real angle ',degrees(temp_ep_approx_gj.Angle(temp_ep.Vect())),temp_ep_approx_gj.Vect().Unit().Px(),temp_ep_approx_gj.Vect().Unit().Py(),temp_ep_approx_gj.Vect().Unit().Pz(),temp_ep.Vect().Unit().Px(),temp_ep.Vect().Unit().Py(),temp_ep.Vect().Unit().Pz(),temp_ep_approx_gj.E(),temp_ep.E(),temp_tot_Event_gj.E(),temp_tot_Event.E()
                    temp_gj1=TLorentzVector(0,0,0,0)
                    temp_gj1.SetPxPyPzE(v_genjet1_sj1_Px[0]+v_genjet1_sj2_Px[0],v_genjet1_sj1_Py[0]+v_genjet1_sj2_Py[0],v_genjet1_sj1_Pz[0]+v_genjet1_sj2_Pz[0],v_genjet1_sj1_E[0]+v_genjet1_sj2_E[0])
                    temp_gj1.Boost(boostE_tot_COM_gj)
                    plane_gj1_ep_approx_gj_boost=TVector3(0,0,0)
                    plane_gj1_ep_approx_gj_boost=temp_gj1.Vect().Cross(temp_ep_approx_gj.Vect()).Unit()
                    temp_gj2=TLorentzVector(0,0,0,0)
                    temp_gj2.SetPxPyPzE(v_genjet2_sj1_Px[0]+v_genjet2_sj2_Px[0],v_genjet2_sj1_Py[0]+v_genjet2_sj2_Py[0],v_genjet2_sj1_Pz[0]+v_genjet2_sj2_Pz[0],v_genjet2_sj1_E[0]+v_genjet2_sj2_E[0])
                    temp_gj2.Boost(boostE_tot_COM_gj)
                    temp_gj2_sj1=TLorentzVector(0,0,0,0)
                    temp_gj2_sj1.SetPxPyPzE(v_genjet2_sj1_Px[0],v_genjet2_sj1_Py[0],v_genjet2_sj1_Pz[0],v_genjet2_sj1_E[0])
                    temp_gj2_sj2=TLorentzVector(0,0,0,0)
                    temp_gj2_sj2.SetPxPyPzE(v_genjet2_sj2_Px[0],v_genjet2_sj2_Py[0],v_genjet2_sj2_Pz[0],v_genjet2_sj2_E[0])
                    temp_gj2_orig=TLorentzVector(0,0,0,0)
                    temp_gj2_orig.SetPxPyPzE(v_genjet2_sj1_Px[0]+v_genjet2_sj2_Px[0],v_genjet2_sj1_Py[0]+v_genjet2_sj2_Py[0],v_genjet2_sj1_Pz[0]+v_genjet2_sj2_Pz[0],v_genjet2_sj1_E[0]+v_genjet2_sj2_E[0])
                    boostgj2_COM=TVector3(0,0,0) 
                    boostgj2_COM=-temp_gj2_orig.BoostVector()
                    temp_gj2_sj1.Boost(boostgj2_COM)
                    temp_gj2_sj2.Boost(boostgj2_COM)

                    h_genjet_theta_gj1_gj2_E_totCOM.Fill(degrees(temp_gj1.Angle(temp_gj2.Vect())),v_eventWeight[0])
                    h_genjet_theta_sj1_sj2_gj2__gj2COM.Fill(degrees(temp_gj2_sj1.Angle(temp_gj2_sj2.Vect())),v_eventWeight[0])
                    h_genjet_theta2_gj1_ep_E_totCOM.Fill(degrees(temp_gj1.Angle(temp_ep_approx_gj.Vect())),v_eventWeight[0])
                    h_deltaAngle_part_gj_theta2_H_ep_E_totCOM.Fill(degrees(temp_gj1.Angle(temp_ep_approx_gj.Vect()))-degrees(temp_H.Angle(temp_ep_approx.Vect())),v_eventWeight[0])
                    #print 'ep approx part',temp_ep_approx.Px(),temp_ep_approx.Py(),temp_ep_approx.Pz(),temp_ep_approx.E()
                    #print 'ep approx gj',temp_ep_approx_gj.Px(),temp_ep_approx_gj.Py(),temp_ep_approx_gj.Pz(),temp_ep_approx_gj.E()
                    ind_gj2_sj_pos_jcE=-1
                    if abs(v_genjet2_sj1_jetChargeE_kappa_0_30[0])> abs(v_genjet2_sj2_jetChargeE_kappa_0_30[0]):
                        if v_genjet2_sj1_jetChargeE_kappa_0_30[0]>0:
                            ind_gj2_sj_pos_jcE=1
                        else:
                            ind_gj2_sj_pos_jcE=2
                    else:
                        if v_genjet2_sj2_jetChargeE_kappa_0_30[0]>0:
                            ind_gj2_sj_pos_jcE=2
                        else:
                            ind_gj2_sj_pos_jcE=1
                    if ind_gj2_sj_pos_jcE==-1:
                        print 'sth wrong in jcE subjet index'
                    temp_gj2_sj_jcE_pos=TLorentzVector(0,0,0,0)
                    if ind_gj2_sj_pos_jcE==1:
                        temp_gj2_sj_jcE_pos=temp_gj2_sj1
                    else:
                        temp_gj2_sj_jcE_pos=temp_gj2_sj2

                    temp_gj2_sj_jcE_neg=TLorentzVector(0,0,0,0)
                    if ind_gj2_sj_pos_jcE==1:
                        temp_gj2_sj_jcE_neg=temp_gj2_sj2
                    else:
                        temp_gj2_sj_jcE_neg=temp_gj2_sj1
                    plane_gj2_gj2_sj_jcE_pos_boost=TVector3(0,0,0)
                    plane_gj2_gj2_sj_jcE_pos=temp_gj2.Vect().Cross(temp_gj2_sj_jcE_pos.Vect()).Unit()

                    h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_gj2com.Fill(degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                    h_deltaAngle_part_gj_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])



                    #print "angle1,angle1 part, deltaAngle, cosAngle1, cosAngle1Part, product cosines pos",degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))
                    #print "angle1,angle1 part, deltaAngle, cosAngle1, cosAngle1Part, product cosines neg",degrees(temp_gj2_sj_jcE_neg.Angle(temp_gj2.Vect())),degrees(temp_Z_qneg_boostZ_COM.Angle(temp_Z.Vect())),degrees(temp_gj2_sj_jcE_neg.Angle(temp_gj2.Vect()))-degrees(temp_Z_qneg_boostZ_COM.Angle(temp_Z.Vect())),cos(temp_gj2_sj_jcE_neg.Angle(temp_gj2.Vect())),cos(temp_Z_qneg_boostZ_COM.Angle(temp_Z.Vect())),cos(temp_gj2_sj_jcE_neg.Angle(temp_gj2.Vect()))*cos(temp_Z_qneg_boostZ_COM.Angle(temp_Z.Vect()))


                    phi_plane_gj1_ep_gj2_sj_jcE_pos=degrees(acos(plane_gj1_ep_approx_gj_boost.Dot(plane_gj2_gj2_sj_jcE_pos)))
                    #plane is already a TVector3
                    if(degrees(temp_gj2_sj_jcE_pos.Angle(plane_gj1_ep_approx_gj_boost))>90.):
                        phi_plane_gj1_ep_gj2_sj_jcE_pos=degrees(acos(plane_gj1_ep_approx_gj_boost.Dot(plane_gj2_gj2_sj_jcE_pos)))+180.
                    h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep.Fill(phi_plane_gj1_ep_gj2_sj_jcE_pos,v_eventWeight[0])
                    h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_gj1_ep_gj2_sj_jcE_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                    if (v_genjet2_sj1_jetChargeE_kappa_0_30[0]*v_genjet2_sj2_jetChargeE_kappa_0_30[0])<0 :
                        h_deltaAngle_part_gj_jcE_0_30_theta1_Z_q_pos_Zcom_oppJC.Fill(degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep_oppJC.Fill(DeltaPhiDirAngles(phi_plane_gj1_ep_gj2_sj_jcE_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                    if abs(v_genjet2_sj1_jetChargeE_kappa_0_30[0]-v_genjet2_sj2_jetChargeE_kappa_0_30[0])>0.20:
                        h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_dm_0_20_gj2com.Fill(degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                        h_deltaAngle_part_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_dm_0_20_vs_plane_gj1_ep.Fill(phi_plane_gj1_ep_gj2_sj_jcE_pos,v_eventWeight[0])
                        h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_dm_0_20_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_gj1_ep_gj2_sj_jcE_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                        if abs(v_genjet2_sj1_jetChargeE_kappa_0_30[0]-v_genjet2_sj2_jetChargeE_kappa_0_30[0])>0.30:
                            h_genjet_theta1_gj2_pos_sj_jetChargeE_0_30_dm_0_30_gj2com.Fill(degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                            h_deltaAngle_part_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                            h_genjet_phi_plane_gj2_pos_sj_jetChargeE_0_30_dm_0_30_vs_plane_gj1_ep.Fill(phi_plane_gj1_ep_gj2_sj_jcE_pos,v_eventWeight[0])
                            h_deltaAngle_part_gj_phi_plane_Z_qpos_jcE_0_30_dm_0_30_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_gj1_ep_gj2_sj_jcE_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                    h_cosAngle_Product_gj_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    if cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0:
                        h_cosAngle_normCosPart_gj_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                    if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.20:
                        h_cosAngle_Product_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        if cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0:
                            h_cosAngle_normCosPart_gj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                            if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.30:
                                h_cosAngle_Product_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                                if cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0:
                                    h_cosAngle_normCosPart_gj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                    #print 'degrees phi plane ',degrees(phi_plane_gj1_ep_gj2_sj_jcE_pos)
                    phi_plane_gj1_ep_gj2_sj_jcE_pos_rad = radians(phi_plane_gj1_ep_gj2_sj_jcE_pos)
                    if sin(phi_plane_gj1_ep_gj2_sj_jcE_pos_rad)!=0:
                        h_genjet_sgnsinphi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep.Fill(sin(phi_plane_gj1_ep_gj2_sj_jcE_pos_rad)/abs(sin(phi_plane_gj1_ep_gj2_sj_jcE_pos_rad))*phi_plane_gj1_ep_gj2_sj_jcE_pos,v_eventWeight[0])
                    if sin(2.*phi_plane_gj1_ep_gj2_sj_jcE_pos_rad)!=0:
                        h_genjet_sgnsin2phi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep.Fill(sin(2.*phi_plane_gj1_ep_gj2_sj_jcE_pos_rad)/abs(sin(2.*phi_plane_gj1_ep_gj2_sj_jcE_pos_rad))*phi_plane_gj1_ep_gj2_sj_jcE_pos,v_eventWeight[0])
                    if cos(phi_plane_gj1_ep_gj2_sj_jcE_pos_rad)!=0:
                        h_genjet_sgncosphi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep.Fill(cos(phi_plane_gj1_ep_gj2_sj_jcE_pos_rad)/abs(cos(phi_plane_gj1_ep_gj2_sj_jcE_pos_rad))*phi_plane_gj1_ep_gj2_sj_jcE_pos,v_eventWeight[0])
                    if cos(2.*phi_plane_gj1_ep_gj2_sj_jcE_pos_rad)!=0:
                        h_genjet_sgncos2phi_phi_plane_gj2_pos_sj_jetChargeE_0_30_vs_plane_gj1_ep.Fill(cos(2.*phi_plane_gj1_ep_gj2_sj_jcE_pos_rad)/abs(cos(2.*phi_plane_gj1_ep_gj2_sj_jcE_pos_rad))*phi_plane_gj1_ep_gj2_sj_jcE_pos,v_eventWeight[0])

                    ind_gj2_sj_pos_jcPt=-1
                    if abs(v_genjet2_sj1_jetChargePt_kappa_0_30[0])> abs(v_genjet2_sj2_jetChargePt_kappa_0_30[0]):
                        if v_genjet2_sj1_jetChargePt_kappa_0_30[0]>0:
                            ind_gj2_sj_pos_jcPt=1
                        else:
                            ind_gj2_sj_pos_jcPt=2
                    else:
                        if v_genjet2_sj2_jetChargePt_kappa_0_30[0]>0:
                            ind_gj2_sj_pos_jcPt=2
                        else:
                            ind_gj2_sj_pos_jcPt=1
                    if ind_gj2_sj_pos_jcPt==-1:
                        print 'sth wrong in jcPt subjet index'
                    temp_gj2_sj_jcPt_pos=TLorentzVector(0,0,0,0)
                    if ind_gj2_sj_pos_jcPt==1:
                        temp_gj2_sj_jcPt_pos=temp_gj2_sj1
                    else:
                        temp_gj2_sj_jcPt_pos=temp_gj2_sj2
                    plane_gj2_gj2_sj_jcPt_pos_boost=TVector3(0,0,0)
                    plane_gj2_gj2_sj_jcPt_pos=temp_gj2.Vect().Cross(temp_gj2_sj_jcPt_pos.Vect()).Unit()
                    h_genjet_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com.Fill(cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                    h_genjet_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com.Fill(cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                    if cos(2.*temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))!=0:
                        h_genjet_sgncos2theta1_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com.Fill(cos(2.*temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))/abs(cos(2.*temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())))*cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                    if cos(2.*temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))!=0:
                        h_genjet_sgncos2theta1_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com.Fill(cos(2.*temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))/abs(cos(2.*temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())))*cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])


                    h_2D_genjet_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com_vs_costheta2_gj1_ep_E_totCOM.Fill(cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),cos(temp_gj1.Angle(temp_ep_approx_gj.Vect())),v_eventWeight[0])
                    h_2D_genjet_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com_vs_costheta2_gj1_ep_E_totCOM.Fill(cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),cos(temp_gj1.Angle(temp_ep_approx_gj.Vect())),v_eventWeight[0])
                    if cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))!=0 and cos(temp_gj1.Angle(temp_ep_approx_gj.Vect()))!=0:
                        h_2D_genjet_sgncostheta1_costheta1_gj2_pos_sj_jetChargeE_0_30_gj2com_vs_sgncostheta2_costheta2_gj1_ep_E_totCOM.Fill(cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect()))/abs(cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())))*cos(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),cos(temp_gj1.Angle(temp_ep_approx_gj.Vect()))/abs(cos(temp_gj1.Angle(temp_ep_approx_gj.Vect())))*cos(temp_gj1.Angle(temp_ep_approx_gj.Vect())),v_eventWeight[0])
                    if cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))!=0 and cos(temp_gj1.Angle(temp_ep_approx_gj.Vect()))!=0:
                        h_2D_genjet_sgncostheta1_costheta1_gj2_pos_sj_jetChargePt_0_30_gj2com_vs_sgncostheta2_costheta2_gj1_ep_E_totCOM.Fill(cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))/abs(cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())))*cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),cos(temp_gj1.Angle(temp_ep_approx_gj.Vect()))/abs(cos(temp_gj1.Angle(temp_ep_approx_gj.Vect())))*cos(temp_gj1.Angle(temp_ep_approx_gj.Vect())),v_eventWeight[0])


                    h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_gj2com.Fill(degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                    h_deltaAngle_part_gj_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])

                    phi_plane_gj1_ep_gj2_sj_jcPt_pos=degrees(acos(plane_gj1_ep_approx_gj_boost.Dot(plane_gj2_gj2_sj_jcPt_pos)))
                    #plane is already a TVector3
                    if(degrees(temp_gj2_sj_jcPt_pos.Angle(plane_gj1_ep_approx_gj_boost))>90.):
                        phi_plane_gj1_ep_gj2_sj_jcPt_pos=degrees(acos(plane_gj1_ep_approx_gj_boost.Dot(plane_gj2_gj2_sj_jcPt_pos)))+180.

                    #print 'degrees phi plane ',degrees(phi_plane_gj1_ep_gj2_sj_jcPt_pos)
                    phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad = radians(phi_plane_gj1_ep_gj2_sj_jcPt_pos)
                    if sin(phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad)!=0:
                        h_genjet_sgnsinphi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep.Fill(sin(phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad)/abs(sin(phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad))*phi_plane_gj1_ep_gj2_sj_jcPt_pos,v_eventWeight[0])
                    if sin(2.*phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad)!=0:
                        h_genjet_sgnsin2phi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep.Fill(sin(2.*phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad)/abs(sin(2.*phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad))*phi_plane_gj1_ep_gj2_sj_jcPt_pos,v_eventWeight[0])
                    if cos(phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad)!=0:
                        h_genjet_sgncosphi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep.Fill(cos(phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad)/abs(cos(phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad))*phi_plane_gj1_ep_gj2_sj_jcPt_pos,v_eventWeight[0])
                    if cos(2.*phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad)!=0:
                        h_genjet_sgncos2phi_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep.Fill(cos(2.*phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad)/abs(cos(2.*phi_plane_gj1_ep_gj2_sj_jcPt_pos_rad))*phi_plane_gj1_ep_gj2_sj_jcPt_pos,v_eventWeight[0])


                    h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_vs_plane_gj1_ep.Fill(phi_plane_gj1_ep_gj2_sj_jcPt_pos,v_eventWeight[0])
                    h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_gj1_ep_gj2_sj_jcPt_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                    if (v_genjet2_sj1_jetChargePt_kappa_0_30[0]*v_genjet2_sj2_jetChargePt_kappa_0_30[0])<0 :
                        h_deltaAngle_part_gj_jcPt_0_30_theta1_Z_q_pos_Zcom_oppJC.Fill(degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep_oppJC.Fill(DeltaPhiDirAngles(phi_plane_gj1_ep_gj2_sj_jcPt_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                    if abs(v_genjet2_sj1_jetChargePt_kappa_0_30[0]-v_genjet2_sj2_jetChargePt_kappa_0_30[0])>0.20:
                        h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_dm_0_20_gj2com.Fill(degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                        h_deltaAngle_part_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_dm_0_20_vs_plane_gj1_ep.Fill(phi_plane_gj1_ep_gj2_sj_jcPt_pos,v_eventWeight[0])
                        h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_dm_0_20_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_gj1_ep_gj2_sj_jcPt_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                        if abs(v_genjet2_sj1_jetChargePt_kappa_0_30[0]-v_genjet2_sj2_jetChargePt_kappa_0_30[0])>0.30:
                            h_genjet_theta1_gj2_pos_sj_jetChargePt_0_30_dm_0_30_gj2com.Fill(degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                            h_deltaAngle_part_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                            h_genjet_phi_plane_gj2_pos_sj_jetChargePt_0_30_dm_0_30_vs_plane_gj1_ep.Fill(phi_plane_gj1_ep_gj2_sj_jcPt_pos,v_eventWeight[0])
                            h_deltaAngle_part_gj_phi_plane_Z_qpos_jcPt_0_30_dm_0_30_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_gj1_ep_gj2_sj_jcPt_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                    h_cosAngle_Product_gj_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    if cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0:
                        h_cosAngle_normCosPart_gj_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                    if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.20:
                        h_cosAngle_Product_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        if cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0:
                            h_cosAngle_normCosPart_gj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                            if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.30:
                                h_cosAngle_Product_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                                if cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0:
                                    h_cosAngle_normCosPart_gj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                    h_deltaAngle_part_rj_theta2_H_ep_E_totCOM.Fill(degrees(temp_rj1.Angle(temp_ep_approx_rj.Vect()))-degrees(temp_H.Angle(temp_ep_approx.Vect())),v_eventWeight[0])
                    h_deltaAngle_part_rj_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_cosAngle_Product_rj_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_cosAngle_Product_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_cosAngle_Product_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_nTracksMax.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_cosAngle_Product_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_chFrac.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_cosAngle_Product_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_chEnergy.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_cosAngle_Product_rj_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_cosAngle_Product_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_cosAngle_Product_rj_nTrack_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_nTracksMax.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_cosAngle_Product_rj_chFrac_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_chFrac.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_cosAngle_Product_rj_chEnergy_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_chEnergy.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    if cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0:
                        h_cosAngle_normCosPart_rj_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                        h_cosAngle_normCosPart_rj_E1_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                        h_cosAngle_normCosPart_rj_nTrack_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_nTracksMax.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                        h_cosAngle_normCosPart_rj_chFrac_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_chFrac.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                        h_cosAngle_normCosPart_rj_chEnergy_jcE_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_chEnergy.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                        h_cosAngle_normCosPart_rj_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                        h_cosAngle_normCosPart_rj_E1_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                        h_cosAngle_normCosPart_rj_nTrack_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_nTracksMax.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                        h_cosAngle_normCosPart_rj_chFrac_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_chFrac.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                        h_cosAngle_normCosPart_rj_chEnergy_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_chEnergy.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                    if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.20:
                        h_cosAngle_Product_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_cosAngle_Product_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_cosAngle_Product_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_cosAngle_Product_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        if cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0:
                            h_cosAngle_normCosPart_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                            h_cosAngle_normCosPart_rj_E1_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                            h_cosAngle_normCosPart_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                            h_cosAngle_normCosPart_rj_E1_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                        if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.30:
                            h_cosAngle_Product_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                            h_cosAngle_Product_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                            h_cosAngle_Product_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                            h_cosAngle_Product_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                            if cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))!=0:
                                h_cosAngle_normCosPart_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                                h_cosAngle_normCosPart_rj_E1_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcE_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                                h_cosAngle_normCosPart_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                                h_cosAngle_normCosPart_rj_E1_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(cos(temp_rj2_sj_jcPt_pos_sj_E1.Angle(temp_rj2.Vect()))*cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))/abs(cos(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect()))),v_eventWeight[0])
                    h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcE_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                    if (abs(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcE_pos,phi_plane_parton_H_ep_Z_Z_qpos)))>180.:
                        print 'deltaphi of angles (dir in angles) bigger than 180',abs(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcE_pos,phi_plane_parton_H_ep_Z_Z_qpos))

                    if (v_jet2_sj1_jetChargeE_kappa_0_30[0]*v_jet2_sj2_jetChargeE_kappa_0_30[0])<0 :
                        h_deltaAngle_part_rj_jcE_0_30_theta1_Z_q_pos_Zcom_oppJC.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_vs_plane_H_ep_oppJC.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcE_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                    if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.20:
                        h_deltaAngle_part_rj_jcE_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_dm_0_20_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcE_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                        if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.30:
                            h_deltaAngle_part_rj_jcE_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                            h_deltaAngle_part_rj_phi_plane_Z_qpos_jcE_0_30_dm_0_30_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcE_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])

                    h_deltaAngle_part_rj_jcPt_0_30_theta1_Z_q_pos_Zcom.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                    h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                    if (v_jet2_sj1_jetChargePt_kappa_0_30[0]*v_jet2_sj2_jetChargePt_kappa_0_30[0])<0 :
                        h_deltaAngle_part_rj_jcPt_0_30_theta1_Z_q_pos_Zcom_oppJC.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_vs_plane_H_ep_oppJC.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                    if abs(v_jet2_sj1_jetChargePt_kappa_0_30[0]-v_jet2_sj2_jetChargePt_kappa_0_30[0])>0.20:
                        h_deltaAngle_part_rj_jcPt_0_30_dm_0_20_theta1_Z_q_pos_Zcom.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                        h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_dm_0_20_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])
                        if abs(v_jet2_sj1_jetChargePt_kappa_0_30[0]-v_jet2_sj2_jetChargePt_kappa_0_30[0])>0.30:
                            h_deltaAngle_part_rj_jcPt_0_30_dm_0_30_theta1_Z_q_pos_Zcom.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_Z_qpos_boostZ_COM.Angle(temp_Z.Vect())),v_eventWeight[0])
                            h_deltaAngle_part_rj_phi_plane_Z_qpos_jcPt_0_30_dm_0_30_vs_plane_H_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_parton_H_ep_Z_Z_qpos),v_eventWeight[0])

                    h_deltaAngle_gj_rj_phi_gj1_ep_E_totCOM.Fill(degrees(temp_rj1.Angle(temp_ep_approx_rj.Vect()))-degrees(temp_gj1.Angle(temp_ep_approx_gj.Vect())),v_eventWeight[0])
                    h_deltaAngle_gj_rj_jcE_0_30_theta_gj_jcE_0_30.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                    h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_vs_plane_gj1_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcE_pos,phi_plane_gj1_ep_gj2_sj_jcE_pos),v_eventWeight[0])
                    if (v_jet2_sj1_jetChargeE_kappa_0_30[0]*v_jet2_sj2_jetChargeE_kappa_0_30[0])<0 :
                        h_deltaAngle_gj_rj_jcE_0_30_theta_gj_jcE_0_30_oppJC.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                        h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_vs_plane_gj1_ep_oppJC.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcE_pos,phi_plane_gj1_ep_gj2_sj_jcE_pos),v_eventWeight[0])
                    if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.20:
                        h_deltaAngle_gj_rj_jcE_0_30_dm_0_20_theta_gj_jcE_0_30.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                        h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_dm_0_20_vs_plane_gj1_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcE_pos,phi_plane_gj1_ep_gj2_sj_jcE_pos),v_eventWeight[0])
                        if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.30:
                            h_deltaAngle_gj_rj_jcE_0_30_dm_0_30_theta_gj_jcE_0_30.Fill(degrees(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcE_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                            h_deltaAngle_gj_rj_phi_plane_gj2_jcE_0_30_dm_0_30_vs_plane_gj1_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcE_pos,phi_plane_gj1_ep_gj2_sj_jcE_pos),v_eventWeight[0])

                    h_deltaAngle_gj_rj_jcPt_0_30_theta_gj_jcPt_0_30.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                    h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_vs_plane_gj1_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_gj1_ep_gj2_sj_jcPt_pos),v_eventWeight[0])
                    if (v_jet2_sj1_jetChargePt_kappa_0_30[0]*v_jet2_sj2_jetChargePt_kappa_0_30[0])<0 :
                        h_deltaAngle_gj_rj_jcPt_0_30_theta_gj_jcPt_0_30_oppJC.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                        h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_vs_plane_gj1_ep_oppJC.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_gj1_ep_gj2_sj_jcPt_pos),v_eventWeight[0])
                    if abs(v_jet2_sj1_jetChargePt_kappa_0_30[0]-v_jet2_sj2_jetChargePt_kappa_0_30[0])>0.20:
                        h_deltaAngle_gj_rj_jcPt_0_30_dm_0_20_theta_gj_jcPt_0_30.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                        h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_dm_0_20_vs_plane_gj1_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_gj1_ep_gj2_sj_jcPt_pos),v_eventWeight[0])
                        if abs(v_jet2_sj1_jetChargePt_kappa_0_30[0]-v_jet2_sj2_jetChargePt_kappa_0_30[0])>0.30:
                            h_deltaAngle_gj_rj_jcPt_0_30_dm_0_30_theta_gj_jcPt_0_30.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                            h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_dm_0_30_vs_plane_gj1_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_gj1_ep_gj2_sj_jcPt_pos),v_eventWeight[0])

                    h_deltaAngle_gj_rj_jcPt_0_30_theta_gj_jcPt_0_30.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                    h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_vs_plane_gj1_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_gj1_ep_gj2_sj_jcPt_pos),v_eventWeight[0])
                    if (v_jet2_sj1_jetChargePt_kappa_0_30[0]*v_jet2_sj2_jetChargePt_kappa_0_30[0])<0 :
                        h_deltaAngle_gj_rj_jcPt_0_30_theta_gj_jcPt_0_30_oppJC.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                        h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_vs_plane_gj1_ep_oppJC.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_gj1_ep_gj2_sj_jcPt_pos),v_eventWeight[0])
                    if abs(v_jet2_sj1_jetChargePt_kappa_0_30[0]-v_jet2_sj2_jetChargePt_kappa_0_30[0])>0.20:
                        h_deltaAngle_gj_rj_jcPt_0_30_dm_0_20_theta_gj_jcPt_0_30.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                        h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_dm_0_20_vs_plane_gj1_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_gj1_ep_gj2_sj_jcPt_pos),v_eventWeight[0])
                        if abs(v_jet2_sj1_jetChargePt_kappa_0_30[0]-v_jet2_sj2_jetChargePt_kappa_0_30[0])>0.30:
                            h_deltaAngle_gj_rj_jcPt_0_30_dm_0_30_theta_gj_jcPt_0_30.Fill(degrees(temp_rj2_sj_jcPt_pos.Angle(temp_rj2.Vect()))-degrees(temp_gj2_sj_jcPt_pos.Angle(temp_gj2.Vect())),v_eventWeight[0])
                            h_deltaAngle_gj_rj_phi_plane_gj2_jcPt_0_30_dm_0_30_vs_plane_gj1_ep.Fill(DeltaPhiDirAngles(phi_plane_rj1_ep_rj2_sj_jcPt_pos,phi_plane_gj1_ep_gj2_sj_jcPt_pos),v_eventWeight[0])



                    if v_jet2_sj1_closestMatch[0]<2:
                        h_jet2_mass_H_matched.Fill(v_jet2_mass[0],v_eventWeight[0])
                    elif v_jet2_sj1_closestMatch[0]!=-1 :
                        h_jet2_mass_Z_matched.Fill(v_jet2_mass[0],v_eventWeight[0])
                    if v_jet2_sj1_closestMatch[0]==0:
                        h_jet2_sj_b_cM_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    elif v_jet2_sj1_closestMatch[0]==1:
                        h_jet2_sj_bbar_cM_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    if v_jet2_sj2_closestMatch[0]==0:
                        h_jet2_sj_b_cM_jetChargeE_kappa_0_30.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    elif v_jet2_sj2_closestMatch[0]==1:
                        h_jet2_sj_bbar_cM_jetChargeE_kappa_0_30.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])

                    if v_jet2_sj1_decMatch[0]==0:
                        h_jet2_sj_b_dM_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    elif v_jet2_sj1_decMatch[0]==1:
                        h_jet2_sj_bbar_dM_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    if v_jet2_sj2_decMatch[0]==0:
                        h_jet2_sj_b_dM_jetChargeE_kappa_0_30.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    elif v_jet2_sj2_decMatch[0]==1:
                        h_jet2_sj_bbar_dM_jetChargeE_kappa_0_30.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    if v_jet2_sj1_closestMatch[0]==2 and v_jet2_sj2_closestMatch[0]==3 :
                        h_jet2_sj_qneg_dM_vs_qpos_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        if (v_jet2_sj1_jetChargeE_kappa_0_30[0]*v_jet2_sj2_jetChargeE_kappa_0_30[0])<0:
                            h_jet2_sj_qneg_dM_vs_qpos_jetChargeE_kappa_0_30_oppCharge.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    elif v_jet2_sj1_closestMatch[0]==3 and v_jet2_sj2_closestMatch[0]==2 :
                        h_jet2_sj_qneg_dM_vs_qpos_jetChargeE_kappa_0_30.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        if (v_jet2_sj1_jetChargeE_kappa_0_30[0]*v_jet2_sj2_jetChargeE_kappa_0_30[0])<0:
                            h_jet2_sj_qneg_dM_vs_qpos_jetChargeE_kappa_0_30_oppCharge.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])

                    if v_jet2_sj1_closestMatch[0]==2:
                        h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    elif v_jet2_sj1_closestMatch[0]==3:
                        h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    if v_jet2_sj2_closestMatch[0]==2:
                        h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    elif v_jet2_sj2_closestMatch[0]==3:
                        h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])


                    if (v_jet2_sj1_jetChargeE_kappa_0_30[0]*v_jet2_sj2_jetChargeE_kappa_0_30[0])<0:
                        if v_jet2_sj1_closestMatch[0]==2:
                            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_oppCharge.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        elif v_jet2_sj1_closestMatch[0]==3:
                            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_oppCharge.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        if v_jet2_sj2_closestMatch[0]==2:
                            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_oppCharge.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        elif v_jet2_sj2_closestMatch[0]==3:
                            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_oppCharge.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])


                    if v_jet2_sj1_decMatch[0]==2:
                        h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.20:
                            if(abs(v_jet2_sj1_jetChargeE_kappa_0_30[0])>abs(v_jet2_sj2_jetChargeE_kappa_0_30[0])):
                                if v_jet2_sj1_jetChargeE_kappa_0_30[0]>0 :
                                    h_jet2_sj_qneg_dM_ChargeHemisphere_kappa_0_30_dm_jetCharge_0_20.Fill(-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                else: 
                                    h_jet2_sj_qneg_dM_ChargeHemisphere_kappa_0_30_dm_jetCharge_0_20.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            else:
                                if v_jet2_sj2_jetChargeE_kappa_0_30[0]>0 :
                                    h_jet2_sj_qneg_dM_ChargeHemisphere_kappa_0_30_dm_jetCharge_0_20.Fill(-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                else: 
                                    h_jet2_sj_qneg_dM_ChargeHemisphere_kappa_0_30_dm_jetCharge_0_20.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            h_jet2_sj_qneg_dM_jetChargeE_vs_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.25:
                                h_jet2_sj_qneg_dM_jetChargeE_vs_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_25.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    elif v_jet2_sj1_decMatch[0]==3:
                        h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.20:
                            h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            h_jet2_sj_qneg_dM_jetChargeE_vs_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_20.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            if abs(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0])>0.25:
                                h_jet2_sj_qneg_dM_jetChargeE_vs_qpos_dM_jetChargeE_kappa_0_30_dm_jetCharge_0_25.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    if v_jet2_sj2_decMatch[0]==2:
                        h_jet2_sj_qneg_dM_jetChargeE_kappa_0_30.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    elif v_jet2_sj2_decMatch[0]==3:
                        h_jet2_sj_qpos_dM_jetChargeE_kappa_0_30.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])

                    if v_jet2_sj1_closestMatch[0]==v_jet2_sj1_decMatch[0] and v_jet2_sj2_closestMatch[0]==v_jet2_sj2_decMatch[0]:
                        h_jet2_E_sj2_over_E_sum_dec_also_close.Fill(v_jet2_sj2_E[0]/(v_jet2_sj1_E[0]+v_jet2_sj2_E[0]),v_eventWeight[0])
                        h_jet2_dAlpha_sj1sj2_dec_also_close.Fill(v_jet2_dAlpha_sj1sj2[0],v_eventWeight[0])
                        if(  ((v_jet2_sj1_closestMatch[0]==2 and v_jet2_sj1_jetChargeE_kappa_0_30[0]<0) and (v_jet2_sj2_closestMatch[0]==3 and v_jet2_sj2_jetChargeE_kappa_0_30[0]>0)) or  ((v_jet2_sj1_closestMatch[0]==3 and v_jet2_sj1_jetChargeE_kappa_0_30[0]>0) and (v_jet2_sj2_closestMatch[0]==2 and v_jet2_sj2_jetChargeE_kappa_0_30[0]<0))):
                            h_jet2_E_sj2_over_E_sum_jetChargeE_kappa_0_30_dec_also_close_rightCharges.Fill(v_jet2_sj2_E[0]/(v_jet2_sj1_E[0]+v_jet2_sj2_E[0]),v_eventWeight[0])
                            h_jet2_dAlpha_sj1sj2_jetChargeE_kappa_0_30_dec_also_close_rightCharges.Fill(v_jet2_dAlpha_sj1sj2[0],v_eventWeight[0])
                        elif( (v_jet2_sj1_closestMatch[0]==2 and v_jet2_sj1_jetChargeE_kappa_0_30[0]>0) or (v_jet2_sj1_closestMatch[0]==3 and v_jet2_sj1_jetChargeE_kappa_0_30[0]<0) or (v_jet2_sj1_closestMatch[0]==3 and v_jet2_sj1_jetChargeE_kappa_0_30[0]<0) or (v_jet2_sj1_closestMatch[0]==2 and v_jet2_sj1_jetChargeE_kappa_0_30[0]>0)):
                            h_jet2_E_sj2_over_E_sum_jetChargeE_kappa_0_30_dec_also_close_wrongCharges.Fill(v_jet2_sj2_E[0]/(v_jet2_sj1_E[0]+v_jet2_sj2_E[0]),v_eventWeight[0])      
                            h_jet2_dAlpha_sj1sj2_jetChargeE_kappa_0_30_dec_also_close_wrongCharges.Fill(v_jet2_dAlpha_sj1sj2[0],v_eventWeight[0])
                        if v_jet2_sj1_closestMatch[0]==0:
                            h_jet2_sj_b_cM_jetChargeE_kappa_0_30_dec_also_close.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        elif v_jet2_sj1_closestMatch[0]==1:
                            h_jet2_sj_bbar_cM_jetChargeE_kappa_0_30_dec_also_close.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        if v_jet2_sj2_closestMatch[0]==0:
                            h_jet2_sj_b_cM_jetChargeE_kappa_0_30_dec_also_close.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        elif v_jet2_sj2_closestMatch[0]==1:
                            h_jet2_sj_bbar_cM_jetChargeE_kappa_0_30_dec_also_close.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        if v_jet2_sj1_closestMatch[0]==2:
                            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            if v_jet2_sj1_jetChargeE_kappa_0_30[0]<0:
                                h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_chFrac.Fill(v_jet2_sj1_chFrac[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_nTrack.Fill(v_jet2_sj1_nTracks[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qneg.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qneg.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            else:
                                h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_chFrac.Fill(v_jet2_sj1_chFrac[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_nTrack.Fill(v_jet2_sj1_nTracks[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                if v_jet2_sj2_jetChargeE_kappa_0_30[0]<0:
                                    h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                    h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        elif v_jet2_sj1_closestMatch[0]==3:
                            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            if v_jet2_sj1_jetChargeE_kappa_0_30[0]>0:
                                h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_chFrac.Fill(v_jet2_sj1_chFrac[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_nTrack.Fill(v_jet2_sj1_nTracks[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qpos.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qpos.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            else:
                                h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_chFrac.Fill(v_jet2_sj1_chFrac[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_nTrack.Fill(v_jet2_sj1_nTracks[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qpos.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qpos.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                if v_jet2_sj2_jetChargeE_kappa_0_30[0]>0:
                                    h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                    h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        if v_jet2_sj2_closestMatch[0]==2:
                            h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            if v_jet2_sj2_jetChargeE_kappa_0_30[0]<0:
                                h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_chFrac.Fill(v_jet2_sj2_chFrac[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_nTrack.Fill(v_jet2_sj2_nTracks[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qneg.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qneg.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            else:
                                h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_chFrac.Fill(v_jet2_sj2_chFrac[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_nTrack.Fill(v_jet2_sj2_nTracks[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                if v_jet2_sj1_jetChargeE_kappa_0_30[0]<0:
                                    h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                    h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                        elif v_jet2_sj2_closestMatch[0]==3:
                            h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            if v_jet2_sj2_jetChargeE_kappa_0_30[0]>0:
                                h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_chFrac.Fill(v_jet2_sj2_chFrac[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_nTrack.Fill(v_jet2_sj2_nTracks[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qpos.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_rightCharge_qpos.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                            else:
                                h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_chFrac.Fill(v_jet2_sj2_chFrac[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_nTrack.Fill(v_jet2_sj2_nTracks[0],v_eventWeight[0])
                                h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qpos.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qpos.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                if v_jet2_sj1_jetChargeE_kappa_0_30[0]>0:
                                    h_jet2_sj_qneg_cM_min_qpos_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos.Fill(v_jet2_sj1_jetChargeE_kappa_0_30[0]-v_jet2_sj2_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                                    h_jet2_sj_qpos_cM_min_qneg_cM_jetChargeE_kappa_0_30_dec_also_close_wrongCharge_qneg_and_qpos.Fill(v_jet2_sj2_jetChargeE_kappa_0_30[0]-v_jet2_sj1_jetChargeE_kappa_0_30[0],v_eventWeight[0])
                    #print 'in signal histos'




        print 'at end of file processing of ',i_final_histo_name_,bdt_value,h_jet1_mass.Integral()
                
        fileout.Write()
    fileout.Close()

 
    return None

def process_files():



    #files_weights_='/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/dataset/filesPolp80GiniIndexNoNormSkipNormalization/dataset/weights/TMVAClassification_BDT.weights.xml'
    #filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage0_75NCuts20
    #GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees175NCuts20
    #files_weights_='/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees175NCuts20_qqqq2TeV_allVar/dataset/weights/TMVAClassification_BDT.weights.xml'
    files_weights_='/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees250AdaBoostBeta020NCuts20_qqqq2TeV_allVar/dataset/weights/TMVAClassification_BDT.weights.xml'

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj.root"  
    process_event(final_histo_name_,input_file_,files_weights_,isSignalData_)

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj_AllEvents.root"  
    process_event(final_histo_name_,input_file_,files_weights_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj.root" 
    process_event(final_histo_name_,input_file_,files_weights_,isSignalData_)

 
    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/test_ee_qqqq_mqqqq_2TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj.root" 
    process_event(final_histo_name_,input_file_,files_weights_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/test_ee_qqqqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj.root" 
    process_event(final_histo_name_,input_file_,files_weights_,isSignalData_)

                     
    #files_weights_='/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_qqqq2TeV_allVar/dataset/weights/TMVAClassification_BDT.weights.xml'           /eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_June24/dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees250AdaBoostBeta020NCuts20_qqqq2TeV_allVar/dataset/weights/
    files_weights_='/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees250AdaBoostBeta020NCuts20_qqqq2TeV_allVar/dataset/weights/TMVAClassification_BDT.weights.xml'

    isSignalData_=True
    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj.root"  
    process_event(final_histo_name_,input_file_,files_weights_,isSignalData_)

    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_AllEvents.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj_AllEvents.root"  
    isSignalData_=True
    process_event(final_histo_name_,input_file_,files_weights_,isSignalData_)

    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj.root" 
    isSignalData_=False
    process_event(final_histo_name_,input_file_,files_weights_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/test_ee_qqqq_mqqqq_2TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj.root" 
    process_event(final_histo_name_,input_file_,files_weights_,isSignalData_)

    isSignalData_=False
    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/test_ee_qqqqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/MVATrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020Trees250NCuts20_AnglesMETProj.root" 
    process_event(final_histo_name_,input_file_,files_weights_,isSignalData_)






    return None

process_files()



