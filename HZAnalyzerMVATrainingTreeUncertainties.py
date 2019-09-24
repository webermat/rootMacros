from ROOT import gROOT, TCanvas, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor
import ROOT as root
from math import cos, sin, pi, degrees, radians, pow, sqrt,acos, log,tan
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

def fill_background_histograms(file,mytree,xsec,usePartonInfo,lumi,performMassCuts,performThetaCuts,fillGenLevel,performJESUnc,JESUncValue,BTagUncValue,norm_up,norm_down):
    print "do something"

    use_EMissNeutrinoProjection=True
    #here use total 4 vector, or 4 vector sum of jet 1 and 2 (see flag below)- isolated photon four vector plus correction with EMiss on both jet axes
    #mass cuts are then also done after projecting the EMiss
    use_sqrtJets=True #in this case use j1+j2-isolated photons and with upper flag still decide if EMiss projection on jets is performed

    fCut_mass_1_min=0.
    fCut_mass1_center=126.
    fCut_mass2_center=92.5
    fCut_mass1_radius=35.
    fCut_mass2_radius=35.


    fCut_thetaWindow=70.
    fCut_thetaRef=90.
    fCut_delta_theta = 100.
 


    fCut_thetaWindow=70.
    fCut_thetaRef=90.


    fCut_JESFactor=1.0
    if performJESUnc:
        fCut_JESFactor+=JESUncValue
        print "value of fCut_JESFactor", fCut_JESFactor
    
    t_var_sqrtS  = array('f',[0])
    t_var_sqrtS_orig  = array('f',[0])
    t_var_sqrtS_gen  = array('f',[0])
    t_var_sqrtS_orig_gen  = array('f',[0])
    t_var_sqrtS_parton  = array('f',[0])
    t_var_nLeptons = array('i',[0])
    t_var_nLeptons_gen = array('i',[0])
    t_var_MET = array('f',[0])

    t_var_reco_y21 = array('f',[0])
    t_var_reco_y32 = array('f',[0])
    t_var_reco_y43 = array('f',[0])

    t_var_jet1_mass = array('f',[0])
    t_var_jet2_mass = array('f',[0])
    t_var_jet1_min_jet2_mass = array('f',[0])
    t_var_jet1_BTag_rfj_BTagMax = array('f',[0])
    t_var_jet1_CTag_rfj_BTagMax = array('f',[0])
    t_var_jet1_LTag_rfj_BTagMax = array('f',[0])
    t_var_jet1_BTag_rfj_CTagMax = array('f',[0])
    t_var_jet1_CTag_rfj_CTagMax = array('f',[0])
    t_var_jet1_LTag_rfj_CTagMax = array('f',[0])
    t_var_jet1_E = array('f',[0])
    t_var_jet2_E = array('f',[0])
    t_var_jet1_Pt = array('f',[0])
    t_var_jet2_Pt = array('f',[0])
    t_var_jet1_theta = array('f',[0])
    t_var_jet2_theta = array('f',[0])
    t_var_jet1_min_jet2_theta = array('f',[0])
    t_var_jet1_phi = array('f',[0])
    t_var_jet2_phi = array('f',[0])
    t_var_dphi_j1j2 = array('f',[0])
    t_var_angle_j1j2 = array('f',[0])
    t_var_jet1_D2_beta1 = array('f',[0])
    t_var_jet2_D2_beta1 = array('f',[0])
    t_var_jet1_D2_beta0_5 = array('f',[0])
    t_var_jet2_D2_beta0_5 = array('f',[0])
    t_var_jet1_C2_beta1 = array('f',[0])
    t_var_jet2_C2_beta1 = array('f',[0])
    t_var_jet1_C2_beta0_5 = array('f',[0])
    t_var_jet2_C2_beta0_5 = array('f',[0])
    t_var_jet1_tau21 = array('f',[0])
    t_var_jet2_tau21 = array('f',[0])
    t_var_jet1_tau32 = array('f',[0])
    t_var_jet2_tau32 = array('f',[0])
    t_var_jet1_C3_beta1 = array('f',[0])
    t_var_jet2_C3_beta1 = array('f',[0])
    t_var_jet1_N2_beta1 = array('f',[0])
    t_var_jet2_N2_beta1 = array('f',[0])
    t_var_jet1_N3_beta1 = array('f',[0])
    t_var_jet2_N3_beta1 = array('f',[0])

    t_var_jet1_D2_beta1_E_theta = array('f',[0])
    t_var_jet2_D2_beta1_E_theta = array('f',[0])
    t_var_jet1_C2_beta1_E_theta = array('f',[0])
    t_var_jet2_C2_beta1_E_theta = array('f',[0])
    t_var_jet1_tau21_E_theta = array('f',[0])
    t_var_jet2_tau21_E_theta = array('f',[0])
    t_var_jet1_tau32_E_theta = array('f',[0])
    t_var_jet2_tau32_E_theta = array('f',[0])
    t_var_jet1_C3_beta1_E_theta = array('f',[0])
    t_var_jet2_C3_beta1_E_theta = array('f',[0])
    t_var_jet1_N2_beta1_E_theta = array('f',[0])
    t_var_jet2_N2_beta1_E_theta = array('f',[0])
    t_var_jet1_N3_beta1_E_theta = array('f',[0])
    t_var_jet2_N3_beta1_E_theta = array('f',[0])

    t_var_jet1_d21 = array('f',[0])
    t_var_jet1_d32 = array('f',[0])
    t_var_jet1_d43 = array('f',[0])
    t_var_jet2_d21 = array('f',[0])
    t_var_jet2_d32 = array('f',[0])
    t_var_jet2_d43 = array('f',[0])

    t_var_jet1_sj1_E = array('f',[0])
    t_var_jet1_sj1_Px = array('f',[0])
    t_var_jet1_sj1_Py = array('f',[0])
    t_var_jet1_sj1_Pz = array('f',[0])
    t_var_jet1_sj1_jetChargePt_kappa_0_30 = array('f',[0])
    t_var_jet1_sj1_jetChargeE_kappa_0_30 = array('f',[0])
    t_var_jet1_sj1_jetChargePt_kappa_0_25 = array('f',[0])
    t_var_jet1_sj1_jetChargeE_kappa_0_25 = array('f',[0])
    t_var_jet1_sj1_nTracks = array('i',[0])
    #charged fraction, sum of charged hadrons, electrons and muons
    t_var_jet1_sj1_chFrac = array('f',[0])
    #closest parton 0 for q- (e.g. b)  from H, 1 for q+ (e.g. bbar) from H, 2 for q- from Z, 3 for q+ from Z, only for signal, else set to -1 
    #it can happen one parton is matched to TWO subjets
    t_var_jet1_sj1_closestMatch = array('i',[0])
    #decision for parton 0 for q- from H, 1 for q+ from H, 2 for q- from Z, 3 for q+ from Z, here one parton is matched to one subjet 
    t_var_jet1_sj1_decMatch = array('i',[0])
    t_var_jet1_sj1_Angle_closestMatch = array('f',[0])
    t_var_jet1_sj1_Angle_decMatch = array('f',[0])

    t_var_jet1_sj2_E = array('f',[0])
    t_var_jet1_sj2_Px = array('f',[0])
    t_var_jet1_sj2_Py = array('f',[0])
    t_var_jet1_sj2_Pz = array('f',[0])
    t_var_jet1_sj2_jetChargePt_kappa_0_30 = array('f',[0])
    t_var_jet1_sj2_jetChargeE_kappa_0_30 = array('f',[0])
    t_var_jet1_sj2_jetChargePt_kappa_0_25 = array('f',[0])
    t_var_jet1_sj2_jetChargeE_kappa_0_25 = array('f',[0])
    t_var_jet1_sj2_nTracks = array('i',[0])
    t_var_jet1_sj2_chFrac = array('f',[0])
    #closest parton 0 for q- from H, 1 for q+ from H, 2 for q- from Z, 3 for q+ from Z, only for signal, else set to -1  
    #it can happen one parton is matched to TWO subjets
    t_var_jet1_sj2_closestMatch = array('i',[0])
    #decision for parton 0 for q- from H, 1 for q+ from H, 2 for q- from Z, 3 for q+ from Z, here one parton is matched to one subjet 
    t_var_jet1_sj2_decMatch = array('i',[0])
    t_var_jet1_sj2_Angle_closestMatch = array('f',[0])
    t_var_jet1_sj2_Angle_decMatch = array('f',[0])
    t_var_jet1_dAlpha_sj1sj2 = array('f',[0])

    t_var_jet2_sj1_E = array('f',[0])
    t_var_jet2_sj1_Px = array('f',[0])
    t_var_jet2_sj1_Py = array('f',[0])
    t_var_jet2_sj1_Pz = array('f',[0])
    t_var_jet2_sj1_jetChargePt_kappa_0_30 = array('f',[0])
    t_var_jet2_sj1_jetChargeE_kappa_0_30 = array('f',[0])
    t_var_jet2_sj1_jetChargePt_kappa_0_25 = array('f',[0])
    t_var_jet2_sj1_jetChargeE_kappa_0_25 = array('f',[0])
    t_var_jet2_sj1_nTracks = array('i',[0])
    t_var_jet2_sj1_chFrac = array('f',[0])
    #closest parton 0 for q- from H, 1 for q+ from H, 2 for q- from Z, 3 for q+ from Z, only for signal, else set to -1  
    #it can happen one parton is matched to TWO subjets
    t_var_jet2_sj1_closestMatch = array('i',[0])
    #decision for parton 0 for q- from H, 1 for q+ from H, 2 for q- from Z, 3 for q+ from Z, here one parton is matched to one subjet 
    t_var_jet2_sj1_decMatch = array('i',[0])
    t_var_jet2_sj1_Angle_closestMatch = array('f',[0])
    t_var_jet2_sj1_Angle_decMatch = array('f',[0])

    t_var_jet2_sj2_E = array('f',[0])
    t_var_jet2_sj2_Px = array('f',[0])
    t_var_jet2_sj2_Py = array('f',[0])
    t_var_jet2_sj2_Pz = array('f',[0])
    t_var_jet2_sj2_jetChargePt_kappa_0_30 = array('f',[0])
    t_var_jet2_sj2_jetChargeE_kappa_0_30 = array('f',[0])
    t_var_jet2_sj2_jetChargePt_kappa_0_25 = array('f',[0])
    t_var_jet2_sj2_jetChargeE_kappa_0_25 = array('f',[0])
    t_var_jet2_sj2_nTracks = array('i',[0])
    t_var_jet2_sj2_chFrac = array('f',[0])
    #closest parton 0 for q-(e.g. b) from H, 1 for q+ (bbar) from H, 2 for q- from Z, 3 for q+ from Z, only for signal, else set to -1  
    #it can happen one parton is matched to TWO subjets
    t_var_jet2_sj2_closestMatch = array('i',[0])
    #decision for parton 0 for q- from H, 1 for q+ from H, 2 for q- from Z, 3 for q+ from Z, here one parton is matched to one subjet 
    t_var_jet2_sj2_decMatch = array('i',[0])
    t_var_jet2_sj2_Angle_closestMatch = array('f',[0])
    t_var_jet2_sj2_Angle_decMatch = array('f',[0])
    t_var_jet2_dAlpha_sj1sj2 = array('f',[0])

    t_var_costheta1_for_Atheta1 = array('f',[0])
    t_var_costheta2_for_Atheta1theta2 = array('f',[0])
    t_var_phi_for_Aphis = array('f',[0])


    t_var_genInv_E = array('f',[0])
    t_var_genInv_Px = array('f',[0])
    t_var_genInv_Py = array('f',[0])
    t_var_genInv_Pz = array('f',[0])

    t_var_gen_y21 = array('f',[0])
    t_var_gen_y32 = array('f',[0])
    t_var_gen_y43 = array('f',[0])
  
    t_var_genjet1_d21 = array('f',[0])
    t_var_genjet1_d32 = array('f',[0])
    t_var_genjet1_d43 = array('f',[0])
    t_var_genjet2_d21 = array('f',[0])
    t_var_genjet2_d32 = array('f',[0])
    t_var_genjet2_d43 = array('f',[0])

    t_var_genjet1_sj1_E = array('f',[0])
    t_var_genjet1_sj1_Px = array('f',[0])
    t_var_genjet1_sj1_Py = array('f',[0])
    t_var_genjet1_sj1_Pz = array('f',[0])
    t_var_genjet1_sj1_jetChargePt_kappa_0_30 = array('f',[0])
    t_var_genjet1_sj1_jetChargeE_kappa_0_30 = array('f',[0])
    t_var_genjet1_sj1_jetChargePt_kappa_0_25 = array('f',[0])
    t_var_genjet1_sj1_jetChargeE_kappa_0_25 = array('f',[0])
    t_var_genjet1_sj1_nTracks = array('i',[0])
    t_var_genjet1_sj1_chFrac = array('f',[0])

    t_var_genjet1_sj2_E = array('f',[0])
    t_var_genjet1_sj2_Px = array('f',[0])
    t_var_genjet1_sj2_Py = array('f',[0])
    t_var_genjet1_sj2_Pz = array('f',[0])
    t_var_genjet1_sj2_jetChargePt_kappa_0_30 = array('f',[0])
    t_var_genjet1_sj2_jetChargeE_kappa_0_30 = array('f',[0])
    t_var_genjet1_sj2_jetChargePt_kappa_0_25 = array('f',[0])
    t_var_genjet1_sj2_jetChargeE_kappa_0_25 = array('f',[0])
    t_var_genjet1_sj2_nTracks = array('i',[0])
    t_var_genjet1_sj2_chFrac = array('f',[0])
    t_var_genjet1_dAlpha_sj1sj2 = array('f',[0])

    t_var_genjet1_dAlpha_sj1_rjsj1 = array('f',[0])
    t_var_genjet1_dAlpha_sj1_rjsj2 = array('f',[0])
    t_var_genjet1_dAlpha_sj1_qmin = array('f',[0])
    t_var_genjet1_dAlpha_sj1_qplus = array('f',[0])
    t_var_genjet1_dAlpha_H = array('f',[0])

    t_var_parton_ep_E = array('f',[0])
    t_var_parton_ep_Px = array('f',[0])
    t_var_parton_ep_Pz = array('f',[0])
    t_var_parton_ep_Py = array('f',[0])

    t_var_parton_em_E = array('f',[0])
    t_var_parton_em_Px = array('f',[0])
    t_var_parton_em_Pz = array('f',[0])
    t_var_parton_em_Py = array('f',[0])

    t_var_parton_H_E = array('f',[0])
    t_var_parton_H_Px = array('f',[0])
    t_var_parton_H_Pz = array('f',[0])
    t_var_parton_H_Py = array('f',[0])
    t_var_parton_H_PDG_Daughter0 = array('i',[0])

    t_var_parton_Z_qpos_E = array('f',[0])
    t_var_parton_Z_qpos_Px = array('f',[0])
    t_var_parton_Z_qpos_Pz = array('f',[0])
    t_var_parton_Z_qpos_Py = array('f',[0])
    t_var_parton_Z_qpos_PDGID = array('i',[0])

    t_var_parton_Z_qneg_E = array('f',[0])
    t_var_parton_Z_qneg_Px = array('f',[0])
    t_var_parton_Z_qneg_Pz = array('f',[0])
    t_var_parton_Z_qneg_Py = array('f',[0])
    t_var_parton_Z_qneg_PDGID = array('i',[0])

    t_var_genjet2_sj1_E = array('f',[0])
    t_var_genjet2_sj1_Px = array('f',[0])
    t_var_genjet2_sj1_Py = array('f',[0])
    t_var_genjet2_sj1_Pz = array('f',[0])
    t_var_genjet2_sj1_jetChargePt_kappa_0_30 = array('f',[0])
    t_var_genjet2_sj1_jetChargeE_kappa_0_30 = array('f',[0])
    t_var_genjet2_sj1_jetChargePt_kappa_0_25 = array('f',[0])
    t_var_genjet2_sj1_jetChargeE_kappa_0_25 = array('f',[0])
    t_var_genjet2_sj1_nTracks = array('i',[0])
    t_var_genjet2_sj1_chFrac = array('f',[0])


    t_var_genjet2_sj2_E = array('f',[0])
    t_var_genjet2_sj2_Px = array('f',[0])
    t_var_genjet2_sj2_Py = array('f',[0])
    t_var_genjet2_sj2_Pz = array('f',[0])
    t_var_genjet2_sj2_jetChargePt_kappa_0_30 = array('f',[0])
    t_var_genjet2_sj2_jetChargeE_kappa_0_30 = array('f',[0])
    t_var_genjet2_sj2_jetChargePt_kappa_0_25 = array('f',[0])
    t_var_genjet2_sj2_jetChargeE_kappa_0_25 = array('f',[0])
    t_var_genjet2_sj2_nTracks = array('i',[0])
    t_var_genjet2_sj2_chFrac = array('f',[0])
    t_var_genjet2_dAlpha_sj1sj2 = array('f',[0])

    t_var_genjet2_dAlpha_sj1_rjsj1 = array('f',[0])
    t_var_genjet2_dAlpha_sj1_rjsj2 = array('f',[0])
    t_var_genjet2_dAlpha_sj1_qmin = array('f',[0])
    t_var_genjet2_dAlpha_sj1_qplus = array('f',[0])
    t_var_genjet2_dAlpha_H = array('f',[0])
    t_var_GenLevelFilled = array('i',[0])
    t_var_eventWeight = array('f',[0])

    mytree.Branch('sqrtS_j1_j2_EMiss', t_var_sqrtS , 'sqrtS_j1_j2_EMiss/F')
    mytree.Branch('sqrtS_j1_j2_orig', t_var_sqrtS_orig , 'sqrtS_j1_j2_orig/F')
    mytree.Branch('sqrtS_j1_j2_EMiss_gen', t_var_sqrtS_gen , 'sqrtS_j1_j2_EMiss_gen/F')
    mytree.Branch('sqrtS_j1_j2_orig_gen', t_var_sqrtS_orig_gen , 'sqrtS_j1_j2_orig_gen/F')
    mytree.Branch('sqrtS_parton', t_var_sqrtS_parton , 'sqrtS_parton/F')
    mytree.Branch('nLeptons', t_var_nLeptons , 'nLeptons/I')
    mytree.Branch('nLeptons_gen', t_var_nLeptons_gen , 'nLeptons_gen/I')

    mytree.Branch('MET', t_var_MET , 'MET/F')
    mytree.Branch('reco_y21', t_var_reco_y21 , 'reco_y21/F')
    mytree.Branch('reco_y32', t_var_reco_y32 , 'reco_y32/F')
    mytree.Branch('reco_y43', t_var_reco_y43 , 'reco_y43/F')
    mytree.Branch('jet1_mass', t_var_jet1_mass , 'jet1_mass/F')
    mytree.Branch('jet2_mass', t_var_jet2_mass , 'jet2_mass/F')
    mytree.Branch('jet1_min_jet2_mass', t_var_jet1_min_jet2_mass , 'jet1_min_jet2_mass/F')
    mytree.Branch('jet1_BTag_rfj_BTagMax', t_var_jet1_BTag_rfj_BTagMax , 'jet1_BTag_rfj_BTagMax/F')
    mytree.Branch('jet1_CTag_rfj_BTagMax', t_var_jet1_CTag_rfj_BTagMax , 'jet1_BTag_rfj_CTagMax/F')
    mytree.Branch('jet1_LTag_rfj_BTagMax', t_var_jet1_LTag_rfj_BTagMax , 'jet1_BTag_rfj_LTagMax/F')
    mytree.Branch('jet1_BTag_rfj_CTagMax', t_var_jet1_BTag_rfj_CTagMax , 'jet1_BTag_rfj_CTagMax/F')
    mytree.Branch('jet1_CTag_rfj_CTagMax', t_var_jet1_CTag_rfj_CTagMax , 'jet1_BTag_rfj_CTagMax/F')
    mytree.Branch('jet1_LTag_rfj_CTagMax', t_var_jet1_LTag_rfj_CTagMax , 'jet1_BTag_rfj_LTagMax/F')
    mytree.Branch('jet1_E', t_var_jet1_E , 'jet1_E/F')
    mytree.Branch('jet2_E', t_var_jet2_E , 'jet2_E/F')
    mytree.Branch('jet1_Pt', t_var_jet1_Pt , 'jet1_Pt/F')
    mytree.Branch('jet2_Pt', t_var_jet2_Pt , 'jet2_Pt/F')
    mytree.Branch('jet1_theta', t_var_jet1_theta , 'jet1_theta/F')
    mytree.Branch('jet2_theta', t_var_jet2_theta , 'jet2_theta/F')
    mytree.Branch('jet1_min_jet2_theta', t_var_jet1_min_jet2_theta , 'jet1_min_jet2_theta/F')
    mytree.Branch('jet1_phi', t_var_jet1_phi , 'jet1_phi/F')
    mytree.Branch('jet2_phi', t_var_jet2_phi , 'jet2_phi/F')
    mytree.Branch('dphi_j1j2', t_var_dphi_j1j2 , 'dphi_j1j2/F')
    mytree.Branch('angle_j1j2', t_var_angle_j1j2 , 'angle_j1j2/F')
    mytree.Branch('jet1_D2_beta1', t_var_jet1_D2_beta1 , 'jet1_D2_beta1/F')
    mytree.Branch('jet2_D2_beta1', t_var_jet2_D2_beta1 , 'jet2_D2_beta1/F')
    mytree.Branch('jet1_D2_beta0_5', t_var_jet1_D2_beta0_5 , 'jet1_D2_beta0_5/F')
    mytree.Branch('jet2_D2_beta0_5', t_var_jet2_D2_beta0_5 , 'jet2_D2_beta0_5/F')
    mytree.Branch('jet1_C2_beta1', t_var_jet1_C2_beta1 , 'jet1_C2_beta1/F')
    mytree.Branch('jet2_C2_beta1', t_var_jet2_C2_beta1 , 'jet2_C2_beta1/F')
    mytree.Branch('jet1_C2_beta0_5', t_var_jet1_C2_beta0_5 , 'jet1_C2_beta0_5/F')
    mytree.Branch('jet2_C2_beta0_5', t_var_jet2_C2_beta0_5 , 'jet2_C2_beta0_5/F')
    mytree.Branch('jet1_tau21', t_var_jet1_tau21 , 'jet1_tau21/F')
    mytree.Branch('jet2_tau21', t_var_jet2_tau21 , 'jet2_tau21/F')
    mytree.Branch('jet1_tau32', t_var_jet1_tau32 , 'jet1_tau32/F')
    mytree.Branch('jet2_tau32', t_var_jet2_tau32 , 'jet2_tau32/F')

    mytree.Branch('jet1_C3_beta1', t_var_jet1_C3_beta1 , 'jet1_C3_beta1/F')
    mytree.Branch('jet2_C3_beta1', t_var_jet2_C3_beta1 , 'jet2_C3_beta1/F')
    mytree.Branch('jet1_N2_beta1', t_var_jet1_N2_beta1 , 'jet1_N2_beta1/F')
    mytree.Branch('jet2_N2_beta1', t_var_jet2_N2_beta1 , 'jet2_N2_beta1/F')
    mytree.Branch('jet1_N3_beta1', t_var_jet1_N3_beta1 , 'jet1_N3_beta1/F')
    mytree.Branch('jet2_N3_beta1', t_var_jet2_N3_beta1 , 'jet2_N3_beta1/F')

    mytree.Branch('jet1_D2_beta1_E_theta', t_var_jet1_D2_beta1_E_theta , 'jet1_D2_beta1_E_theta/F')
    mytree.Branch('jet2_D2_beta1_E_theta', t_var_jet2_D2_beta1_E_theta , 'jet2_D2_beta1_E_theta/F')
    mytree.Branch('jet1_C2_beta1_E_theta', t_var_jet1_C2_beta1_E_theta , 'jet1_C2_beta1_E_theta/F')
    mytree.Branch('jet2_C2_beta1_E_theta', t_var_jet2_C2_beta1_E_theta , 'jet2_C2_beta1_E_theta/F')
    mytree.Branch('jet1_tau21_E_theta', t_var_jet1_tau21_E_theta , 'jet1_tau21/F')
    mytree.Branch('jet2_tau21_E_theta', t_var_jet2_tau21_E_theta , 'jet2_tau21/F')
    mytree.Branch('jet1_tau32_E_theta', t_var_jet1_tau32_E_theta , 'jet1_tau32_E_theta/F')
    mytree.Branch('jet2_tau32_E_theta', t_var_jet2_tau32_E_theta , 'jet2_tau32_E_theta/F')
    mytree.Branch('jet1_C3_beta1_E_theta', t_var_jet1_C3_beta1_E_theta , 'jet1_C3_beta1_E_theta/F')
    mytree.Branch('jet2_C3_beta1_E_theta', t_var_jet2_C3_beta1_E_theta , 'jet2_C3_beta1_E_theta/F')
    mytree.Branch('jet1_N2_beta1_E_theta', t_var_jet1_N2_beta1_E_theta , 'jet1_N2_beta1_E_theta/F')
    mytree.Branch('jet2_N2_beta1_E_theta', t_var_jet2_N2_beta1_E_theta , 'jet2_N2_beta1_E_theta/F')
    mytree.Branch('jet1_N3_beta1_E_theta', t_var_jet1_N3_beta1_E_theta , 'jet1_N3_beta1_E_theta/F')
    mytree.Branch('jet2_N3_beta1_E_theta', t_var_jet2_N3_beta1_E_theta , 'jet2_N3_beta1_E_theta/F')

    #order subjets by energy
    mytree.Branch('jet1_sj1_E', t_var_jet1_sj1_E , 'jet1_sj1_E/F')
    mytree.Branch('jet1_sj1_Px', t_var_jet1_sj1_Px , 'jet1_sj1_Px/F')
    mytree.Branch('jet1_sj1_Py', t_var_jet1_sj1_Py , 'jet1_sj1_Py/F')
    mytree.Branch('jet1_sj1_Pz', t_var_jet1_sj1_Pz , 'jet1_sj1_Pz/F')
    mytree.Branch('jet1_sj1_jetChargePt_kappa_0_30', t_var_jet1_sj1_jetChargePt_kappa_0_30 , 'jet1_sj1_jetChargePt_kappa_0_30/F')
    mytree.Branch('jet1_sj1_jetChargeE_kappa_0_30', t_var_jet1_sj1_jetChargeE_kappa_0_30 , 'jet1_sj1_jetChargeE_kappa_0_30/F')
    mytree.Branch('jet1_sj1_jetChargePt_kappa_0_25', t_var_jet1_sj1_jetChargePt_kappa_0_25 , 'jet1_sj1_jetChargePt_kappa_0_25/F')
    mytree.Branch('jet1_sj1_jetChargeE_kappa_0_25', t_var_jet1_sj1_jetChargeE_kappa_0_25 , 'jet1_sj1_jetChargeE_kappa_0_25/F')
    mytree.Branch('jet1_sj1_chFrac', t_var_jet1_sj1_chFrac, 'jet1_sj1_chFrac/F')
    mytree.Branch('jet1_sj1_nTracks', t_var_jet1_sj1_nTracks , 'jet1_sj1_nTracks/I')
    mytree.Branch('jet1_sj1_closestMatch', t_var_jet1_sj1_closestMatch, 'jet1_sj1_closestMatch/I')
    mytree.Branch('jet1_sj1_decMatch', t_var_jet1_sj1_decMatch, 'jet1_sj1_decMatch/I')
    mytree.Branch('jet1_sj1_Angle_closestMatch', t_var_jet1_sj1_Angle_closestMatch, 'jet1_sj1_Angle_closestMatch/F')
    mytree.Branch('jet1_sj1_Angle_decMatch', t_var_jet1_sj1_Angle_decMatch, 'jet1_sj1_Angle_decMatch/F')
    mytree.Branch('jet1_sj2_E', t_var_jet1_sj2_E , 'jet1_sj2_E/F')
    mytree.Branch('jet1_sj2_Px', t_var_jet1_sj2_Px , 'jet1_sj2_Px/F')
    mytree.Branch('jet1_sj2_Py', t_var_jet1_sj2_Py , 'jet1_sj2_Py/F')
    mytree.Branch('jet1_sj2_Pz', t_var_jet1_sj2_Pz , 'jet1_sj2_Pz/F')
    mytree.Branch('jet1_sj2_jetChargePt_kappa_0_30', t_var_jet1_sj2_jetChargePt_kappa_0_30 , 'jet1_sj2_jetChargePt_kappa_0_30/F')
    mytree.Branch('jet1_sj2_jetChargeE_kappa_0_30', t_var_jet1_sj2_jetChargeE_kappa_0_30 , 'jet1_sj2_jetChargeE_kappa_0_30/F')
    mytree.Branch('jet1_sj2_jetChargePt_kappa_0_25', t_var_jet1_sj2_jetChargePt_kappa_0_25 , 'jet1_sj2_jetChargePt_kappa_0_25/F')
    mytree.Branch('jet1_sj2_jetChargeE_kappa_0_25', t_var_jet1_sj2_jetChargeE_kappa_0_25 , 'jet1_sj2_jetChargeE_kappa_0_25/F')
    mytree.Branch('jet1_sj2_nTracks', t_var_jet1_sj2_nTracks , 'jet1_sj2_nTracks/I')
    mytree.Branch('jet1_sj2_chFrac', t_var_jet1_sj2_chFrac, 'jet1_sj2_chFrac/F')
    mytree.Branch('jet1_sj2_closestMatch', t_var_jet1_sj2_closestMatch, 'jet1_sj2_closestMatch/I')
    mytree.Branch('jet1_sj2_decMatch', t_var_jet1_sj2_decMatch, 'jet1_sj2_decMatch/I')
    mytree.Branch('jet1_sj2_Angle_closestMatch', t_var_jet1_sj2_Angle_closestMatch, 'jet1_sj2_Angle_closestMatch/F')
    mytree.Branch('jet1_sj2_Angle_decMatch', t_var_jet1_sj2_Angle_decMatch, 'jet1_sj2_Angle_decMatch/F')
    mytree.Branch('jet1_dAlpha_sj1sj2',t_var_jet1_dAlpha_sj1sj2,'jet1_dAlpha_sj1sj2/F')
    mytree.Branch('jet1_d21', t_var_jet1_d21, 'jet1_d21/F');
    mytree.Branch('jet1_d32', t_var_jet1_d32, 'jet1_d32/F');    
    mytree.Branch('jet1_d43', t_var_jet1_d43, 'jet1_d43/F');
    mytree.Branch('jet2_sj1_E', t_var_jet2_sj1_E , 'jet2_sj1_E/F')
    mytree.Branch('jet2_sj1_Px', t_var_jet2_sj1_Px , 'jet2_sj1_Px/F')
    mytree.Branch('jet2_sj1_Py', t_var_jet2_sj1_Py , 'jet2_sj1_Py/F')
    mytree.Branch('jet2_sj1_Pz', t_var_jet2_sj1_Pz , 'jet2_sj1_Pz/F')
    mytree.Branch('jet2_sj1_jetChargePt_kappa_0_30', t_var_jet2_sj1_jetChargePt_kappa_0_30 , 'jet2_sj1_jetChargePt_kappa_0_30/F')
    mytree.Branch('jet2_sj1_jetChargeE_kappa_0_30', t_var_jet2_sj1_jetChargeE_kappa_0_30 , 'jet2_sj1_jetChargeE_kappa_0_30/F')
    mytree.Branch('jet2_sj1_jetChargePt_kappa_0_25', t_var_jet2_sj1_jetChargePt_kappa_0_25 , 'jet2_sj1_jetChargePt_kappa_0_25/F')
    mytree.Branch('jet2_sj1_jetChargeE_kappa_0_25', t_var_jet2_sj1_jetChargeE_kappa_0_25 , 'jet2_sj1_jetChargeE_kappa_0_25/F')
    mytree.Branch('jet2_sj1_nTracks', t_var_jet2_sj1_nTracks , 'jet2_sj1_nTracks/I')
    mytree.Branch('jet2_sj1_chFrac', t_var_jet2_sj1_chFrac, 'jet2_sj1_chFrac/F')
    mytree.Branch('jet2_sj1_closestMatch', t_var_jet2_sj1_closestMatch, 'jet2_sj1_closestMatch/I')
    mytree.Branch('jet2_sj1_decMatch', t_var_jet2_sj1_decMatch, 'jet2_sj1_decMatch/I')
    mytree.Branch('jet2_sj1_Angle_closestMatch', t_var_jet2_sj1_Angle_closestMatch, 'jet2_sj1_Angle_closestMatch/F')
    mytree.Branch('jet2_sj1_Angle_decMatch', t_var_jet2_sj1_Angle_decMatch, 'jet2_sj1_Angle_decMatch/F')
    mytree.Branch('jet2_sj2_E', t_var_jet2_sj2_E , 'jet2_sj2_E/F')
    mytree.Branch('jet2_sj2_Px', t_var_jet2_sj2_Px , 'jet2_sj2_Px/F')
    mytree.Branch('jet2_sj2_Py', t_var_jet2_sj2_Py , 'jet2_sj2_Py/F')
    mytree.Branch('jet2_sj2_Pz', t_var_jet2_sj2_Pz , 'jet2_sj2_Pz/F')
    mytree.Branch('jet2_sj2_jetChargePt_kappa_0_30', t_var_jet2_sj2_jetChargePt_kappa_0_30 , 'jet2_sj2_jetChargePt_kappa_0_30/F')
    mytree.Branch('jet2_sj2_jetChargeE_kappa_0_30', t_var_jet2_sj2_jetChargeE_kappa_0_30 , 'jet2_sj2_jetChargeE_kappa_0_30/F')
    mytree.Branch('jet2_sj2_jetChargePt_kappa_0_25', t_var_jet2_sj2_jetChargePt_kappa_0_25 , 'jet2_sj2_jetChargePt_kappa_0_25/F')
    mytree.Branch('jet2_sj2_jetChargeE_kappa_0_25', t_var_jet2_sj2_jetChargeE_kappa_0_25 , 'jet2_sj2_jetChargeE_kappa_0_25/F')
    mytree.Branch('jet2_sj2_nTracks', t_var_jet2_sj2_nTracks , 'jet2_sj2_nTracks/I')
    mytree.Branch('jet2_sj2_chFrac', t_var_jet2_sj2_chFrac, 'jet2_sj2_chFrac/F')
    mytree.Branch('jet2_sj2_closestMatch', t_var_jet2_sj2_closestMatch, 'jet2_sj2_closestMatch/I')
    mytree.Branch('jet2_sj2_decMatch', t_var_jet2_sj2_decMatch, 'jet2_sj2_decMatch/I')
    mytree.Branch('jet2_sj2_Angle_closestMatch', t_var_jet2_sj2_Angle_closestMatch, 'jet2_sj2_Angle_closestMatch/F')
    mytree.Branch('jet2_sj2_Angle_decMatch', t_var_jet2_sj2_Angle_decMatch, 'jet2_sj2_Angle_decMatch/F')
    mytree.Branch('jet2_dAlpha_sj1sj2',t_var_jet2_dAlpha_sj1sj2,'jet2_dAlpha_sj1sj2/F')
    mytree.Branch('jet2_d21', t_var_jet2_d21, 'jet2_d21/F');
    mytree.Branch('jet2_d32', t_var_jet2_d32, 'jet2_d32/F');    
    mytree.Branch('jet2_d43', t_var_jet2_d43, 'jet2_d43/F');
    mytree.Branch('costheta1_for_Atheta1',t_var_costheta1_for_Atheta1,'costheta1_for_Atheta1/F')
    mytree.Branch('costheta2_for_Atheta1theta2',t_var_costheta2_for_Atheta1theta2,'costheta2_for_Atheta1theta2/F')
    mytree.Branch('phi_for_Aphis',t_var_phi_for_Aphis,'phi_for_Aphis/F')

    mytree.Branch('GenLevelFilled',t_var_GenLevelFilled,'GenLevelFilled/I')
    mytree.Branch('eventWeight', t_var_eventWeight , 'eventWeight/F')

    if fillGenLevel:
        mytree.Branch('gen_y21', t_var_gen_y21 , 'gen_y21/F')
        mytree.Branch('gen_y32', t_var_gen_y32 , 'gen_y32/F')
        mytree.Branch('gen_y43', t_var_gen_y43 , 'gen_y43/F')
        mytree.Branch('genInv_E', t_var_genInv_E , 'genInv_E/F')
        mytree.Branch('genInv_Px', t_var_genInv_Px , 'genInv_Px/F')
        mytree.Branch('genInv_Py', t_var_genInv_Py , 'genInv_Py/F')
        mytree.Branch('genInv_Pz', t_var_genInv_Pz , 'genInv_Pz/F')
        #order here genjets with respect to recojet --> aka genjet 1 is the one matched to recojet 1
        mytree.Branch('genjet1_sj1_E', t_var_genjet1_sj1_E , 'genjet1_sj1_E/F')
        mytree.Branch('genjet1_sj1_Px', t_var_genjet1_sj1_Px , 'genjet1_sj1_Px/F')
        mytree.Branch('genjet1_sj1_Py', t_var_genjet1_sj1_Py , 'genjet1_sj1_Py/F')
        mytree.Branch('genjet1_sj1_Pz', t_var_genjet1_sj1_Pz , 'genjet1_sj1_Pz/F')
        mytree.Branch('genjet1_sj1_jetChargePt_kappa_0_30', t_var_genjet1_sj1_jetChargePt_kappa_0_30 , 'genjet1_sj1_jetChargePt_kappa_0_30/F')
        mytree.Branch('genjet1_sj1_jetChargeE_kappa_0_30', t_var_genjet1_sj1_jetChargeE_kappa_0_30 , 'genjet1_sj1_jetChargeE_kappa_0_30/F')
        mytree.Branch('genjet1_sj1_jetChargePt_kappa_0_25', t_var_genjet1_sj1_jetChargePt_kappa_0_25 , 'genjet1_sj1_jetChargePt_kappa_0_25/F')
        mytree.Branch('genjet1_sj1_jetChargeE_kappa_0_25', t_var_genjet1_sj1_jetChargeE_kappa_0_25 , 'genjet1_sj1_jetChargeE_kappa_0_25/F')
        mytree.Branch('genjet1_sj1_nTracks', t_var_genjet1_sj1_nTracks , 'genjet1_sj1_nTracks/I')
        mytree.Branch('genjet1_sj1_chFrac', t_var_genjet1_sj1_chFrac, 'genjet1_sj1_chFrac/F')
        #rjsj1 sj1 (in E) of closest rj to gj
        mytree.Branch('genjet1_dAlpha_sj1_rjsj1',t_var_genjet1_dAlpha_sj1_rjsj1,'genjet1_dAlpha_sj1_rjsj1/F')
        mytree.Branch('genjet1_dAlpha_sj1_rjsj2',t_var_genjet1_dAlpha_sj1_rjsj2,'genjet1_dAlpha_sj1_rjsj2/F')
        #closest negative parton to to gj-sj (independent if H or Z jet)
        mytree.Branch('genjet1_dAlpha_sj1_qmin',t_var_genjet1_dAlpha_sj1_qmin,'genjet1_dAlpha_sj1_qmin/F')
        mytree.Branch('genjet1_dAlpha_sj1_qplus',t_var_genjet1_dAlpha_sj1_qplus,'genjet1_dAlpha_sj1_qplus/F')
        mytree.Branch('genjet1_dAlpha_H',t_var_genjet1_dAlpha_H,'genjet1_dAlpha_H/F')
        mytree.Branch('genjet1_sj2_E', t_var_genjet1_sj2_E , 'genjet1_sj2_E/F')
        mytree.Branch('genjet1_sj2_Px', t_var_genjet1_sj2_Px , 'genjet1_sj2_Px/F')
        mytree.Branch('genjet1_sj2_Py', t_var_genjet1_sj2_Py , 'genjet1_sj2_Py/F')
        mytree.Branch('genjet1_sj2_Pz', t_var_genjet1_sj2_Pz , 'genjet1_sj2_Pz/F')
        mytree.Branch('genjet1_sj2_jetChargePt_kappa_0_30', t_var_genjet1_sj2_jetChargePt_kappa_0_30 , 'genjet1_sj2_jetChargePt_kappa_0_30/F')
        mytree.Branch('genjet1_sj2_jetChargeE_kappa_0_30', t_var_genjet1_sj2_jetChargeE_kappa_0_30 , 'genjet1_sj2_jetChargeE_kappa_0_30/F')
        mytree.Branch('genjet1_sj2_jetChargePt_kappa_0_25', t_var_genjet1_sj2_jetChargePt_kappa_0_25 , 'genjet1_sj2_jetChargePt_kappa_0_25/F')
        mytree.Branch('genjet1_sj2_jetChargeE_kappa_0_25', t_var_genjet1_sj2_jetChargeE_kappa_0_25 , 'genjet1_sj2_jetChargeE_kappa_0_25/F')
        mytree.Branch('genjet1_sj2_nTracks', t_var_genjet1_sj2_nTracks , 'genjet1_sj2_nTracks/I')
        mytree.Branch('genjet1_sj2_chFrac', t_var_genjet1_sj2_chFrac, 'genjet1_sj2_chFrac/F')
        mytree.Branch('genjet1_dAlpha_sj1sj2',t_var_genjet1_dAlpha_sj1sj2,'genjet1_dAlpha_sj1sj2/F')
        mytree.Branch('genjet1_d21', t_var_genjet1_d21, 'genjet1_d21/F');
        mytree.Branch('genjet1_d32', t_var_genjet1_d32, 'genjet1_d32/F');    
        mytree.Branch('genjet1_d43', t_var_genjet1_d43, 'genjet1_d43/F');
        mytree.Branch('genjet2_sj1_E', t_var_genjet2_sj1_E , 'genjet2_sj1_E/F')
        mytree.Branch('genjet2_sj1_Px', t_var_genjet2_sj1_Px , 'genjet2_sj1_Px/F')
        mytree.Branch('genjet2_sj1_Py', t_var_genjet2_sj1_Py , 'genjet2_sj1_Py/F')
        mytree.Branch('genjet2_sj1_Pz', t_var_genjet2_sj1_Pz , 'genjet2_sj1_Pz/F')
        mytree.Branch('genjet2_sj1_jetChargePt_kappa_0_30', t_var_genjet2_sj1_jetChargePt_kappa_0_30 , 'genjet2_sj1_jetChargePt_kappa_0_30/F')
        mytree.Branch('genjet2_sj1_jetChargeE_kappa_0_30', t_var_genjet2_sj1_jetChargeE_kappa_0_30 , 'genjet2_sj1_jetChargeE_kappa_0_30/F')
        mytree.Branch('genjet2_sj1_jetChargePt_kappa_0_25', t_var_genjet2_sj1_jetChargePt_kappa_0_25 , 'genjet2_sj1_jetChargePt_kappa_0_25/F')
        mytree.Branch('genjet2_sj1_jetChargePt_kappa_0_25', t_var_genjet2_sj1_jetChargeE_kappa_0_25 , 'genjet2_sj1_jetChargeE_kappa_0_25/F')
        mytree.Branch('genjet2_sj1_nTracks', t_var_genjet2_sj1_nTracks , 'genjet2_sj1_nTracks/I')
        mytree.Branch('genjet2_sj1_chFrac', t_var_genjet2_sj1_chFrac, 'genjet2_sj1_chFrac/F')
        mytree.Branch('genjet2_sj2_E', t_var_genjet2_sj2_E , 'genjet2_sj2_E/F')
        mytree.Branch('genjet2_sj2_Px', t_var_genjet2_sj2_Px , 'genjet2_sj2_Px/F')
        mytree.Branch('genjet2_sj2_Py', t_var_genjet2_sj2_Py , 'genjet2_sj2_Py/F')
        mytree.Branch('genjet2_sj2_Pz', t_var_genjet2_sj2_Pz , 'genjet2_sj2_Pz/F')
        mytree.Branch('genjet2_sj2_jetChargePt_kappa_0_30', t_var_genjet2_sj2_jetChargePt_kappa_0_30 , 'genjet2_sj2_jetChargePt_kappa_0_30/F')
        mytree.Branch('genjet2_sj2_jetChargeE_kappa_0_30', t_var_genjet2_sj2_jetChargeE_kappa_0_30 , 'genjet2_sj2_jetChargeE_kappa_0_30/F')
        mytree.Branch('genjet2_sj2_jetChargePt_kappa_0_25', t_var_genjet2_sj2_jetChargePt_kappa_0_25 , 'genjet2_sj2_jetChargePt_kappa_0_25/F')
        mytree.Branch('genjet2_sj2_jetChargeE_kappa_0_25', t_var_genjet2_sj2_jetChargeE_kappa_0_25 , 'genjet2_sj2_jetChargeE_kappa_0_25/F')
        mytree.Branch('genjet2_sj2_nTracks', t_var_genjet2_sj2_nTracks , 'genjet2_sj2_nTracks/I')
        mytree.Branch('genjet2_sj2_chFrac', t_var_genjet2_sj2_chFrac, 'genjet2_sj2_chFrac/F')
        mytree.Branch('genjet2_dAlpha_sj1sj2',t_var_genjet2_dAlpha_sj1sj2,'genjet2_dAlpha_sj1sj2/F')
        #rjsj1 sj1 (in E) of closest rj to gj
        mytree.Branch('genjet2_dAlpha_sj1_rjsj1',t_var_genjet2_dAlpha_sj1_rjsj1,'genjet2_dAlpha_sj1_rjsj1/F')
        mytree.Branch('genjet2_dAlpha_sj1_rjsj2',t_var_genjet2_dAlpha_sj1_rjsj2,'genjet2_dAlpha_sj1_rjsj2/F')
        #closest negative parton to to gj-sj (independent if H or Z jet)
        mytree.Branch('genjet2_dAlpha_sj1_qmin',t_var_genjet2_dAlpha_sj1_qmin,'genjet2_dAlpha_sj1_qmin/F')
        mytree.Branch('genjet2_dAlpha_sj1_qplus',t_var_genjet2_dAlpha_sj1_qplus,'genjet2_dAlpha_sj1_qplus/F')
        mytree.Branch('genjet2_dAlpha_H',t_var_genjet2_dAlpha_H,'genjet2_dAlpha_H/F')
        mytree.Branch('genjet2_d21', t_var_genjet2_d21, 'genjet2_d21/F');
        mytree.Branch('genjet2_d32', t_var_genjet2_d32, 'genjet2_d32/F');    
        mytree.Branch('genjet2_d43', t_var_genjet2_d43, 'genjet2_d43/F');

        mytree.Branch('parton_ep_E', t_var_parton_ep_E , 'parton_ep_E/F')
        mytree.Branch('parton_ep_Px', t_var_parton_ep_Px , 'parton_ep_Px/F')
        mytree.Branch('parton_ep_Py', t_var_parton_ep_Py , 'parton_ep_Py/F')
        mytree.Branch('parton_ep_Pz', t_var_parton_ep_Pz , 'parton_ep_Pz/F')
 
        mytree.Branch('parton_em_E', t_var_parton_em_E , 'parton_em_E/F')
        mytree.Branch('parton_em_Px', t_var_parton_em_Px , 'parton_em_Px/F')
        mytree.Branch('parton_em_Py', t_var_parton_em_Py , 'parton_em_Py/F')
        mytree.Branch('parton_em_Pz', t_var_parton_em_Pz , 'parton_em_Pz/F')

        mytree.Branch('parton_H_E', t_var_parton_H_E , 'parton_H_E/F')
        mytree.Branch('parton_H_Px', t_var_parton_H_Px , 'parton_H_Px/F')
        mytree.Branch('parton_H_Py', t_var_parton_H_Py , 'parton_H_Py/F')
        mytree.Branch('parton_H_Pz', t_var_parton_H_Pz , 'parton_H_Pz/F')
        mytree.Branch('parton_H_PDG_Daughter0', t_var_parton_H_PDG_Daughter0 , 'parton_H_PDG_Daughter0/I')

        mytree.Branch('parton_Z_qpos_E', t_var_parton_Z_qpos_E , 'parton_Z_qpos_E/F')
        mytree.Branch('parton_Z_qpos_Px', t_var_parton_Z_qpos_Px , 'parton_Z_qpos_Px/F')
        mytree.Branch('parton_Z_qpos_Py', t_var_parton_Z_qpos_Py , 'parton_Z_qpos_Py/F')
        mytree.Branch('parton_Z_qpos_Pz', t_var_parton_Z_qpos_Pz , 'parton_Z_qpos_Pz/F')
        mytree.Branch('parton_Z_qpos_PDGID', t_var_parton_Z_qpos_PDGID , 'parton_Z_qpos_PDGID/I')

        mytree.Branch('parton_Z_qneg_E', t_var_parton_Z_qneg_E , 'parton_Z_qneg_E/F')
        mytree.Branch('parton_Z_qneg_Px', t_var_parton_Z_qneg_Px , 'parton_Z_qneg_Px/F')
        mytree.Branch('parton_Z_qneg_Py', t_var_parton_Z_qneg_Py , 'parton_Z_qneg_Py/F')
        mytree.Branch('parton_Z_qneg_Pz', t_var_parton_Z_qneg_Pz , 'parton_Z_qneg_Pz/F')
        mytree.Branch('parton_Z_qneg_PDGID', t_var_parton_Z_qneg_PDGID , 'parton_Z_qneg_PDGID/I')

  


    #hist = file2.Get("h_runstatistics")

    #if tree.GetEntries()!=hist.GetBinContent(1) :
    #    print "tree_entries not hist content/diff", tree.GetEntries(),hist.GetBinContent(1), tree.GetEntries()-hist.GetBinContent(1)

    n_bins_high=100
    lim_theta_low=0.
    lim_theta_high=180.

    lim_costheta_low=-1.
    lim_costheta_high=1.

    lim_cosProd_low=-1.05
    lim_cosProd_high=1.05

    n_bins_high_2D=40
    
    lim_phi_low=-0.0
    lim_phi_high=360.0

    h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom = TH1F( "h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)
    h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom = TH1F( "h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)


    h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
    h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
    h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
    h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
    
    h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
    h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
    h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
    h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, lim_phi_low, lim_phi_high)

    h_2D_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM = TH2F( "h_2D_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)
    h_2D_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM = TH2F( "h_2D_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)

    h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom.Sumw2()
    h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom.Sumw2()

    h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Sumw2()
    h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Sumw2()
    h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Sumw2()
    h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Sumw2()
    
    h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Sumw2()
    h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Sumw2()
    h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Sumw2()
    h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Sumw2()

    h_2D_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Sumw2()
    h_2D_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Sumw2()

    h_parton_costheta1_Z_qpos_Zcom = TH1F( "h_parton_costheta1_Z_qpos_Zcom", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)
    h_parton_costheta2_H_ep_HZ_COM  = TH1F( "h_parton_costheta2_H_ep_HZ_COM", "", n_bins_high, lim_cosProd_low,lim_cosProd_high)
    h_parton_phi_plane_Z_qpos_vs_plane_H_ep= TH1F( "h_parton_phi_plane_Z_qpos_vs_plane_H_ep", "", n_bins_high, lim_phi_low, lim_phi_high)
    h_2D_parton_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM = TH2F( "h_2D_parton_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM", "", n_bins_high_2D, lim_cosProd_low,lim_cosProd_high,n_bins_high_2D, lim_cosProd_low,lim_cosProd_high)

    h_parton_costheta1_Z_qpos_Zcom.Sumw2()
    h_parton_costheta2_H_ep_HZ_COM.Sumw2()
    h_parton_phi_plane_Z_qpos_vs_plane_H_ep.Sumw2()
    h_2D_parton_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Sumw2()


    h_parton_theta1_Z_qpos_vs_Z = TH1F( "h_parton_theta1_Z_qpos_vs_Z", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_Z_qpos_vs_Z = TH1F( "h_parton_costheta1_Z_qpos_vs_Z", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_Z_qpos_vs_Z_sqrtS_750 = TH1F( "h_parton_theta1_Z_qpos_vs_Z_sqrtS_750", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_Z_qpos_vs_Z_sqrtS_750 = TH1F( "h_parton_costheta1_Z_qpos_vs_Z_sqrtS_750", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_Z_qpos_vs_Z_sqrtS_750_2500 = TH1F( "h_parton_theta1_Z_qpos_vs_Z_sqrtS_750_2500", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_Z_qpos_vs_Z_sqrtS_750_2500 = TH1F( "h_parton_costheta1_Z_qpos_vs_Z_sqrtS_750_2500", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_Z_qpos_vs_Z_sqrtS_2500 = TH1F( "h_parton_theta1_Z_qpos_vs_Z_sqrtS_2500", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_Z_qpos_vs_Z_sqrtS_2500 = TH1F( "h_parton_costheta1_Z_qpos_vs_Z_sqrtS_2500", "", n_bins_high, lim_costheta_low, lim_costheta_high)

    h_parton_theta1_Z_qneg_vs_Z = TH1F( "h_parton_theta1_Z_qneg_vs_Z", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_Z_qneg_vs_Z = TH1F( "h_parton_costheta1_Z_qneg_vs_Z", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_Z_qneg_vs_Z_sqrtS_750 = TH1F( "h_parton_theta1_Z_qneg_vs_Z_sqrtS_750", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_Z_qneg_vs_Z_sqrtS_750 = TH1F( "h_parton_costheta1_Z_qneg_vs_Z_sqrtS_750", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_Z_qneg_vs_Z_sqrtS_750_2500 = TH1F( "h_parton_theta1_Z_qneg_vs_Z_sqrtS_750_2500", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_Z_qneg_vs_Z_sqrtS_750_2500 = TH1F( "h_parton_costheta1_Z_qneg_vs_Z_sqrtS_750_2500", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_Z_qneg_vs_Z_sqrtS_2500 = TH1F( "h_parton_theta1_Z_qneg_vs_Z_sqrtS_2500", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_Z_qneg_vs_Z_sqrtS_2500 = TH1F( "h_parton_costheta1_Z_qneg_vs_Z_sqrtS_2500", "", n_bins_high, lim_costheta_low, lim_costheta_high)

    h_parton_theta1_H_bbar_vs_H = TH1F( "h_parton_theta1_H_bbar_vs_H", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_H_bbar_vs_H = TH1F( "h_parton_costheta1_H_bbar_vs_H", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_H_bbar_vs_H_sqrtS_750 = TH1F( "h_parton_theta1_H_bbar_vs_H_sqrtS_750", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_H_bbar_vs_H_sqrtS_750 = TH1F( "h_parton_costheta1_H_bbar_vs_H_sqrtS_750", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_H_bbar_vs_H_sqrtS_750_2500 = TH1F( "h_parton_theta1_H_bbar_vs_H_sqrtS_750_2500", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_H_bbar_vs_H_sqrtS_750_2500 = TH1F( "h_parton_costheta1_H_bbar_vs_H_sqrtS_750_2500", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_H_bbar_vs_H_sqrtS_2500 = TH1F( "h_parton_theta1_H_bbar_vs_H_sqrtS_2500", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_H_bbar_vs_H_sqrtS_2500 = TH1F( "h_parton_costheta1_H_bbar_vs_H_sqrtS_2500", "", n_bins_high, lim_costheta_low, lim_costheta_high)

    h_parton_theta1_H_b_vs_H = TH1F( "h_parton_theta1_H_b_vs_H", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_H_b_vs_H = TH1F( "h_parton_costheta1_H_b_vs_H", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_H_b_vs_H_sqrtS_750 = TH1F( "h_parton_theta1_H_b_vs_H_sqrtS_750", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_H_b_vs_H_sqrtS_750 = TH1F( "h_parton_costheta1_H_b_vs_H_sqrtS_750", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_H_b_vs_H_sqrtS_750_2500 = TH1F( "h_parton_theta1_H_b_vs_H_sqrtS_750_2500", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_H_b_vs_H_sqrtS_750_2500 = TH1F( "h_parton_costheta1_H_b_vs_H_sqrtS_750_2500", "", n_bins_high, lim_costheta_low, lim_costheta_high)
    h_parton_theta1_H_b_vs_H_sqrtS_2500 = TH1F( "h_parton_theta1_H_b_vs_H_sqrtS_2500", "", n_bins_high, lim_theta_low, lim_theta_high)
    h_parton_costheta1_H_b_vs_H_sqrtS_2500 = TH1F( "h_parton_costheta1_H_b_vs_H_sqrtS_2500", "", n_bins_high, lim_costheta_low, lim_costheta_high)

    sqrtS_high = 2500.0
    #number of leptons need to be smaller than this number
    m_cut_nLeptons = 1
 
    weight_total=0.0
    weight_total2=0.0
    tree = file.Get("showerData")

    weight = xsec*lumi/tree.GetEntries()
    print "tree-entries ",tree.GetEntries(), " weight ",weight, "xsec",xsec,"lumi",lumi,"total original ",weight*tree.GetEntries(), 'to', xsec*lumi

    num_entry=-1
    num_total_exception=0
    if fillGenLevel:
        t_var_GenLevelFilled[0]=1
    else: 
        t_var_GenLevelFilled[0]=0

    for ientry in tree:
        pass_reco_selection=True
        pass_parton_selection = False
        num_entry+=1
        #if num_entry!=-1:
        #    continue

        t_var_costheta1_for_Atheta1[0]=-10.
        t_var_costheta2_for_Atheta1theta2[0]=-10.
        t_var_phi_for_Aphis[0]=-10.

        t_var_sqrtS[0]=-10
        t_var_sqrtS_orig[0]=-10
        t_var_sqrtS_gen[0]=-10
        t_var_sqrtS_orig_gen[0]=-10
        t_var_sqrtS_parton[0]=-10
        t_var_nLeptons[0]=-10
        t_var_nLeptons_gen[0]=-10
        t_var_MET[0]=-10

        t_var_reco_y21[0]=-1
        t_var_reco_y32[0]=-1
        t_var_reco_y43[0]=-1

        t_var_jet1_mass[0]=-10
        t_var_jet2_mass[0]=-10
        t_var_jet1_min_jet2_mass[0]=-10
        t_var_jet1_BTag_rfj_BTagMax[0]=-10
        t_var_jet1_CTag_rfj_BTagMax[0]=-10
        t_var_jet1_LTag_rfj_BTagMax[0]=-10
        t_var_jet1_BTag_rfj_CTagMax[0]=-10
        t_var_jet1_CTag_rfj_CTagMax[0]=-10
        t_var_jet1_LTag_rfj_CTagMax[0]=-10
        t_var_jet1_E[0]=-10
        t_var_jet2_E[0]=-10
        t_var_jet1_Pt[0]=-10
        t_var_jet2_Pt[0]=-10
        t_var_jet1_theta[0]=-10
        t_var_jet2_theta[0]=-10
        t_var_jet1_min_jet2_theta[0]=-10
        t_var_jet1_phi[0]=-10
        t_var_jet2_phi[0]=-10
        t_var_dphi_j1j2[0]=-10
        t_var_angle_j1j2[0]=-10
        t_var_jet1_D2_beta1[0]=-10
        t_var_jet2_D2_beta1[0]=-10
        t_var_jet1_D2_beta0_5[0]=-10
        t_var_jet2_D2_beta0_5[0]=-10
        t_var_jet1_C2_beta1[0]=-10
        t_var_jet2_C2_beta1[0]=-10
        t_var_jet1_C2_beta0_5[0]=-10
        t_var_jet2_C2_beta0_5[0]=-10
        t_var_jet1_tau21[0]=-10
        t_var_jet2_tau21[0]=-10
        t_var_jet1_tau32[0]=-10
        t_var_jet2_tau32[0]=-10
        t_var_jet1_C3_beta1[0]=-10
        t_var_jet2_C3_beta1[0]=-10
        t_var_jet1_N2_beta1[0]=-10
        t_var_jet2_N2_beta1[0]=-10        
        t_var_jet1_N3_beta1[0]=-10
        t_var_jet2_N3_beta1[0]=-10

        t_var_jet1_D2_beta1_E_theta[0]=-10
        t_var_jet2_D2_beta1_E_theta[0]=-10
        t_var_jet1_C2_beta1_E_theta[0]=-10
        t_var_jet2_C2_beta1_E_theta[0]=-10
        t_var_jet1_tau21_E_theta[0]=-10
        t_var_jet2_tau21_E_theta[0]=-10
        t_var_jet1_tau32_E_theta[0]=-10
        t_var_jet2_tau32_E_theta[0]=-10
        t_var_jet1_C3_beta1_E_theta[0]=-10
        t_var_jet2_C3_beta1_E_theta[0]=-10
        t_var_jet1_N2_beta1_E_theta[0]=-10
        t_var_jet2_N2_beta1_E_theta[0]=-10        
        t_var_jet1_N3_beta1_E_theta[0]=-10
        t_var_jet2_N3_beta1_E_theta[0]=-10

        t_var_jet1_sj1_E[0]=-10
        t_var_jet1_sj1_Px[0]=-10
        t_var_jet1_sj1_Py[0]=-10
        t_var_jet1_sj1_Pz[0]=-10
        t_var_jet1_sj1_jetChargePt_kappa_0_30[0]=-10
        t_var_jet1_sj1_jetChargeE_kappa_0_30[0]=-10
        t_var_jet1_sj1_jetChargePt_kappa_0_25[0]=-10
        t_var_jet1_sj1_jetChargeE_kappa_0_25[0]=-10
        t_var_jet1_sj1_nTracks[0]=-10
        t_var_jet1_sj1_chFrac[0]=-10
        t_var_jet1_sj1_closestMatch[0]=-10
        t_var_jet1_sj1_decMatch[0]=-10
        t_var_jet1_sj1_Angle_closestMatch[0]=-10
        t_var_jet1_sj1_Angle_decMatch[0]=-10
        
        t_var_jet1_sj2_E[0]=-10
        t_var_jet1_sj2_Px[0]=-10
        t_var_jet1_sj2_Py[0]=-10
        t_var_jet1_sj2_Pz[0]=-10
        t_var_jet1_sj2_jetChargePt_kappa_0_30[0]=-10
        t_var_jet1_sj2_jetChargeE_kappa_0_30[0]=-10
        t_var_jet1_sj2_jetChargePt_kappa_0_25[0]=-10
        t_var_jet1_sj2_jetChargeE_kappa_0_25[0]=-10
        t_var_jet1_sj2_nTracks[0]=-10
        t_var_jet1_sj2_chFrac[0]=-10

        t_var_jet1_d21[0]=-1
        t_var_jet1_d32[0]=-1
        t_var_jet1_d43[0]=-1
        t_var_jet2_d21[0]=-1
        t_var_jet2_d32[0]=-1
        t_var_jet2_d43[0]=-1


        t_var_jet1_sj2_closestMatch[0]=-10
    
        t_var_jet1_sj2_decMatch[0]=-10
        t_var_jet1_sj2_Angle_closestMatch[0]=-10
        t_var_jet1_sj2_Angle_decMatch[0]=-10
        t_var_jet1_dAlpha_sj1sj2[0]=-10
        
        t_var_jet2_sj1_E[0]=-10
        t_var_jet2_sj1_Px[0]=-10
        t_var_jet2_sj1_Py[0]=-10
        t_var_jet2_sj1_Pz[0]=-10
        t_var_jet2_sj1_jetChargePt_kappa_0_30[0]=-10
        t_var_jet2_sj1_jetChargeE_kappa_0_30[0]=-10
        t_var_jet2_sj1_jetChargePt_kappa_0_25[0]=-10
        t_var_jet2_sj1_jetChargeE_kappa_0_25[0]=-10
        t_var_jet2_sj1_nTracks[0]=-10
        t_var_jet2_sj1_chFrac[0]=-10
        t_var_jet2_sj1_closestMatch[0]=-10
        t_var_jet2_sj1_decMatch[0]=-10
        t_var_jet2_sj1_Angle_closestMatch[0]=-10
        t_var_jet2_sj1_Angle_decMatch[0]=-10
        
        t_var_jet2_sj2_E[0]=-10
        t_var_jet2_sj2_Px[0]=-10
        t_var_jet2_sj2_Py[0]=-10
        t_var_jet2_sj2_Pz[0]=-10
        t_var_jet2_sj2_jetChargePt_kappa_0_30[0]=-10
        t_var_jet2_sj2_jetChargeE_kappa_0_30[0]=-10
        t_var_jet2_sj2_jetChargePt_kappa_0_25[0]=-10
        t_var_jet2_sj2_jetChargeE_kappa_0_25[0]=-10
        t_var_jet2_sj2_nTracks[0]=-10
        t_var_jet2_sj2_chFrac[0]=-10
        t_var_jet2_sj2_closestMatch[0]=-10
        t_var_jet2_sj2_decMatch[0]=-10
        t_var_jet2_sj2_Angle_closestMatch[0]=-10
        t_var_jet2_sj2_Angle_decMatch[0]=-10
        t_var_jet2_dAlpha_sj1sj2[0]=-10
        
        t_var_genInv_E[0]=-10
        t_var_genInv_Px[0]=-10
        t_var_genInv_Py[0]=-10
        t_var_genInv_Pz[0]=-10
        
        t_var_gen_y21[0]=-1
        t_var_gen_y32[0]=-1
        t_var_gen_y43[0]=-1
        
        t_var_genjet1_sj1_E[0]=-10
        t_var_genjet1_sj1_Px[0]=-10
        t_var_genjet1_sj1_Py[0]=-10
        t_var_genjet1_sj1_Pz[0]=-10
        t_var_genjet1_sj1_jetChargePt_kappa_0_30[0]=-10
        t_var_genjet1_sj1_jetChargeE_kappa_0_30[0]=-10
        t_var_genjet1_sj1_jetChargePt_kappa_0_25[0]=-10
        t_var_genjet1_sj1_jetChargeE_kappa_0_25[0]=-10
        t_var_genjet1_sj1_nTracks[0]=-10
        t_var_genjet1_sj1_chFrac[0]=-10
        
        t_var_genjet1_sj2_E[0]=-10
        t_var_genjet1_sj2_Px[0]=-10
        t_var_genjet1_sj2_Py[0]=-10
        t_var_genjet1_sj2_Pz[0]=-10
        t_var_genjet1_sj2_jetChargePt_kappa_0_30[0]=-10
        t_var_genjet1_sj2_jetChargeE_kappa_0_30[0]=-10
        t_var_genjet1_sj2_jetChargePt_kappa_0_25[0]=-10
        t_var_genjet1_sj2_jetChargeE_kappa_0_25[0]=-10
        t_var_genjet1_sj2_nTracks[0]=-10
        t_var_genjet1_sj2_chFrac[0]=-10
        t_var_genjet1_dAlpha_sj1sj2[0]=-10
        t_var_genjet1_d21[0]=-1
        t_var_genjet1_d32[0]=-1
        t_var_genjet1_d43[0]=-1
        t_var_genjet2_d21[0]=-1
        t_var_genjet2_d32[0]=-1
        t_var_genjet2_d43[0]=-1
        
        t_var_genjet1_dAlpha_sj1_rjsj1[0]=-10
        t_var_genjet1_dAlpha_sj1_rjsj2[0]=-10
        t_var_genjet1_dAlpha_sj1_qmin[0]=-10
        t_var_genjet1_dAlpha_sj1_qplus[0]=-10
        t_var_genjet1_dAlpha_H[0]=-10
        
        t_var_parton_ep_E[0]=-10
        t_var_parton_ep_Px[0]=-10
        t_var_parton_ep_Pz[0]=-10
        t_var_parton_ep_Py[0]=-10
        
        t_var_parton_em_E[0]=-10
        t_var_parton_em_Px[0]=-10
        t_var_parton_em_Pz[0]=-10
        t_var_parton_em_Py[0]=-10
        
        t_var_parton_H_E[0]=-10
        t_var_parton_H_Px[0]=-10
        t_var_parton_H_Pz[0]=-10
        t_var_parton_H_Py[0]=-10
        t_var_parton_H_PDG_Daughter0[0]=-10
        
        t_var_parton_Z_qpos_E[0]=-10
        t_var_parton_Z_qpos_Px[0]=-10
        t_var_parton_Z_qpos_Pz[0]=-10
        t_var_parton_Z_qpos_Py[0]=-10
        t_var_parton_Z_qpos_PDGID[0]=-10
        
        t_var_parton_Z_qneg_E[0]=-10
        t_var_parton_Z_qneg_Px[0]=-10
        t_var_parton_Z_qneg_Pz[0]=-10
        t_var_parton_Z_qneg_Py[0]=-10
        t_var_parton_Z_qneg_PDGID[0]=-10
        
        t_var_genjet2_sj1_E[0]=-10
        t_var_genjet2_sj1_Px[0]=-10
        t_var_genjet2_sj1_Py[0]=-10
        t_var_genjet2_sj1_Pz[0]=-10
        t_var_genjet2_sj1_jetChargePt_kappa_0_30[0]=-10
        t_var_genjet2_sj1_jetChargeE_kappa_0_30[0]=-10
        t_var_genjet2_sj1_jetChargePt_kappa_0_25[0]=-10
        t_var_genjet2_sj1_jetChargeE_kappa_0_25[0]=-10
        t_var_genjet2_sj1_nTracks[0]=-10
        t_var_genjet2_sj1_chFrac[0]=-10
        
        
        t_var_genjet2_sj2_E[0]=-10
        t_var_genjet2_sj2_Px[0]=-10
        t_var_genjet2_sj2_Py[0]=-10
        t_var_genjet2_sj2_Pz[0]=-10
        t_var_genjet2_sj2_jetChargePt_kappa_0_30[0]=-10
        t_var_genjet2_sj2_jetChargeE_kappa_0_30[0]=-10
        t_var_genjet2_sj2_jetChargePt_kappa_0_25[0]=-10
        t_var_genjet2_sj2_jetChargeE_kappa_0_25[0]=-10
        t_var_genjet2_sj2_nTracks[0]=-10
        t_var_genjet2_sj2_chFrac[0]=-10
        t_var_genjet2_dAlpha_sj1sj2[0]=-10
        
        t_var_genjet2_dAlpha_sj1_rjsj1[0]=-10
        t_var_genjet2_dAlpha_sj1_rjsj2[0]=-10
        t_var_genjet2_dAlpha_sj1_qmin[0]=-10
        t_var_genjet2_dAlpha_sj1_qplus[0]=-10
        t_var_genjet2_dAlpha_H[0]=-10
        t_var_GenLevelFilled[0]=-10
        t_var_eventWeight[0]=-10






        t_var_eventWeight[0]=weight


        if num_entry%(int(tree.GetEntries()/5.)) == 0:
            print "sig BG in entry ",num_entry


        reco_pass_mass_cuts=False

        #first reco jet on number of identified subjets follows below --> in order to check parton level quantities now
        quark_Z_decays = 0
        H_decays_bbar = True

        temp_ep=TLorentzVector(0,0,0,0);
        temp_em=TLorentzVector(0,0,0,0);

        tempHP4=TLorentzVector(0,0,0,0);
        tempH_b=TLorentzVector(0,0,0,0);
        tempH_bbar=TLorentzVector(0,0,0,0);
        tempH_p1=TLorentzVector(0,0,0,0);
        tempH_p2=TLorentzVector(0,0,0,0);

        tempZP4=TLorentzVector(0,0,0,0);
        tempZ_q_pos=TLorentzVector(0,0,0,0);
        tempZ_q_neg=TLorentzVector(0,0,0,0);
      
        t_var_sqrtS_parton[0]=-1
        t_var_sqrtS_gen[0]=-1
        t_var_sqrtS_orig_gen[0]=-1
        t_var_nLeptons_gen[0]=-1

        t_var_sqrtS[0]=-1
        t_var_sqrtS_orig[0]=-1
        t_var_nLeptons[0] =-1
        t_var_parton_H_PDG_Daughter0[0]=0

        index_firstH=-1
        index_firstZ=-1

        pass_gen_parton_cuts=False

        if usePartonInfo :
            trueME_E=ientry.trueME_E
            trueME_Px=ientry.trueME_Px
            trueME_Py=ientry.trueME_Py
            trueME_Pz=ientry.trueME_Pz
            trueME_PDGID=ientry.trueME_PDGID

            for i in range(len(trueME_E)):
                if trueME_PDGID[i]==25:
                   index_firstH=i
                   break
            for i in range(len(trueME_E)):
                if trueME_PDGID[i]==23:
                   index_firstZ=i
                   break
            #then 4 is Z,5 is H typically, but stupid enough we have exceptions for whatever reason
            if range(len(trueME_E)) > 7:
                if trueME_PDGID[4]!=23 or trueME_PDGID[5]!=25 :
                    num_total_exception+=1
                quark_Z_decays=abs(trueME_PDGID[index_firstH+1])

                if(abs(trueME_PDGID[index_firstH+3])!=5 or abs(trueME_PDGID[index_firstH+4])!=5) :
                   H_decays_bbar=False

            t_var_parton_H_PDG_Daughter0[0]=trueME_PDGID[index_firstH+3]

            ind_ep=-1
            ind_em=-1
            if trueME_PDGID[0]==-11 and trueME_PDGID[1]==11 : 
                ind_ep=0
                ind_em=1
            elif trueME_PDGID[0]==11 and trueME_PDGID[1]==-11 :
                ind_ep=1
                ind_em=0
            else:
                print 'indices should be consistent with outgoing e+e-',trueME_PDGID[0],trueME_PDGID[1]

            t_var_parton_em_Px[0]=trueME_Px[ind_em]
            t_var_parton_em_Py[0]=trueME_Py[ind_em]
            t_var_parton_em_Pz[0]=trueME_Pz[ind_em]
            t_var_parton_em_E[0]=trueME_E[ind_em]
            t_var_parton_ep_Px[0]=trueME_Px[ind_ep]
            t_var_parton_ep_Py[0]=trueME_Py[ind_ep]
            t_var_parton_ep_Pz[0]=trueME_Pz[ind_ep]
            t_var_parton_ep_E[0]=trueME_E[ind_ep]
            tempep=TLorentzVector(0,0,0,0);
            tempem=TLorentzVector(0,0,0,0);
            tempep.SetPxPyPzE(trueME_Px[ind_ep],trueME_Py[ind_ep],trueME_Pz[ind_ep],trueME_E[ind_ep])
            tempem.SetPxPyPzE(trueME_Px[ind_em],trueME_Py[ind_em],trueME_Pz[ind_em],trueME_E[ind_em])
            t_var_sqrtS_parton[0]=(tempep+tempem).M()


            if trueME_PDGID[index_firstH+3]==5 :
                tempH_b.SetPxPyPzE(trueME_Px[index_firstH+3],trueME_Py[index_firstH+3],trueME_Pz[index_firstH+3],trueME_E[index_firstH+3])
                tempH_bbar.SetPxPyPzE(trueME_Px[index_firstH+4],trueME_Py[index_firstH+4],trueME_Pz[index_firstH+4],trueME_E[index_firstH+4])
                tempH_p1=tempH_b
                tempH_p2=tempH_bbar
                if trueME_PDGID[index_firstH+3] != -trueME_PDGID[index_firstH+4]:
                    print 'pdg id should for H be the same, what is wrong ',trueME_PDGID[index_firstH+3],-trueME_PDGID[index_firstH+4]
            elif trueME_PDGID[index_firstH+3]==-5 :
                tempH_bbar.SetPxPyPzE(trueME_Px[index_firstH+3],trueME_Py[index_firstH+3],trueME_Pz[index_firstH+3],trueME_E[index_firstH+3])
                tempH_b.SetPxPyPzE(trueME_Px[index_firstH+4],trueME_Py[index_firstH+4],trueME_Pz[index_firstH+4],trueME_E[index_firstH+4])
                tempH_p1=tempH_b
                tempH_p2=tempH_bbar
                if trueME_PDGID[index_firstH+3] != -trueME_PDGID[index_firstH+4]:
                    print 'pdg id should for H be the same, what is wrong ',trueME_PDGID[index_firstH+3],-trueME_PDGID[index_firstH+4]
            else:
                tempH_p1.SetPxPyPzE(trueME_Px[index_firstH+3],trueME_Py[index_firstH+3],trueME_Pz[index_firstH+3],trueME_E[index_firstH+3])
                tempH_p1.SetPxPyPzE(trueME_Px[index_firstH+4],trueME_Py[index_firstH+4],trueME_Pz[index_firstH+4],trueME_E[index_firstH+4])
                if trueME_PDGID[index_firstH+3] != -trueME_PDGID[index_firstH+4]:
                    if H_decays_bbar:
                        print 'pdg id should for H be the same, what is wrong ',trueME_PDGID[index_firstH+3],-trueME_PDGID[index_firstH+4]

            #positive charge, aka up type quarks or down bar type quarks
            if(trueME_PDGID[index_firstH+1]==-1 or trueME_PDGID[index_firstH+1]==2 or trueME_PDGID[index_firstH+1]==-3 or trueME_PDGID[index_firstH+1]==4 or trueME_PDGID[index_firstH+1]==-5 or trueME_PDGID[index_firstH+1]==6):
                if trueME_PDGID[index_firstH+1] != -trueME_PDGID[index_firstH+2]:
                    print 'pdg id should be the same, what is wrong ',trueME_PDGID[index_firstH+1],-trueME_PDGID[index_firstH+2]
                tempZ_q_pos.SetPxPyPzE(trueME_Px[index_firstH+1],trueME_Py[index_firstH+1],trueME_Pz[index_firstH+1],trueME_E[index_firstH+1])
                tempZ_q_neg.SetPxPyPzE(trueME_Px[index_firstH+2],trueME_Py[index_firstH+2],trueME_Pz[index_firstH+2],trueME_E[index_firstH+2])
                t_var_parton_Z_qpos_PDGID[0]=trueME_PDGID[index_firstH+1]
                t_var_parton_Z_qneg_PDGID[0]=trueME_PDGID[index_firstH+2]
            else:
                #quark index 6 is negatively charged
                if trueME_PDGID[index_firstH+1] != -trueME_PDGID[index_firstH+2]:
                    print 'pdg id should be the same, what is wrong ',trueME_PDGID[index_firstH+1],-trueME_PDGID[index_firstH+2]
                tempZ_q_neg.SetPxPyPzE(trueME_Px[index_firstH+1],trueME_Py[index_firstH+1],trueME_Pz[index_firstH+1],trueME_E[index_firstH+1])
                tempZ_q_pos.SetPxPyPzE(trueME_Px[index_firstH+2],trueME_Py[index_firstH+2],trueME_Pz[index_firstH+2],trueME_E[index_firstH+2])
                t_var_parton_Z_qneg_PDGID[0]=trueME_PDGID[index_firstH+1]
                t_var_parton_Z_qpos_PDGID[0]=trueME_PDGID[index_firstH+2]

            tempZ_first=TLorentzVector(0,0,0,0)
            tempZ_first.SetPxPyPzE(trueME_Px[index_firstZ],trueME_Py[index_firstZ],trueME_Pz[index_firstZ],trueME_E[index_firstZ])

            tempHP4=tempH_p1+tempH_p2
            t_var_parton_H_E[0]=tempHP4.E()
            t_var_parton_H_Px[0]=tempHP4.Px()
            t_var_parton_H_Py[0]=tempHP4.Py()
            t_var_parton_H_Pz[0]=tempHP4.Pz()
                      
            tempZP4=tempZ_q_pos+tempZ_q_neg
            #print 'tempZp4 mass, delta wrt to H',tempZP4.M(),degrees(tempZP4.Angle(tempHP4.Vect())),(tempep+tempem).M(),tempZ_first.M()
            #if(tempZP4.M()-91.0)>20 :
            #    print 'mass far off',tempZP4.M(),degrees(tempZP4.Angle(tempHP4.Vect())),(tempep+tempem).M(),trueME_PDGID[index_firstH+1],trueME_PDGID[index_firstH+2],tempZ_first.M(),index_firstH,index_firstZ

            boostZ_COM_parton=TVector3(0,0,0) 
            boostZ_COM_parton=-tempZP4.BoostVector()



            t_var_parton_Z_qneg_E[0]=tempZ_q_neg.E()
            t_var_parton_Z_qneg_Px[0]=tempZ_q_neg.Px()
            t_var_parton_Z_qneg_Py[0]=tempZ_q_neg.Py()
            t_var_parton_Z_qneg_Pz[0]=tempZ_q_neg.Pz()
            t_var_parton_Z_qpos_E[0]=tempZ_q_pos.E()
            t_var_parton_Z_qpos_Px[0]=tempZ_q_pos.Px()
            t_var_parton_Z_qpos_Py[0]=tempZ_q_pos.Py()
            t_var_parton_Z_qpos_Pz[0]=tempZ_q_pos.Pz()
 
            boostH_COM_parton=TVector3(0,0,0) 
            boostH_COM_parton=-tempHP4.BoostVector()
            boostZ_COM_parton=TVector3(0,0,0) 
            boostZ_COM_parton=-tempZP4.BoostVector()
            boostE_tot_COM_parton=TVector3(0,0,0) 
            boostE_tot_COM_parton=-(tempep+tempem).BoostVector()

            tempep_in_E_tot_COM=TLorentzVector(0,0,0,0);            
            tempep_in_E_tot_COM.SetPxPyPzE(tempep.Px(),tempep.Py(),tempep.Pz(),tempep.E())
            tempep_in_E_tot_COM.Boost(boostE_tot_COM_parton)
           
            tempZP4_in_E_tot_COM_parton=TLorentzVector(0,0,0,0)
            tempZP4_in_E_tot_COM_parton.SetPxPyPzE(tempZP4.Px(),tempZP4.Py(),tempZP4.Pz(),tempZP4.E())
            tempZP4_in_E_tot_COM_parton.Boost(boostE_tot_COM_parton)
            tempHP4_in_E_tot_COM_parton=TLorentzVector(0,0,0,0)
            tempHP4_in_E_tot_COM_parton.SetPxPyPzE(tempHP4.Px(),tempHP4.Py(),tempHP4.Pz(),tempHP4.E())
            tempHP4_in_E_tot_COM_parton.Boost(boostE_tot_COM_parton)


            if (abs(tempZP4_in_E_tot_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect()) -pi) >1.e-4 and H_decays_bbar ):
                print 'H and Z in E tot COM not back to back',degrees(tempZP4_in_E_tot_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect()))

            tempZ_q_pos_Z_COM_parton=TLorentzVector(0,0,0,0)
            tempZ_q_pos_Z_COM_parton.SetPxPyPzE(tempZ_q_pos.Px(),tempZ_q_pos.Py(),tempZ_q_pos.Pz(),tempZ_q_pos.E())
            tempZ_q_pos_Z_COM_parton.Boost(boostZ_COM_parton)
            tempZ_q_neg_Z_COM_parton=TLorentzVector(0,0,0,0)
            tempZ_q_neg_Z_COM_parton.SetPxPyPzE(tempZ_q_neg.Px(),tempZ_q_neg.Py(),tempZ_q_neg.Pz(),tempZ_q_neg.E())
            tempZ_q_neg_Z_COM_parton.Boost(boostZ_COM_parton)
           

            #angular histograms for parton level at HIGH sqrtS
            #if (tempep+tempem).M()>sqrtS_high and abs(degrees(tempZP4.Theta())-90.)<55. and abs(degrees(tempHP4.Theta())-90.)<55. and abs(degrees(tempHP4.Theta()-tempZP4.Theta()))<100 :
            if (tempep+tempem).M()>sqrtS_high :
                #pass_parton_selection = True
                #print 'parton Z/H/delta ',degrees(tempZP4.Theta()),degrees(tempHP4.Theta()),degrees(tempHP4.Theta()-tempZP4.Theta())
                h_parton_costheta1_Z_qpos_vs_Z.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)

                #angle theta1 --> between tempZ in E_tot_COM system and tempZ_q_pos_Z_COM_parton direction
                h_parton_costheta1_Z_qpos_Zcom.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                if( cos(2.*tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())) /abs( cos(2.*tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())) ))>0:
                    h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                elif( cos(2.*tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())) /abs( cos(2.*tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())) ))<0:
                    h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),-weight)

                h_parton_costheta2_H_ep_HZ_COM.Fill(cos(tempHP4_in_E_tot_COM_parton.Angle(tempep_in_E_tot_COM.Vect())),weight)
                h_2D_parton_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),cos(tempHP4_in_E_tot_COM_parton.Angle(tempep_in_E_tot_COM.Vect())),weight)
                #angle theta2 --> between tempH in E_tot_COM system and temp_ep in E_tot_COM system
                if(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect()))/abs(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())))*cos(tempHP4_in_E_tot_COM_parton.Angle(tempep_in_E_tot_COM.Vect()))/abs(cos(tempHP4_in_E_tot_COM_parton.Angle(tempep_in_E_tot_COM.Vect()))))>0:
                    h_2D_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),cos(tempHP4_in_E_tot_COM_parton.Angle(tempep_in_E_tot_COM.Vect())),weight)
                elif(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect()))/abs(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())))*cos(tempHP4_in_E_tot_COM_parton.Angle(tempep_in_E_tot_COM.Vect()))/abs(cos(tempHP4_in_E_tot_COM_parton.Angle(tempep_in_E_tot_COM.Vect()))))<0:
                    h_2D_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),cos(tempHP4_in_E_tot_COM_parton.Angle(tempep_in_E_tot_COM.Vect())),-weight)
                plane_Z_Z_qpos_boost=TVector3(0,0,0)
                plane_Z_Z_qpos_boost=tempZP4_in_E_tot_COM_parton.Vect().Cross(tempZ_q_pos_Z_COM_parton.Vect()).Unit()          
                plane_H_ep_boost=TVector3(0,0,0)
                plane_H_ep_boost=tempHP4_in_E_tot_COM_parton.Vect().Cross(tempep_in_E_tot_COM.Vect()).Unit()
                phi_plane_parton_H_ep_Z_Z_qpos = degrees(acos(plane_H_ep_boost.Dot(plane_Z_Z_qpos_boost)))

                if degrees(tempZ_q_pos_Z_COM_parton.Angle(plane_H_ep_boost))>90.:
                    phi_plane_parton_H_ep_Z_Z_qpos=degrees(acos(plane_H_ep_boost.Dot(plane_Z_Z_qpos_boost)))+180.

                h_parton_phi_plane_Z_qpos_vs_plane_H_ep.Fill(phi_plane_parton_H_ep_Z_Z_qpos,weight)
                phi_plane_parton_H_ep_Z_Z_qpos_rad = radians(phi_plane_parton_H_ep_Z_Z_qpos)
                if sin(phi_plane_parton_H_ep_Z_Z_qpos_rad)!=0:
                    if(sin(phi_plane_parton_H_ep_Z_Z_qpos_rad)/abs(sin(phi_plane_parton_H_ep_Z_Z_qpos_rad)))>0:
                        h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(phi_plane_parton_H_ep_Z_Z_qpos,weight)
                    else:
                        h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(phi_plane_parton_H_ep_Z_Z_qpos,-weight)
                if sin(2.*phi_plane_parton_H_ep_Z_Z_qpos_rad)!=0:
                    if(sin(2*phi_plane_parton_H_ep_Z_Z_qpos_rad)/abs(sin(2*phi_plane_parton_H_ep_Z_Z_qpos_rad)))>0:
                        h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(phi_plane_parton_H_ep_Z_Z_qpos,weight)
                    else:
                        h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(phi_plane_parton_H_ep_Z_Z_qpos,-weight)
                if cos(phi_plane_parton_H_ep_Z_Z_qpos_rad)!=0:
                    if(cos(phi_plane_parton_H_ep_Z_Z_qpos_rad)/abs(cos(phi_plane_parton_H_ep_Z_Z_qpos_rad)))>0:
                        h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(phi_plane_parton_H_ep_Z_Z_qpos,weight)
                    else:
                        h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(phi_plane_parton_H_ep_Z_Z_qpos,-weight)
                if cos(2.*phi_plane_parton_H_ep_Z_Z_qpos_rad)!=0:
                    if(cos(2*phi_plane_parton_H_ep_Z_Z_qpos_rad)/abs(cos(2*phi_plane_parton_H_ep_Z_Z_qpos_rad)))>0:
                        h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(phi_plane_parton_H_ep_Z_Z_qpos,weight)
                    else:
                        h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Fill(phi_plane_parton_H_ep_Z_Z_qpos,-weight)


                if (abs(tempZ_q_pos_Z_COM_parton.Angle(tempZ_q_neg_Z_COM_parton.Vect()) -pi) >1.e-4):
                    print 'qpos and qneg in Z COM not back to back',degrees(tempZ_q_pos_Z_COM_parton.Angle(tempZ_q_neg_Z_COM_parton.Vect()))

            tempH_bbar_H_COM_parton=TLorentzVector(0,0,0,0)
            tempH_bbar_H_COM_parton=tempH_bbar
            tempH_bbar_H_COM_parton.Boost(boostH_COM_parton)
            tempH_b_H_COM_parton=TLorentzVector(0,0,0,0)
            tempH_b_H_COM_parton=tempH_b
            tempH_b_H_COM_parton.Boost(boostH_COM_parton)
                
            if (abs(tempH_bbar_H_COM_parton.Angle(tempH_b_H_COM_parton.Vect())-pi)>1.e-4 and H_decays_bbar):
                print 'b and bbar in H COM not back to back now ',degrees(tempH_bbar_H_COM_parton.Angle(tempH_b_H_COM_parton.Vect())),H_decays_bbar,tempH_b_H_COM_parton.P(),tempH_bbar_H_COM_parton.P()


 
            h_parton_theta1_Z_qpos_vs_Z.Fill(degrees(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
            h_parton_costheta1_Z_qpos_vs_Z.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
            h_parton_theta1_H_bbar_vs_H.Fill(degrees(tempH_bbar_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
            h_parton_costheta1_H_bbar_vs_H.Fill(cos(tempH_bbar_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
            if (tempep+tempem).M()<750. :
                h_parton_theta1_Z_qpos_vs_Z_sqrtS_750.Fill(degrees(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_Z_qpos_vs_Z_sqrtS_750.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_theta1_H_bbar_vs_H_sqrtS_750.Fill(degrees(tempH_bbar_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_H_bbar_vs_H_sqrtS_750.Fill(cos(tempH_bbar_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
            elif (tempep+tempem).M()<2500. :
                h_parton_theta1_Z_qpos_vs_Z_sqrtS_750_2500.Fill(degrees(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_Z_qpos_vs_Z_sqrtS_750_2500.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_theta1_H_bbar_vs_H_sqrtS_750_2500.Fill(degrees(tempH_bbar_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_H_bbar_vs_H_sqrtS_750_2500.Fill(cos(tempH_bbar_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
            else :
                #histograms at highest energy sqrtS>2500 GeV
                h_parton_theta1_Z_qpos_vs_Z_sqrtS_2500.Fill(degrees(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_Z_qpos_vs_Z_sqrtS_2500.Fill(cos(tempZ_q_pos_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_theta1_H_bbar_vs_H_sqrtS_2500.Fill(degrees(tempH_bbar_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_H_bbar_vs_H_sqrtS_2500.Fill(cos(tempH_bbar_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)


 


            h_parton_theta1_Z_qneg_vs_Z.Fill(degrees(tempZ_q_neg_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
            h_parton_costheta1_Z_qneg_vs_Z.Fill(cos(tempZ_q_neg_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
            h_parton_theta1_H_b_vs_H.Fill(degrees(tempH_b_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
            h_parton_costheta1_H_b_vs_H.Fill(cos(tempH_b_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
            if (tempep+tempem).M()<750. :
                h_parton_theta1_Z_qneg_vs_Z_sqrtS_750.Fill(degrees(tempZ_q_neg_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_Z_qneg_vs_Z_sqrtS_750.Fill(cos(tempZ_q_neg_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_theta1_H_b_vs_H_sqrtS_750.Fill(degrees(tempH_b_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_H_b_vs_H_sqrtS_750.Fill(cos(tempH_b_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
            elif (tempep+tempem).M()<2500. :
                h_parton_theta1_Z_qneg_vs_Z_sqrtS_750_2500.Fill(degrees(tempZ_q_neg_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_Z_qneg_vs_Z_sqrtS_750_2500.Fill(cos(tempZ_q_neg_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_theta1_H_b_vs_H_sqrtS_750_2500.Fill(degrees(tempH_b_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_H_b_vs_H_sqrtS_750_2500.Fill(cos(tempH_b_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
            else :
                h_parton_theta1_Z_qneg_vs_Z_sqrtS_2500.Fill(degrees(tempZ_q_neg_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_Z_qneg_vs_Z_sqrtS_2500.Fill(cos(tempZ_q_neg_Z_COM_parton.Angle(tempZP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_theta1_H_b_vs_H_sqrtS_2500.Fill(degrees(tempH_b_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)
                h_parton_costheta1_H_b_vs_H_sqrtS_2500.Fill(cos(tempH_b_H_COM_parton.Angle(tempHP4_in_E_tot_COM_parton.Vect())),weight)


            #print outs below show, once we move away from decays of H to bbar, the filling of H and Z quantities are far off
            #determination by H from saved decay products doesn't work, since we would need to save the structure differently
            #bbar decay structure is though still fine, once we select H-->bbar things are as expected and both ways determine sqrtS the same
            #if (((tempep+tempem).M()-(tempHP4+tempZP4).M())/(tempep+tempem).M())>1.e-4:
            #    print 'epem M, H Z M',(tempep+tempem).M(),(tempHP4+tempZP4).M(),'delta',((tempep+tempem).M()-(tempHP4+tempZP4).M())/(tempep+tempem).M()
            #else:
            #    print 'things are fine?'

        #if(H_decays_bbar==False):
        #    continue

        tempTotRecoP4=TLorentzVector(0,0,0,0);
        tempTotRecoP4.SetPxPyPzE(ientry.totPFO_Px,ientry.totPFO_Py,ientry.totPFO_Pz,ientry.totPFO_E)


        recojet_subjet_rfj_j_E=ientry.recojet_subjet_rfj_j_E
        recojet_subjet_E=ientry.recojet_subjet_E
 
        if len(recojet_subjet_rfj_j_E)<4 or len(recojet_subjet_E)<4 :
            #print "too few subjets or refined jets ",len(recojet_subjet_E),len(recojet_subjet_rfj_j_E),num_entry,tempTotRecoP4.E()
            pass_reco_selection=False


        n_IsoPh_gen=0
        n_IsoLep_gen=0
        sqrtS_eff_gen=-1
        tempGenIsoPhP4=TLorentzVector(0,0,0,0)
        tempGenIsoLepP4=TLorentzVector(0,0,0,0)
        isoPartGenDR10_E=ientry.isoPartGenDR10_E
        isoPartGenDR10_Px=ientry.isoPartGenDR10_Px
        isoPartGenDR10_Py=ientry.isoPartGenDR10_Py
        isoPartGenDR10_Pz=ientry.isoPartGenDR10_Pz
        isoPartGenDR10_relIso=ientry.isoPartGenDR10_relIso
        isoPartGenDR10_PDGID=ientry.isoPartGenDR10_PDGID
        
        for i in range(len(isoPartGenDR10_E)):
            if isoPartGenDR10_PDGID[i]==22:
                if isoPartGenDR10_relIso[i]<0.10:
                    temp=TLorentzVector(0,0,0,0)
                    temp.SetPxPyPzE(isoPartGenDR10_Px[i],isoPartGenDR10_Py[i],isoPartGenDR10_Pz[i],isoPartGenDR10_E[i])
                    tempGenIsoPhP4+=temp
                    n_IsoPh_gen+=1
            elif (abs(isoPartGenDR10_PDGID[i])==11 or abs(isoPartGenDR10_PDGID[i])==13) :
                if isoPartGenDR10_relIso[i]<0.10:
                    temp=TLorentzVector(0,0,0,0)
                    temp.SetPxPyPzE(isoPartGenDR10_Px[i],isoPartGenDR10_Py[i],isoPartGenDR10_Pz[i],isoPartGenDR10_E[i])
                    tempGenIsoLepP4+=temp
                    n_IsoLep_gen+=1

        genjet_E=ientry.genjet_E
        #if len(genjet_E)<2 : 
        #    print "less than two jets in event, skip"
        #    continue;
        genjet_Px=ientry.genjet_Px
        genjet_Py=ientry.genjet_Py
        genjet_Pz=ientry.genjet_Pz

        t_var_gen_y21[0]=ientry.gen_y21_max
        t_var_gen_y32[0]=ientry.gen_y32_max
        t_var_gen_y43[0]=ientry.gen_y43_max


        gj_m1_orig=TLorentzVector(0,0,0,0);
        gj_m2_orig=TLorentzVector(0,0,0,0);
        if len(genjet_E)>=2:
            gj_m1_orig.SetPxPyPzE(genjet_Px[0],genjet_Py[0],genjet_Pz[0],genjet_E[0]);
            gj_m2_orig.SetPxPyPzE(genjet_Px[1],genjet_Py[1],genjet_Pz[1],genjet_E[1]);

            tempGenEMissP4=TLorentzVector(0,0,0,0);
            #true sums over all visible particles, leptons and photons not removed, thus gives us the full total MET
            tempGenEMissP4.SetPxPyPzE(-(ientry.true_Px),-(ientry.true_Py),-(ientry.true_Pz),sqrt(pow(ientry.true_Px,2)+pow(ientry.true_Py,2)+pow(ientry.true_Pz,2)))
            tempGenMETP4=TLorentzVector(0,0,0,0)
            tempGenMETP4.SetPxPyPzE(tempGenEMissP4.Px(),tempGenEMissP4.Py(),0,tempGenEMissP4.Pt())
            
            tempGenEMissCorrP4=TLorentzVector(0,0,0,0)
            
            gj1_EMissProjVecProp=TLorentzVector(0,0,0,0)
            gj2_EMissProjVecProp=TLorentzVector(0,0,0,0)
            
            if(tempGenMETP4.Pt()>0):
                gj1_METProjProp=(tempGenMETP4.Vect().Dot(gj_m1_orig.Vect().Unit()))*gj_m1_orig.Vect().Unit();
                #check if gj1 and MET in same hemisphere
                if(tempGenMETP4.Vect().Unit().Dot(gj_m1_orig.Vect().Unit())>0):
                    gj1_EMissProjVecProp.SetPxPyPzE(gj1_METProjProp.Px(),gj1_METProjProp.Py(),gj1_METProjProp.Pt()*gj_m1_orig.Pz()/gj_m1_orig.Pt(),gj1_METProjProp.Pt()*gj_m1_orig.P()/gj_m1_orig.Pt());
                gj2_METProjProp=(tempGenMETP4.Vect().Dot(gj_m2_orig.Vect().Unit()))*gj_m2_orig.Vect().Unit();
                #check if gj2 and MET in same hemisphere
                if(tempGenMETP4.Vect().Unit().Dot(gj_m2_orig.Vect().Unit())>0):
                    gj2_EMissProjVecProp.SetPxPyPzE(gj2_METProjProp.Px(),gj2_METProjProp.Py(),gj2_METProjProp.Pt()*gj_m2_orig.Pz()/gj_m2_orig.Pt(),gj2_METProjProp.Pt()*gj_m2_orig.P()/gj_m2_orig.Pt());
            gj1_EMiss=TLorentzVector(0,0,0,0)
            gj1_EMiss=gj_m1_orig+gj1_EMissProjVecProp
            gj2_EMiss=TLorentzVector(0,0,0,0)
            gj2_EMiss=gj_m2_orig+gj2_EMissProjVecProp
            tempGenEMissCorrP4=gj1_EMiss+gj2_EMiss-gj_m1_orig-gj_m2_orig
            
            gj_m1=gj_m1_orig
            gj_m2=gj_m2_orig

            if use_sqrtJets :
                #for jet clustering isolated photons and leptons are removed beforehand, but they are included in sum of gennstructed particles
                sqrtS_eff_gen=(gj_m1+gj_m2).M()
            else :
                sqrtS_eff_gen=(tempTotGenP4-tempGenIsoPhP4).M()

            if use_EMissNeutrinoProjection :
                gj_m1=gj1_EMiss
                gj_m2=gj2_EMiss
                if use_sqrtJets :
                #in jet clustering isolated photons and leptons are removed beforehand, but they are included in sum of gennstructed particles
                    sqrtS_eff_gen=(gj_m1+gj_m2+tempGenEMissCorrP4).M()
                else :
                    sqrtS_eff_gen=(tempTotGenP4-tempGenIsoPhP4+tempGenEMissCorrP4).M()

            ind_genjetM1=0  
            ind_genjetM2=1           
            if(gj_m1.M()<gj_m2.M()) :
                ind_genjetM2=0  
                ind_genjetM1=1  
                temp=TLorentzVector(0,0,0,0)
                temp=gj_m1
                gj_m1=gj_m2
                gj_m2=temp

            fCut_pass_ellipse_genmass_cut=False
            if (gj_m1.M()>fCut_mass_1_min and (pow((gj_m2.M()-fCut_mass2_center)/fCut_mass2_radius,2)+pow((gj_m1.M()-fCut_mass1_center)/fCut_mass1_radius,2))<1.) :
                fCut_pass_ellipse_genmass_cut=True 


            pass_gen_mass_cut = True
            if performMassCuts and not fCut_pass_ellipse_genmass_cut :
                pass_gen_mass_cut=False
           
        t_var_nLeptons_gen[0] = n_IsoLep_gen

        t_var_sqrtS_gen[0]=sqrtS_eff_gen
        t_var_sqrtS_orig_gen[0]=(gj_m1_orig+gj_m2_orig).M()

        n_IsoPh_reco=0
        n_IsoLep_reco=0
        tempRecoIsoPhP4=TLorentzVector(0,0,0,0)
        tempRecoIsoLepP4=TLorentzVector(0,0,0,0)
        isoPartRecoDR10_E=ientry.isoPartRecoDR10_E
        isoPartRecoDR10_Px=ientry.isoPartRecoDR10_Px
        isoPartRecoDR10_Py=ientry.isoPartRecoDR10_Py
        isoPartRecoDR10_Pz=ientry.isoPartRecoDR10_Pz
        isoPartRecoDR10_relIso=ientry.isoPartRecoDR10_relIso
        isoPartRecoDR10_PDGID=ientry.isoPartRecoDR10_PDGID
        
        for i in range(len(isoPartRecoDR10_E)):
            if isoPartRecoDR10_PDGID[i]==22:
                if isoPartRecoDR10_relIso[i]<0.10:
                    temp=TLorentzVector(0,0,0,0)
                    temp.SetPxPyPzE(isoPartRecoDR10_Px[i],isoPartRecoDR10_Py[i],isoPartRecoDR10_Pz[i],isoPartRecoDR10_E[i])
                    tempRecoIsoPhP4+=temp
                    n_IsoPh_reco+=1
            elif (abs(isoPartRecoDR10_PDGID[i])==11 or abs(isoPartRecoDR10_PDGID[i])==13) :
                if isoPartRecoDR10_relIso[i]<0.10:
                    temp=TLorentzVector(0,0,0,0)
                    temp.SetPxPyPzE(isoPartRecoDR10_Px[i],isoPartRecoDR10_Py[i],isoPartRecoDR10_Pz[i],isoPartRecoDR10_E[i])
                    tempRecoIsoLepP4+=temp
                    n_IsoLep_reco+=1

        recojet_E=ientry.recojet_E
        if len(recojet_E)<2 : 
            #print "less than two jets in event, skip"
            pass_reco_selection=False
            continue;
        recojet_Px=ientry.recojet_Px
        recojet_Py=ientry.recojet_Py
        recojet_Pz=ientry.recojet_Pz
        t_var_reco_y21[0]=ientry.reco_y21_max
        t_var_reco_y32[0]=ientry.reco_y32_max
        t_var_reco_y43[0]=ientry.reco_y43_max

        rj_m1_true=TLorentzVector(0,0,0,0);
        rj_m2_true=TLorentzVector(0,0,0,0);
        rj_m1_true.SetPxPyPzE(recojet_Px[0],recojet_Py[0],recojet_Pz[0],recojet_E[0]);
        rj_m2_true.SetPxPyPzE(recojet_Px[1],recojet_Py[1],recojet_Pz[1],recojet_E[1]);

        rj_m1_orig=TLorentzVector(0,0,0,0);
        rj_m2_orig=TLorentzVector(0,0,0,0);
        rj_m1_orig.SetPxPyPzE(fCut_JESFactor*recojet_Px[0],fCut_JESFactor*recojet_Py[0],fCut_JESFactor*recojet_Pz[0],fCut_JESFactor*recojet_E[0]);
        rj_m2_orig.SetPxPyPzE(fCut_JESFactor*recojet_Px[1],fCut_JESFactor*recojet_Py[1],fCut_JESFactor*recojet_Pz[1],fCut_JESFactor*recojet_E[1]);

        tempRecoEMissP4=TLorentzVector(0,0,0,0);
        #totPFO is ALL PFOs, leptons and photons not removed, thus gives us the full total MET
        tempRecoEMissP4.SetPxPyPzE(-(ientry.totPFO_Px),-(ientry.totPFO_Py),-(ientry.totPFO_Pz),sqrt(pow(ientry.totPFO_Px,2)+pow(ientry.totPFO_Py,2)+pow(ientry.totPFO_Pz,2)))

        tempRecoMETP4=TLorentzVector(0,0,0,0)
        tempRecoMETP4.SetPxPyPzE(tempRecoEMissP4.Px(),tempRecoEMissP4.Py(),0,tempRecoEMissP4.Pt())

        tempRecoEMissCorrP4=TLorentzVector(0,0,0,0)

        rj1_EMissProjVecProp=TLorentzVector(0,0,0,0)
        rj2_EMissProjVecProp=TLorentzVector(0,0,0,0)
        rj1_EMissProjVecProp_true=TLorentzVector(0,0,0,0)
        rj2_EMissProjVecProp_true=TLorentzVector(0,0,0,0)

        if(tempRecoMETP4.Pt()>0):
            rj1_METProjProp=(tempRecoMETP4.Vect().Dot(rj_m1_orig.Vect().Unit()))*rj_m1_orig.Vect().Unit();
            #check if rj1 and MET in same hemisphere
            if(tempRecoMETP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit())>0):
                rj1_EMissProjVecProp.SetPxPyPzE(rj1_METProjProp.Px(),rj1_METProjProp.Py(),rj1_METProjProp.Pt()*rj_m1_orig.Pz()/rj_m1_orig.Pt(),rj1_METProjProp.Pt()*rj_m1_orig.P()/rj_m1_orig.Pt());
            rj2_METProjProp=(tempRecoMETP4.Vect().Dot(rj_m2_orig.Vect().Unit()))*rj_m2_orig.Vect().Unit();
            #check if rj2 and MET in same hemisphere
            if(tempRecoMETP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit())>0):
                rj2_EMissProjVecProp.SetPxPyPzE(rj2_METProjProp.Px(),rj2_METProjProp.Py(),rj2_METProjProp.Pt()*rj_m2_orig.Pz()/rj_m2_orig.Pt(),rj2_METProjProp.Pt()*rj_m2_orig.P()/rj_m2_orig.Pt());
            rj1_METProjProp_true=(tempRecoMETP4.Vect().Dot(rj_m1_true.Vect().Unit()))*rj_m1_true.Vect().Unit();
            #check if rj1 and MET in same hemisphere
            if(tempRecoMETP4.Vect().Unit().Dot(rj_m1_true.Vect().Unit())>0):
                rj1_EMissProjVecProp_true.SetPxPyPzE(rj1_METProjProp_true.Px(),rj1_METProjProp_true.Py(),rj1_METProjProp_true.Pt()*rj_m1_true.Pz()/rj_m1_true.Pt(),rj1_METProjProp_true.Pt()*rj_m1_true.P()/rj_m1_true.Pt());
            rj2_METProjProp_true=(tempRecoMETP4.Vect().Dot(rj_m2_true.Vect().Unit()))*rj_m2_true.Vect().Unit();
            #check if rj2 and MET in same hemisphere
            if(tempRecoMETP4.Vect().Unit().Dot(rj_m2_true.Vect().Unit())>0):
                rj2_EMissProjVecProp_true.SetPxPyPzE(rj2_METProjProp_true.Px(),rj2_METProjProp_true.Py(),rj2_METProjProp_true.Pt()*rj_m2_true.Pz()/rj_m2_true.Pt(),rj2_METProjProp_true.Pt()*rj_m2_true.P()/rj_m2_true.Pt());
        rj1_EMiss=TLorentzVector(0,0,0,0)
        rj1_EMiss=rj_m1_orig+rj1_EMissProjVecProp
        rj2_EMiss=TLorentzVector(0,0,0,0)
        rj2_EMiss=rj_m2_orig+rj2_EMissProjVecProp
        tempRecoEMissCorrP4=rj1_EMiss+rj2_EMiss-rj_m1_orig-rj_m2_orig
        rj1_EMiss_true=TLorentzVector(0,0,0,0)
        rj1_EMiss_true=rj_m1_true+rj1_EMissProjVecProp_true
        rj2_EMiss_true=TLorentzVector(0,0,0,0)
        rj2_EMiss_true=rj_m2_true+rj2_EMissProjVecProp_true
        tempRecoEMissCorrP4_true=rj1_EMiss_true+rj2_EMiss_true-rj_m1_true-rj_m2_true

        rj_m1=rj_m1_orig
        rj_m2=rj_m2_orig
        rj_m1_true=rj_m1_true
        rj_m2_true=rj_m2_true

        if use_sqrtJets :
            #for jet clustering isolated photons and leptons are removed beforehand, but they are included in sum of reconstructed particles
            sqrtS_eff_reco=(rj_m1+rj_m2).M()
        else :
            sqrtS_eff_reco=(tempTotRecoP4-tempRecoIsoPhP4).M()

        if use_EMissNeutrinoProjection :
            rj_m1=rj1_EMiss
            rj_m2=rj2_EMiss
            if use_sqrtJets :
                #in jet clustering isolated photons and leptons are removed beforehand, but they are included in sum of reconstructed particles
                sqrtS_eff_reco=(rj_m1+rj_m2+tempRecoEMissCorrP4).M()
            else :
                sqrtS_eff_reco=(tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M()

        sqrtS_eff_reco_true=sqrtS_eff_reco
        if use_EMissNeutrinoProjection :
            rj_m1_true=rj1_EMiss_true
            rj_m2_true=rj2_EMiss_true
            if use_sqrtJets :
                #in jet clustering isolated photons and leptons are removed beforehand, but they are included in sum of reconstructed particles
                sqrtS_eff_reco_true=(rj_m1_true+rj_m2_true+tempRecoEMissCorrP4_true).M()

        ind_jetM1=0  
        ind_jetM2=1           
        if(rj_m1.M()<rj_m2.M()) :
            ind_jetM2=0  
            ind_jetM1=1  
            temp=TLorentzVector(0,0,0,0)
            temp=rj_m1
            rj_m1=rj_m2
            rj_m2=temp

        if(rj_m1_true.M()<rj_m2_true.M()):
            temp=TLorentzVector(0,0,0,0)
            temp=rj_m1_true
            rj_m1_true=rj_m2_true
            rj_m2_true=temp

        t_var_sqrtS[0]=sqrtS_eff_reco
        t_var_sqrtS_orig[0]=(rj_m1_orig+rj_m2_orig).M()
        t_var_nLeptons[0] = n_IsoLep_reco
        

        #if abs(t_var_sqrtS[0]-t_var_sqrtS_parton[0])<1.0:
        #    print 'values damned close ',t_var_sqrtS[0],t_var_sqrtS_parton[0],(rj_m1+rj_m2+tempRecoEMissCorrP4).M(),(tempep+tempem).M(),(rj_m1_orig+rj_m2_orig).M(),(tempTotRecoP4-tempRecoIsoPhP4+tempRecoEMissCorrP4).M()

        #if t_var_sqrtS_parton[0]<t_var_sqrtS[0]:
        #    print 'for once the value is smaller',t_var_sqrtS_parton[0],t_var_sqrtS[0]

        pass_mass_cut=True

        fCut_pass_ellipse_mass_cut=False
        if (rj_m1.M()>fCut_mass_1_min and (pow((rj_m2.M()-fCut_mass2_center)/fCut_mass2_radius,2)+pow((rj_m1.M()-fCut_mass1_center)/fCut_mass1_radius,2))<1.) :
            fCut_pass_ellipse_mass_cut=True 
            #print "pass mass cuts",pass_mass_cut

 

        if performMassCuts and not fCut_pass_ellipse_mass_cut :
            #print "fail mass cuts",pass_mass_cut
            pass_mass_cut=False
            pass_reco_selection=False

        pass_theta_cut=True

        if performThetaCuts and ((abs(degrees(rj_m1.Theta())-fCut_thetaRef))>fCut_thetaWindow or abs(degrees(rj_m1.Theta()-rj_m2.Theta()))>fCut_delta_theta) :
            print "fail theta cuts",pass_mass_cut
            pass_theta_cut=False
            pass_reco_selection=False

        #pass_gen_parton_cuts = False

        #if ((len(genjet_E)==2 and n_IsoLep_gen<m_cut_nLeptons and sqrtS_eff_gen>sqrtS_high and pass_gen_mass_cut) or t_var_sqrtS_parton[0]>sqrtS_high):
        #    pass_gen_parton_cuts=True
      
        if n_IsoLep_reco>=m_cut_nLeptons or sqrtS_eff_reco<sqrtS_high:
            #print "fail lepton and sqrtS cut "
            pass_reco_selection=False

        #if not pass_reco_selection:
        #    print 'fail sqrtS part,sqrtReco,nlep,masscut/thetacut/pass_parton ',(tempep+tempem).M(),sqrtS_eff_reco,n_IsoLep_reco,pass_mass_cut,pass_theta_cut,pass_parton_selection,pass_reco_selection

        #if (len(genjet_E)==2 and n_IsoLep_gen<m_cut_nLeptons and sqrtS_eff_gen>sqrtS_high and pass_gen_mass_cut):
        #if not pass_reco_selection:
        #    if n_IsoLep_reco<m_cut_nLeptons and sqrtS_eff_reco>sqrtS_high:
        #        print 'basic cuts survived, but still fail ',fCut_pass_ellipse_mass_cut,len(recojet_subjet_rfj_j_E),len(recojet_subjet_E),len(recojet_E)

        #if(sqrtS_eff_reco_true>sqrtS_high and sqrtS_eff_reco<sqrtS_high and (num_entry%100==0)):
        #    print "fail now orig sqrtS, after JES sqrtS",sqrtS_eff_reco_true,sqrtS_eff_reco,sqrtS_eff_reco/sqrtS_eff_reco_true,num_entry
        #    print "fail now orig jetE, after JES jetE",rj_m1_true.E(),rj_m1.E(),rj_m1.E()/rj_m1_true.E(),rj_m2_true.E(),rj_m2.E(),rj_m2.E()/rj_m2_true.E(),num_entry
        #    print "fail now orig jetMass, after JES jetMass",rj_m1_true.M(),rj_m1.M(),rj_m1.M()/rj_m1_true.M(),rj_m2_true.M(),rj_m2.M(),rj_m2.M()/rj_m2_true.M(),num_entry

        if pass_reco_selection:
            #print 'pass sqrtS part,sqrtReco,nlep,masscut/thetacut/pass_parton ',(tempep+tempem).M(),sqrtS_eff_reco,n_IsoLep_reco,pass_mass_cut,pass_theta_cut,pass_parton_selection
        #if t_var_sqrtS_parton[0]>sqrtS_high:
            if sqrtS_eff_reco>sqrtS_high : 
                #if (sqrtS_eff_reco>sqrtS_high and sqrtS_eff_reco_true<sqrtS_high and (num_entry%100==0)):
                #    print "pass now, orig sqrtS, after JES sqrtS",sqrtS_eff_reco_true,sqrtS_eff_reco,sqrtS_eff_reco/sqrtS_eff_reco_true,num_entry
                #    print "pass now, orig jetE, after JES jetE",rj_m1_true.E(),rj_m1.E(),rj_m1.E()/rj_m1_true.E(),rj_m2_true.E(),rj_m2.E(),rj_m2.E()/rj_m2_true.E(),num_entry
                #    print "pass now, orig jetMass, after JES jetMass",rj_m1_true.M(),rj_m1.M(),rj_m1.M()/rj_m1_true.M(),rj_m2_true.M(),rj_m2.M(),rj_m2.M()/rj_m2_true.M(),num_entry
                recojet_subjet_rfj_j_E=ientry.recojet_subjet_rfj_j_E
                #length check done at the beginning
                #if len(recojet_subjet_rfj_j_E)!=4 : 
                #print "not four subjets in event, skip"
                #continue;
                recojet_subjet_rfj_j_Px                    =ientry.recojet_subjet_rfj_j_Px
                recojet_subjet_rfj_j_Py                    =ientry.recojet_subjet_rfj_j_Py
                recojet_subjet_rfj_j_Pz                    =ientry.recojet_subjet_rfj_j_Pz
                #recojet_subjet_rfj_j_jetChargePt_kappa_0_30=ientry.recojet_subjet_rfj_j_jetChargePt_kappa_0_30
                recojet_subjet_rfj_j_jetindex              =ientry.recojet_subjet_rfj_j_jetindex
                recojet_subjet_rfj_j_subjetindex           =ientry.recojet_subjet_rfj_j_subjetindex
                recojet_subjet_rfj_j_BTag                  =ientry.recojet_subjet_rfj_j_BTag
                recojet_subjet_rfj_j_CTag                  =ientry.recojet_subjet_rfj_j_CTag
                recojet_subjet_rfj_j_OTag                  =ientry.recojet_subjet_rfj_j_OTag

                ind_j_BTagSumMax=recojet_subjet_rfj_j_jetindex[0]
                ind_j_CTagSumMax=recojet_subjet_rfj_j_jetindex[0]
                ind_j_LTagSumMax=recojet_subjet_rfj_j_jetindex[0]

                if(recojet_subjet_rfj_j_BTag[0]+recojet_subjet_rfj_j_BTag[1]) < (recojet_subjet_rfj_j_BTag[2]+recojet_subjet_rfj_j_BTag[3]) :
                    ind_j_BTagSumMax=recojet_subjet_rfj_j_jetindex[2]
                if(recojet_subjet_rfj_j_CTag[0]+recojet_subjet_rfj_j_CTag[1]) < (recojet_subjet_rfj_j_CTag[2]+recojet_subjet_rfj_j_CTag[3]) :
                    ind_j_CTagSumMax=recojet_subjet_rfj_j_jetindex[2]
                if(recojet_subjet_rfj_j_OTag[0]+recojet_subjet_rfj_j_OTag[1]) < (recojet_subjet_rfj_j_OTag[2]+recojet_subjet_rfj_j_OTag[3]) :
                    ind_j_LTagSumMax=recojet_subjet_rfj_j_jetindex[2]

                if(recojet_subjet_rfj_j_jetindex[0]!=recojet_subjet_rfj_j_jetindex[1]) or (recojet_subjet_rfj_j_jetindex[2]!=recojet_subjet_rfj_j_jetindex[3]):
                    print "from first principle these rfj indices should be the same, but why aren't they "<<recojet_subjet_rfj_j_jetindex[0],recojet_subjet_rfj_j_jetindex[1],recojet_subjet_rfj_j_jetindex[2],recojet_subjet_rfj_j_jetindex[3]


                #redefine index of BTagMax, which is more suited to reject backgrounds,
                #BTagMax is leading btag of subjets of jet 1
                ind_sj_BTagMax=-1
                ind_sj_CTagMax=-1
                ind_sj_LTagMax=-1

                BTagMax=0
                CTagMax=0
                LTagMax=0

                if recojet_subjet_rfj_j_jetindex[0]==ind_jetM1: 
                    BTagMax=recojet_subjet_rfj_j_BTag[0]
                    ind_sj_BTagMax=0
                    if recojet_subjet_rfj_j_BTag[1]>recojet_subjet_rfj_j_BTag[0]:
                        BTagMax=recojet_subjet_rfj_j_BTag[1]
                        ind_sj_BTagMax=1
                else:
                    BTagMax=recojet_subjet_rfj_j_BTag[2]
                    ind_sj_BTagMax=2
                    if recojet_subjet_rfj_j_BTag[3]>recojet_subjet_rfj_j_BTag[2]:
                        BTagMax=recojet_subjet_rfj_j_BTag[3]
                        ind_sj_BTagMax=3

                if recojet_subjet_rfj_j_jetindex[0]==ind_jetM1: 
                    CTagMax=recojet_subjet_rfj_j_BTag[0]
                    ind_sj_CTagMax=0
                    if recojet_subjet_rfj_j_BTag[1]>recojet_subjet_rfj_j_BTag[0]:
                        CTagMax=recojet_subjet_rfj_j_BTag[1]
                        ind_sj_CTagMax=1
                else:
                    CTagMax=recojet_subjet_rfj_j_BTag[2]
                    ind_sj_CTagMax=2
                    if recojet_subjet_rfj_j_BTag[3]>recojet_subjet_rfj_j_BTag[2]:
                        CTagMax=recojet_subjet_rfj_j_BTag[3]
                        ind_sj_CTagMax=3

                recojet_subjet_Px                    =ientry.recojet_subjet_Px
                recojet_subjet_Py                    =ientry.recojet_subjet_Py
                recojet_subjet_Pz                    =ientry.recojet_subjet_Pz
                recojet_subjet_jetChargePt_kappa_0_30=ientry.recojet_subjet_jetChargePt_kappa_0_30
                recojet_subjet_jetChargeE_kappa_0_30=ientry.recojet_subjet_jetChargeE_kappa_0_30
                recojet_subjet_jetChargePt_kappa_0_25=ientry.recojet_subjet_jetChargePt_kappa_0_25
                recojet_subjet_jetChargeE_kappa_0_25=ientry.recojet_subjet_jetChargeE_kappa_0_25
                recojet_subjet_jetindex              =ientry.recojet_subjet_jetindex
                recojet_subjet_NCH                  =ientry.recojet_subjet_NCH
                recojet_subjet_CHFraction                 =ientry.recojet_subjet_CHFraction
                recojet_subjet_ElFraction                 =ientry.recojet_subjet_ElFraction
                recojet_subjet_MuFraction                 =ientry.recojet_subjet_MuFraction

 


                ind_j1_sj1=-1
                ind_j1_sj2=-1
                ind_j2_sj1=-1
                ind_j2_sj2=-1
                #jets always ordered as pair, check again to be sure
                if(recojet_subjet_jetindex[0] != recojet_subjet_jetindex[1] or  recojet_subjet_jetindex[2] != recojet_subjet_jetindex[3]):
                    print 'indices not in the expected order',recojet_subjet_jetindex[0],recojet_subjet_jetindex[1],recojet_subjet_jetindex[2],recojet_subjet_jetindex[3]

                if(recojet_subjet_jetindex[0]==ind_jetM1):
                    if recojet_subjet_E[0]>recojet_subjet_E[1]:
                        ind_j1_sj1=0
                        ind_j1_sj2=1
                    else:
                        ind_j1_sj1=1
                        ind_j1_sj2=0
                    if recojet_subjet_E[2]>recojet_subjet_E[3]:
                        ind_j2_sj1=2
                        ind_j2_sj2=3
                    else:
                        ind_j2_sj1=3
                        ind_j2_sj2=2
                else:
                    if recojet_subjet_E[0]>recojet_subjet_E[1]:
                        ind_j2_sj1=0
                        ind_j2_sj2=1
                    else:
                        ind_j2_sj1=1
                        ind_j2_sj2=0
                    if recojet_subjet_E[2]>recojet_subjet_E[3]:
                        ind_j1_sj1=2
                        ind_j1_sj2=3
                    else:
                        ind_j1_sj1=3
                        ind_j1_sj2=2   
                if  (ind_j1_sj1==-1 or ind_j1_sj2==-1 or ind_j2_sj1==-1 or ind_j2_sj2==-1):
                    print'indices not assigned properly',ind_j1_sj1,ind_j1_sj2,ind_j2_sj1,ind_j2_sj2

                #print 'is there anything correct ',recojet_subjet_E[0],recojet_subjet_Px[0],recojet_subjet_Py[0],recojet_subjet_Pz[0]

                t_var_jet1_sj1_E[0]=recojet_subjet_E[ind_j1_sj1]
                t_var_jet1_sj1_Px[0]=recojet_subjet_Px[ind_j1_sj1]
                t_var_jet1_sj1_Py[0]=recojet_subjet_Py[ind_j1_sj1]
                t_var_jet1_sj1_Pz[0]=recojet_subjet_Pz[ind_j1_sj1]
                t_var_jet1_sj1_jetChargePt_kappa_0_30[0]=recojet_subjet_jetChargePt_kappa_0_30[ind_j1_sj1]
                t_var_jet1_sj1_jetChargeE_kappa_0_30[0]=recojet_subjet_jetChargeE_kappa_0_30[ind_j1_sj1]
                t_var_jet1_sj1_jetChargePt_kappa_0_25[0]=recojet_subjet_jetChargePt_kappa_0_30[ind_j1_sj1]
                t_var_jet1_sj1_jetChargeE_kappa_0_25[0]=recojet_subjet_jetChargeE_kappa_0_25[ind_j1_sj1]
                t_var_jet1_sj1_chFrac[0]=recojet_subjet_CHFraction[ind_j1_sj1]+recojet_subjet_ElFraction[ind_j1_sj1]+recojet_subjet_MuFraction[ind_j1_sj1]
                t_var_jet1_sj1_nTracks[0]=recojet_subjet_NCH[ind_j1_sj1]
                t_var_jet1_sj2_E[0]=recojet_subjet_E[ind_j1_sj2]
                t_var_jet1_sj2_Px[0]=recojet_subjet_Px[ind_j1_sj2]
                t_var_jet1_sj2_Py[0]=recojet_subjet_Py[ind_j1_sj2]
                t_var_jet1_sj2_Pz[0]=recojet_subjet_Pz[ind_j1_sj2]
                t_var_jet1_sj2_jetChargePt_kappa_0_30[0]=recojet_subjet_jetChargePt_kappa_0_30[ind_j1_sj2]
                t_var_jet1_sj2_jetChargeE_kappa_0_30[0]=recojet_subjet_jetChargeE_kappa_0_30[ind_j1_sj2]
                t_var_jet1_sj2_jetChargePt_kappa_0_25[0]=recojet_subjet_jetChargePt_kappa_0_30[ind_j1_sj2]
                t_var_jet1_sj2_jetChargeE_kappa_0_25[0]=recojet_subjet_jetChargeE_kappa_0_25[ind_j1_sj2]
                t_var_jet1_sj2_chFrac[0]=recojet_subjet_CHFraction[ind_j1_sj2]+recojet_subjet_ElFraction[ind_j1_sj2]+recojet_subjet_MuFraction[ind_j1_sj2]
                t_var_jet1_sj2_nTracks[0]=recojet_subjet_NCH[ind_j2_sj2]
                rj1_sj1=TLorentzVector(0,0,0,0)
                rj1_sj2=TLorentzVector(0,0,0,0)
                rj1_sj1.SetPxPyPzE(t_var_jet1_sj1_Px[0],t_var_jet1_sj1_Py[0],t_var_jet1_sj1_Pz[0],t_var_jet1_sj1_E[0])
                rj1_sj2.SetPxPyPzE(t_var_jet1_sj2_Px[0],t_var_jet1_sj2_Py[0],t_var_jet1_sj2_Pz[0],t_var_jet1_sj2_E[0])
                t_var_jet1_dAlpha_sj1sj2[0]=degrees(rj1_sj1.Angle(rj1_sj2.Vect()))
                t_var_jet2_sj1_E[0]=recojet_subjet_E[ind_j2_sj1]
                t_var_jet2_sj1_Px[0]=recojet_subjet_Px[ind_j2_sj1]
                t_var_jet2_sj1_Py[0]=recojet_subjet_Py[ind_j2_sj1]
                t_var_jet2_sj1_Pz[0]=recojet_subjet_Pz[ind_j2_sj1]
                t_var_jet2_sj1_jetChargePt_kappa_0_30[0]=recojet_subjet_jetChargePt_kappa_0_30[ind_j2_sj1]
                t_var_jet2_sj1_jetChargeE_kappa_0_30[0]=recojet_subjet_jetChargeE_kappa_0_30[ind_j2_sj1]
                t_var_jet2_sj1_jetChargePt_kappa_0_25[0]=recojet_subjet_jetChargePt_kappa_0_30[ind_j2_sj1]
                t_var_jet2_sj1_jetChargeE_kappa_0_25[0]=recojet_subjet_jetChargeE_kappa_0_25[ind_j2_sj1]
                t_var_jet2_sj1_chFrac[0]=recojet_subjet_CHFraction[ind_j2_sj1]+recojet_subjet_ElFraction[ind_j2_sj1]+recojet_subjet_MuFraction[ind_j2_sj1]
                t_var_jet2_sj1_nTracks[0]=recojet_subjet_NCH[ind_j2_sj1]
                t_var_jet2_sj2_E[0]=recojet_subjet_E[ind_j2_sj2]
                t_var_jet2_sj2_Px[0]=recojet_subjet_Px[ind_j2_sj2]
                t_var_jet2_sj2_Py[0]=recojet_subjet_Py[ind_j2_sj2]
                t_var_jet2_sj2_Pz[0]=recojet_subjet_Pz[ind_j2_sj2]
                t_var_jet2_sj2_jetChargePt_kappa_0_30[0]=recojet_subjet_jetChargePt_kappa_0_30[ind_j2_sj2]
                t_var_jet2_sj2_jetChargeE_kappa_0_30[0]=recojet_subjet_jetChargeE_kappa_0_30[ind_j2_sj2]
                t_var_jet2_sj2_jetChargePt_kappa_0_25[0]=recojet_subjet_jetChargePt_kappa_0_30[ind_j2_sj2]
                t_var_jet2_sj2_jetChargeE_kappa_0_25[0]=recojet_subjet_jetChargeE_kappa_0_25[ind_j2_sj2]
                t_var_jet2_sj2_chFrac[0]=recojet_subjet_CHFraction[ind_j2_sj2]+recojet_subjet_ElFraction[ind_j2_sj2]+recojet_subjet_MuFraction[ind_j2_sj2]
                t_var_jet2_sj2_nTracks[0]=recojet_subjet_NCH[ind_j2_sj2]
                rj2_sj1=TLorentzVector(0,0,0,0)
                rj2_sj2=TLorentzVector(0,0,0,0)
                rj2_sj1.SetPxPyPzE(t_var_jet2_sj1_Px[0],t_var_jet2_sj1_Py[0],t_var_jet2_sj1_Pz[0],t_var_jet2_sj1_E[0])
                rj2_sj2.SetPxPyPzE(t_var_jet2_sj2_Px[0],t_var_jet2_sj2_Py[0],t_var_jet2_sj2_Pz[0],t_var_jet2_sj2_E[0])
                t_var_jet2_dAlpha_sj1sj2[0]=degrees(rj2_sj1.Angle(rj2_sj2.Vect()))



                #scale subjets
                sf_rj1=rj_m1.Px()/(rj1_sj1.Px()+rj1_sj2.Px())
                if abs(sf_rj1-1.)<1.e-6 :
                    sf_rj1=1.
                #if sf_rj1>1:
                #    print 'sf1 large ',sf_rj1
                #take recojets from above to boost subjets in com system
                temp_rj1_sj1=TLorentzVector(0,0,0,0)
                temp_rj1_sj2=TLorentzVector(0,0,0,0)
                if sf_rj1==1.:
                    temp_rj1_sj1=rj1_sj1
                    temp_rj1_sj2=rj1_sj2
                else:
                    temp_rj1_sj1.SetPxPyPzE(sf_rj1*rj1_sj1.Px(),sf_rj1*rj1_sj1.Py(),sf_rj1*rj1_sj1.Pz(),rj1_sj1.E()+rj1_sj1.E()/(rj1_sj1.E()+rj1_sj2.E())*(rj_m1.E()-(rj1_sj1.E()+rj1_sj2.E())))
                    temp_rj1_sj1.SetPxPyPzE(sf_rj1*rj1_sj2.Px(),sf_rj1*rj1_sj2.Py(),sf_rj1*rj1_sj2.Pz(),rj1_sj2.E()+rj1_sj2.E()/(rj1_sj1.E()+rj1_sj2.E())*(rj_m1.E()-(rj1_sj1.E()+rj1_sj2.E())))
                sf_rj2=rj_m2.Px()/(rj2_sj1.Px()+rj2_sj2.Px())
                if abs(sf_rj2-1.)<1.e-6 :
                    sf_rj2=1.
                #if sf_rj2>1:
                #    print 'sf2 large ',sf_rj2
                #take recojets from above to boost subjets in com system
                temp_rj2_sj1=TLorentzVector(0,0,0,0)
                temp_rj2_sj2=TLorentzVector(0,0,0,0)
                if sf_rj2==1.:
                    temp_rj2_sj1=rj2_sj1
                    temp_rj2_sj2=rj2_sj2
                else:
                    temp_rj2_sj1.SetPxPyPzE(sf_rj2*rj2_sj1.Px(),sf_rj2*rj2_sj1.Py(),sf_rj2*rj2_sj1.Pz(),rj2_sj1.E()+rj2_sj1.E()/(rj2_sj1.E()+rj2_sj2.E())*(rj_m2.E()-(rj2_sj1.E()+rj2_sj2.E())))
                    temp_rj2_sj2.SetPxPyPzE(sf_rj2*rj2_sj2.Px(),sf_rj2*rj2_sj2.Py(),sf_rj2*rj2_sj2.Pz(),rj2_sj2.E()+rj2_sj2.E()/(rj2_sj1.E()+rj2_sj2.E())*(rj_m2.E()-(rj2_sj1.E()+rj2_sj2.E())))


                temp_tot_Event_rj=TLorentzVector(0,0,0,0)
                temp_tot_Event_rj=rj_m1+rj_m2
                temp_ep_real_rj=TLorentzVector(0,0,0,0)
                temp_pz_ep_reco=temp_tot_Event_rj.Pz()*0.5-0.5*sqrt(temp_tot_Event_rj.E()*temp_tot_Event_rj.E()-(temp_tot_Event_rj.Pz()*temp_tot_Event_rj.Pz()+4*pow(root.TDatabasePDG.Instance().GetParticle(11).Mass(),2)))
                temp_ep_real_E_rj=sqrt(pow(temp_tot_Event_rj.Pz()*0.5-0.5*sqrt(temp_tot_Event_rj.E()*temp_tot_Event_rj.E()-(temp_tot_Event_rj.Pz()*temp_tot_Event_rj.Pz()+4*pow(root.TDatabasePDG.Instance().GetParticle(11).Mass(),2))),2)+pow(root.TDatabasePDG.Instance().GetParticle(11).Mass(),2))
                temp_ep_real_rj.SetPxPyPzE(0,0,temp_pz_ep_reco,temp_ep_real_E_rj)
                temp_ep_real_rj_boost=TLorentzVector(0,0,0,0)
                temp_ep_real_rj_boost.SetPxPyPzE(0,0,temp_pz_ep_reco,temp_ep_real_E_rj)
                boostE_tot_COM_rj=TVector3(0,0,0) 
                boostE_tot_COM_rj=-temp_tot_Event_rj.BoostVector()
                temp_ep_real_rj_boost.Boost(boostE_tot_COM_rj)


                temp_rj1=TLorentzVector(0,0,0,0)
                temp_rj1.SetPxPyPzE(rj_m1.Px(),rj_m1.Py(),rj_m1.Pz(),rj_m1.E())
                temp_rj1.Boost(boostE_tot_COM_rj)
                temp_rj2=TLorentzVector(0,0,0,0)
                temp_rj2.SetPxPyPzE(rj_m2.Px(),rj_m2.Py(),rj_m2.Pz(),rj_m2.E())
                temp_rj2.Boost(boostE_tot_COM_rj)
                plane_rj1_ep_approx_rj_boost=TVector3(0,0,0)
                plane_rj1_ep_approx_rj_boost=temp_rj1.Vect().Cross(temp_ep_real_rj_boost.Vect()).Unit()
                #if temp_rj2_sj1.M()<=0 or temp_rj2_sj2.M()<=0 :
                #    print 'wtf with masses rj2 ',temp_rj2_sj1.M(),temp_rj2_sj1.P(),temp_rj2_sj1.E(),temp_rj2_sj2.M(),temp_rj2_sj2.P(),temp_rj2_sj2.E(),sf_rj2,temp_tot_Event_rj.M()/temp_tot_Event_rj_orig.M(),temp_tot_Event_rj.M()
                #from now on ALL rj2 quantities are in COM regime
                temp_rj2_orig_corr=TLorentzVector(0,0,0,0)
                temp_rj2_orig_corr.SetPxPyPzE(rj_m2.Px(),rj_m2.Py(),rj_m2.Pz(),rj_m2.E())
                boostrj2_COM=TVector3(0,0,0) 
                boostrj2_COM=-temp_rj2_orig_corr.BoostVector()
                #boost subjets into rj2 COM
                temp_rj2_sj1.Boost(boostrj2_COM)
                temp_rj2_sj2.Boost(boostrj2_COM)
         

                #decide about positive quark here
                #print 'ep approx rj',temp_ep_approx_rj.Px(),temp_ep_approx_rj.Py(),temp_ep_approx_rj.Pz(),temp_ep_approx_rj.E()
                ind_rj2_sj_pos_jcE=-1
                if abs(t_var_jet2_sj1_jetChargeE_kappa_0_30[0])> abs(t_var_jet2_sj2_jetChargeE_kappa_0_30[0]):
                    if t_var_jet2_sj1_jetChargeE_kappa_0_30[0]>0:
                        ind_rj2_sj_pos_jcE=1
                    else:
                        ind_rj2_sj_pos_jcE=2
                else:
                    if t_var_jet2_sj2_jetChargeE_kappa_0_30[0]>0:
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

          
                plane_rj2_rj2_sj_jcE_pos_boost=TVector3(0,0,0)
                plane_rj2_rj2_sj_jcE_pos=temp_rj2.Vect().Cross(temp_rj2_sj_jcE_pos.Vect()).Unit()

                t_var_costheta1_for_Atheta1[0]=cos(temp_rj2_sj_jcE_pos.Angle(temp_rj2.Vect()))
                t_var_costheta2_for_Atheta1theta2[0]=cos(temp_rj1.Angle(temp_ep_real_rj_boost.Vect()))
 
                phi_plane_rj1_ep_rj2_sj_jcE_pos=degrees(acos(plane_rj1_ep_approx_rj_boost.Dot(plane_rj2_rj2_sj_jcE_pos)))
                #plane is already a TVector3
                if(degrees(temp_rj2_sj_jcE_pos.Angle(plane_rj1_ep_approx_rj_boost))>90.):
                    phi_plane_rj1_ep_rj2_sj_jcE_pos=degrees(acos(plane_rj1_ep_approx_rj_boost.Dot(plane_rj2_rj2_sj_jcE_pos)))+180.
                t_var_phi_for_Aphis[0]=phi_plane_rj1_ep_rj2_sj_jcE_pos


                t_var_MET[0]=tempRecoMETP4.Pt()
                t_var_jet1_mass[0]=rj_m1.M()
                t_var_jet2_mass[0]=rj_m2.M()
                t_var_jet1_min_jet2_mass[0]=rj_m1.M()-rj_m2.M()
                t_var_jet1_BTag_rfj_BTagMax[0]=recojet_subjet_rfj_j_BTag[ind_sj_BTagMax]
                t_var_jet1_CTag_rfj_BTagMax[0]=recojet_subjet_rfj_j_CTag[ind_sj_BTagMax]
                t_var_jet1_LTag_rfj_BTagMax[0]=recojet_subjet_rfj_j_OTag[ind_sj_BTagMax]
                t_var_jet1_BTag_rfj_CTagMax[0]=recojet_subjet_rfj_j_BTag[ind_sj_CTagMax]
                t_var_jet1_CTag_rfj_CTagMax[0]=recojet_subjet_rfj_j_CTag[ind_sj_CTagMax]
                t_var_jet1_LTag_rfj_CTagMax[0]=recojet_subjet_rfj_j_OTag[ind_sj_CTagMax]
                t_var_jet1_E[0]=rj_m1.E()
                t_var_jet2_E[0]=rj_m2.E()
                t_var_jet1_Pt[0]=rj_m1.Pt()
                t_var_jet2_Pt[0]=rj_m2.Pt()
                t_var_jet1_theta[0]=degrees(rj_m1.Theta())
                t_var_jet2_theta[0]=degrees(rj_m2.Theta())
                t_var_jet1_min_jet2_theta[0]=degrees(rj_m1.Theta()-rj_m2.Theta())
                t_var_jet1_phi[0]=degrees(rj_m1.Phi())
                t_var_jet2_phi[0]=degrees(rj_m2.Phi())
                t_var_dphi_j1j2[0]=degrees(rj_m1.DeltaPhi(rj_m2))
                t_var_angle_j1j2[0]=degrees(rj_m1.Angle(rj_m2.Vect()))

                t_var_jet1_d21[0]=ientry.recojet_dij_21_max[ind_jetM1]
                t_var_jet1_d32[0]=ientry.recojet_dij_32_max[ind_jetM1]
                t_var_jet1_d43[0]=ientry.recojet_dij_43_max[ind_jetM1]
                t_var_jet2_d21[0]=ientry.recojet_dij_21_max[ind_jetM2]
                t_var_jet2_d32[0]=ientry.recojet_dij_32_max[ind_jetM2]
                t_var_jet2_d43[0]=ientry.recojet_dij_43_max[ind_jetM2]

                t_var_jet1_D2_beta1[0]=ientry.recojet_beta1_D2[ind_jetM1]
                t_var_jet2_D2_beta1[0]=ientry.recojet_beta1_D2[ind_jetM2]
                t_var_jet1_D2_beta0_5[0]=ientry.recojet_beta0_5_D2[ind_jetM1]
                t_var_jet2_D2_beta0_5[0]=ientry.recojet_beta0_5_D2[ind_jetM2]
                t_var_jet1_C2_beta1[0]=ientry.recojet_beta1_C2[ind_jetM1]
                t_var_jet2_C2_beta1[0]=ientry.recojet_beta1_C2[ind_jetM2]
                t_var_jet1_C2_beta0_5[0]=ientry.recojet_beta0_5_C2[ind_jetM1]
                t_var_jet2_C2_beta0_5[0]=ientry.recojet_beta0_5_C2[ind_jetM2]

                t_var_jet1_C3_beta1[0]=ientry.recojet_beta1_C3[ind_jetM1]
                t_var_jet2_C3_beta1[0]=ientry.recojet_beta1_C3[ind_jetM2]
                t_var_jet1_N2_beta1[0]=ientry.recojet_beta1_N2[ind_jetM1]
                t_var_jet2_N2_beta1[0]=ientry.recojet_beta1_N2[ind_jetM2]
                t_var_jet1_N3_beta1[0]=ientry.recojet_beta1_N3[ind_jetM1]
                t_var_jet2_N3_beta1[0]=ientry.recojet_beta1_N3[ind_jetM2]


                t_var_jet1_D2_beta1_E_theta[0]=ientry.recojet_beta1_D2_E_theta[ind_jetM1]
                t_var_jet2_D2_beta1_E_theta[0]=ientry.recojet_beta1_D2_E_theta[ind_jetM2]
                t_var_jet1_C2_beta1_E_theta[0]=ientry.recojet_beta1_C2_E_theta[ind_jetM1]
                t_var_jet2_C2_beta1_E_theta[0]=ientry.recojet_beta1_C2_E_theta[ind_jetM2]
                t_var_jet1_C3_beta1_E_theta[0]=ientry.recojet_beta1_C3_E_theta[ind_jetM1]
                t_var_jet2_C3_beta1_E_theta[0]=ientry.recojet_beta1_C3_E_theta[ind_jetM2]
                t_var_jet1_N2_beta1_E_theta[0]=ientry.recojet_beta1_N2_E_theta[ind_jetM1]
                t_var_jet2_N2_beta1_E_theta[0]=ientry.recojet_beta1_N2_E_theta[ind_jetM2]
                t_var_jet1_N3_beta1_E_theta[0]=ientry.recojet_beta1_N3_E_theta[ind_jetM1]
                t_var_jet2_N3_beta1_E_theta[0]=ientry.recojet_beta1_N3_E_theta[ind_jetM2]

                eta1=-log(tan(radians(t_var_jet1_theta[0])/2.))
                eta2=-log(tan(radians(t_var_jet2_theta[0])/2.))
                temp_test1 = TLorentzVector(0,0,0,0)
                temp_test1.SetPtEtaPhiE(t_var_jet1_Pt[0],eta1,radians(t_var_jet1_phi[0]),t_var_jet1_E[0])
                temp_test2 = TLorentzVector(0,0,0,0)
                temp_test2.SetPtEtaPhiE(t_var_jet2_Pt[0],eta2,radians(t_var_jet2_phi[0]),t_var_jet2_E[0])
                test_mass1=(temp_test1+temp_test2).M()
                temp_test1.SetPtEtaPhiM(t_var_jet1_Pt[0],eta1,radians(t_var_jet1_phi[0]),t_var_jet1_mass[0])
                temp_test2.SetPtEtaPhiM(t_var_jet2_Pt[0],eta2,radians(t_var_jet2_phi[0]),t_var_jet2_mass[0])
                test_mass2=(temp_test1+temp_test2).M()
                #print 'sqrtS tests with jets ',test_mass1, test_mass2,(rj_m1+rj_m2).M()
                #print 'jet test',rj_m1.M(),t_var_jet1_mass[0],degrees(rj_m1.Theta()),t_var_jet1_theta[0],rj_m1.Pt(),t_var_jet1_Pt[0],rj_m1.E(),t_var_jet1_E[0],rj_m1.M(),t_var_jet1_mass[0]
                if(ientry.recojet_nsubjettiness1[ind_jetM1]!=0):
                    t_var_jet1_tau21[0]=ientry.recojet_nsubjettiness2[ind_jetM1]/ientry.recojet_nsubjettiness1[ind_jetM1]
                else:
                    t_var_jet1_tau21[0]=-1
                if(ientry.recojet_nsubjettiness1[ind_jetM2]!=0):
                    t_var_jet2_tau21[0]=ientry.recojet_nsubjettiness2[ind_jetM2]/ientry.recojet_nsubjettiness1[ind_jetM2]
                else:
                    t_var_jet2_tau21[0]=-1
                if(ientry.recojet_nsubjettiness2[ind_jetM1]!=0):
                    t_var_jet1_tau32[0]=ientry.recojet_nsubjettiness3[ind_jetM1]/ientry.recojet_nsubjettiness2[ind_jetM1]
                else:
                    t_var_jet1_tau32[0]=-1
                if(ientry.recojet_nsubjettiness2[ind_jetM2]!=0):
                    t_var_jet2_tau32[0]=ientry.recojet_nsubjettiness3[ind_jetM2]/ientry.recojet_nsubjettiness2[ind_jetM2]
                else:
                    t_var_jet2_tau32[0]=-1
                if ientry.recojet_nsubjettiness1_lrz[ind_jetM1]!=0:
                    t_var_jet1_tau21_E_theta[0]=ientry.recojet_nsubjettiness2_lrz[ind_jetM1]/ientry.recojet_nsubjettiness1_lrz[ind_jetM1]
                else:
                    t_var_jet1_tau21_E_theta[0]=-1
                if ientry.recojet_nsubjettiness1_lrz[ind_jetM2]!=0:
                    t_var_jet2_tau21_E_theta[0]=ientry.recojet_nsubjettiness2_lrz[ind_jetM2]/ientry.recojet_nsubjettiness1_lrz[ind_jetM2]
                else:
                    t_var_jet2_tau21_E_theta[0]=-1
                if(ientry.recojet_nsubjettiness2_lrz[ind_jetM1]!=0):
                    t_var_jet1_tau32_E_theta[0]=ientry.recojet_nsubjettiness3[ind_jetM1]/ientry.recojet_nsubjettiness2_lrz[ind_jetM1]
                else:
                    t_var_jet1_tau32_E_theta[0]=-1
                if(ientry.recojet_nsubjettiness2_lrz[ind_jetM2]!=0):
                    t_var_jet2_tau32_E_theta[0]=ientry.recojet_nsubjettiness3[ind_jetM2]/ientry.recojet_nsubjettiness2_lrz[ind_jetM2]
                else:
                    t_var_jet2_tau32_E_theta[0]=-1

                #print 'get to filling of weight ',t_var_eventWeight

                t_var_jet1_sj1_closestMatch[0]=-1
                t_var_jet1_sj2_closestMatch[0]=-1
                t_var_jet2_sj1_closestMatch[0]=-1
                t_var_jet2_sj2_closestMatch[0]=-1
                t_var_jet1_sj1_Angle_closestMatch[0]=-1
                t_var_jet1_sj2_Angle_closestMatch[0]=-1
                t_var_jet2_sj1_Angle_closestMatch[0]=-1
                t_var_jet2_sj2_Angle_closestMatch[0]=-1

                t_var_jet1_sj1_decMatch[0]=-1
                t_var_jet1_sj2_decMatch[0]=-1
                t_var_jet2_sj1_decMatch[0]=-1
                t_var_jet2_sj2_decMatch[0]=-1
                t_var_jet1_sj1_Angle_decMatch[0]=-1
                t_var_jet1_sj2_Angle_decMatch[0]=-1
                t_var_jet2_sj1_Angle_decMatch[0]=-1
                t_var_jet2_sj2_Angle_decMatch[0]=-1


                #check if MC info should be used, aka if it is available for signal sample
                if usePartonInfo:
                    #print 'angle leading jet H/Z ',rj_m1.Angle(tempHP4.Vect()),rj_m1.Angle(tempZP4.Vect()),rj1_sj1.Angle(tempH_b.Vect()),rj1_sj1.Angle(tempH_bbar.Vect()),rj1_sj1.Angle(tempZ_q_neg.Vect()),rj1_sj1.Angle(tempZ_q_pos.Vect())
                    if(rj_m1.Angle(tempHP4.Vect())< rj_m1.Angle(tempZP4.Vect())):
                        #jet1 is matched to H, jet2 to Z
                        if rj1_sj1.Angle(tempH_b.Vect())<rj1_sj1.Angle(tempH_bbar.Vect()) :
                            t_var_jet1_sj1_closestMatch[0]=0
                            t_var_jet1_sj1_Angle_closestMatch[0]=degrees(rj1_sj1.Angle(tempH_b.Vect()))
                        else:
                            t_var_jet1_sj1_closestMatch[0]=1
                            t_var_jet1_sj1_Angle_closestMatch[0]=degrees(rj1_sj1.Angle(tempH_bbar.Vect()))
                        if rj1_sj2.Angle(tempH_b.Vect())<rj1_sj2.Angle(tempH_bbar.Vect()) :
                            t_var_jet1_sj2_closestMatch[0]=0
                            t_var_jet1_sj2_Angle_closestMatch[0]=degrees(rj1_sj2.Angle(tempH_b.Vect()))
                        else:
                            t_var_jet1_sj2_closestMatch[0]=1
                            t_var_jet1_sj2_Angle_closestMatch[0]=degrees(rj1_sj2.Angle(tempH_bbar.Vect()))
                        #check now if both subjets matched to same quark
                        if t_var_jet1_sj1_closestMatch[0]==t_var_jet1_sj2_closestMatch[0]:
                          if t_var_jet1_sj1_closestMatch[0]==0 :
                              if rj1_sj1.Angle(tempH_b.Vect())<rj1_sj2.Angle(tempH_b.Vect()):
                                  t_var_jet1_sj1_decMatch[0]=0
                                  t_var_jet1_sj2_decMatch[0]=1
                                  t_var_jet1_sj1_Angle_decMatch[0]=degrees(rj1_sj1.Angle(tempH_b.Vect()))
                                  t_var_jet1_sj2_Angle_decMatch[0]=degrees(rj1_sj2.Angle(tempH_bbar.Vect()))
                              else:
                                  t_var_jet1_sj1_decMatch[0]=1
                                  t_var_jet1_sj2_decMatch[0]=0
                                  t_var_jet1_sj1_Angle_decMatch[0]=degrees(rj1_sj1.Angle(tempH_bbar.Vect()))
                                  t_var_jet1_sj2_Angle_decMatch[0]=degrees(rj1_sj2.Angle(tempH_b.Vect()))
                          else:
                              if rj1_sj1.Angle(tempH_bbar.Vect())<rj1_sj2.Angle(tempH_bbar.Vect()):
                                  t_var_jet1_sj1_decMatch[0]=1
                                  t_var_jet1_sj2_decMatch[0]=0
                                  t_var_jet1_sj1_Angle_decMatch[0]=degrees(rj1_sj1.Angle(tempH_bbar.Vect()))
                                  t_var_jet1_sj2_Angle_decMatch[0]=degrees(rj1_sj2.Angle(tempH_b.Vect()))
                              else:
                                  t_var_jet1_sj1_decMatch[0]=0
                                  t_var_jet1_sj2_decMatch[0]=1
                                  t_var_jet1_sj1_Angle_decMatch[0]=degrees(rj1_sj1.Angle(tempH_b.Vect()))
                                  t_var_jet1_sj2_Angle_decMatch[0]=degrees(rj1_sj2.Angle(tempH_bbar.Vect()))
                        else:
                           t_var_jet1_sj1_decMatch[0]= t_var_jet1_sj1_closestMatch[0]
                           t_var_jet1_sj2_decMatch[0]= t_var_jet1_sj2_closestMatch[0]
                           t_var_jet1_sj1_Angle_decMatch[0]=t_var_jet1_sj1_Angle_closestMatch[0]
                           t_var_jet1_sj2_Angle_decMatch[0]=t_var_jet1_sj2_Angle_closestMatch[0]
                        #jet1 is matched to H, jet2 to Z
                        if rj2_sj1.Angle(tempZ_q_neg.Vect())<rj2_sj1.Angle(tempZ_q_pos.Vect()) :
                            t_var_jet2_sj1_closestMatch[0]=2
                            t_var_jet2_sj1_Angle_closestMatch[0]=degrees(rj2_sj1.Angle(tempZ_q_neg.Vect()))
                        else:
                            t_var_jet2_sj1_closestMatch[0]=3
                            t_var_jet2_sj1_Angle_closestMatch[0]=degrees(rj2_sj1.Angle(tempZ_q_pos.Vect()))
                        if rj2_sj2.Angle(tempZ_q_neg.Vect())<rj2_sj2.Angle(tempZ_q_pos.Vect()) :
                            t_var_jet2_sj2_closestMatch[0]=2
                            t_var_jet2_sj2_Angle_closestMatch[0]=degrees(rj2_sj2.Angle(tempZ_q_neg.Vect()))
                        else:
                            t_var_jet2_sj2_closestMatch[0]=3
                            t_var_jet2_sj2_Angle_closestMatch[0]=degrees(rj2_sj2.Angle(tempZ_q_pos.Vect()))
                        #check now if both subjets matched to same quark
                        if t_var_jet2_sj1_closestMatch[0]==t_var_jet2_sj2_closestMatch[0]:
                          if t_var_jet2_sj1_closestMatch[0]==2 :
                              if rj2_sj1.Angle(tempZ_q_neg.Vect())<rj2_sj2.Angle(tempZ_q_neg.Vect()):
                                  t_var_jet2_sj1_decMatch[0]=2
                                  t_var_jet2_sj2_decMatch[0]=3
                                  t_var_jet2_sj1_Angle_decMatch[0]=degrees(rj2_sj1.Angle(tempZ_q_neg.Vect()))
                                  t_var_jet2_sj2_Angle_decMatch[0]=degrees(rj2_sj2.Angle(tempZ_q_pos.Vect()))
                              else:
                                  t_var_jet2_sj1_decMatch[0]=3
                                  t_var_jet2_sj2_decMatch[0]=2
                                  t_var_jet2_sj1_Angle_decMatch[0]=degrees(rj2_sj1.Angle(tempZ_q_pos.Vect()))
                                  t_var_jet2_sj2_Angle_decMatch[0]=degrees(rj2_sj2.Angle(tempZ_q_neg.Vect()))
                          else:
                              if rj2_sj1.Angle(tempZ_q_pos.Vect())<rj2_sj2.Angle(tempZ_q_pos.Vect()):
                                  t_var_jet2_sj1_decMatch[0]=3
                                  t_var_jet2_sj2_decMatch[0]=2
                                  t_var_jet2_sj1_Angle_decMatch[0]=degrees(rj2_sj1.Angle(tempZ_q_pos.Vect()))
                                  t_var_jet2_sj2_Angle_decMatch[0]=degrees(rj2_sj2.Angle(tempZ_q_neg.Vect()))
                              else:
                                  t_var_jet2_sj1_decMatch[0]=2
                                  t_var_jet2_sj2_decMatch[0]=3
                                  t_var_jet2_sj1_Angle_decMatch[0]=degrees(rj2_sj1.Angle(tempZ_q_neg.Vect()))
                                  t_var_jet2_sj2_Angle_decMatch[0]=degrees(rj2_sj2.Angle(tempZ_q_pos.Vect()))
                        else:
                           t_var_jet2_sj1_decMatch[0]= t_var_jet2_sj1_closestMatch[0]
                           t_var_jet2_sj2_decMatch[0]= t_var_jet2_sj2_closestMatch[0]
                           t_var_jet2_sj1_Angle_decMatch[0]=t_var_jet2_sj1_Angle_closestMatch[0]
                           t_var_jet2_sj2_Angle_decMatch[0]=t_var_jet2_sj2_Angle_closestMatch[0]
                    else:
                        #discuss cases where leading reconstructed jet matched to Z
                        #jet1 is matched to Z, jet2 to H
                        if rj1_sj1.Angle(tempZ_q_neg.Vect())<rj1_sj1.Angle(tempZ_q_pos.Vect()) :
                            t_var_jet1_sj1_closestMatch[0]=2
                            t_var_jet1_sj1_Angle_closestMatch[0]=degrees(rj1_sj1.Angle(tempZ_q_neg.Vect()))
                        else:
                            t_var_jet1_sj1_closestMatch[0]=3
                            t_var_jet1_sj1_Angle_closestMatch[0]=degrees(rj1_sj1.Angle(tempZ_q_pos.Vect()))
                        if rj1_sj2.Angle(tempZ_q_neg.Vect())<rj1_sj2.Angle(tempZ_q_pos.Vect()) :
                            t_var_jet1_sj2_closestMatch[0]=2
                            t_var_jet1_sj2_Angle_closestMatch[0]=degrees(rj1_sj2.Angle(tempZ_q_neg.Vect()))
                        else:
                            t_var_jet1_sj2_closestMatch[0]=3
                            t_var_jet1_sj2_Angle_closestMatch[0]=degrees(rj1_sj2.Angle(tempZ_q_pos.Vect()))
                        #check now if both subjets matched to same quark
                        if t_var_jet1_sj1_closestMatch[0]==t_var_jet1_sj2_closestMatch[0]:
                          if t_var_jet1_sj1_closestMatch[0]==2 :
                              if rj1_sj1.Angle(tempZ_q_neg.Vect())<rj1_sj2.Angle(tempZ_q_neg.Vect()):
                                  t_var_jet1_sj1_decMatch[0]=2
                                  t_var_jet1_sj2_decMatch[0]=3
                                  t_var_jet1_sj1_Angle_decMatch[0]=degrees(rj1_sj1.Angle(tempZ_q_neg.Vect()))
                                  t_var_jet1_sj2_Angle_decMatch[0]=degrees(rj1_sj2.Angle(tempZ_q_pos.Vect()))
                              else:
                                  t_var_jet1_sj1_decMatch[0]=3
                                  t_var_jet1_sj2_decMatch[0]=2
                                  t_var_jet1_sj1_Angle_decMatch[0]=degrees(rj1_sj1.Angle(tempZ_q_pos.Vect()))
                                  t_var_jet1_sj2_Angle_decMatch[0]=degrees(rj1_sj2.Angle(tempZ_q_neg.Vect()))
                          else:
                              if rj1_sj1.Angle(tempZ_q_pos.Vect())<rj1_sj2.Angle(tempZ_q_pos.Vect()):
                                  t_var_jet1_sj1_decMatch[0]=3
                                  t_var_jet1_sj2_decMatch[0]=2
                                  t_var_jet1_sj1_Angle_decMatch[0]=degrees(rj1_sj1.Angle(tempZ_q_pos.Vect()))
                                  t_var_jet1_sj2_Angle_decMatch[0]=degrees(rj1_sj2.Angle(tempZ_q_neg.Vect()))
                              else:
                                  t_var_jet1_sj1_decMatch[0]=2
                                  t_var_jet1_sj2_decMatch[0]=3
                                  t_var_jet1_sj1_Angle_decMatch[0]=degrees(rj1_sj1.Angle(tempZ_q_neg.Vect()))
                                  t_var_jet1_sj2_Angle_decMatch[0]=degrees(rj1_sj2.Angle(tempZ_q_pos.Vect()))
                        else:
                           t_var_jet1_sj1_decMatch[0]= t_var_jet1_sj1_closestMatch[0]
                           t_var_jet1_sj2_decMatch[0]= t_var_jet1_sj2_closestMatch[0]
                           t_var_jet1_sj1_Angle_decMatch[0]=t_var_jet1_sj1_Angle_closestMatch[0]
                           t_var_jet1_sj2_Angle_decMatch[0]=t_var_jet1_sj2_Angle_closestMatch[0]
                        #jet1 is matched to Z, jet2 to Q
                        if rj2_sj1.Angle(tempH_b.Vect())<rj2_sj1.Angle(tempH_bbar.Vect()) :
                            t_var_jet2_sj1_closestMatch[0]=0
                            t_var_jet2_sj1_Angle_closestMatch[0]=degrees(rj2_sj1.Angle(tempH_b.Vect()))
                        else:
                            t_var_jet2_sj1_closestMatch[0]=1
                            t_var_jet2_sj1_Angle_closestMatch[0]=degrees(rj2_sj1.Angle(tempH_bbar.Vect()))
                        if rj2_sj2.Angle(tempH_b.Vect())<rj2_sj2.Angle(tempH_bbar.Vect()) :
                            t_var_jet2_sj2_closestMatch[0]=0
                            t_var_jet2_sj2_Angle_closestMatch[0]=degrees(rj2_sj2.Angle(tempH_b.Vect()))
                        else:
                            t_var_jet2_sj2_closestMatch[0]=1
                            t_var_jet2_sj2_Angle_closestMatch[0]=degrees(rj2_sj2.Angle(tempH_bbar.Vect()))
                        #check now if both subjets matched to same quark
                        if t_var_jet2_sj1_closestMatch[0]==t_var_jet2_sj2_closestMatch[0]:
                          if t_var_jet2_sj1_closestMatch[0]==0 :
                              if rj2_sj1.Angle(tempH_b.Vect())<rj2_sj2.Angle(tempH_b.Vect()):
                                  t_var_jet2_sj1_decMatch[0]=0
                                  t_var_jet2_sj2_decMatch[0]=1
                                  t_var_jet2_sj1_Angle_decMatch[0]=degrees(rj2_sj1.Angle(tempH_b.Vect()))
                                  t_var_jet2_sj2_Angle_decMatch[0]=degrees(rj2_sj2.Angle(tempH_bbar.Vect()))
                              else:
                                  t_var_jet2_sj1_decMatch[0]=1
                                  t_var_jet2_sj2_decMatch[0]=0
                                  t_var_jet2_sj1_Angle_decMatch[0]=degrees(rj2_sj1.Angle(tempH_bbar.Vect()))
                                  t_var_jet2_sj2_Angle_decMatch[0]=degrees(rj2_sj2.Angle(tempH_b.Vect()))
                          else:
                              if rj2_sj1.Angle(tempH_bbar.Vect())<rj2_sj2.Angle(tempH_bbar.Vect()):
                                  t_var_jet2_sj1_decMatch[0]=1
                                  t_var_jet2_sj2_decMatch[0]=0
                                  t_var_jet2_sj1_Angle_decMatch[0]=degrees(rj2_sj1.Angle(tempH_bbar.Vect()))
                                  t_var_jet2_sj2_Angle_decMatch[0]=degrees(rj2_sj2.Angle(tempH_b.Vect()))
                              else:
                                  t_var_jet2_sj1_decMatch[0]=0
                                  t_var_jet2_sj2_decMatch[0]=1
                                  t_var_jet2_sj1_Angle_decMatch[0]=degrees(rj2_sj1.Angle(tempH_b.Vect()))
                                  t_var_jet2_sj2_Angle_decMatch[0]=degrees(rj2_sj2.Angle(tempH_bbar.Vect()))
                        else:
                           t_var_jet2_sj1_decMatch[0]= t_var_jet2_sj1_closestMatch[0]
                           t_var_jet2_sj2_decMatch[0]= t_var_jet2_sj2_closestMatch[0]
                           t_var_jet2_sj1_Angle_decMatch[0]=t_var_jet2_sj1_Angle_closestMatch[0]
                           t_var_jet2_sj2_Angle_decMatch[0]=t_var_jet2_sj2_Angle_closestMatch[0]
                    #if t_var_jet1_sj1_Angle_decMatch[0]>90 or t_var_jet1_sj2_Angle_decMatch[0]>90 or t_var_jet2_sj1_Angle_decMatch[0]>90 or t_var_jet2_sj2_Angle_decMatch[0]>90:
                    #    print 'what went all wrong for the indices ',t_var_jet1_sj1_Angle_decMatch[0],t_var_jet1_sj2_Angle_decMatch[0],t_var_jet2_sj1_Angle_decMatch[0],t_var_jet2_sj2_Angle_decMatch[0],t_var_jet1_sj1_decMatch[0],t_var_jet1_sj2_decMatch[0],t_var_jet2_sj1_decMatch[0],t_var_jet2_sj2_decMatch[0]
                    #    print 'j1_sj1 wrt H_b/H_bbar/Z_q_neg/Z_q_pos',degrees(rj1_sj1.Angle(tempH_b.Vect())),degrees(rj1_sj1.Angle(tempH_bbar.Vect())),degrees(rj1_sj1.Angle(tempZ_q_neg.Vect())),degrees(rj1_sj1.Angle(tempZ_q_pos.Vect())),' H b-b, Z qq ',degrees(tempH_bbar.Angle(tempH_b.Vect())),degrees(tempZ_q_pos.Angle(tempZ_q_neg.Vect())),degrees(tempH_bbar.Angle(tempHP4.Vect())),degrees(tempH_b.Angle(tempHP4.Vect())),degrees(tempZ_q_neg.Angle(tempZP4.Vect())),degrees(tempZ_q_pos.Angle(tempZP4.Vect())),' ratios', tempH_bbar.E()/tempHP4.E(),tempH_b.E()/tempHP4.E(), tempZ_q_neg.E()/tempZP4.E(),tempZ_q_pos.E()/tempZP4.E(),tempHP4.E(),tempZP4.E(),' true values?', trueME_PDGID[4],trueME_E[4], trueME_PDGID[5],trueME_E[5]
                    #    print 'j1_sj2 wrt H_b/H_bbar/Z_q_neg/Z_q_pos',degrees(rj1_sj2.Angle(tempH_b.Vect())),degrees(rj1_sj2.Angle(tempH_bbar.Vect())),degrees(rj1_sj2.Angle(tempZ_q_neg.Vect())),degrees(rj1_sj2.Angle(tempZ_q_pos.Vect()))
                    #    print 'j2_sj1 wrt H_b/H_bbar/Z_q_neg/Z_q_pos',degrees(rj2_sj1.Angle(tempH_b.Vect())),degrees(rj2_sj1.Angle(tempH_bbar.Vect())),degrees(rj2_sj1.Angle(tempZ_q_neg.Vect())),degrees(rj2_sj1.Angle(tempZ_q_pos.Vect()))
                    #    print 'j2_sj2 wrt H_b/H_bbar/Z_q_neg/Z_q_pos',degrees(rj2_sj2.Angle(tempH_b.Vect())),degrees(rj2_sj2.Angle(tempH_bbar.Vect())),degrees(rj2_sj2.Angle(tempZ_q_neg.Vect())),degrees(rj2_sj2.Angle(tempZ_q_pos.Vect()))
                if fillGenLevel :
                    t_var_genInv_E[0] =ientry.true_inv_E
                    t_var_genInv_Px[0] =ientry.true_inv_Px
                    t_var_genInv_Py[0] =ientry.true_inv_Py
                    t_var_genInv_Pz[0] =ientry.true_inv_Pz
                    genjet_E  =ientry.genjet_E
                    genjet_Px  =ientry.genjet_Px
                    genjet_Py  =ientry.genjet_Py
                    genjet_Pz  =ientry.genjet_Pz
                    if len(genjet_E)==2 :
                        temp_gen=TLorentzVector(0,0,0,0)
                        temp_gen.SetPxPyPzE(genjet_Px[0],genjet_Py[0],genjet_Pz[0],genjet_E[0])
                        ind_genjetM1=0
                        ind_genjetM2=1
                        if(temp_gen.Angle(rj_m1.Vect())>temp_gen.Angle(rj_m2.Vect())):
                            ind_genjetM1=1
                            ind_genjetM2=0
                        genjet_subjet_E                  =ientry.genjet_subjet_E
                        if len(genjet_subjet_E)==4 : 
                            genjet_subjet_Px                    =ientry.genjet_subjet_Px
                            genjet_subjet_Px                    =ientry.genjet_subjet_Px
                            genjet_subjet_Py                    =ientry.genjet_subjet_Py
                            genjet_subjet_Pz                    =ientry.genjet_subjet_Pz
                            genjet_subjet_jetChargePt_kappa_0_30=ientry.genjet_subjet_jetChargePt_kappa_0_30
                            genjet_subjet_jetChargeE_kappa_0_30=ientry.genjet_subjet_jetChargeE_kappa_0_30
                            genjet_subjet_jetChargePt_kappa_0_25=ientry.genjet_subjet_jetChargePt_kappa_0_25
                            genjet_subjet_jetChargeE_kappa_0_25=ientry.genjet_subjet_jetChargeE_kappa_0_25
                            genjet_subjet_jetindex              =ientry.genjet_subjet_jetindex
                            genjet_subjet_NCH                  =ientry.genjet_subjet_NCH
                            genjet_subjet_CHFraction                 =ientry.genjet_subjet_CHFraction
                            genjet_subjet_ElFraction                 =ientry.genjet_subjet_ElFraction
                            genjet_subjet_MuFraction                 =ientry.genjet_subjet_MuFraction
                            ind_gj1_sj1=-1
                            ind_gj1_sj2=-1
                            ind_gj2_sj1=-1
                            ind_gj2_sj2=-1
                            if(genjet_subjet_jetindex[0]==ind_genjetM1):
                                if genjet_subjet_E[0]>genjet_subjet_E[1]:
                                    ind_gj1_sj1=0
                                    ind_gj1_sj2=1
                                else:
                                    ind_gj1_sj1=1
                                    ind_gj1_sj2=0
                                if genjet_subjet_E[2]>genjet_subjet_E[3]:
                                    ind_gj2_sj1=2
                                    ind_gj2_sj2=3
                                else:
                                    ind_gj2_sj1=3
                                    ind_gj2_sj2=2
                            else:
                                if genjet_subjet_E[0]>genjet_subjet_E[1]:
                                    ind_gj2_sj1=0
                                    ind_gj2_sj2=1
                                else:
                                    ind_gj2_sj1=1
                                    ind_gj2_sj2=0
                                if genjet_subjet_E[2]>genjet_subjet_E[3]:
                                    ind_gj1_sj1=2
                                    ind_gj1_sj2=3
                                else:
                                    ind_gj1_sj1=3
                                    ind_gj1_sj2=2   
                            if  (ind_gj1_sj1==-1 or ind_gj1_sj2==-1 or ind_gj2_sj1==-1 or ind_gj2_sj2==-1):
                                print 'genjet indices not assigned properly',ind_gj1_sj1,ind_gj1_sj2,ind_gj2_sj1,ind_gj2_sj2
                            t_var_genjet1_sj1_E[0]=genjet_subjet_E[ind_gj1_sj1]
                            t_var_genjet1_sj1_Px[0]=genjet_subjet_Px[ind_gj1_sj1]
                            t_var_genjet1_sj1_Py[0]=genjet_subjet_Py[ind_gj1_sj1]
                            t_var_genjet1_sj1_Pz[0]=genjet_subjet_Pz[ind_gj1_sj1]
                            t_var_genjet1_sj1_jetChargePt_kappa_0_30[0]=genjet_subjet_jetChargePt_kappa_0_30[ind_gj1_sj1]
                            t_var_genjet1_sj1_jetChargeE_kappa_0_30[0]=genjet_subjet_jetChargeE_kappa_0_30[ind_gj1_sj1]
                            t_var_genjet1_sj1_jetChargePt_kappa_0_25[0]=genjet_subjet_jetChargePt_kappa_0_30[ind_gj1_sj1]
                            t_var_genjet1_sj1_jetChargeE_kappa_0_25[0]=genjet_subjet_jetChargeE_kappa_0_25[ind_gj1_sj1]
                            t_var_genjet1_sj1_chFrac[0]=genjet_subjet_CHFraction[ind_gj1_sj1]+genjet_subjet_ElFraction[ind_gj1_sj1]+genjet_subjet_MuFraction[ind_gj1_sj1]
                            t_var_genjet1_sj1_nTracks[0]=genjet_subjet_NCH[ind_gj1_sj1]
                            t_var_genjet1_sj2_E[0]=genjet_subjet_E[ind_gj1_sj2]
                            t_var_genjet1_sj2_Px[0]=genjet_subjet_Px[ind_gj1_sj2]
                            t_var_genjet1_sj2_Py[0]=genjet_subjet_Py[ind_gj1_sj2]
                            t_var_genjet1_sj2_Pz[0]=genjet_subjet_Pz[ind_gj1_sj2]
                            t_var_genjet1_sj2_jetChargePt_kappa_0_30[0]=genjet_subjet_jetChargePt_kappa_0_30[ind_gj1_sj2]
                            t_var_genjet1_sj2_jetChargeE_kappa_0_30[0]=genjet_subjet_jetChargeE_kappa_0_30[ind_gj1_sj2]
                            t_var_genjet1_sj2_jetChargePt_kappa_0_25[0]=genjet_subjet_jetChargePt_kappa_0_30[ind_gj1_sj2]
                            t_var_genjet1_sj2_jetChargeE_kappa_0_25[0]=genjet_subjet_jetChargeE_kappa_0_25[ind_gj1_sj2]
                            t_var_genjet1_sj2_chFrac[0]=genjet_subjet_CHFraction[ind_gj1_sj2]+genjet_subjet_ElFraction[ind_gj1_sj2]+genjet_subjet_MuFraction[ind_gj1_sj2]
                            t_var_genjet1_sj2_nTracks[0]=genjet_subjet_NCH[ind_gj2_sj2]


                            t_var_genjet1_d21[0]=ientry.genjet_dij_21_max[ind_genjetM1]
                            t_var_genjet1_d32[0]=ientry.genjet_dij_32_max[ind_genjetM1]
                            t_var_genjet1_d43[0]=ientry.genjet_dij_43_max[ind_genjetM1]
                            t_var_genjet2_d21[0]=ientry.genjet_dij_21_max[ind_genjetM2]
                            t_var_genjet2_d32[0]=ientry.genjet_dij_32_max[ind_genjetM2]
                            t_var_genjet2_d43[0]=ientry.genjet_dij_43_max[ind_genjetM2]

                            gj1_sj1=TLorentzVector(0,0,0,0)
                            gj1_sj2=TLorentzVector(0,0,0,0)
                            gj1_sj1.SetPxPyPzE(t_var_genjet1_sj1_Px[0],t_var_genjet1_sj1_Py[0],t_var_genjet1_sj1_Pz[0],t_var_genjet1_sj1_E[0])
                            gj1_sj2.SetPxPyPzE(t_var_genjet1_sj2_Px[0],t_var_genjet1_sj2_Py[0],t_var_genjet1_sj2_Pz[0],t_var_genjet1_sj2_E[0])
                            t_var_genjet1_dAlpha_sj1sj2[0]=degrees(gj1_sj1.Angle(gj1_sj2.Vect()))
                            t_var_genjet1_dAlpha_sj1_rjsj1[0]=degrees(gj1_sj1.Angle(rj1_sj1.Vect()))
                            t_var_genjet1_dAlpha_sj1_rjsj2[0]=degrees(gj1_sj1.Angle(rj1_sj2.Vect()))
                            t_var_genjet1_dAlpha_H[0]=-10
                            t_var_genjet1_dAlpha_sj1_qmin[0]=-10
                            t_var_genjet1_dAlpha_sj1_qplus[0]=-10
                            if usePartonInfo: 
                                gj1=TLorentzVector(0,0,0,0)
                                gj1=gj1_sj1+gj1_sj2
                                t_var_genjet1_dAlpha_H[0]=degrees(gj1.Angle(tempHP4.Vect()))
                                if gj1.Angle(tempHP4.Vect())< gj1.Angle(tempZP4.Vect()):
                                    #gj2 matched to H, compare now to b(neg) and bbar(pos charge)
                                    t_var_genjet1_dAlpha_sj1_qmin[0]=degrees(gj1_sj1.Angle(tempH_b.Vect()))
                                    t_var_genjet1_dAlpha_sj1_qplus[0]=degrees(gj1_sj1.Angle(tempH_bbar.Vect()))
                                else:
                                    t_var_genjet1_dAlpha_sj1_qmin[0]=degrees(gj1_sj1.Angle(tempZ_q_neg.Vect()))
                                    t_var_genjet1_dAlpha_sj1_qplus[0]=degrees(gj1_sj1.Angle(tempZ_q_pos.Vect()))
                            t_var_genjet2_sj1_E[0]=genjet_subjet_E[ind_gj2_sj1]
                            t_var_genjet2_sj1_Px[0]=genjet_subjet_Px[ind_gj2_sj1]
                            t_var_genjet2_sj1_Py[0]=genjet_subjet_Py[ind_gj2_sj1]
                            t_var_genjet2_sj1_Pz[0]=genjet_subjet_Pz[ind_gj2_sj1]
                            t_var_genjet2_sj1_jetChargePt_kappa_0_30[0]=genjet_subjet_jetChargePt_kappa_0_30[ind_gj2_sj1]
                            t_var_genjet2_sj1_jetChargeE_kappa_0_30[0]=genjet_subjet_jetChargeE_kappa_0_30[ind_gj2_sj1]
                            t_var_genjet2_sj1_jetChargePt_kappa_0_25[0]=genjet_subjet_jetChargePt_kappa_0_30[ind_gj2_sj1]
                            t_var_genjet2_sj1_jetChargeE_kappa_0_25[0]=genjet_subjet_jetChargeE_kappa_0_25[ind_gj2_sj1]
                            t_var_genjet2_sj1_chFrac[0]=genjet_subjet_CHFraction[ind_gj2_sj1]+genjet_subjet_ElFraction[ind_gj2_sj1]+genjet_subjet_MuFraction[ind_gj2_sj1]
                            t_var_genjet2_sj1_nTracks[0]=genjet_subjet_NCH[ind_gj2_sj1]
                            t_var_genjet2_sj2_E[0]=genjet_subjet_E[ind_gj2_sj2]
                            t_var_genjet2_sj2_Px[0]=genjet_subjet_Px[ind_gj2_sj2]
                            t_var_genjet2_sj2_Py[0]=genjet_subjet_Py[ind_gj2_sj2]
                            t_var_genjet2_sj2_Pz[0]=genjet_subjet_Pz[ind_gj2_sj2]
                            t_var_genjet2_sj2_jetChargePt_kappa_0_30[0]=genjet_subjet_jetChargePt_kappa_0_30[ind_gj2_sj2]
                            t_var_genjet2_sj2_jetChargeE_kappa_0_30[0]=genjet_subjet_jetChargeE_kappa_0_30[ind_gj2_sj2]
                            t_var_genjet2_sj2_jetChargePt_kappa_0_25[0]=genjet_subjet_jetChargePt_kappa_0_30[ind_gj2_sj2]
                            t_var_genjet2_sj2_jetChargeE_kappa_0_25[0]=genjet_subjet_jetChargeE_kappa_0_25[ind_gj2_sj2]
                            t_var_genjet2_sj2_chFrac[0]=genjet_subjet_CHFraction[ind_gj2_sj2]+genjet_subjet_ElFraction[ind_gj2_sj2]+genjet_subjet_MuFraction[ind_gj2_sj2]
                            t_var_genjet2_sj2_nTracks[0]=genjet_subjet_NCH[ind_gj2_sj2]
                            gj2_sj1=TLorentzVector(0,0,0,0)
                            gj2_sj2=TLorentzVector(0,0,0,0)
                            gj2_sj1.SetPxPyPzE(t_var_genjet2_sj1_Px[0],t_var_genjet2_sj1_Py[0],t_var_genjet2_sj1_Pz[0],t_var_genjet2_sj1_E[0])
                            gj2_sj2.SetPxPyPzE(t_var_genjet2_sj2_Px[0],t_var_genjet2_sj2_Py[0],t_var_genjet2_sj2_Pz[0],t_var_genjet2_sj2_E[0])
                            t_var_genjet2_dAlpha_sj1sj2[0]=degrees(gj2_sj1.Angle(gj2_sj2.Vect()))
                            t_var_genjet2_dAlpha_sj1_rjsj1[0]=degrees(gj2_sj1.Angle(rj2_sj1.Vect()))
                            t_var_genjet2_dAlpha_sj1_rjsj2[0]=degrees(gj2_sj1.Angle(rj2_sj2.Vect()))
                            t_var_genjet2_dAlpha_H[0]=-10
                            t_var_genjet2_dAlpha_sj1_qmin[0]=-10
                            t_var_genjet2_dAlpha_sj1_qplus[0]=-10
                            if usePartonInfo: 
                                gj2=TLorentzVector(0,0,0,0)
                                gj2=gj2_sj1+gj2_sj2
                                t_var_genjet2_dAlpha_H[0]=degrees(gj2.Angle(tempHP4.Vect()))
                                if gj2.Angle(tempHP4.Vect())< gj2.Angle(tempZP4.Vect()):
                                    #gj2 matched to H, compare now to b(neg) and bbar(pos charge)
                                    t_var_genjet2_dAlpha_sj1_qmin[0]=degrees(gj2_sj1.Angle(tempH_b.Vect()))
                                    t_var_genjet2_dAlpha_sj1_qplus[0]=degrees(gj2_sj1.Angle(tempH_bbar.Vect()))
                                else:
                                    t_var_genjet2_dAlpha_sj1_qmin[0]=degrees(gj2_sj1.Angle(tempZ_q_neg.Vect()))
                                    t_var_genjet2_dAlpha_sj1_qplus[0]=degrees(gj2_sj1.Angle(tempZ_q_pos.Vect()))
                                    
                
            weight_total+=weight
        if(pass_reco_selection or pass_parton_selection) :
            weight_total2+=weight;
            mytree.Fill()

    print 'total events after all running ', weight_total," tree filling weight ", weight_total2,"weight times tree entries ",weight*(float(mytree.GetEntries())),"masscuts/thetacuts?",performMassCuts,performThetaCuts
    h_parton_theta1_Z_qpos_vs_Z.Write()
    h_parton_costheta1_Z_qpos_vs_Z.Write()
    h_parton_theta1_Z_qpos_vs_Z_sqrtS_750.Write()
    h_parton_costheta1_Z_qpos_vs_Z_sqrtS_750.Write()
    h_parton_theta1_Z_qpos_vs_Z_sqrtS_750_2500.Write()
    h_parton_costheta1_Z_qpos_vs_Z_sqrtS_750_2500.Write()
    h_parton_theta1_Z_qpos_vs_Z_sqrtS_2500.Write()
    h_parton_costheta1_Z_qpos_vs_Z_sqrtS_2500.Write()

    h_parton_theta1_Z_qneg_vs_Z.Write()
    h_parton_costheta1_Z_qneg_vs_Z.Write()
    h_parton_theta1_Z_qneg_vs_Z_sqrtS_750.Write()
    h_parton_costheta1_Z_qneg_vs_Z_sqrtS_750.Write()
    h_parton_theta1_Z_qneg_vs_Z_sqrtS_750_2500.Write()
    h_parton_costheta1_Z_qneg_vs_Z_sqrtS_750_2500.Write()
    h_parton_theta1_Z_qneg_vs_Z_sqrtS_2500.Write()
    h_parton_costheta1_Z_qneg_vs_Z_sqrtS_2500.Write()

    h_parton_theta1_H_bbar_vs_H.Write()
    h_parton_costheta1_H_bbar_vs_H.Write()
    h_parton_theta1_H_bbar_vs_H_sqrtS_750.Write()
    h_parton_costheta1_H_bbar_vs_H_sqrtS_750.Write()
    h_parton_theta1_H_bbar_vs_H_sqrtS_750_2500.Write()
    h_parton_costheta1_H_bbar_vs_H_sqrtS_750_2500.Write()
    h_parton_theta1_H_bbar_vs_H_sqrtS_2500.Write()
    h_parton_costheta1_H_bbar_vs_H_sqrtS_2500.Write()

    h_parton_theta1_H_b_vs_H.Write()
    h_parton_costheta1_H_b_vs_H.Write()
    h_parton_theta1_H_b_vs_H_sqrtS_750.Write()
    h_parton_costheta1_H_b_vs_H_sqrtS_750.Write()
    h_parton_theta1_H_b_vs_H_sqrtS_750_2500.Write()
    h_parton_costheta1_H_b_vs_H_sqrtS_750_2500.Write()
    h_parton_theta1_H_b_vs_H_sqrtS_2500.Write()
    h_parton_costheta1_H_b_vs_H_sqrtS_2500.Write()

    h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom.Write()
    h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom.Write()

    h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Write()
    h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Write()
    h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Write()
    h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Write()
    
    h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Write()
    h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Write()
    h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Write()
    h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Write()

    h_2D_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Write()
    h_2D_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Write()
    h_parton_costheta1_Z_qpos_Zcom.Write()
    h_parton_costheta2_H_ep_HZ_COM.Write()
    h_parton_phi_plane_Z_qpos_vs_plane_H_ep.Write()
    h_2D_parton_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Write()


    if usePartonInfo:
        A_pos_theta1_parton=h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral()
        B_neg_theta1_parton=h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral()
        alpha_theta1_HZ_error=sqrt(4.*abs(A_pos_theta1_parton)*abs(B_neg_theta1_parton)/pow((abs(A_pos_theta1_parton)+abs(B_neg_theta1_parton)),3))
        print 'A_theta1',(h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral()+h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral())/(abs(h_parton_pos_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral())+abs(h_parton_neg_sgncos2theta1_costheta1_Z_qpos_Zcom.Integral())), 'theta1_error',alpha_theta1_HZ_error
        
        A_pos_phi1_parton=h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
        B_neg_phi1_parton=h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
        alpha_phi1_HZ_error=sqrt(4.*abs(A_pos_phi1_parton)*abs(B_neg_phi1_parton)/pow((abs(A_pos_phi1_parton)+abs(B_neg_phi1_parton)),3))
        print 'A_phi1',(h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()+h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())/(abs(h_parton_pos_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())+abs(h_parton_neg_sgnsinphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())), 'phi1_error',alpha_phi1_HZ_error
        
        A_pos_phi2_parton=h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
        B_neg_phi2_parton=h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
        alpha_phi2_HZ_error=sqrt(4.*abs(A_pos_phi2_parton)*abs(B_neg_phi2_parton)/pow((abs(A_pos_phi2_parton)+abs(B_neg_phi2_parton)),3))
        print 'A_phi2',(h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()+h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())/(abs(h_parton_pos_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())+abs(h_parton_neg_sgnsin2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())), 'phi2_error',alpha_phi2_HZ_error
        
        A_pos_phi3_parton=h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
        B_neg_phi3_parton=h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
        alpha_phi3_HZ_error=sqrt(4.*abs(A_pos_phi3_parton)*abs(B_neg_phi3_parton)/pow((abs(A_pos_phi3_parton)+abs(B_neg_phi3_parton)),3))
        print 'A_phi3',(h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()+h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())/(abs(h_parton_pos_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())+abs(h_parton_neg_sgncosphi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())), 'phi3_error',alpha_phi3_HZ_error
        
        A_pos_phi4_parton=h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
        B_neg_phi4_parton=h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()
        alpha_phi4_HZ_error=sqrt(4.*abs(A_pos_phi4_parton)*abs(B_neg_phi4_parton)/pow((abs(A_pos_phi4_parton)+abs(B_neg_phi4_parton)),3))
        print 'A_phi4',(h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral()+h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())/(abs(h_parton_pos_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())+abs(h_parton_neg_sgncos2phi_phi_plane_Z_qpos_vs_plane_H_ep.Integral())), 'phi4_error',alpha_phi4_HZ_error
        
        A_pos_theta12_parton=h_2D_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Integral()
        B_neg_theta12_parton=h_2D_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Integral()
        alpha_theta12_HZ_error=sqrt(4.*abs(A_pos_theta12_parton)*abs(B_neg_theta12_parton)/pow((abs(A_pos_theta12_parton)+abs(B_neg_theta12_parton)),3))
        print 'A_theta1theta2',(h_2D_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Integral()+h_2D_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Integral())/(abs(h_2D_parton_pos_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Integral())+abs(h_2D_parton_neg_sgncostheta1_times_sgncostheta2_costheta1_Z_qpos_Zcom_vs_costheta2_H_ep_HZ_COM.Integral())), 'theta12_error',alpha_theta12_HZ_error

    return None


    

def process_event(i_final_histo_name_,i_input_file_name_,i_xsec_,i_lumi_, i_use_partonInfo_,i_bool_applyMassCuts_,i_bool_applyThetaCuts_,i_bool_fillGenLevel_,i_bool_performJESUnc_,i_JESUncValue_):
    print "at start of process event"         
    input_file_=root.TFile.Open(i_input_file_name_)
    #input_file2_=root.TFile.Open(i_input_file_name2_)
    lumi=i_lumi_
    xsec_=i_xsec_
    use_partonInfo=i_use_partonInfo_
    bool_performMassCut=i_bool_applyMassCuts_
    bool_performthetaCut=i_bool_applyThetaCuts_
    bool_fillGenLevel=i_bool_fillGenLevel_
    bool_performJESUnc=i_bool_performJESUnc_
    JESUncValue=i_JESUncValue_
 
    print 'process event, parton/signal/mass/massrect/theta/BTagCut/C2/D2 cut ',use_partonInfo,bool_performMassCut,bool_performthetaCut,bool_performJESUnc,JESUncValue
          
                            
    file_histogram = root.TFile(i_final_histo_name_, "RECREATE")

    _mytree = TTree('MVATrainingVariables', 'MVATrainingVariables')
 

    fill_background_histograms(input_file_,_mytree,xsec_,use_partonInfo,lumi,bool_performMassCut,bool_performthetaCut,bool_fillGenLevel,bool_performJESUnc,JESUncValue)
    #file_histogram.cd()
    file_histogram.Write()
    file_histogram.Close()


    CLICdpStyle()



  
    return None

def process_files():

    lumi_=4000.

    use_partonInfo_=True
    #performMassCuts_=False
    performMassCuts_=True
    performJESUnc_=True
    performThetaCuts_=False
    fillGenInfo_=True

    JESUncValue_=0.005

    cross_section_= 3.83
    #final_histo_name="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/HZqq_Analyzer_noMassCuts_190709_hzqq_13391_polm80_sqrtS_j1_j2_EMissProj_H_bb_Zquark_VLC7VtxRFJVLC7.root"
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_hzqq_13391_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/hzqq_13391_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_hzqq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"  
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005_withPartonHistos_AllEvents.root"  
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_hzqq_noMass_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005_withPartonHistos.root"  
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_hzqq_noMassCut_noThetaCut_MVATrainingTree_Test_OnlyPartonSelection.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_,fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_13394_to_13397_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qqqq_13394_to_13397_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qqqq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 902.
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_mqqqq_2_TeV_13696_to_13699_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_mqq_1TeV_13425_to_13428_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qqqq_mqqqq_2TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=  369.8 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_13399_to_13402_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_13399_to_13402_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 1269.
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_


    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_mqq_1_TeV_13425_to_13428_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_mqq_1TeV_13425_to_13428_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 170.8
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbcbbc_13094_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_bbcbbc_13094_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_bbcbbc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 9.2271753e-3 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbubbu_13095_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_bbubbu_13095_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_bbubbu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=  9.1731760e-3 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ddcyyc_13096_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ddcyyc_13096_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_ddcyyc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=  1.3757137
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_dduyyu_13097_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_dduyyu_13097_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_dduyyu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=  14.498909
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscbbc_13098_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_sscbbc_13098_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_sscbbc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=  12.499614
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscssc_13099_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_sscssc_13099_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_sscssc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=  1.1651315
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssussu_13123_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ssussu_13123_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_ssussu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=  1.2615661e-2
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssubbu_13292_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ssubbu_13292_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_ssubbu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=   5.4145233e-2 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycbbu_13318_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycbbu_13318_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yycbbu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=  13.394883
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycddu_13326_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycddu_13326_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yycddu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=   	2.0054737
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycssu_13323_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycssu_13323_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yycssu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 2.0248353
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyubbc_13320_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyubbc_13320_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yyubbc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=13.330064
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyuddc_13328_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyuddc_13328_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yyuddc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=2.0034170
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyussc_13325_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyussc_13325_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yyussc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=2.0189010
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_


    #NOW WE GO TO THE POSITIVE POLARIZATION CASES

    lumi_=1000.

    use_partonInfo_=True
    fillGenInfo_=True
    #use_partonInfo_=False
    cross_section_= 2.67
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_hzqq_13392_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/hzqq_13392_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_hzqq_noMassCut_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005_OnlyPartonSelection.root"  
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005_withPartonHistos_AllEvents.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_13393_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qqqq_13393_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qqqq_noCuts_noSqrtS.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root"
    cross_section_= 120.
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_mqqqq_2_TeV_13700_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7ee_qq_mqq_1_TeV_13429_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qqqq_mqqqq_2TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root"
    cross_section_=   49.2
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_13398_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_13398_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root"
    cross_section_= 786.
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_mqq_1_TeV_13429_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7ee_qq_mqq_1_TeV_13429_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root"
    cross_section_=  73.5 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_



    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbcbbc_13071_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_bbcbbc_13071_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_bbcbbc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 2.9986901e-3 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbubbu_13072_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_bbubbu_13072_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_bbubbu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 2.9825397e-3 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ddcyyc_13073_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ddcyyc_13073_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_ddcyyc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=1.7824610e-1
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_dduyyu_13074_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_dduyyu_13074_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_dduyyu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 5.0109474
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscbbc_13075_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_sscbbc_13075_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_sscbbc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=4.8938333
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscssc_13076_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_sscssc_13076_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_sscssc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 1.3677677e-1 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssussu_13077_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ssussu_13077_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_ssussu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 3.3776171e-3 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssubbu_13293_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ssubbu_13293_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_ssubbu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 2.3216638e-2 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycbbu_13322_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycbbu_13322_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yycbbu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 5.2101109
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycddu_13319_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycddu_13319_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yycddu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 4.0984879e-1 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycssu_13327_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycssu_13327_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yycssu_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=  4.1853929e-1 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyubbc_13324_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyubbc_13324_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yyubbc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 5.2070149
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyuddc_13321_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyuddc_13321_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yyuddc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_= 4.1203686e-1 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyussc_13329_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyussc_13329_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yyussc_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_JESUnc1_005.root" 
    cross_section_=4.2245034e-1
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_, fillGenInfo_,performJESUnc_,JESUncValue_)
    print 'finished file', final_histo_name_


    return None

process_files()



