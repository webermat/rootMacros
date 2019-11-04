from ROOT import gROOT, TCanvas, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor
import ROOT as root
from math import cos, sin, pi, degrees, radians, pow, sqrt,acos

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

def fill_HZ_histograms(file,xsec,hist_vec_reco_1D,lumi):
#, h_hist_parton, h_hist_vec_gen, h_hist_vec_reco_1D, usePartonInfo, x_sec, fill_partonInfo, fill_genInfo):
#    print "size of histogram vectors ",len(h_hist_parton),len(h_hist_vec_gen),len(h_hist_vec_reco_1D), " booleans ",usePartonInfo,fill_partonInfo,fill_genInfo, " xsec ",xsec 
    
    use_EMissNeutrinoProjection=True
    #here use total 4 vector, or 4 vector sum of jet 1 and 2 (see flag below)- isolated photon four vector plus correction with EMiss on both jet axes
    #mass cuts are then also done after projecting the EMiss
    use_sqrtJets=True #in this case use j1+j2-isolated photons and with upper flag still decide if EMiss projection on jets is performed
    use_MHMiss_over_PFOMiss = False #if recoil jet missing energy or PFO missing energy is used in the missing energy projection

    fCut_mass_1_min=0.
    fCut_mass1_center=126.
    fCut_mass2_center=92.5
    fCut_mass1_radius=35.
    fCut_mass2_radius=35.

    fCut_thetaWindow=70.
    fCut_thetaRef=90.
    fCut_delta_theta = 100.
    
    tree = file.Get("showerData")

    #hist = file2.Get("h_runstatistics")

    #if tree.GetEntries()!=hist.GetBinContent(1):
    #    print "tree_entries not hist content/diff", tree.GetEntries(),hist.GetBinContent(1), tree.GetEntries()-hist.GetBinContent(1)

    sqrtS_low = 750.0
    #sqrtS_high = 2500.0
    sqrtS_high = 2500.0
    #number of leptons need to be smaller than this number
    m_cut_nLeptons = 1
 
    weight = xsec*lumi/tree.GetEntries()
    print "tree-entries ",tree.GetEntries(), " weight ",weight, "xsec",xsec,"lumi",lumi,"total event",weight*tree.GetEntries(), 'to', xsec*lumi

    num_entry=-1
    num_total_exception=0
    num_count=0

    for ientry in tree:
        num_entry+=1

        if num_entry%(int(tree.GetEntries()/5.)) == 0:
            print "in entry ",num_entry

        tempTotEventP4=TLorentzVector(0,0,0,0)
        tempTotEventP4HHZ=TLorentzVector(0,0,0,0)
        tempZP4=TLorentzVector(0,0,0,0)
        tempH1P4=TLorentzVector(0,0,0,0)
        tempH2P4=TLorentzVector(0,0,0,0)
        tempHHP4=TLorentzVector(0,0,0,0)
        
        quark_Z_decays=0
        
        tempH1_b=TLorentzVector(0,0,0,0);
        tempH1_bbar=TLorentzVector(0,0,0,0);
        tempH2_b=TLorentzVector(0,0,0,0);
        tempH2_bbar=TLorentzVector(0,0,0,0);

        tempZ_q_pos=TLorentzVector(0,0,0,0)
        tempZ_q_neg=TLorentzVector(0,0,0,0)

        temp_H1_q1=TLorentzVector(0,0,0,0)
        temp_H1_q2=TLorentzVector(0,0,0,0)

        temp_H2_q1=TLorentzVector(0,0,0,0)
        temp_H2_q2=TLorentzVector(0,0,0,0)

        temp_H1_V1_q1=TLorentzVector(0,0,0,0)
        temp_H1_V1_q2=TLorentzVector(0,0,0,0)

        temp_H1_V2_q1=TLorentzVector(0,0,0,0)
        temp_H1_V2_q2=TLorentzVector(0,0,0,0)

        temp_H2_V1_q1=TLorentzVector(0,0,0,0)
        temp_H2_V1_q2=TLorentzVector(0,0,0,0)

        temp_H2_V2_q1=TLorentzVector(0,0,0,0)
        temp_H2_V2_q2=TLorentzVector(0,0,0,0)


        H1H2_decays_bbar = True

        trueME_E=ientry.trueME_E
        trueME_Px=ientry.trueME_Px
        trueME_Py=ientry.trueME_Py
        trueME_Pz=ientry.trueME_Pz
        trueME_PDGID=ientry.trueME_PDGID
 
        
        index_firstH=-1
        index_secondH=-1
        index_firstH=-1
        
        for i in range(len(trueME_E)):
            if trueME_PDGID[i]==25 and index_firstH==-1:
                index_firstH=i
            if trueME_PDGID[i]==25 and index_firstH!=-1 and index_secondH==-1 :
                index_secondH=i
            if trueME_PDGID[i]==23 and index_firstZ==-1:
                index_firstZ=i
            if index_firstZ!=1 and index_firstH!=-1 and index_secondH!=-1 :
                break

            #then 4 is Z,5 is H typically, but stupid enough we have exceptions for whatever reason
            if range(len(trueME_E)) > 7:
                if trueME_PDGID[4]!=23 or trueME_PDGID[5]!=25 :
                    num_total_exception+=1
                quark_Z_decays=abs(trueME_PDGID[index_firstH+1])

                if(abs(trueME_PDGID[index_firstH+3])!=5 or abs(trueME_PDGID[index_firstH+4])!=5) :
                   H_decays_bbar=False

            if(H_decays_bbar==False):
                continue;
        
            #indices 0-1 outgoing electrons after ISR and beam strahlung (i.e. the "real collision")
            #index 2,3 ISR photons of incoming electrons
            #index 4 is Z, index 5 is H
            #index 6 and 7 are quarks from Z
            #index 8 and 9 are decay products of H
            for i in range(len(trueME_E)):
                temp=TLorentzVector(0,0,0,0)
                temp.SetPxPyPzE(trueME_Px[i],trueME_Py[i],trueME_Pz[i],trueME_E[i])
                #print 'particle ',i,' pdg ',trueME_PDGID[i]
                if i<2:
                    tempTotEventP4+=temp
                #for vec_ind in recojet_subjet_rfj_j_E:
                #print "vec_entry ",vec_ind

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
            else:
                #quark index 6 is negatively charged
                if trueME_PDGID[index_firstH+1] != -trueME_PDGID[index_firstH+2]:
                    print 'pdg id should be the same, what is wrong ',trueME_PDGID[index_firstH+1],-trueME_PDGID[index_firstH+2]
                tempZ_q_neg.SetPxPyPzE(trueME_Px[index_firstH+1],trueME_Py[index_firstH+1],trueME_Pz[index_firstH+1],trueME_E[index_firstH+1])
                tempZ_q_pos.SetPxPyPzE(trueME_Px[index_firstH+2],trueME_Py[index_firstH+2],trueME_Pz[index_firstH+2],trueME_E[index_firstH+2])
            
            tempZ_first=TLorentzVector(0,0,0,0)
            tempZ_first.SetPxPyPzE(trueME_Px[index_firstZ],trueME_Py[index_firstZ],trueME_Pz[index_firstZ],trueME_E[index_firstZ])

            tempHP4=tempH_p1+tempH_p2
            tempZP4=tempZ_q_pos+tempZ_q_neg
            hist_sqrtS_1D[0].Fill(tempTotEventP4.M(),weight)

        if(H_decays_bbar==False):
            continue;

    print 'total events after all running signal histos', hist_vec_reco_1D[0].Integral(0,hist_vec_reco_1D[0].GetNbinsX()+1),"masscuts/thetacuts?/betacuts/C2/D2",performMassCuts,performThetaCuts,performBTagCuts,performC2Cuts,performD2Cuts,num_count
    return None


    

def process_event(i_final_histo_name_,i_input_file_name_,i_xsec_,i_lumi_):
    print "at start of process event"         
    input_file_=root.TFile.Open(i_input_file_name_)
    #input_file2_=root.TFile.Open(i_input_file_name2_)
    lumi=i_lumi_
    xsec_=i_xsec_
    print 'process event'
          
                            
    file_histogram = root.TFile(i_final_histo_name_, "RECREATE")

    n_bins_high=200;
    
    lim_energy_low=0.;
    lim_energy_high=3050.;
    
    h_sqrtS_e1_e2_effective = new TH1F("h_sqrtS_e1_e2_effective","", n_bins_high, lim_energy_low,lim_energy_high);
    h_H1_H2_E_sum = new TH1F("h_H_H2_E_sum","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_H1_E = new TH1F("h_H1_E","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_H2_E = new TH1F("h_H1_E","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_Z_E = new TH1F("h_Z_E_sum","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    
    lim_mass_low=0;
    lim_mass_high=500;
    
    h_H1_H2_mass = new TH1F("h_H_H2_E_sum","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    
    double lim_dalpha_low=0.;
    double lim_dalpha_high=180.;
    
    h_dalpha_H1_H2_comb_vs_Z = TH1F("h_dalpha_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dphi_H1_H2_comb_vs_Z = TH1F("h_dphi_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dtheta_H1_H2_comb_vs_Z = TH1F("h_dtheta_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    
    h_dalpha_H1_H2 = TH1F("h_dalpha_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dphi_H1_H2 = TH1F("h_dphi_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dtheta_H1_H2 = TH1F("h_dtheta_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    
    lim_dalpha_qqbar_low=0.;
    lim_dalpha_qqbar_high=30.;
    
    h_dalpha_H1_bbar = TH1F("h_dalpha_H1_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dphi_H1_bbar = TH1F("h_dphi_H1_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dtheta_H1_bbar = TH1F("h_dtheta_H1_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    
    h_dalpha_H2_bbar = TH1F("h_dalpha_H2_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dphi_H2_bbar = TH1F("h_dphi_H2_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dtheta_H2_bbar = TH1F("h_dtheta_H2_bbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    
    h_dalpha_H1_qqbar = TH1F("h_dalpha_H1_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dphi_H1_qqbar = TH1F("h_dphi_H1_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dtheta_H1_qqbar = TH1F("h_dtheta_H1_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    
    h_dalpha_H2_qqbar = TH1F("h_dalpha_H2_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dphi_H2_qqbar = TH1F("h_dphi_H2_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dtheta_H2_qqbar = TH1F("h_dtheta_H2_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    
    h_dalpha_Z_qqbar = TH1F("h_dalpha_Z_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dphi_Z_qqbar = TH1F("h_dphi_Z_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dtheta_Z_qqbar = TH1F("h_dtheta_Z_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    
    
    lim_dalpha_qqbar_allH_low=0.;
    lim_dalpha_qqbar_allH_high=100.;
    
    h_dalpha_max_allH_qqbar = TH1F("h_dalpha_max_allH_qqbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
    h_dphi_max_allH_qqbar = TH1F("h_dphi_max_allH_qqbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
    h_dtheta_max_allH_qqbar = TH1F("h_dtheta_max_allH_qqbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
    
    h_dalpha_max_allH_bbar = TH1F("h_dalpha_max_allH_bbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
    h_dphi_max_allH_bbar = TH1F("h_dphi_max_allH_bbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
    h_dtheta_max_allH_bbar = TH1F("h_dtheta_max_allH_bbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
    
    h_dalpha_allH_diboson_qqbar = TH1F("h_dalpha_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dphi_allH_diboson_qqbar = TH1F("h_dphi_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dtheta_allH_diboson_qqbar = TH1F("h_dtheta_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    
    lim_H1_E_over_allH_E_low=0.5;
    lim_H1_E_over_allH_allH_high=1.0;
    
    h_E_H1_over_E_allH = TH1F("h_E_H1_over_E_allH","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
    h_E_H1_over_E_allH_bbar = TH1F("h_E_H1_over_E_allH_bbar","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
    h_E_H1_over_E_allH_qqbar = TH1F("h_E_H1_over_E_allH_qqbar","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
    
    hist_vec_HZ_parton_list=[]
    
    hist_vec_HZ_parton.append(h_sqrtS_e1_e2_effective)
    hist_vec_HZ_parton.append(h_H1_H2_E_sum)
    hist_vec_HZ_parton.append(h_H1_E)
    hist_vec_HZ_parton.append(h_H2_E)
    hist_vec_HZ_parton.append(h_Z_E)
    #5 done here
    hist_vec_HZ_parton.append(h_H1_H2_mass)
    hist_vec_HZ_parton.append(h_dalpha_H1_H2_comb_vs_Z)
    hist_vec_HZ_parton.append(h_dphi_H1_H2_comb_vs_Z)
    hist_vec_HZ_parton.append(h_dtheta_H1_H2_comb_vs_Z)
    hist_vec_HZ_parton.append(h_dalpha_H1_H2)
    #10 done here
    hist_vec_HZ_parton.append(h_dphi_H1_H2)
    hist_vec_HZ_parton.append(h_dtheta_H1_H2)
    hist_vec_HZ_parton.append(h_dalpha_H1_bbar)
    hist_vec_HZ_parton.append(h_dphi_H1_bbar)
    hist_vec_HZ_parton.append(h_dtheta_H1_bbar)
    #15 done here
    hist_vec_HZ_parton.append(h_dalpha_H2_bbar)
    hist_vec_HZ_parton.append(h_dphi_H2_bbar)
    hist_vec_HZ_parton.append(h_dtheta_H2_bbar)
    hist_vec_HZ_parton.append(h_dalpha_H1_qqbar)
    hist_vec_HZ_parton.append(h_dphi_H1_qqbar)
    #20 done here
    hist_vec_HZ_parton.append(h_dtheta_H1_qqbar)
    hist_vec_HZ_parton.append(h_dalpha_H2_qqbar)
    hist_vec_HZ_parton.append(h_dphi_H2_qqbar)
    hist_vec_HZ_parton.append(h_dtheta_H2_qqbar)
    hist_vec_HZ_parton.append(h_dalpha_Z_qqbar)
    #25 done here
    hist_vec_HZ_parton.append(h_dphi_Z_qqbar)
    hist_vec_HZ_parton.append(h_dtheta_Z_qqbar)
    hist_vec_HZ_parton.append(h_dalpha_max_allH_qqbar)
    hist_vec_HZ_parton.append(h_dphi_max_allH_qqbar)
    hist_vec_HZ_parton.append(h_dtheta_max_allH_qqbar)
    #30 done here
    hist_vec_HZ_parton.append(h_dalpha_max_allH_bbar)
    hist_vec_HZ_parton.append(h_dphi_max_allH_bbar)
    hist_vec_HZ_parton.append(h_dtheta_max_allH_bbar)
    hist_vec_HZ_parton.append(h_dalpha_allH_diboson_qqbar)
    hist_vec_HZ_parton.append(h_dphi_allH_diboson_qqbar)
    #35 done here
    hist_vec_HZ_parton.append(h_dtheta_allH_diboson_qqbar)
    hist_vec_HZ_parton.append(h_E_H1_over_E_allH)
    hist_vec_HZ_parton.append(h_E_H1_over_E_allH_bbar)
    hist_vec_HZ_parton.append(h_E_H1_over_E_allH_qqbar)
    #39 in total


    for hist in hist_vec_HZ_parton:
        hist.Sumw2()


    print 'length of hist list ',len(hist_vec_HZ_parton)
    fill_signal_parton_histograms(input_file_,xsec_,hist_vec_HZ_parton,lumi)
    
    CLICdpStyle()


    file_histogram.Write()
    file_histogram.Close()
  
    return None

def process_files():

    lumi_=4000.
    print 'start processing of files'

    cross_section_= 6.06e-02
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/HHZStudy_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190904Prod/VLC7VtxRFJVLC7/polm80/test_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14_DR7_partonlevelOnly.root"  
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_)
    print 'finished file', final_histo_name_

    lumi_=1000.
    cross_section_= 4.23e-02
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/HHZStudy_hhz_14344_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190904Prod/VLC7VtxRFJVLC7/polp80/test_hhz_14344_polp80_3TeV_wO_CLIC_o3_v14_DR7_partonlevelOnly.root"  
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_)
    print 'finished file', final_histo_name_

    return None

process_files()



