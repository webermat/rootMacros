from ROOT import gROOT, TCanvas, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor
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


def fill_background_histograms(file,mytree,xsec,usePartonInfo,lumi,performMassCuts,performThetaCuts):
    print "do something"

    t_var_sqrtS  = array('f',[0])
    t_var_sqrtS_orig  = array('f',[0])
    t_var_MET = array('f',[0])
    t_var_jet1_mass = array('f',[0])
    t_var_jet2_mass = array('f',[0])
    t_var_jet1_min_jet2_mass = array('f',[0])
    t_var_jet1_BTag_rfj_BTagMax = array('f',[0])
    t_var_jet1_CTag_rfj_BTagMax = array('f',[0])
    t_var_jet1_LTag_rfj_BTagMax = array('f',[0])
    t_var_jet1_E = array('f',[0])
    t_var_jet2_E = array('f',[0])
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
    t_var_eventWeight = array('f',[0])
    mytree.Branch('sqrtS_j1_j2_EMiss', t_var_sqrtS , 'sqrtS_j1_j2_EMiss/F')
    mytree.Branch('sqrtS_j1_j2_orig', t_var_sqrtS_orig , 'sqrtS_j1_j2_orig/F')
    mytree.Branch('MET', t_var_MET , 'MET/F')
    mytree.Branch('jet1_mass', t_var_jet1_mass , 'jet1_mass/F')
    mytree.Branch('jet2_mass', t_var_jet2_mass , 'jet2_mass/F')
    mytree.Branch('jet1_min_jet2_mass', t_var_jet1_min_jet2_mass , 'jet1_min_jet2_mass/F')
    mytree.Branch('jet1_BTag_rfj_BTagMax', t_var_jet1_BTag_rfj_BTagMax , 'jet1_BTag_rfj_BTagMax/F')
    mytree.Branch('jet1_CTag_rfj_BTagMax', t_var_jet1_CTag_rfj_BTagMax , 'jet1_BTag_rfj_CTagMax/F')
    mytree.Branch('jet1_LTag_rfj_BTagMax', t_var_jet1_LTag_rfj_BTagMax , 'jet1_BTag_rfj_LTagMax/F')
    mytree.Branch('jet1_E', t_var_jet1_E , 'jet1_E/F')
    mytree.Branch('jet2_E', t_var_jet2_E , 'jet2_E/F')
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
    mytree.Branch('jet1_tau21', t_var_jet1_tau21 , 'jet1_tau12/F')
    mytree.Branch('jet2_tau21', t_var_jet2_tau21 , 'jet2_tau12/F')
    mytree.Branch('jet1_tau32', t_var_jet1_tau32 , 'jet1_tau32/F')
    mytree.Branch('jet2_tau32', t_var_jet2_tau32 , 'jet2_tau32/F')
    mytree.Branch('eventWeight', t_var_eventWeight , 'eventWeight/F')
  
    use_EMissNeutrinoProjection=True
    #here use total 4 vector, or 4 vector sum of jet 1 and 2 (see flag below)- isolated photon four vector plus correction with EMiss on both jet axes
    #mass cuts are then also done after projecting the EMiss
    use_sqrtJets=True #in this case use j1+j2-isolated photons and with upper flag still decide if EMiss projection on jets is performed
  
    fCut_mass_1_min=0.
    fCut_mass1_center=126.
    fCut_mass2_center=92.5
    fCut_mass1_radius=40.
    fCut_mass2_radius=40.


    fCut_thetaWindow=70.
    fCut_thetaRef=90.
    fCut_delta_theta = 100.
 


    fCut_thetaWindow=70.
    fCut_thetaRef=90.

    #hist = file2.Get("h_runstatistics")

    #if tree.GetEntries()!=hist.GetBinContent(1) :
    #    print "tree_entries not hist content/diff", tree.GetEntries(),hist.GetBinContent(1), tree.GetEntries()-hist.GetBinContent(1)

    sqrtS_high = 2500.0
    #number of leptons need to be smaller than this number
    m_cut_nLeptons = 1
 
    weight_total=0
    tree = file.Get("showerData")

    weight = xsec*lumi/tree.GetEntries()
    print "tree-entries ",tree.GetEntries(), " weight ",weight, "xsec",xsec,"lumi",lumi,"total original ",weight*tree.GetEntries(), 'to', xsec*lumi

    num_entry=-1

    for ientry in tree:
        num_entry+=1

        if num_entry%(int(tree.GetEntries()/5.)) == 0:
            print "sig BG in entry ",num_entry

        tempTotRecoP4=TLorentzVector(0,0,0,0);
        tempTotRecoP4.SetPxPyPzE(ientry.totPFO_Px,ientry.totPFO_Py,ientry.totPFO_Pz,ientry.totPFO_E)


        recojet_subjet_rfj_j_E=ientry.recojet_subjet_rfj_j_E
        recojet_subjet_E=ientry.recojet_subjet_E
 
        if len(recojet_subjet_rfj_j_E)<4 or len(recojet_subjet_E)<4 :
            #print "too few subjets or refined jets ",len(recojet_subjet_E),len(recojet_subjet_rfj_j_E),num_entry,tempTotRecoP4.E()
            continue

        reco_pass_mass_cuts=False

        quark_Z_decays = 0
        H_decays_bbar = True

        if usePartonInfo :
            trueME_E=ientry.trueME_E
            trueME_Px=ientry.trueME_Px
            trueME_Py=ientry.trueME_Py
            trueME_Pz=ientry.trueME_Pz
            trueME_PDGID=ientry.trueME_PDGID
            
            if range(len(trueME_E)) > 7:
                quark_Z_decays=abs(trueME_PDGID[6])
                if(abs(trueME_PDGID[8])!=5 or abs(trueME_PDGID[9])!=5) :
                    H_decays_bbar=False
            if(H_decays_bbar==False):
                continue;
        
 

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
            print "less than two jets in event, skip"
            continue;
        recojet_Px=ientry.recojet_Px
        recojet_Py=ientry.recojet_Py
        recojet_Pz=ientry.recojet_Pz

        rj_m1_orig=TLorentzVector(0,0,0,0);
        rj_m2_orig=TLorentzVector(0,0,0,0);
        rj_m1_orig.SetPxPyPzE(recojet_Px[0],recojet_Py[0],recojet_Pz[0],recojet_E[0]);
        rj_m2_orig.SetPxPyPzE(recojet_Px[1],recojet_Py[1],recojet_Pz[1],recojet_E[1]);

        tempRecoEMissP4=TLorentzVector(0,0,0,0);
        #totPFO is ALL PFOs, leptons and photons not removed, thus gives us the full total MET
        tempRecoEMissP4.SetPxPyPzE(-(ientry.totPFO_Px),-(ientry.totPFO_Py),-(ientry.totPFO_Pz),sqrt(pow(ientry.totPFO_Px,2)+pow(ientry.totPFO_Py,2)+pow(ientry.totPFO_Pz,2)))

        tempRecoMETP4=TLorentzVector(0,0,0,0)
        tempRecoMETP4.SetPxPyPzE(tempRecoEMissP4.Px(),tempRecoEMissP4.Py(),0,tempRecoEMissP4.Pt())

        tempRecoEMissCorrP4=TLorentzVector(0,0,0,0)

        rj1_EMissProjVecProp=TLorentzVector(0,0,0,0)
        rj2_EMissProjVecProp=TLorentzVector(0,0,0,0)

        if(tempRecoMETP4.Pt()>0):
            rj1_METProjProp=(tempRecoMETP4.Vect().Dot(rj_m1_orig.Vect().Unit()))*rj_m1_orig.Vect().Unit();
            #check if rj1 and MET in same hemisphere
            if(tempRecoMETP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit())>0):
                rj1_EMissProjVecProp.SetPxPyPzE(rj1_METProjProp.Px(),rj1_METProjProp.Py(),rj1_METProjProp.Pt()*rj_m1_orig.Pz()/rj_m1_orig.Pt(),rj1_METProjProp.Pt()*rj_m1_orig.P()/rj_m1_orig.Pt());
            rj2_METProjProp=(tempRecoMETP4.Vect().Dot(rj_m2_orig.Vect().Unit()))*rj_m2_orig.Vect().Unit();
            #check if rj2 and MET in same hemisphere
            if(tempRecoMETP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit())>0):
                rj2_EMissProjVecProp.SetPxPyPzE(rj2_METProjProp.Px(),rj2_METProjProp.Py(),rj2_METProjProp.Pt()*rj_m2_orig.Pz()/rj_m2_orig.Pt(),rj2_METProjProp.Pt()*rj_m2_orig.P()/rj_m2_orig.Pt());
        rj1_EMiss=TLorentzVector(0,0,0,0)
        rj1_EMiss=rj_m1_orig+rj1_EMissProjVecProp
        rj2_EMiss=TLorentzVector(0,0,0,0)
        rj2_EMiss=rj_m2_orig+rj2_EMissProjVecProp
        tempRecoEMissCorrP4=rj1_EMiss+rj2_EMiss-rj_m1_orig-rj_m2_orig

        rj_m1=rj_m1_orig
        rj_m2=rj_m2_orig

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

        ind_jetM1=0  
        ind_jetM2=1           
        if(rj_m1.M()<rj_m2.M()) :
            ind_jetM2=0  
            ind_jetM1=1  
            temp=TLorentzVector(0,0,0,0)
            temp=rj_m1
            rj_m1=rj_m2
            rj_m2=temp

        fCut_pass_ellipse_mass_cut=False
        if (rj_m1.M()>fCut_mass_1_min and (pow((rj_m2.M()-fCut_mass2_center)/fCut_mass2_radius,2)+pow((rj_m1.M()-fCut_mass1_center)/fCut_mass1_radius,2))<1.) :
            fCut_pass_ellipse_mass_cut=True 

        if performMassCuts and not fCut_pass_ellipse_mass_cut :
            continue
        #if performMassCuts and ((rj_m1.M()<fCut_mass_1_min or rj_m1.M()>fCut_mass_1_max) or (rj_m2.M()<fCut_mass_2_min or rj_m2.M()>fCut_mass_2_max)) :
        #     continue

        if performThetaCuts and ((abs(degrees(rj_m1.Theta())-fCut_thetaRef))>fCut_thetaWindow or abs(degrees(rj_m1.Theta()-rj_m2.Theta()))>fCut_delta_theta) :
             continue
        
 

	if(len(recojet_E)==2 and n_IsoLep_reco<m_cut_nLeptons):
    
            if sqrtS_eff_reco>sqrtS_high :
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




                t_var_sqrtS[0]=sqrtS_eff_reco
                t_var_sqrtS_orig[0]=(rj_m1_orig+rj_m2_orig).M()
                t_var_MET[0]=tempRecoMETP4.Pt()
                t_var_jet1_mass[0]=rj_m1.M()
                t_var_jet2_mass[0]=rj_m2.M()
                t_var_jet1_min_jet2_mass[0]=rj_m1.M()-rj_m2.M()
                t_var_jet1_BTag_rfj_BTagMax[0]=recojet_subjet_rfj_j_BTag[ind_sj_BTagMax]
                t_var_jet1_CTag_rfj_BTagMax[0]=recojet_subjet_rfj_j_CTag[ind_sj_BTagMax]
                t_var_jet1_LTag_rfj_BTagMax[0]=recojet_subjet_rfj_j_OTag[ind_sj_BTagMax]
                t_var_jet1_E[0]=rj_m1.E()
                t_var_jet2_E[0]=rj_m2.E()
                t_var_jet1_theta[0]=degrees(rj_m1.Theta())
                t_var_jet2_theta[0]=degrees(rj_m2.Theta())
                t_var_jet1_min_jet2_theta[0]=degrees(rj_m1.Theta()-rj_m2.Theta())
                t_var_jet1_phi[0]=degrees(rj_m1.Phi())
                t_var_jet2_phi[0]=degrees(rj_m2.Phi())
                t_var_dphi_j1j2[0]=degrees(rj_m1.DeltaPhi(rj_m2))
                t_var_angle_j1j2[0]=degrees(rj_m1.Angle(rj_m2.Vect()))
                t_var_jet1_D2_beta1[0]=ientry.recojet_beta1_D2[ind_jetM1]
                t_var_jet2_D2_beta1[0]=ientry.recojet_beta1_D2[ind_jetM2]
                t_var_jet1_D2_beta0_5[0]=ientry.recojet_beta0_5_D2[ind_jetM1]
                t_var_jet2_D2_beta0_5[0]=ientry.recojet_beta0_5_D2[ind_jetM2]
                t_var_jet1_C2_beta1[0]=ientry.recojet_beta1_C2[ind_jetM1]
                t_var_jet2_C2_beta1[0]=ientry.recojet_beta1_C2[ind_jetM2]
                t_var_jet1_C2_beta0_5[0]=ientry.recojet_beta0_5_C2[ind_jetM1]
                t_var_jet2_C2_beta0_5[0]=ientry.recojet_beta0_5_C2[ind_jetM2]
                t_var_jet1_tau21[0]=ientry.recojet_nsubjettiness2[ind_jetM1]/ientry.recojet_nsubjettiness1[ind_jetM1]
                t_var_jet2_tau21[0]=ientry.recojet_nsubjettiness2[ind_jetM2]/ientry.recojet_nsubjettiness1[ind_jetM2]
                if(ientry.recojet_nsubjettiness2[ind_jetM1]!=0):
                    t_var_jet1_tau32[0]=ientry.recojet_nsubjettiness3[ind_jetM1]/ientry.recojet_nsubjettiness2[ind_jetM1]
                else:
                    t_var_jet1_tau32[0]=-1
                if(ientry.recojet_nsubjettiness2[ind_jetM2]!=0):
                    t_var_jet2_tau32[0]=ientry.recojet_nsubjettiness3[ind_jetM2]/ientry.recojet_nsubjettiness2[ind_jetM2]
                else:
                    t_var_jet2_tau32[0]=-1
                t_var_eventWeight[0]=weight
                #print 'get to filling of weight ',t_var_eventWeight

                mytree.Fill()
                weight_total+=weight
    print 'total events after all running ', weight_total," masscuts/thetacuts?",performMassCuts,performThetaCuts
    return None


    

def process_event(i_final_histo_name_,i_input_file_name_,i_xsec_,i_lumi_, i_use_partonInfo_,i_bool_applyMassCuts_,i_bool_applyThetaCuts_):
    print "at start of process event"         
    input_file_=root.TFile.Open(i_input_file_name_)
    #input_file2_=root.TFile.Open(i_input_file_name2_)
    lumi=i_lumi_
    xsec_=i_xsec_
    use_partonInfo=i_use_partonInfo_
    bool_performMassCut=i_bool_applyMassCuts_
    bool_performthetaCut=i_bool_applyThetaCuts_
 
    print 'process event, parton/signal/mass/massrect/theta/BTagCut/C2/D2 cut ',use_partonInfo,bool_performMassCut,bool_performthetaCut
          
                            
    file_histogram = root.TFile(i_final_histo_name_, "RECREATE")

    _mytree = TTree('MVATrainingVariables', 'MVATrainingVariables')
 

    fill_background_histograms(input_file_,_mytree,xsec_,use_partonInfo,lumi,bool_performMassCut,bool_performthetaCut)
    #file_histogram.cd()
    file_histogram.Write()
    file_histogram.Close()


    CLICdpStyle()



  
    return None

def process_files():

    lumi_=4000.
    use_partonInfo_=True
    performMassCuts_=True
    performThetaCuts_=False

    cross_section_= 3.83
    #final_histo_name="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/HZqq_Analyzer_noMassCuts_190417_hzqq_13391_polm80_sqrtS_j1_j2_EMissProj_H_bb_Zquark_VLC7VtxRFJVLC7.root"
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_hzqq_13391_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/hzqq_13391_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_hzqq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"  
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_hzqq_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root"  
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_13394_to_13397_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    ##input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_qqqq_13394_to_13397_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qqqq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qqqq_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 902.
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_mqqqq_2_TeV_13696_to_13699_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_mqq_1TeV_13425_to_13428_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qqqq_mqqqq_2TeV_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=  369.8 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_13399_to_13402_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_13399_to_13402_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 1269.
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_


    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_mqq_1_TeV_13425_to_13428_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_mqq_1TeV_13425_to_13428_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 170.8
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbcbbc_13094_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_bbcbbc_13094_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_bbcbbc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 9.2271753e-3 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbubbu_13095_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_bbubbu_13095_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_bbubbu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=  9.1731760e-3 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ddcyyc_13096_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_ddcyyc_13096_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_ddcyyc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=  1.3757137
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_dduyyu_13097_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_dduyyu_13097_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_dduyyu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=  14.498909
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscbbc_13098_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_sscbbc_13098_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_sscbbc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=  12.499614
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscssc_13099_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_sscssc_13099_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_sscssc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=  1.1651315
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssussu_13123_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_ssussu_13123_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_ssussu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=  1.2615661e-2
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssubbu_13292_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_ssubbu_13292_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_ssubbu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=   5.4145233e-2 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycbbu_13318_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yycbbu_13318_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_yycbbu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=  13.394883
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycddu_13326_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yycddu_13326_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_yycddu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=   	2.0054737
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycssu_13323_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yycssu_13323_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_yycssu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 2.0248353
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyubbc_13320_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yyubbc_13320_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_yyubbc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=13.330064
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyuddc_13328_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yyuddc_13328_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_yyuddc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=2.0034170
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyussc_13325_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yyussc_13325_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/test_ee_yyussc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=2.0189010
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_


    lumi_=1000.

    use_partonInfo_=True
    #use_partonInfo_=False
    cross_section_= 2.67
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_hzqq_13392_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/hzqq_13392_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_hzqq_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root"  
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_13393_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_qqqq_13393_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qqqq_noCuts_noSqrtS.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qqqq_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root"
    cross_section_= 120.
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_mqqqq_2_TeV_13700_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7ee_qq_mqq_1_TeV_13429_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qqqq_mqqqq_2TeV_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root"
    cross_section_=   49.2
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_13398_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_13398_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root"
    cross_section_= 786.
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_mqq_1_TeV_13429_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7ee_qq_mqq_1_TeV_13429_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root"
    cross_section_=  73.5 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_



    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbcbbc_13071_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_bbcbbc_13071_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_bbcbbc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 2.9986901e-3 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbubbu_13072_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_bbubbu_13072_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_bbubbu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 2.9825397e-3 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ddcyyc_13073_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_ddcyyc_13073_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_ddcyyc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=1.7824610e-1
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_dduyyu_13074_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_dduyyu_13074_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_dduyyu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 5.0109474
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscbbc_13075_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_sscbbc_13075_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_sscbbc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=4.8938333
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscssc_13076_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_sscssc_13076_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_sscssc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 1.3677677e-1 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssussu_13077_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_ssussu_13077_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_ssussu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 3.3776171e-3 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssubbu_13293_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_ssubbu_13293_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_ssubbu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 2.3216638e-2 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycbbu_13322_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yycbbu_13322_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_yycbbu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 5.2101109
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycddu_13319_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yycddu_13319_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_yycddu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 4.0984879e-1 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycssu_13327_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yycssu_13327_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_yycssu_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=  4.1853929e-1 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyubbc_13324_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yyubbc_13324_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_yyubbc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 5.2070149
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyuddc_13321_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yyuddc_13321_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_yyuddc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_= 4.1203686e-1 
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyussc_13329_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190417_gcc62_CT/VtxVLC7RFJVLC7/ee_yyussc_13329_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/test_ee_yyussc_ellipse_m1_126_dm_40_m2_92_5_dm_40_noThetaCut_MVATrainingTree.root" 
    cross_section_=4.2245034e-1
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,performMassCuts_,performThetaCuts_)
    print 'finished file', final_histo_name_





    return None

process_files()



