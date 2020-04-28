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

def fill_HZ_histograms(file,xsec, usePartonInfo,hist_vec_reco_1D,hist_sqrtS_1D,hist_vec_part_vs_gen_or_reco_2D,lumi,performMassCuts,useMassRectangleCuts,performThetaCuts,performBTagCuts,performC2Cuts,performD2Cuts):
#, h_hist_parton, h_hist_vec_gen, h_hist_vec_reco_1D, usePartonInfo, x_sec, fill_partonInfo, fill_genInfo):
#    print "size of histogram vectors ",len(h_hist_parton),len(h_hist_vec_gen),len(h_hist_vec_reco_1D), " booleans ",usePartonInfo,fill_partonInfo,fill_genInfo, " xsec ",xsec 
    
    use_EMissNeutrinoProjection=True
    #here use total 4 vector, or 4 vector sum of jet 1 and 2 (see flag below)- isolated photon four vector plus correction with EMiss on both jet axes
    #mass cuts are then also done after projecting the EMiss
    use_sqrtJets=True #in this case use j1+j2-isolated photons and with upper flag still decide if EMiss projection on jets is performed
    use_MHMiss_over_PFOMiss = False #if recoil jet missing energy or PFO missing energy is used in the missing energy projection

    fCut_sj_BTagMax_min=0.9



    fCut_mass_1_min=100.
    fCut_mass_1_max=160.
    fCut_mass_2_min=50.
    fCut_mass_2_max=135.

    #if rectangle or ellipse should be used
    #fCut_use_rectangle=True
    fCut_use_rectangle=useMassRectangleCuts
    #fCut_mass_1_min=100.
    #fCut_mass1_center=126.
    #fCut_mass1_radius=20.
    #fCut_mass2_center=92.5
    #fCut_mass2_radius=20.
    fCut_mass_1_min=0.
    fCut_mass1_center=126.
    fCut_mass2_center=92.5
    fCut_mass1_radius=35.
    fCut_mass2_radius=35.

    fCut_C2_j1_center=0.075
    fCut_C2_j1_radius=0.075
    fCut_C2_j2_center=0.055
    fCut_C2_j2_radius=0.070

    fCut_D2_j1_center=1.75
    fCut_D2_j1_radius=1.25
    fCut_D2_j2_center=1.75
    fCut_D2_j2_radius=1.25

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
        recojet_subjet_rfj_j_E=ientry.recojet_subjet_rfj_j_E
        if len(recojet_subjet_rfj_j_E)<4 :
            #print "too few subjets ",len(recojet_subjet_rfj_j_E)
            continue
        recojet_subjet_rfj_j_Px=ientry.recojet_subjet_rfj_j_Px
        recojet_subjet_rfj_j_Py=ientry.recojet_subjet_rfj_j_Py
        recojet_subjet_rfj_j_Pz=ientry.recojet_subjet_rfj_j_Pz
        #print "in ientry ", ientry

        if num_entry%(int(tree.GetEntries()/5.)) == 0:
            print "in entry ",num_entry
            
        gen_pass_mass_cuts=False
        reco_pass_mass_cuts=False

        tempTotEventP4=TLorentzVector(0,0,0,0)
        tempTotEventP4HZ=TLorentzVector(0,0,0,0)
        tempZP4=TLorentzVector(0,0,0,0)
        tempHP4=TLorentzVector(0,0,0,0)
        
        quark_Z_decays=0
        
        tempH_b=TLorentzVector(0,0,0,0);
        tempH_bbar=TLorentzVector(0,0,0,0);
        tempH_p1=TLorentzVector(0,0,0,0);
        tempH_p2=TLorentzVector(0,0,0,0);
        tempZ_q_pos=TLorentzVector(0,0,0,0)
        tempZ_q_neg=TLorentzVector(0,0,0,0)

        angle_min_sj_rj_Z_q_pos=pi
        ind_Z_q_pos_sj_rj_angle=-1
        angle_min_sj_rj_Z_q_neg=pi
        ind_Z_q_neg_sj_rj_angle=-1

        angle_min_sj_rj_H_b=pi
        ind_H_b_sj_rj_angle=-1
        angle_min_sj_rj_H_bbar_neg=pi
        ind_H_bbar_neg_sj_rj_angle=-1


        #now fill the subjets matched to b and bbar of the H matched large jet
        #use here MC truth matching to find out good parameters
        temp_sj_b_rj_H=TLorentzVector(0,0,0,0)
        temp_sj_bbar_rj_H=TLorentzVector(0,0,0,0)

        temp_sj_q_pos_rj_Z=TLorentzVector(0,0,0,0)
        temp_sj_q_neg_rj_Z=TLorentzVector(0,0,0,0)

        H_decays_bbar = True

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

        tempTotGenP4=TLorentzVector(0,0,0,0)
        tempTotGenP4.SetPxPyPzE(ientry.true_Px,ientry.true_Py,ientry.true_Pz,ientry.true_E);
        n_IsoPh_gen=0
        n_IsoLep_gen=0
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
                        
        sqrtS_eff_gen=(tempTotGenP4-tempGenIsoPhP4).M();

        tempGenTrueEMissP4=TLorentzVector(0,0,0,0);
        tempGenTrueEMissP4.SetPxPyPzE(ientry.true_inv_Px,ientry.true_inv_Py,ientry.true_inv_Pz,ientry.true_inv_E)

        tempGenEMissP4=TLorentzVector(0,0,0,0);
        tempGenEMissP4.SetPxPyPzE(-ientry.true_Px,-ientry.true_Py,-ientry.true_Pz,sqrt(pow(ientry.true_Px,2)+pow(ientry.true_Py,2)+pow(ientry.true_Pz,2)))

        #print "entries in px/py/pz/E ",ientry.true_inv_Px,-ientry.true_Px,ientry.true_inv_Py,-ientry.true_Py,ientry.true_inv_Pz,-ientry.true_Pz


        tempGenMETP4=TLorentzVector(0,0,0,0)
        tempGenMETP4.SetPxPyPzE(tempGenEMissP4.Px(),tempGenEMissP4.Py(),0,tempGenEMissP4.Pt())

        tempGenTrueMETP4=TLorentzVector(0,0,0,0)
        tempGenTrueMETP4.SetPxPyPzE(tempGenTrueEMissP4.Px(),tempGenTrueEMissP4.Py(),0,tempGenTrueEMissP4.Pt())

        #print "difference between gm and gmt ",tempGenMETP4.Px()-tempGenTrueMETP4.Px(),tempGenMETP4.Py()-tempGenTrueMETP4.Py(),tempGenMETP4.Pz()-tempGenTrueMETP4.Pz(),tempGenMETP4.E()-tempGenTrueMETP4.E()

        tempGenTrueEMissCorrP4=TLorentzVector(0,0,0,0)
        tempGenEMissCorrP4=TLorentzVector(0,0,0,0)



        genjet_E=ientry.genjet_E
        if len(genjet_E)==2 : 
            genjet_Px=ientry.genjet_Px
            genjet_Py=ientry.genjet_Py
            genjet_Pz=ientry.genjet_Pz
            
            gj_m1_orig=TLorentzVector(0,0,0,0);
            gj_m2_orig=TLorentzVector(0,0,0,0);
            gj_m1_orig.SetPxPyPzE(genjet_Px[0],genjet_Py[0],genjet_Pz[0],genjet_E[0]);
            gj_m2_orig.SetPxPyPzE(genjet_Px[1],genjet_Py[1],genjet_Pz[1],genjet_E[1]);
            
            gj1_EMissProjVecProp=TLorentzVector(0,0,0,0)
            gj2_EMissProjVecProp=TLorentzVector(0,0,0,0)
            gj1_TrueEMissProjVecProp=TLorentzVector(0,0,0,0)
            gj2_TrueEMissProjVecProp=TLorentzVector(0,0,0,0)
            
            if(tempGenMETP4.Pt()>0):
                gj1_METProjProp=(tempGenMETP4.Vect().Dot(gj_m1_orig.Vect().Unit()))*gj_m1_orig.Vect().Unit();
                #check if gj1 and MET in same hemisphere
                if(tempGenMETP4.Vect().Unit().Dot(gj_m1_orig.Vect().Unit())>0):
                    gj1_EMissProjVecProp.SetPxPyPzE(gj1_METProjProp.Px(),gj1_METProjProp.Py(),gj1_METProjProp.Pt()*gj_m1_orig.Pz()/gj_m1_orig.Pt(),gj1_METProjProp.Pt()*gj_m1_orig.P()/gj_m1_orig.Pt());
                gj2_METProjProp=(tempGenMETP4.Vect().Dot(gj_m2_orig.Vect().Unit()))*gj_m2_orig.Vect().Unit();
                #check if gj1 and MET in same hemisphere
                if(tempGenMETP4.Vect().Unit().Dot(gj_m2_orig.Vect().Unit())>0):
                   gj2_EMissProjVecProp.SetPxPyPzE(gj2_METProjProp.Px(),gj2_METProjProp.Py(),gj2_METProjProp.Pt()*gj_m2_orig.Pz()/gj_m2_orig.Pt(),gj2_METProjProp.Pt()*gj_m2_orig.P()/gj_m2_orig.Pt());
                        
                        
            gj1_EMiss=TLorentzVector(0,0,0,0)
            gj1_EMiss=gj_m1_orig+gj1_EMissProjVecProp
            gj2_EMiss=TLorentzVector(0,0,0,0)
            gj2_EMiss=gj_m2_orig+gj2_EMissProjVecProp
            tempGenEMissCorrP4=gj1_EMiss+gj2_EMiss-gj_m1_orig-gj_m2_orig
                        
            if(tempGenTrueMETP4.Pt()>0):
                gj1_TrueMETProjProp=(tempGenTrueMETP4.Vect().Dot(gj_m1_orig.Vect().Unit()))*gj_m1_orig.Vect().Unit();
               #check if gj1 and TrueMET in same hemisphere
                if(tempGenTrueMETP4.Vect().Unit().Dot(gj_m1_orig.Vect().Unit())>0):
                    gj1_TrueEMissProjVecProp.SetPxPyPzE(gj1_TrueMETProjProp.Px(),gj1_TrueMETProjProp.Py(),gj1_TrueMETProjProp.Pt()*gj_m1_orig.Pz()/gj_m1_orig.Pt(),gj1_TrueMETProjProp.Pt()*gj_m1_orig.P()/gj_m1_orig.Pt());
                gj2_TrueMETProjProp=(tempGenTrueMETP4.Vect().Dot(gj_m2_orig.Vect().Unit()))*gj_m2_orig.Vect().Unit();
                #check if gj2 and TrueMET in same hemisphere
                if(tempGenTrueMETP4.Vect().Unit().Dot(gj_m2_orig.Vect().Unit())>0):
                    gj2_TrueEMissProjVecProp.SetPxPyPzE(gj2_TrueMETProjProp.Px(),gj2_TrueMETProjProp.Py(),gj2_TrueMETProjProp.Pt()*gj_m2_orig.Pz()/gj_m2_orig.Pt(),gj2_TrueMETProjProp.Pt()*gj_m2_orig.P()/gj_m2_orig.Pt());


            gj1_TrueEMiss=TLorentzVector(0,0,0,0)
            gj1_TrueEMiss=gj_m1_orig+gj1_TrueEMissProjVecProp
            gj2_TrueEMiss=TLorentzVector(0,0,0,0)
            gj2_TrueEMiss=gj_m2_orig+gj2_TrueEMissProjVecProp
            tempGenTrueEMissCorrP4=gj1_TrueEMiss+gj2_TrueEMiss-gj_m1_orig-gj_m2_orig

            gj_H_orig=TLorentzVector(0,0,0,0);
            gj_Z_orig=TLorentzVector(0,0,0,0);
            #fill now 2D histograms for sqrtS parton vs sqrtS gen
            if(usePartonInfo):
                hist_vec_part_vs_gen_or_reco_2D[0].Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4).M(),weight)
                hist_vec_part_vs_gen_or_reco_2D[1].Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4+(gj1_EMiss+gj2_EMiss-gj_m1_orig-gj_m2_orig)).M(),weight)
                hist_vec_part_vs_gen_or_reco_2D[2].Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4+(gj1_TrueEMiss+gj2_TrueEMiss-gj_m1_orig-gj_m2_orig)).M(),weight)
                hist_vec_part_vs_gen_or_reco_2D[3].Fill(tempTotEventP4.M(),(tempTotGenP4-tempGenIsoPhP4+tempGenTrueEMissP4).M(),weight)
                hist_vec_part_vs_gen_or_reco_2D[4].Fill(tempTotEventP4.M(),(gj_m1_orig+gj_m2_orig).M(),weight)
                hist_vec_part_vs_gen_or_reco_2D[5].Fill(tempTotEventP4.M(),(gj1_EMiss+gj2_EMiss).M(),weight)
                hist_vec_part_vs_gen_or_reco_2D[6].Fill(tempTotEventP4.M(),(gj1_TrueEMiss+gj2_TrueEMiss).M(),weight)

            hist_sqrtS_1D[1].Fill((tempTotGenP4-tempGenIsoPhP4).M(),weight)
            hist_sqrtS_1D[2].Fill((tempTotGenP4-tempGenIsoPhP4+(gj1_EMiss+gj2_EMiss-gj_m1_orig-gj_m2_orig)).M(),weight)
            hist_sqrtS_1D[3].Fill((tempTotGenP4-tempGenIsoPhP4+(gj1_TrueEMiss+gj2_TrueEMiss-gj_m1_orig-gj_m2_orig)).M(),weight)
            hist_sqrtS_1D[4].Fill((tempTotGenP4-tempGenIsoPhP4+tempGenTrueEMissP4).M(),weight)
            hist_sqrtS_1D[5].Fill((gj_m1_orig+gj_m2_orig).M(),weight)
            hist_sqrtS_1D[6].Fill((gj1_EMiss+gj2_EMiss).M(),weight)
            hist_sqrtS_1D[7].Fill((gj1_TrueEMiss+gj2_TrueEMiss).M(),weight)


        tempTotRecoP4=TLorentzVector(0,0,0,0);
        tempTotRecoP4.SetPxPyPzE(ientry.totPFO_Px,ientry.totPFO_Py,ientry.totPFO_Pz,ientry.totPFO_E)

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
            #if(abs(degrees(tempRecoMETP4.DeltaPhi(rj_m1_orig)))<90.):
                #if(tempRecoMETP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit())<0):
                #print "DeltaPhi small j1, but vector large? ",degrees(tempRecoMETP4.DeltaPhi(rj_m1_orig)),tempRecoMETP4.DeltaPhi(rj_m1_orig),tempRecoMETP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit()),degrees(acos(tempRecoMETP4.Vect().Unit().Dot(rj_m1_orig.Vect().Unit())))
                rj1_EMissProjVecProp.SetPxPyPzE(rj1_METProjProp.Px(),rj1_METProjProp.Py(),rj1_METProjProp.Pt()*rj_m1_orig.Pz()/rj_m1_orig.Pt(),rj1_METProjProp.Pt()*rj_m1_orig.P()/rj_m1_orig.Pt());
            rj2_METProjProp=(tempRecoMETP4.Vect().Dot(rj_m2_orig.Vect().Unit()))*rj_m2_orig.Vect().Unit();
            #check if rj1 and MET in same hemisphere
            if(tempRecoMETP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit())>0):
            #if(abs(degrees(tempRecoMETP4.DeltaPhi(rj_m2_orig)))<90.):
                #if(tempRecoMETP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit())<0):
                #print "DeltaPhi small j2, but vector large? ",degrees(tempRecoMETP4.DeltaPhi(rj_m1_orig)),tempRecoMETP4.DeltaPhi(rj_m1_orig),tempRecoMETP4.Vect().Unit().Dot(rj_m2_orig.Vect().Unit())
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
           
          #fill now 2D histograms for sqrtS parton vs sqrtS gen
        if(usePartonInfo):
            hist_vec_part_vs_gen_or_reco_2D[7].Fill(tempTotEventP4.M(),(tempTotRecoP4-tempRecoIsoPhP4).M(),weight)
            hist_vec_part_vs_gen_or_reco_2D[8].Fill(tempTotEventP4.M(),(tempTotRecoP4-tempRecoIsoPhP4+(rj1_EMiss+rj2_EMiss-rj_m1_orig-rj_m2_orig)).M(),weight)
            hist_vec_part_vs_gen_or_reco_2D[9].Fill(tempTotEventP4.M(),(rj_m1_orig+rj_m2_orig).M(),weight)
            hist_vec_part_vs_gen_or_reco_2D[10].Fill(tempTotEventP4.M(),(rj1_EMiss+rj2_EMiss).M(),weight)

        hist_sqrtS_1D[8].Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight)
        hist_sqrtS_1D[9].Fill((tempTotRecoP4-tempRecoIsoPhP4+(rj1_EMiss+rj2_EMiss-rj_m1_orig-rj_m2_orig)).M(),weight)
        hist_sqrtS_1D[10].Fill((rj_m1_orig+rj_m2_orig).M(),weight)
        hist_sqrtS_1D[11].Fill((rj1_EMiss+rj2_EMiss).M(),weight)
        if (tempTotRecoP4-tempRecoIsoPhP4).M()> sqrtS_high :
            hist_sqrtS_1D[12].Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight)
            if tempTotEventP4.M()> sqrtS_high : 
                 hist_sqrtS_1D[16].Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight)
        if (tempTotRecoP4-tempRecoIsoPhP4+(rj1_EMiss+rj2_EMiss-rj_m1_orig-rj_m2_orig)).M()> sqrtS_high :
            hist_sqrtS_1D[13].Fill((tempTotRecoP4-tempRecoIsoPhP4+(rj1_EMiss+rj2_EMiss-rj_m1_orig-rj_m2_orig)).M(),weight)
            if tempTotEventP4.M()> sqrtS_high : 
                 hist_sqrtS_1D[17].Fill((tempTotRecoP4-tempRecoIsoPhP4+(rj1_EMiss+rj2_EMiss-rj_m1_orig-rj_m2_orig)).M(),weight)
        if (rj_m1_orig+rj_m2_orig).M()> sqrtS_high :
            hist_sqrtS_1D[14].Fill((rj_m1_orig+rj_m2_orig).M(),weight)
            if tempTotEventP4.M()> sqrtS_high : 
                 hist_sqrtS_1D[18].Fill((rj_m1_orig+rj_m2_orig).M(),weight)
        if (rj1_EMiss+rj2_EMiss).M()> sqrtS_high :
            hist_sqrtS_1D[15].Fill((rj1_EMiss+rj2_EMiss).M(),weight)
            if tempTotEventP4.M()> sqrtS_high : 
                 hist_sqrtS_1D[19].Fill((rj1_EMiss+rj2_EMiss).M(),weight)
        if tempTotEventP4.M()> sqrtS_high :
           hist_sqrtS_1D[20].Fill((tempTotRecoP4-tempRecoIsoPhP4).M(),weight)
           hist_sqrtS_1D[21].Fill((tempTotRecoP4-tempRecoIsoPhP4+(rj1_EMiss+rj2_EMiss-rj_m1_orig-rj_m2_orig)).M(),weight)
           hist_sqrtS_1D[22].Fill((rj_m1_orig+rj_m2_orig).M(),weight)
           hist_sqrtS_1D[23].Fill((rj1_EMiss+rj2_EMiss).M(),weight)

           hist_sqrtS_1D[24].Fill(((tempTotRecoP4-tempRecoIsoPhP4).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight)
           hist_sqrtS_1D[25].Fill(((tempTotRecoP4-tempRecoIsoPhP4+(rj1_EMiss+rj2_EMiss-rj_m1_orig-rj_m2_orig)).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight)
           hist_sqrtS_1D[26].Fill(((rj_m1_orig+rj_m2_orig).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight)
           hist_sqrtS_1D[27].Fill(((rj1_EMiss+rj2_EMiss).M()-tempTotEventP4.M())/tempTotEventP4.M(),weight)


        ind_jetM1=0  
        ind_jetM2=1           
        if(rj_m1.M()<rj_m2.M()) :
            ind_jetM2=0  
            ind_jetM1=1  
            temp_M=rj_m1
            rj_m1=rj_m2;
            rj_m2=temp_M

        rj_H=TLorentzVector(0,0,0,0)
        ind_rj_H=-1;
        rj_Z=TLorentzVector(0,0,0,0)
        ind_rj_Z=-1;
        
        
        #float cut_m_j1_c=126;
        #float cut_delta_m_j1_e = 20;
        #float cut_m_j2_c=95;
        #float cut_delta_m_j2_e = 25;
        #float cut_m_j1_c_veto=86;
        #float cut_delta_m_j1_e_veto = 40;=
        #float cut_m_j2_c_veto=84;
        #float cut_delta_m_j2_e_veto_plus = 40;
        #float cut_delta_m_j2_e_veto_minus = 40;

        #bool veto_qqqq_mass_ellipse_reco = false;
        #if( (rj_m2.M()<cut_m_j2_c_veto && ((pow((rj_m2.M()-cut_m_j2_c_veto)/cut_delta_m_j2_e_veto_minus,2)+pow((rj_m1.M()-cut_m_j1_c_veto)/cut_delta_m_j1_e_veto,2))<1.)) || (rj_m2.M()>cut_m_j2_c_veto && ((pow((rj_m2.M()-cut_m_j2_c_veto)/cut_delta_m_j2_e_veto_plus,2)+pow((rj_m1.M()-cut_m_j1_c_veto)/cut_delta_m_j1_e_veto,2))<1.))){
	#veto_qqqq_mass_ellipse_reco =true;
	#//std::cout<<"veto reco is hit"<<std::endl;
        #}

        #if( (pow((rj_m2.M()-fCut_mass2_center)/fCut_mass2_radius,2)+pow((rj_m1.M()-fCut_mass2_center)/fCut_mass1_radius,2))>1.)
	#  reco_pass_mass_cuts=true;
	#}
        fCut_pass_rect_mass_cut=False
        if (rj_m1.M()>fCut_mass_1_min and rj_m1.M()<fCut_mass_1_max) and (rj_m2.M()>fCut_mass_2_min and rj_m2.M()<fCut_mass_2_max):
            fCut_pass_rect_mass_cut=True 
            
        fCut_pass_ellipse_mass_cut=False
        if (rj_m1.M()>fCut_mass_1_min and  (pow((rj_m2.M()-fCut_mass2_center)/fCut_mass2_radius,2)+pow((rj_m1.M()-fCut_mass1_center)/fCut_mass1_radius,2))<1.) :
            fCut_pass_ellipse_mass_cut=True 




        if(len(recojet_E)==2 and n_IsoLep_reco<m_cut_nLeptons) and  sqrtS_eff_reco>sqrtS_high :
            num_count+=weight


        if performMassCuts and ((fCut_use_rectangle and not fCut_pass_rect_mass_cut) or (not fCut_use_rectangle and not fCut_pass_ellipse_mass_cut)):
            continue

        if performThetaCuts and ((abs(degrees(rj_m1.Theta())-fCut_thetaRef))>fCut_thetaWindow or abs(degrees(rj_m1.Theta()-rj_m2.Theta()))>fCut_delta_theta) :
            continue


        if(len(recojet_E)==2 and n_IsoLep_reco<m_cut_nLeptons):
            #tie now H and Z to one of the jets
            if rj_m1.Angle(tempHP4.Vect())<rj_m1.Angle(tempZP4.Vect()) :
                ind_rj_H=ind_jetM1
                rj_H=rj_m1
                ind_rj_Z=ind_jetM2
                rj_Z=rj_m2
                if rj_m2.Angle(tempHP4.Vect())<rj_m2.Angle(tempZP4.Vect()) :
                    #print "seems also m2 is closer to H than to Z ",rj_m1.Angle(tempHP4.Vect()),rj_m1.Angle(tempZP4.Vect()),rj_m2.Angle(tempHP4.Vect()),rj_m2.Angle(tempZP4.Vect()),sqrtS_eff_reco,sqrtS_eff_gen,tempTotEventP4.M()
                    if rj_m2.Angle(tempHP4.Vect())<rj_m1.Angle(tempHP4.Vect()) :  
                        ind_rj_H=ind_jetM2
                        rj_H=rj_m2
                        ind_rj_Z=ind_jetM1
                        rj_Z=rj_m1
            else:
                ind_rj_Z=ind_jetM1
                rj_Z=rj_m1
                ind_rj_H=ind_jetM2
                rj_H=rj_m2
                if rj_m2.Angle(tempHP4.Vect())>rj_m2.Angle(tempZP4.Vect()) :
                    #print "seems also m2 is closer to Z than to H ",rj_m1.Angle(tempHP4.Vect()),rj_m1.Angle(tempZP4.Vect()),rj_m2.Angle(tempHP4.Vect()),rj_m2.Angle(tempZP4.Vect()),sqrtS_eff_reco,sqrtS_eff_gen,tempTotEventP4.M()
                    if rj_m2.Angle(tempZP4.Vect())<rj_m1.Angle(tempZP4.Vect()) :  
                        ind_rj_Z=ind_jetM2
                        rj_Z=rj_m2
                        ind_rj_H=ind_jetM1
                        rj_H=rj_m1
            
            if (ind_rj_H==-1 or ind_rj_Z==-1) or (ind_rj_H==ind_rj_Z) :
                print "problem with jet indices association ",ind_rj_H,ind_rj_Z
            

            #performC2Cuts

            #fCut_C2_j1_center=0.075
            #fCut_C2_j1_radius=0.075
            #fCut_C2_j2_center=0.055
            #fCut_C2_j2_radius=0.070

            #fCut_D2_j1_center=1.25
            #fCut_D2_j1_radius=1.25
            #fCut_D2_j2_center=1.25
            #fCut_D2_j2_radius=1.25


            fCut_pass_upper_ellipse_C2_cut=False
            if((ientry.recojet_beta1_C2[ind_jetM1]<fCut_C2_j1_center and ientry.recojet_beta1_C2[ind_jetM2]<fCut_C2_j2_center) or  (pow((ientry.recojet_beta1_C2[ind_jetM2]-fCut_C2_j2_center)/fCut_C2_j2_radius,2)+pow((ientry.recojet_beta1_C2[ind_jetM1]-fCut_C2_j1_center)/fCut_C2_j1_radius,2))<1.) : 
                 fCut_pass_upper_ellipse_C2_cut=True


            if performC2Cuts and not fCut_pass_upper_ellipse_C2_cut :
                 continue


            fCut_pass_upper_ellipse_D2_cut=False
            if((ientry.recojet_beta1_D2[ind_jetM1]<fCut_D2_j1_center and ientry.recojet_beta1_D2[ind_jetM2]<fCut_D2_j2_center) or  (pow((ientry.recojet_beta1_D2[ind_jetM2]-fCut_D2_j2_center)/fCut_D2_j2_radius,2)+pow((ientry.recojet_beta1_D2[ind_jetM1]-fCut_D2_j1_center)/fCut_D2_j1_radius,2))<1.) : 
                 fCut_pass_upper_ellipse_D2_cut=True


            if performD2Cuts and not fCut_pass_upper_ellipse_D2_cut :
                 continue



            

            if sqrtS_eff_reco>sqrtS_high :

                recojet_subjet_rfj_j_E=ientry.recojet_subjet_rfj_j_E
                if len(recojet_subjet_rfj_j_E)!=4 : 
                    print "not four subjets in event, skip"
                    continue;
                recojet_subjet_rfj_j_Px                    =ientry.recojet_subjet_rfj_j_Px
                recojet_subjet_rfj_j_Py                    =ientry.recojet_subjet_rfj_j_Py
                recojet_subjet_rfj_j_Pz                    =ientry.recojet_subjet_rfj_j_Pz
                #recojet_subjet_rfj_j_jetChargePt_kappa_0_30=ientry.recojet_subjet_rfj_j_jetChargePt_kappa_0_30
                recojet_subjet_rfj_j_jetindex              =ientry.recojet_subjet_rfj_j_jetindex
                recojet_subjet_rfj_j_subjetindex           =ientry.recojet_subjet_rfj_j_subjetindex
                recojet_subjet_rfj_j_BTag                  =ientry.recojet_subjet_rfj_j_BTag
                recojet_subjet_rfj_j_CTag                  =ientry.recojet_subjet_rfj_j_CTag
                recojet_subjet_rfj_j_OTag                  =ientry.recojet_subjet_rfj_j_OTag
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

                if BTagMax==0:
                    continue;


                if performBTagCuts and BTagMax<fCut_sj_BTagMax_min :
                    continue;

                #histograms [0] mass of matched H jet, [1] mass of matched Z jet, [2] b-tag of matched b-rfjet of H, [3] b-tag of matched bbar rfjet of H
                hist_vec_reco_1D[0].Fill(rj_H.M(),weight);
                hist_vec_reco_1D[1].Fill(rj_Z.M(),weight);

                recojet_subjet_rfj_j_E=ientry.recojet_subjet_rfj_j_E
                if len(recojet_subjet_rfj_j_E)!=4 : 
                    print "not four subjets in event, skip"
                    continue;
                recojet_subjet_rfj_j_Px                    =ientry.recojet_subjet_rfj_j_Px
                recojet_subjet_rfj_j_Py                    =ientry.recojet_subjet_rfj_j_Py
                recojet_subjet_rfj_j_Pz                    =ientry.recojet_subjet_rfj_j_Pz
                #recojet_subjet_rfj_j_jetChargePt_kappa_0_30=ientry.recojet_subjet_rfj_j_jetChargePt_kappa_0_30
                recojet_subjet_rfj_j_jetindex              =ientry.recojet_subjet_rfj_j_jetindex
                recojet_subjet_rfj_j_subjetindex           =ientry.recojet_subjet_rfj_j_subjetindex
                recojet_subjet_rfj_j_BTag                  =ientry.recojet_subjet_rfj_j_BTag
                recojet_subjet_rfj_j_CTag                  =ientry.recojet_subjet_rfj_j_CTag
                recojet_subjet_rfj_j_OTag                  =ientry.recojet_subjet_rfj_j_OTag

                #for matching between the original jet and refined subjet
                recojet_subjet_E                     =ientry.recojet_subjet_E
                recojet_subjet_Px                    =ientry.recojet_subjet_Px
                recojet_subjet_Py                    =ientry.recojet_subjet_Py
                recojet_subjet_Pz                    =ientry.recojet_subjet_Pz
                #recojet_subjet_jetChargeE_kappa_0_15 =ientry.recojet_subjet_jetChargeE_kappa_0_15
                recojet_subjet_jetChargeE_kappa_0_20 =ientry.recojet_subjet_jetChargeE_kappa_0_20
                recojet_subjet_jetChargeE_kappa_0_25 =ientry.recojet_subjet_jetChargeE_kappa_0_25
                recojet_subjet_jetChargeE_kappa_0_30 =ientry.recojet_subjet_jetChargeE_kappa_0_30
                recojet_subjet_jetChargeE_kappa_0_50 =ientry.recojet_subjet_jetChargeE_kappa_0_50
                recojet_subjet_jetindex               =ientry.recojet_subjet_jetindex

                if len(recojet_subjet_rfj_j_E)!=len(recojet_subjet_E):
                    print "subjet and rf subjet length not the same ",len(recojet_subjet_rfj_j_E),len(recojet_subjet_E)
                    continue;
                #now loop over subjets to find out which one is closer to the b-subjet
                #sanity plots --> check that combined rfjet separated 4 jets combine into a jet close to the fat jet again 
                #particularly the original VLC10 2 fat jet, separated into 4 jets with VLC7 by LCFIPlus
                rj_m1_rfj=TLorentzVector(0,0,0,0);
                rj_m2_rfj=TLorentzVector(0,0,0,0);

                ind_sj_to_H=0
                ind_sj_to_Z=0
                ind_sj_qpos_Z=-1
                ind_sj_qneg_Z=-1
                ind_sj_bbar_H=-1
                ind_sj_b_H=-1

                for i in range(len(recojet_subjet_rfj_j_E)):
                    #print "sqrtS angles ",sqrtS_eff_reco,H_decays_bbar
                    #print "loop over subjet ",i,ind_rj_H,ind_rj_Z
                    temp=TLorentzVector(0,0,0,0)
                    temp.SetPxPyPzE(recojet_subjet_rfj_j_Px[i],recojet_subjet_rfj_j_Py[i],recojet_subjet_rfj_j_Pz[i],recojet_subjet_rfj_j_E[i])
                    if (recojet_subjet_jetindex[i]==0):
                        rj_m1_rfj+=temp;
                    else:
                        rj_m2_rfj+=temp; 
                    #check if subjets belong to higgs 
                    if recojet_subjet_rfj_j_jetindex[i]==ind_rj_H:
                        ind_sj_to_H+=1
                        if temp.Angle(tempH_b.Vect())<temp.Angle(tempH_bbar.Vect()):
                            #check if index already assigned to previously, then now redecide
                            if(ind_sj_b_H!=-1):
                                temp_sj2=TLorentzVector(0,0,0,0)
                                temp_sj2.SetPxPyPzE(recojet_subjet_rfj_j_Px[ind_sj_b_H],recojet_subjet_rfj_j_Py[ind_sj_b_H],recojet_subjet_rfj_j_Pz[ind_sj_b_H],recojet_subjet_rfj_j_E[ind_sj_b_H])
                                if temp.Angle(tempH_b.Vect())<temp_sj2.Angle(tempH_b.Vect()):
                                    ind_sj_bbar_H=ind_sj_b_H
                                    ind_sj_b_H=i
                                else:
                                    ind_sj_bbar_H=i
                            else:
                                ind_sj_b_H=i
                        else:
                           #check if H bbar index already assigned to previously, then now redecide
                            if(ind_sj_bbar_H!=-1):
                                temp_sj2=TLorentzVector(0,0,0,0)
                                temp_sj2.SetPxPyPzE(recojet_subjet_rfj_j_Px[ind_sj_bbar_H],recojet_subjet_rfj_j_Py[ind_sj_bbar_H],recojet_subjet_rfj_j_Pz[ind_sj_bbar_H],recojet_subjet_rfj_j_E[ind_sj_bbar_H])
                                if temp.Angle(tempH_bbar.Vect())<temp_sj2.Angle(tempH_bbar.Vect()):
                                    ind_sj_b_H=ind_sj_bbar_H
                                    ind_sj_bbar_H=i
                                else:
                                    ind_sj_b_H=i
                            else:
                                ind_sj_bbar_H=i
                    else:
                        ind_sj_to_Z+=1
                        if temp.Angle(tempZ_q_neg.Vect())<temp.Angle(tempZ_q_pos.Vect()):
                            #check if index already assigned to previously, then now redecide
                            if(ind_sj_qneg_Z!=-1):
                                temp_sj2=TLorentzVector(0,0,0,0)
                                temp_sj2.SetPxPyPzE(recojet_subjet_rfj_j_Px[ind_sj_qneg_Z],recojet_subjet_rfj_j_Py[ind_sj_qneg_Z],recojet_subjet_rfj_j_Pz[ind_sj_qneg_Z],recojet_subjet_rfj_j_E[ind_sj_qneg_Z])
                                if temp.Angle(tempZ_q_neg.Vect())<temp_sj2.Angle(tempZ_q_neg.Vect()):
                                    ind_sj_qpos_Z=ind_sj_qneg_Z
                                    ind_sj_qneg_Z=i
                                else:
                                    ind_sj_qpos_Z=i
                            else:
                                ind_sj_qneg_Z=i
                        else:
                           #check if H bbar index already assigned to previously, then now redecide
                            if(ind_sj_qpos_Z!=-1):
                                temp_sj2=TLorentzVector(0,0,0,0)
                                temp_sj2.SetPxPyPzE(recojet_subjet_rfj_j_Px[ind_sj_qpos_Z],recojet_subjet_rfj_j_Py[ind_sj_qpos_Z],recojet_subjet_rfj_j_Pz[ind_sj_qpos_Z],recojet_subjet_rfj_j_E[ind_sj_qpos_Z])
                                if temp.Angle(tempZ_q_pos.Vect())<temp_sj2.Angle(tempZ_q_pos.Vect()):
                                    ind_sj_qneg_Z=ind_sj_qpos_Z
                                    ind_sj_qpos_Z=i
                                else:
                                    ind_sj_qneg_Z=i
                            else:
                                ind_sj_qpos_Z=i

                if ind_sj_to_Z!=2 or ind_sj_to_H!=2:
                    print "HUGE HUGE HUGE mistake in subjet assignment ",ind_sj_to_Z,ind_sj_to_H,(len(recojet_subjet_rfj_j_E))

                if ind_sj_qpos_Z==-1 or ind_sj_qneg_Z==-1:
                   print "Z subjet index assignment off ",ind_sj_qpos_Z,ind_sj_qneg_Z

                if ind_sj_b_H==-1 or ind_sj_bbar_H==-1:
                    print "H subjet index assignment off ",ind_sj_b_H,ind_sj_bbar_H
                    
                if recojet_subjet_rfj_j_BTag[ind_sj_b_H]>recojet_subjet_rfj_j_BTag[ind_sj_bbar_H] :
                    hist_vec_reco_1D[28].Fill(recojet_subjet_rfj_j_BTag[ind_sj_b_H],weight);
                    hist_vec_reco_1D[29].Fill(recojet_subjet_rfj_j_BTag[ind_sj_bbar_H],weight);
                else: 
                    hist_vec_reco_1D[28].Fill(recojet_subjet_rfj_j_BTag[ind_sj_bbar_H],weight);
                    hist_vec_reco_1D[29].Fill(recojet_subjet_rfj_j_BTag[ind_sj_b_H],weight);
                hist_vec_reco_1D[30].Fill(recojet_subjet_rfj_j_BTag[ind_sj_b_H]+recojet_subjet_rfj_j_BTag[ind_sj_bbar_H],weight);
                hist_vec_reco_1D[31].Fill(recojet_subjet_rfj_j_CTag[ind_sj_b_H]+recojet_subjet_rfj_j_CTag[ind_sj_bbar_H],weight);
                hist_vec_reco_1D[32].Fill(recojet_subjet_rfj_j_OTag[ind_sj_b_H]+recojet_subjet_rfj_j_OTag[ind_sj_bbar_H],weight);
                
                if(quark_Z_decays==5):
                   if recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z]>recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z] :
                      hist_vec_reco_1D[33].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z],weight);
                      hist_vec_reco_1D[34].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z],weight);
                   else:
                      hist_vec_reco_1D[33].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z],weight);
                      hist_vec_reco_1D[34].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z],weight);
                   hist_vec_reco_1D[35].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z]+recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z],weight);
                   hist_vec_reco_1D[36].Fill(recojet_subjet_rfj_j_CTag[ind_sj_qpos_Z]+recojet_subjet_rfj_j_CTag[ind_sj_qneg_Z],weight);
                   hist_vec_reco_1D[37].Fill(recojet_subjet_rfj_j_OTag[ind_sj_qpos_Z]+recojet_subjet_rfj_j_OTag[ind_sj_qneg_Z],weight);
                elif quark_Z_decays==4 :
                   if recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z]>recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z] :
                      hist_vec_reco_1D[38].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z],weight);
                      hist_vec_reco_1D[39].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z],weight);
                   else:
                      hist_vec_reco_1D[38].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z],weight);
                      hist_vec_reco_1D[39].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z],weight);
                   hist_vec_reco_1D[40].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z]+recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z],weight);
                   hist_vec_reco_1D[41].Fill(recojet_subjet_rfj_j_CTag[ind_sj_qpos_Z]+recojet_subjet_rfj_j_CTag[ind_sj_qneg_Z],weight);
                   hist_vec_reco_1D[42].Fill(recojet_subjet_rfj_j_OTag[ind_sj_qpos_Z]+recojet_subjet_rfj_j_OTag[ind_sj_qneg_Z],weight);
                else:
                   if recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z]>recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z] :
                      hist_vec_reco_1D[43].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z],weight);
                      hist_vec_reco_1D[44].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z],weight);
                   else:
                      hist_vec_reco_1D[43].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z],weight);
                      hist_vec_reco_1D[44].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z],weight);
                   hist_vec_reco_1D[45].Fill(recojet_subjet_rfj_j_BTag[ind_sj_qpos_Z]+recojet_subjet_rfj_j_BTag[ind_sj_qneg_Z],weight);
                   hist_vec_reco_1D[46].Fill(recojet_subjet_rfj_j_CTag[ind_sj_qpos_Z]+recojet_subjet_rfj_j_CTag[ind_sj_qneg_Z],weight);
                   hist_vec_reco_1D[47].Fill(recojet_subjet_rfj_j_OTag[ind_sj_qpos_Z]+recojet_subjet_rfj_j_OTag[ind_sj_qneg_Z],weight);

                ind_sj_o_b_H=-1
                ind_sj_o_bbar_H=-1
                ind_sj_o_qpos_Z=-1
                ind_sj_o_qneg_Z=-1
                index_sj_o=0
                ind_sj_o_to_H=0
                ind_sj_o_to_Z=0
                #now loop over subjets to find out which one is closer to the b-subjet
                for i in range(len(recojet_subjet_E)):
                    #print "sqrtS angles ",sqrtS_eff_reco,H_decays_bbar
                    #print "loop over subjet ",i,ind_rj_H,ind_rj_Z
                    temp=TLorentzVector(0,0,0,0)
                    temp.SetPxPyPzE(recojet_subjet_Px[i],recojet_subjet_Py[i],recojet_subjet_Pz[i],recojet_subjet_E[i])
                    #check if subjets belong to higgs 
                    if recojet_subjet_jetindex[i]==ind_rj_H:
                        ind_sj_o_to_H+=1
                        if temp.Angle(tempH_b.Vect())<temp.Angle(tempH_bbar.Vect()):
                            #check if index already assigned to previously, then now redecide
                            if(ind_sj_o_b_H!=-1):
                                temp_sj_o2=TLorentzVector(0,0,0,0)
                                temp_sj_o2.SetPxPyPzE(recojet_subjet_Px[ind_sj_o_b_H],recojet_subjet_Py[ind_sj_o_b_H],recojet_subjet_Pz[ind_sj_o_b_H],recojet_subjet_E[ind_sj_o_b_H])
                                if temp.Angle(tempH_b.Vect())<temp_sj_o2.Angle(tempH_b.Vect()):
                                    ind_sj_o_bbar_H=ind_sj_o_b_H
                                    ind_sj_o_b_H=i
                                else:
                                    ind_sj_o_bbar_H=i
                            else:
                                ind_sj_o_b_H=i
                        else:
                           #check if H bbar index already assigned to previously, then now redecide
                            if(ind_sj_o_bbar_H!=-1):
                                temp_sj_o2=TLorentzVector(0,0,0,0)
                                temp_sj_o2.SetPxPyPzE(recojet_subjet_Px[ind_sj_o_bbar_H],recojet_subjet_Py[ind_sj_o_bbar_H],recojet_subjet_Pz[ind_sj_o_bbar_H],recojet_subjet_E[ind_sj_o_bbar_H])
                                if temp.Angle(tempH_bbar.Vect())<temp_sj_o2.Angle(tempH_bbar.Vect()):
                                    ind_sj_o_b_H=ind_sj_o_bbar_H
                                    ind_sj_o_bbar_H=i
                                else:
                                    ind_sj_o_b_H=i
                            else:
                                ind_sj_o_bbar_H=i
                    else:
                        ind_sj_o_to_Z+=1
                        if temp.Angle(tempZ_q_neg.Vect())<temp.Angle(tempZ_q_pos.Vect()):
                            #check if index already assigned to previously, then now redecide
                            if(ind_sj_o_qneg_Z!=-1):
                                temp_sj_o2=TLorentzVector(0,0,0,0)
                                temp_sj_o2.SetPxPyPzE(recojet_subjet_Px[ind_sj_o_qneg_Z],recojet_subjet_Py[ind_sj_o_qneg_Z],recojet_subjet_Pz[ind_sj_o_qneg_Z],recojet_subjet_E[ind_sj_o_qneg_Z])
                                if temp.Angle(tempZ_q_neg.Vect())<temp_sj_o2.Angle(tempZ_q_neg.Vect()):
                                    ind_sj_o_qpos_Z=ind_sj_o_qneg_Z
                                    ind_sj_o_qneg_Z=i
                                else:
                                    ind_sj_o_qpos_Z=i
                            else:
                                ind_sj_o_qneg_Z=i
                        else:
                           #check if H bbar index already assigned to previously, then now redecide
                            if(ind_sj_o_qpos_Z!=-1):
                                temp_sj_o2=TLorentzVector(0,0,0,0)
                                temp_sj_o2.SetPxPyPzE(recojet_subjet_Px[ind_sj_o_qpos_Z],recojet_subjet_Py[ind_sj_o_qpos_Z],recojet_subjet_Pz[ind_sj_o_qpos_Z],recojet_subjet_E[ind_sj_o_qpos_Z])
                                if temp.Angle(tempZ_q_pos.Vect())<temp_sj_o2.Angle(tempZ_q_pos.Vect()):
                                    ind_sj_o_qneg_Z=ind_sj_o_qpos_Z
                                    ind_sj_o_qpos_Z=i
                                else:
                                    ind_sj_o_qneg_Z=i
                            else:
                                ind_sj_o_qpos_Z=i

                if ind_sj_o_to_Z!=2 or ind_sj_o_to_H!=2:
                    print "HUGE HUGE HUGE mistake in orig subjet assignment ",ind_sj_o_to_Z,ind_sj_o_to_H,(len(recojet_subjet_E))

                if ind_sj_o_qpos_Z==-1 or ind_sj_o_qneg_Z==-1:
                    print "Z orig subjet index assignment off ",ind_sj_o_qpos_Z,ind_sj_o_qneg_Z

                if ind_sj_o_b_H==-1 or ind_sj_o_bbar_H==-1:
                    print "H orig subjet index assignment off ",ind_sj_o_b_H,ind_sj_o_bbar_H
                ind_rfj_H=-1
                ind_rfj_Z=-1

                ind_rfjetM1=0
                ind_rfjetM2=1

                if rj_m1_rfj.Angle(tempHP4.Vect())<rj_m1_rfj.Angle(tempZP4.Vect()) :
                    ind_rfj_H=ind_rfjetM1
                    ind_rfj_Z=ind_rfjetM2
                    if rj_m2_rfj.Angle(tempHP4.Vect())<rj_m2_rfj.Angle(tempZP4.Vect()) :
                        if rj_m2_rfj.Angle(tempHP4.Vect())<rj_m1_rfj.Angle(tempHP4.Vect()) :  
                          ind_rfj_H=ind_rfjetM2
                          ind_rfj_Z=ind_rfjetM1   
                else:
                    ind_rfj_Z=ind_rfjetM1
                    ind_rfj_H=ind_rfjetM2
                    if rj_m2_rfj.Angle(tempHP4.Vect())>rj_m2_rfj.Angle(tempZP4.Vect()) :
                    #print "seems also m2 is closer to Z than to H ",rj_m1.Angle(tempHP4.Vect()),rj_m1.Angle(tempZP4.Vect()),rj_m2.Angle(tempHP4.Vect()),rj_m2.Angle(tempZP4.Vect()),sqrtS_eff_reco,sqrtS_eff_gen,tempTotEventP4.M()
                        if rj_m2_rfj.Angle(tempZP4.Vect())<rj_m1_rfj.Angle(tempZP4.Vect()) :  
                           ind_rfj_Z=ind_rfjetM2
                           ind_rfj_H=ind_rfjetM1
                        
                if (ind_rfj_H==-1 or ind_rfj_Z==-1) or (ind_rfj_H==ind_rfj_Z) :
                   print "problem with rf jet indices association ",ind_rfj_H,ind_rfj_Z
                rfj_H=TLorentzVector(0,0,0,0)
                rfj_Z=TLorentzVector(0,0,0,0)
                if(ind_rfj_H==ind_rfjetM1) :
                    rfj_H=rj_m1_rfj
                    rfj_Z=rj_m2_rfj
                else: 
                    rfj_H=rj_m2_rfj
                    rfj_Z=rj_m1_rfj

                rj_o_H=TLorentzVector(0,0,0,0)
                rj_o_Z=TLorentzVector(0,0,0,0)
                if(ind_rj_H==0) :
                    rj_o_H=rj_m1_orig
                    rj_o_Z=rj_m2_orig
                else: 
                    rj_o_H=rj_m2_orig
                    rj_o_Z=rj_m1_orig


                #hist_vec_reco_1D[12].Fill(degrees(rj_H.Angle(rfj_H.Vect())),weight)
                #hist_vec_reco_1D[13].Fill(degrees(rj_Z.Angle(rfj_Z.Vect())),weight)
                #hist_vec_reco_1D[14].Fill((rfj_H.E()-rj_H.E())/rj_H.E(),weight)
                #hist_vec_reco_1D[15].Fill((rfj_Z.E()-rj_Z.E())/rj_Z.E(),weight)
                hist_vec_reco_1D[12].Fill(degrees(rj_o_H.Angle(rfj_H.Vect())),weight)
                hist_vec_reco_1D[13].Fill(degrees(rj_o_Z.Angle(rfj_Z.Vect())),weight)
                hist_vec_reco_1D[14].Fill((rfj_H.E()-rj_o_H.E())/rj_o_H.E(),weight)
                hist_vec_reco_1D[15].Fill((rfj_Z.E()-rj_o_Z.E())/rj_o_Z.E(),weight)
                #hist_vec_reco_1D[48].Fill(degrees(rfj_H.Angle(rfj_Z.Vect())),weight)
                hist_vec_reco_1D[48].Fill(degrees(rj_m1_rfj.Angle(rj_m2_rfj.Vect())),weight)
                hist_vec_reco_1D[49].Fill(degrees(rj_H.Angle(rj_Z.Vect())),weight)
                hist_vec_reco_1D[50].Fill(rj_o_H.M(),weight)
                hist_vec_reco_1D[51].Fill(rj_o_Z.M(),weight)
                hist_vec_reco_1D[52].Fill(rfj_H.M(),weight)
                hist_vec_reco_1D[53].Fill(rfj_Z.M(),weight)
                hist_vec_reco_1D[54].Fill(rj_o_H.E(),weight)
                hist_vec_reco_1D[55].Fill(rj_o_Z.E(),weight)
                hist_vec_reco_1D[56].Fill(rj_H.E(),weight)
                hist_vec_reco_1D[57].Fill(rj_Z.E(),weight)
                hist_vec_reco_1D[58].Fill((rj_o_H.E()-tempHP4.E())/tempHP4.E(),weight)
                hist_vec_reco_1D[59].Fill((rj_o_Z.E()-tempZP4.E())/tempZP4.E(),weight)
                hist_vec_reco_1D[60].Fill((rj_H.E()-tempHP4.E())/tempHP4.E(),weight)
                hist_vec_reco_1D[61].Fill((rj_Z.E()-tempZP4.E())/tempZP4.E(),weight)
                hist_vec_reco_1D[62].Fill(degrees(rj_m1.Angle(tempHP4.Vect())),weight)
                hist_vec_reco_1D[63].Fill(degrees(rj_m2.Angle(tempZP4.Vect())),weight)
                #hist_vec_reco_1D[64].Fill((recojet_subjet_jetChargeE_kappa_0_15)[ind_sj_o_qpos_Z],weight);
                hist_vec_reco_1D[65].Fill((recojet_subjet_jetChargeE_kappa_0_20)[ind_sj_o_qpos_Z],weight);
                hist_vec_reco_1D[66].Fill((recojet_subjet_jetChargeE_kappa_0_25)[ind_sj_o_qpos_Z],weight);
                hist_vec_reco_1D[67].Fill((recojet_subjet_jetChargeE_kappa_0_30)[ind_sj_o_qpos_Z],weight);
                hist_vec_reco_1D[68].Fill((recojet_subjet_jetChargeE_kappa_0_50)[ind_sj_o_qpos_Z],weight);
                #hist_vec_reco_1D[69].Fill((recojet_subjet_jetChargeE_kappa_0_15)[ind_sj_o_qneg_Z],weight);
                hist_vec_reco_1D[70].Fill((recojet_subjet_jetChargeE_kappa_0_20)[ind_sj_o_qneg_Z],weight);
                hist_vec_reco_1D[71].Fill((recojet_subjet_jetChargeE_kappa_0_25)[ind_sj_o_qneg_Z],weight);
                hist_vec_reco_1D[72].Fill((recojet_subjet_jetChargeE_kappa_0_30)[ind_sj_o_qneg_Z],weight);
                hist_vec_reco_1D[73].Fill((recojet_subjet_jetChargeE_kappa_0_50)[ind_sj_o_qneg_Z],weight);



                #hist[8] angle of rfjet and original subjet for H matched b jet, [9] DR of rfjet and original sj for H matched bbar jet, [10] rel deltaE of rfjet and original subjet for H matched jet, 
                #[11] rel deltaE of rfjet and original subjet for Z matched jet, [12] vtx mass of rfjets for H matched jet, [13] vtx mass for rfjets for Z matched light decays, [14] vtx mass for rfjet of Z matched c decays
                #[15] vtx mass for rfjet of Z matched b decays [16] b-tag of Z b decay rfjets [17] b-tag of Z c decay rfjets [18] b-tag of Z light decays 
                #[19] jet charge b-tag matched H rfjet, [20] jet charge bbar-tagged matched H rfjet, [21] jet charged pos quark matched Z jet, [22] jet charge neg quark matched Z jet 
                #[23] c-tag H matched rfjets [24] O-tag of H matched rfjets [25] c-tag of Z b decay rfjets [26] c-tag of Z c decays [27] c-tag of Z light decays 
                #[28] O-tag of Z b decay rfjets [29] O-tag of Z c decays [30] O-tag of Z light decays 
                #[31] delta subjet charge of rfjet and original subjet for H matched jet, [32] delta subjet charge of rfjet and original sj for Z matched jet
                hist_vec_reco_1D[2].Fill(recojet_subjet_rfj_j_BTag[ind_sj_b_H],weight);
                hist_vec_reco_1D[3].Fill(recojet_subjet_rfj_j_BTag[ind_sj_bbar_H],weight);  
                #print "sqrtS ",sqrtS_eff_reco,H_decays_bbar


                temp_sj1=TLorentzVector(0,0,0,0)
                temp_sj2=TLorentzVector(0,0,0,0)
                #first do difference histos between refined jets and original jets for b's from H
                temp_sj1.SetPxPyPzE(recojet_subjet_rfj_j_Px[ind_sj_b_H],recojet_subjet_rfj_j_Py[ind_sj_b_H],recojet_subjet_rfj_j_Pz[ind_sj_b_H],recojet_subjet_rfj_j_E[ind_sj_b_H])
                temp_sj2.SetPxPyPzE(recojet_subjet_Px[ind_sj_o_b_H],recojet_subjet_Py[ind_sj_o_b_H],recojet_subjet_Pz[ind_sj_o_b_H],recojet_subjet_E[ind_sj_o_b_H])
                hist_vec_reco_1D[8].Fill(degrees(temp_sj1.Angle(temp_sj2.Vect())),weight);
                hist_vec_reco_1D[10].Fill((temp_sj1.E()-temp_sj2.E())/temp_sj2.E(),weight);

                #first do difference histos between refined jets and original jets for bbar's from H
                temp_sj1.SetPxPyPzE(recojet_subjet_rfj_j_Px[ind_sj_bbar_H],recojet_subjet_rfj_j_Py[ind_sj_bbar_H],recojet_subjet_rfj_j_Pz[ind_sj_bbar_H],recojet_subjet_rfj_j_E[ind_sj_bbar_H])
                temp_sj2.SetPxPyPzE(recojet_subjet_Px[ind_sj_o_bbar_H],recojet_subjet_Py[ind_sj_o_bbar_H],recojet_subjet_Pz[ind_sj_o_bbar_H],recojet_subjet_E[ind_sj_o_bbar_H])
                hist_vec_reco_1D[8].Fill(degrees(temp_sj1.Angle(temp_sj2.Vect())),weight);
                hist_vec_reco_1D[10].Fill((temp_sj1.E()-temp_sj2.E())/temp_sj2.E(),weight);

                #now do difference histos between refined jets and original jets for quarks's from Z
                temp_sj1.SetPxPyPzE(recojet_subjet_rfj_j_Px[ind_sj_qneg_Z],recojet_subjet_rfj_j_Py[ind_sj_qneg_Z],recojet_subjet_rfj_j_Pz[ind_sj_qneg_Z],recojet_subjet_rfj_j_E[ind_sj_qneg_Z])
                temp_sj2.SetPxPyPzE(recojet_subjet_Px[ind_sj_o_qneg_Z],recojet_subjet_Py[ind_sj_o_qneg_Z],recojet_subjet_Pz[ind_sj_o_qneg_Z],recojet_subjet_E[ind_sj_o_qneg_Z])
                hist_vec_reco_1D[9].Fill(degrees(temp_sj1.Angle(temp_sj2.Vect())),weight);
                hist_vec_reco_1D[11].Fill((temp_sj1.E()-temp_sj2.E())/temp_sj2.E(),weight);
                #now do difference histos between refined jets and original jets for quarks's from Z
                temp_sj1.SetPxPyPzE(recojet_subjet_rfj_j_Px[ind_sj_qpos_Z],recojet_subjet_rfj_j_Py[ind_sj_qpos_Z],recojet_subjet_rfj_j_Pz[ind_sj_qpos_Z],recojet_subjet_rfj_j_E[ind_sj_qpos_Z])
                temp_sj2.SetPxPyPzE(recojet_subjet_Px[ind_sj_o_qpos_Z],recojet_subjet_Py[ind_sj_o_qpos_Z],recojet_subjet_Pz[ind_sj_o_qpos_Z],recojet_subjet_E[ind_sj_o_qpos_Z])
                hist_vec_reco_1D[9].Fill(degrees(temp_sj1.Angle(temp_sj2.Vect())),weight);
                hist_vec_reco_1D[11].Fill((temp_sj1.E()-temp_sj2.E())/temp_sj2.E(),weight);

                temp_sj1.SetPxPyPzE(recojet_subjet_Px[ind_sj_b_H],recojet_subjet_Py[ind_sj_b_H],recojet_subjet_Pz[ind_sj_b_H],recojet_subjet_E[ind_sj_b_H])
                temp_sj2.SetPxPyPzE(recojet_subjet_Px[ind_sj_bbar_H],recojet_subjet_Py[ind_sj_bbar_H],recojet_subjet_Pz[ind_sj_bbar_H],recojet_subjet_E[ind_sj_bbar_H])
                #hist_vec_reco_1D[12].Fill(degrees(temp_sj1.Angle(temp_sj2.Vect())),weight);

                temp_sj1.SetPxPyPzE(recojet_subjet_rfj_j_Px[ind_sj_qpos_Z],recojet_subjet_rfj_j_Py[ind_sj_qpos_Z],recojet_subjet_rfj_j_Pz[ind_sj_qpos_Z],recojet_subjet_rfj_j_E[ind_sj_qpos_Z])
                temp_sj2.SetPxPyPzE(recojet_subjet_rfj_j_Px[ind_sj_qneg_Z],recojet_subjet_rfj_j_Py[ind_sj_qneg_Z],recojet_subjet_rfj_j_Pz[ind_sj_qneg_Z],recojet_subjet_rfj_j_E[ind_sj_qneg_Z])

                #hist_vec_reco_1D[13].Fill(degrees(temp_sj1.Angle(temp_sj2.Vect())),weight);

                temp_sj1.SetPxPyPzE(recojet_subjet_Px[ind_sj_qpos_Z],recojet_subjet_Py[ind_sj_qpos_Z],recojet_subjet_Pz[ind_sj_qpos_Z],recojet_subjet_E[ind_sj_qpos_Z])
                temp_sj2.SetPxPyPzE(recojet_subjet_Px[ind_sj_qneg_Z],recojet_subjet_Py[ind_sj_qneg_Z],recojet_subjet_Pz[ind_sj_qneg_Z],recojet_subjet_E[ind_sj_qneg_Z])

                #hist_vec_reco_1D[14].Fill(degrees(temp_sj1.Angle(temp_sj2.Vect())),weight);



                recojet_subjet_rfj_j_svtx_Charge           =ientry.recojet_subjet_rfj_j_svtx_Charge
                recojet_subjet_rfj_j_svtx_nTrack           =ientry.recojet_subjet_rfj_j_svtx_nTrack
                recojet_subjet_rfj_j_svtx_Mass             =ientry.recojet_subjet_rfj_j_svtx_Mass
                recojet_subjet_rfj_j_svtx_E                =ientry.recojet_subjet_rfj_j_svtx_E
                recojet_subjet_rfj_j_svtx_jetindex         =ientry.recojet_subjet_rfj_j_svtx_jetindex

                if quark_Z_decays==0:
                    print "Z decays wrong ",quark_Z_decays


                veto_double_vertex=False;
                #histogram [4] vertex-charge of b matched rfjet of H, [5] vtx charge of bbar matched rfjet of H, [6] vtx ntracks of b match, [7] vtx ntracks of bbar match for H jet
                #loop over events to take out events, where any of the vertices is double
                for s in range(len(recojet_subjet_rfj_j_svtx_Charge)):
                    #print "now in loop ",s
                    for s2 in range(len(recojet_subjet_rfj_j_svtx_Charge)):
                        if (s2!=s and abs(recojet_subjet_rfj_j_svtx_E[s]-recojet_subjet_rfj_j_svtx_E[s2])<1.e-6):
                            veto_double_vertex=True;
                            #print "should break ",s,s2,num_entry
                            print "seems vertices are the same ",s2,s,recojet_subjet_rfj_j_svtx_E[s],recojet_subjet_rfj_j_svtx_E[s2],recojet_subjet_rfj_j_svtx_E[s]-recojet_subjet_rfj_j_svtx_E[s2], len(recojet_subjet_rfj_j_svtx_Charge),num_entry
                            break

                    if veto_double_vertex:
                        #print "we should break now ",s,num_entry
                        break

                if (not veto_double_vertex):
                  for s in range(len(recojet_subjet_rfj_j_svtx_Charge)):
                    #if recojet_subjet_rfj_j_svtx_jetindex[s] == ind_rj_H:
                        #print "fill histo ",s,num_entry
                        #hist_vec_reco_1D[12].Fill(recojet_subjet_rfj_j_svtx_Mass[s],weight); 
                    #else:
                        #if quark_Z_decays<4:
                           #hist_vec_reco_1D[13].Fill(recojet_subjet_rfj_j_svtx_Mass[s],weight); 
                        #elif quark_Z_decays==4:
                           #hist_vec_reco_1D[14].Fill(recojet_subjet_rfj_j_svtx_Mass[s],weight); 
                        #else:
                           #hist_vec_reco_1D[15].Fill(recojet_subjet_rfj_j_svtx_Mass[s],weight); 
                           #if quark_Z_decays!=5:
                               #print "at this point should have only b's left",quark_Z_decays
                    if recojet_subjet_rfj_j_svtx_jetindex[s] == ind_rj_H:
                        hist_vec_reco_1D[4].Fill(recojet_subjet_rfj_j_svtx_Charge[s],weight);   
                        hist_vec_reco_1D[6].Fill(recojet_subjet_rfj_j_svtx_nTrack[s],weight);  
                    if recojet_subjet_rfj_j_svtx_jetindex[s] == ind_rj_Z:
                        hist_vec_reco_1D[5].Fill(recojet_subjet_rfj_j_svtx_Charge[s],weight); 
                        hist_vec_reco_1D[7].Fill(recojet_subjet_rfj_j_svtx_nTrack[s],weight);

                # indices rfj1 b-matched or neg quark, rfj2 bbar matched or pos quark
                #2D hist [0] btag rfj1 and rfj2 for H matched jet, [1] btag rfj and rfj2 for Z matched jet, light decay, [2] btag rfj1 and rfj2 for Z matched jet and ccbar decay  [3] btag rfj1 and rfj2 for Z matched jet and bb decay 
                #2D hist [4] n svtx rfj1 and n svtx rfj2 for H matched jet, [5] n svtx rfj1 and n svtx rfj2 for Z matched, but b decays, [6] n svtx rfj1 and n svtx rfj2 for Z matched, cc decays, [7] n svtx rfj1 and n svtx rfj2, Z matched light decays
                #2D [8] Btag vs Ctag for rfjets from H, [9] Btag vs CTag for Z rfjets, bb decays, [10] BTag vs CTag for Z rfjets, cc decays, [11] BTag vs OTag for Z rfjets, light decays
                #2D hist [12] vtx charge rfj1 and rfj2 for H matched jet, [13] vtx charge rfj and rfj2 for Z matched jet, light decay, [14] vtx charge rfj1 and rfj2 for Z matched jet and ccbar decay  [15] vtx charge rfj1 and rfj2 for Z matched jet and bb decay 
                #2D [16] ctag rfj and rfj2 for Z matched jet, bbbar decay, [17] ctag rfj1 and rfj2 for Z matched jet and ccbar decay [18] ctag rfj1 and rfj2 for Z matched jet light decay 
                #2D [19] Otag rfj and rfj2 for Z matched jet, bbbar decay, [20] Otag rfj1 and rfj2 for Z matched jet and ccbar decay [21] Otag rfj1 and rfj2 for Z matched jet light decay 
                #2D [22] subjet charge b vs bbar matched rfjet of H, [23] subjet charge q pos vs q neg matched rfjet of Z
                #2D [24] subjet charge b matched rfjet vs subjet of H, [25] subjet charge bbar matched rfjet vs subjet of H
                #2D [25] subjet charge q pos matched rfjet vs subjet of Z, [26] subjet charge q neg matched rfjet vs subjet of Z
                #2D [27] btag rfjet 1 vs 2 of H, if m_H jet < m_Z jet, [28] vtx charge rfjet 1 vs 2 of H, if m_H jet < m_Z jet, [29] vtx ntracks rfjet 1 vs 2 of H, if m_H jet < m_Z jet, [30] subjet charge rfjet 1 vs rfjet 2 of H jet with mH<mZ
                #2D [31] btag rfjet 1 vs 2 of Z, if m_H jet < m_Z jet, [32] vtx charge rfjet 1 vs 2 of Z, if m_H jet < m_Z jet, [33] vtx ntracks rfjet 1 vs 2 of H, if m_H jet < m_Z jet, [34] subjet charge rfjet 1 vs rfjet 2 of H jet with mH<mZ
                #2D [35] H vs Z matched jet mass, [36] H matched jet b-tag sum vs Z matched jet b-tag sum, Z b decay, [37] H matched jet b-tag sum vs Z matched jet b-tag sum, Z c decay, [38] H matched jet b-tag sum vs Z matched jet b-tag sum, Z light
                #2D [39] H matched jet b-tag product vs Z matched jet b-tag product, Z b decay, [40] H matched jet b-tag product vs Z matched jet b-tag product, Z c decay, [41] H matched jet b-tag product vs Z matched jet b-tag product, Z light
                #2D [42] H matched jet b-tag c-tag sum vs Z matched jet b-tag c-tag sum, Z b decay, [43] H matched jet b-tag c-tag sum vs Z matched jet b-tag c-tag sum, Z c decay, [44] H matched jet b-tag c-tag sum vs Z matched jet b-tag c-tag sum, Z light
                #2D [45] H matched jet b-tag minus O-tag sum vs Z matched jet b-tag minus O-tag sum, Z b decay, [46] H matched jet b-tag minus O-tag sum vs Z matched jet b-tag minus 0 Tag sum, Z c decay, [47] H matched jet b-tag minus O-tag vs Z matched jet b-tag minus O-0tag, Z light
                #2D [48] H matched jet O-tag sum vs Z matched jet O-tag, Z b decay, [49] H matched jet O-tag vs Z matched jet O-tag, Z c decay, [50] H matched jet O-tag sum vs Z matched jet O-tag sum, Z light
                #2D [51] H matched jet c-tag sum vs Z matched jet c-tag, Z b decay, [52] H matched jet c-tag vs Z matched jet c-tag, Z c decay, [53] H matched jet c-tag sum vs Z matched jet c-tag sum, Z light

            #if use_sqrtJets :
            #    sqrtS_eff_reco=(rj_m1+rj_m2-tempGenIsoPhP4).M();
        
 


        #here use total 4 vector, or 4 vector sum of jet 1 and 2 (see flag below)- isolated photon four vector plus correction with EMiss on both jet axes
        #mass cuts are then also done after projecting the EMiss

    print 'total events after all running signal histos', hist_vec_reco_1D[0].Integral(0,hist_vec_reco_1D[0].GetNbinsX()+1),"masscuts/thetacuts?/betacuts/C2/D2",performMassCuts,performThetaCuts,performBTagCuts,performC2Cuts,performD2Cuts,num_count
    return None



def fill_background_histograms(file,xsec,usePartonInfo,hist_vec_reco_1D,hist_vec_reco_2D,lumi,performMassCuts,useMassRectangleCuts,performThetaCuts,performBTagCuts,performC2Cuts,performD2Cuts):
    print "do something"
  
    use_EMissNeutrinoProjection=True
    #here use total 4 vector, or 4 vector sum of jet 1 and 2 (see flag below)- isolated photon four vector plus correction with EMiss on both jet axes
    #mass cuts are then also done after projecting the EMiss
    use_sqrtJets=True #in this case use j1+j2-isolated photons and with upper flag still decide if EMiss projection on jets is performed
   
    fCut_sj_BTagMax_min=0.9

    fCut_mass_1_min=100.
    fCut_mass_1_max=160.
    fCut_mass_2_min=50.
    fCut_mass_2_max=135.
    #if rectangle or ellipse should be used
    #fCut_use_rectangle=True
    fCut_use_rectangle=useMassRectangleCuts
    #fCut_mass_1_min=100.
    #fCut_mass1_center=126.
    #fCut_mass2_center=92.5
    #fCut_mass1_radius=20.
    #fCut_mass2_radius=20.
    fCut_mass_1_min=0.
    fCut_mass1_center=126.
    fCut_mass2_center=92.5
    fCut_mass1_radius=35.
    fCut_mass2_radius=35.

    fCut_C2_j1_center=0.075
    fCut_C2_j1_radius=0.075
    fCut_C2_j2_center=0.055
    fCut_C2_j2_radius=0.070

    fCut_D2_j1_center=1.75
    fCut_D2_j1_radius=1.25
    fCut_D2_j2_center=1.75
    fCut_D2_j2_radius=1.25

    fCut_thetaWindow=70.
    fCut_thetaRef=90.
    fCut_delta_theta = 100.
 


    fCut_thetaWindow=70.
    fCut_thetaRef=90.

    tree = file.Get("showerData")

    #hist = file2.Get("h_runstatistics")

    #if tree.GetEntries()!=hist.GetBinContent(1) :
    #    print "tree_entries not hist content/diff", tree.GetEntries(),hist.GetBinContent(1), tree.GetEntries()-hist.GetBinContent(1)

    sqrtS_high = 2500.0
    #number of leptons need to be smaller than this number
    m_cut_nLeptons = 1
 
    weight = xsec*lumi/tree.GetEntries()
    print "tree-entries ",tree.GetEntries(), " weight ",weight, "xsec",xsec,"lumi",lumi,"total original ",weight*tree.GetEntries(), 'to', xsec*lumi

    num_entry=-1
    num_count=0
    num_total_exception=0

    for ientry in tree:
        num_entry+=1

        if num_entry%(int(tree.GetEntries()/5.)) == 0:
            print "sig BG in entry ",num_entry

        tempTotRecoP4=TLorentzVector(0,0,0,0);
        tempTotRecoP4.SetPxPyPzE(ientry.totPFO_Px,ientry.totPFO_Py,ientry.totPFO_Pz,ientry.totPFO_E)


        recojet_subjet_rfj_j_E=ientry.recojet_subjet_rfj_j_E
        recojet_subjet_E=ientry.recojet_subjet_E
        #if len(recojet_subjet_E)<4 :
        #if len(recojet_subjet_rfj_j_E)==4 :
        #   print "too few subjets, but enough refined subjets ",len(recojet_subjet_rfj_j_E),num_entry,tempTotRecoP4.E(),len(recojet_subjet_E)
        #else: 
        #   print "too few refined and subjets ",len(recojet_subjet_rfj_j_E),num_entry,len(recojet_subjet_E),tempTotRecoP4.E()
        #continue
        #if len(recojet_subjet_rfj_j_E)<4 :
        #if len(recojet_subjet_E)==4 :
        #   print "too few refined subjets, but enough subjets ",len(recojet_subjet_rfj_j_E),num_entry,len(recojet_subjet_E),tempTotRecoP4.E()
        #else: 
        #print "too few refined and subjets ",len(recojet_subjet_rfj_j_E),num_entry,len(recojet_subjet_E),tempTotRecoP4.E()
        #continue
 
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

        if(len(recojet_E)==2 and n_IsoLep_reco<m_cut_nLeptons) and  sqrtS_eff_reco>sqrtS_high :
            num_count+=weight


        ind_jetM1=0  
        ind_jetM2=1           
        if(rj_m1.M()<rj_m2.M()) :
            ind_jetM2=0  
            ind_jetM1=1  
            temp=TLorentzVector(0,0,0,0)
            temp=rj_m1
            rj_m1=rj_m2
            rj_m2=temp


        fCut_pass_rect_mass_cut=False
        if (rj_m1.M()>fCut_mass_1_min and rj_m1.M()<fCut_mass_1_max) and (rj_m2.M()>fCut_mass_2_min and rj_m2.M()<fCut_mass_2_max):
            fCut_pass_rect_mass_cut=True 
        fCut_pass_ellipse_mass_cut=False
        if (rj_m1.M()>fCut_mass_1_min and (pow((rj_m2.M()-fCut_mass2_center)/fCut_mass2_radius,2)+pow((rj_m1.M()-fCut_mass1_center)/fCut_mass1_radius,2))<1.) :
            fCut_pass_ellipse_mass_cut=True 


        if performMassCuts and (  (fCut_use_rectangle and not fCut_pass_rect_mass_cut)  or (not fCut_use_rectangle and not fCut_pass_ellipse_mass_cut)) :
            continue
        #if performMassCuts and ((rj_m1.M()<fCut_mass_1_min or rj_m1.M()>fCut_mass_1_max) or (rj_m2.M()<fCut_mass_2_min or rj_m2.M()>fCut_mass_2_max)) :
        #     continue

        if performThetaCuts and ((abs(degrees(rj_m1.Theta())-fCut_thetaRef))>fCut_thetaWindow or abs(degrees(rj_m1.Theta()-rj_m2.Theta()))>fCut_delta_theta) :
             continue
        
  
        fCut_pass_upper_ellipse_C2_cut=False
        if((ientry.recojet_beta1_C2[ind_jetM1]<fCut_C2_j1_center and ientry.recojet_beta1_C2[ind_jetM2]<fCut_C2_j2_center) or  (pow((ientry.recojet_beta1_C2[ind_jetM2]-fCut_C2_j2_center)/fCut_C2_j2_radius,2)+pow((ientry.recojet_beta1_C2[ind_jetM1]-fCut_C2_j1_center)/fCut_C2_j1_radius,2))<1.) : 
            fCut_pass_upper_ellipse_C2_cut=True


        if performC2Cuts and not fCut_pass_upper_ellipse_C2_cut :
            continue


        fCut_pass_upper_ellipse_D2_cut=False
        if((ientry.recojet_beta1_D2[ind_jetM1]<fCut_D2_j1_center and ientry.recojet_beta1_D2[ind_jetM2]<fCut_D2_j2_center) or  (pow((ientry.recojet_beta1_D2[ind_jetM2]-fCut_D2_j2_center)/fCut_D2_j2_radius,2)+pow((ientry.recojet_beta1_D2[ind_jetM1]-fCut_D2_j1_center)/fCut_D2_j1_radius,2))<1.) : 
            fCut_pass_upper_ellipse_D2_cut=True


        if performD2Cuts and not fCut_pass_upper_ellipse_D2_cut :
            continue

	if(len(recojet_E)==2 and n_IsoLep_reco<m_cut_nLeptons):
            hist_vec_reco_1D[53].Fill((rj_m1_orig+rj_m2_orig).M(),weight)
            hist_vec_reco_1D[52].Fill((rj1_EMiss+rj2_EMiss).M(),weight)
            if(rj1_EMiss+rj2_EMiss).M()>sqrtS_high:
               hist_vec_reco_1D[54].Fill(rj1_EMiss.E()+rj2_EMiss.E(),weight)
               hist_vec_reco_1D[59].Fill(tempRecoMETP4.Pt()/(rj1_EMiss.E()+rj2_EMiss.E()),weight)
               hist_vec_reco_1D[60].Fill(tempRecoMETP4.Pt()/(rj_m1_orig+rj_m2_orig).M(),weight)
               hist_vec_reco_1D[62].Fill(tempRecoMETP4.Pt()/ientry.totPFO_E,weight)
               hist_vec_reco_1D[63].Fill(rj_m1_orig.E()+rj_m2_orig.E(),weight)
               hist_vec_reco_1D[64].Fill((rj_m1_orig+rj_m2_orig).M(),weight)

            if(rj_m1_orig+rj_m2_orig).M()>sqrtS_high:
               hist_vec_reco_1D[55].Fill(rj1_EMiss.E()+rj2_EMiss.E(),weight)
               hist_vec_reco_1D[56].Fill(rj_m1_orig.E()+rj_m2_orig.E(),weight)
               hist_vec_reco_1D[57].Fill(tempRecoMETP4.Pt()/(rj1_EMiss.E()+rj2_EMiss.E()),weight)
               hist_vec_reco_1D[58].Fill(tempRecoMETP4.Pt()/(rj_m1_orig+rj_m2_orig).M(),weight)
               hist_vec_reco_1D[61].Fill(tempRecoMETP4.Pt()/ientry.totPFO_E,weight)


            #if usePartonInfo and ((rj_m1_orig+rj_m2_orig).M()>sqrtS_high or (rj1_EMiss+rj2_EMiss).M()>sqrtS_high):
            if((rj_m1_orig+rj_m2_orig).M()>sqrtS_high or (rj1_EMiss+rj2_EMiss).M()>sqrtS_high):
               trueME_E=ientry.trueME_E
               trueME_Px=ientry.trueME_Px
               trueME_Py=ientry.trueME_Py
               trueME_Pz=ientry.trueME_Pz
               trueME_PDGID=ientry.trueME_PDGID
               tempMEDieeout=TLorentzVector(0,0,0,0);
               tempMEDiQuark4=TLorentzVector(0,0,0,0);
               tempMEFourQuark4=TLorentzVector(0,0,0,0);
               tempMESixQuarkNoTop4=TLorentzVector(0,0,0,0);
               countQuark=0
               countQuarkNoTop=0
               for i in range(len(trueME_E)):
                   #indices 0-1 outgoing electrons after ISR and beam strahlung (i.e. the "real collision")
                   if(i<2):
                     temp=TLorentzVector(0,0,0,0)
                     temp.SetPxPyPzE(trueME_Px[i],trueME_Py[i],trueME_Pz[i],trueME_E[i])
                     tempMEDieeout+=temp
                   if(trueME_PDGID[i]!=0 and abs(trueME_PDGID[i])<7):
                     temp=TLorentzVector(0,0,0,0)
                     temp.SetPxPyPzE(trueME_Px[i],trueME_Py[i],trueME_Pz[i],trueME_E[i])
                     if countQuark<2:
                        tempMEDiQuark4+=temp
                     if countQuark<4:
                         tempMEFourQuark4+=temp
                     countQuark+=1
                     if abs(trueME_PDGID[i])!=6:
                        if countQuarkNoTop<6:
                            tempMESixQuarkNoTop4+=temp
                        countQuarkNoTop+=1
               if(rj_m1_orig+rj_m2_orig).M()>sqrtS_high:
                   hist_vec_reco_1D[111].Fill(tempMEDiQuark4.M(),weight)
                   hist_vec_reco_1D[113].Fill(tempMEFourQuark4.M(),weight)
                   hist_vec_reco_1D[115].Fill(tempMESixQuarkNoTop4.M(),weight)
                   hist_vec_reco_1D[117].Fill(tempMEDieeout.M(),weight)
                   hist_vec_reco_1D[119].Fill(tempMEDiQuark4.E(),weight)
                   hist_vec_reco_1D[121].Fill(tempMEFourQuark4.E(),weight)
                   hist_vec_reco_1D[123].Fill(tempMESixQuarkNoTop4.E(),weight)
                   hist_vec_reco_1D[125].Fill(tempMEDieeout.E(),weight)
               if(rj1_EMiss+rj2_EMiss).M()>sqrtS_high:
                   hist_vec_reco_1D[112].Fill(tempMEDiQuark4.M(),weight)
                   hist_vec_reco_1D[114].Fill(tempMEFourQuark4.M(),weight)
                   hist_vec_reco_1D[116].Fill(tempMESixQuarkNoTop4.M(),weight)
                   hist_vec_reco_1D[118].Fill(tempMEDieeout.M(),weight)
                   hist_vec_reco_1D[120].Fill(tempMEDiQuark4.E(),weight)
                   hist_vec_reco_1D[122].Fill(tempMEFourQuark4.E(),weight)
                   hist_vec_reco_1D[124].Fill(tempMESixQuarkNoTop4.E(),weight)
                   hist_vec_reco_1D[126].Fill(tempMEDieeout.E(),weight)
 
    
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


                #checked at the beginning by requirement of 4 indices
                #if len(recojet_subjet_rfj_j_E)!=len(recojet_subjet_E):
                #print "subjet and rf subjet length not the same ",len(recojet_subjet_rfj_j_E),len(recojet_subjet_E)
                #continue;
                rj_m1_rfj=TLorentzVector(0,0,0,0);
                rj_m2_rfj=TLorentzVector(0,0,0,0);
                fill_jet_mass_1_subjet_histos=True
                fill_jet_mass_2_subjet_histos=True
                for i in range(len(recojet_subjet_rfj_j_E)):
                    if recojet_subjet_rfj_j_CTag[i]> CTagMax :
                        CTagMax=recojet_subjet_rfj_j_CTag[i]
                        ind_sj_CTagMax=i
                    if recojet_subjet_rfj_j_OTag[i]> LTagMax :
                        LTagMax=recojet_subjet_rfj_j_OTag[i]
                        ind_sj_LTagMax=i
                    #subjets filled in order, meaning if sj x is the first one with jetindex in question, sj x+1 is the other one of the pair
                    #check for exactly 4 subjets in collection happened earlier on
                    temp=TLorentzVector(0,0,0,0)
                    temp.SetPxPyPzE(recojet_subjet_rfj_j_Px[i],recojet_subjet_rfj_j_Py[i],recojet_subjet_rfj_j_Pz[i],recojet_subjet_rfj_j_E[i])
                    if(fill_jet_mass_1_subjet_histos and recojet_subjet_rfj_j_jetindex[i]==ind_jetM1):
                        temp2=TLorentzVector(0,0,0,0)
                        temp2.SetPxPyPzE(recojet_subjet_rfj_j_Px[i+1],recojet_subjet_rfj_j_Py[i+1],recojet_subjet_rfj_j_Pz[i],recojet_subjet_rfj_j_E[i+1])
                        if recojet_subjet_rfj_j_jetindex[i]!=recojet_subjet_rfj_j_jetindex[i+1]:
                            print 'sj indices were supposed to be the same, rfj mass1'
                        ind_sj1_j1=i
                        ind_sj2_j1=i+1
                        if recojet_subjet_rfj_j_BTag[i+1]>recojet_subjet_rfj_j_BTag[i]:
                           ind_sj2_j1=i
                           ind_sj1_j1=i+1
                        hist_vec_reco_1D[33].Fill(recojet_subjet_rfj_j_BTag[ind_sj1_j1],weight) 
                        if performBTagCuts and BTagMax>fCut_sj_BTagMax_min :
                            hist_vec_reco_1D[38].Fill(recojet_subjet_rfj_j_E[ind_sj1_j1]/(recojet_subjet_rfj_j_E[ind_sj1_j1]+recojet_subjet_rfj_j_E[ind_sj2_j1]),weight) 
                            hist_vec_reco_1D[44].Fill(degrees(temp.Angle(temp2.Vect())),weight)
                            if recojet_subjet_rfj_j_E[ind_sj1_j1]>recojet_subjet_rfj_j_E[ind_sj2_j1]:
                                hist_vec_reco_1D[39].Fill(recojet_subjet_rfj_j_E[ind_sj1_j1]/(recojet_subjet_rfj_j_E[ind_sj1_j1]+recojet_subjet_rfj_j_E[ind_sj2_j1]),weight) 
                            else:
                                hist_vec_reco_1D[39].Fill(recojet_subjet_rfj_j_E[ind_sj2_j1]/(recojet_subjet_rfj_j_E[ind_sj1_j1]+recojet_subjet_rfj_j_E[ind_sj2_j1]),weight) 
                        fill_jet_mass_1_subjet_histos=False
                    if(fill_jet_mass_2_subjet_histos and recojet_subjet_rfj_j_jetindex[i]==ind_jetM2):
                        temp2=TLorentzVector(0,0,0,0)
                        temp2.SetPxPyPzE(recojet_subjet_rfj_j_Px[i+1],recojet_subjet_rfj_j_Py[i+1],recojet_subjet_rfj_j_Pz[i],recojet_subjet_rfj_j_E[i+1])
                        if recojet_subjet_rfj_j_jetindex[i]!=recojet_subjet_rfj_j_jetindex[i+1]:
                            print 'sj indices were supposed to be the same, rfj mass2'
                        ind_sj1_j2=i
                        ind_sj2_j2=i+1
                        if recojet_subjet_rfj_j_E[i+1]>recojet_subjet_rfj_j_E[i]:
                           ind_sj2_j2=i
                           ind_sj1_j2=i+1
                        hist_vec_reco_1D[34].Fill(recojet_subjet_rfj_j_BTag[ind_sj1_j2],weight) 
                        hist_vec_reco_1D[35].Fill(recojet_subjet_rfj_j_CTag[ind_sj1_j2],weight) 
                        hist_vec_reco_1D[36].Fill(recojet_subjet_rfj_j_OTag[ind_sj1_j2],weight) 
                        if performBTagCuts and BTagMax>fCut_sj_BTagMax_min :
                            hist_vec_reco_1D[47].Fill(recojet_subjet_rfj_j_BTag[ind_sj1_j2],weight) 
                            hist_vec_reco_1D[48].Fill(recojet_subjet_rfj_j_CTag[ind_sj1_j2],weight) 
                            hist_vec_reco_1D[49].Fill(recojet_subjet_rfj_j_OTag[ind_sj1_j2],weight) 
                            hist_vec_reco_1D[41].Fill(recojet_subjet_rfj_j_E[ind_sj1_j2]/(recojet_subjet_rfj_j_E[ind_sj1_j2]+recojet_subjet_rfj_j_E[ind_sj2_j2]),weight) 
                            hist_vec_reco_1D[46].Fill(degrees(temp.Angle(temp2.Vect())),weight)
                        if recojet_subjet_rfj_j_BTag[ind_sj1_j2]>recojet_subjet_rfj_j_BTag[ind_sj2_j2]:
                           hist_vec_reco_1D[37].Fill(recojet_subjet_rfj_j_BTag[ind_sj1_j2],weight) 
                           if performBTagCuts and BTagMax>fCut_sj_BTagMax_min :
                               hist_vec_reco_1D[50].Fill(recojet_subjet_rfj_j_BTag[ind_sj1_j2],weight) 
                        else :
                           hist_vec_reco_1D[37].Fill(recojet_subjet_rfj_j_BTag[ind_sj2_j2],weight) 
                           if performBTagCuts and BTagMax>fCut_sj_BTagMax_min :
                               hist_vec_reco_1D[50].Fill(recojet_subjet_rfj_j_BTag[ind_sj2_j2],weight) 
                        fill_jet_mass_2_subjet_histos=False
                    #print "sqrtS angles ",sqrtS_eff_reco,H_decays_bbar
                    #print "loop over subjet ",i,ind_rj_H,ind_rj_Z
                    if (recojet_subjet_rfj_j_jetindex[i]==ind_jetM1):
                        rj_m1_rfj+=temp;
                    else:
                        rj_m2_rfj+=temp; 
    

                #for matching between the original jet and refined subjet
                recojet_subjet_E                     =ientry.recojet_subjet_E
                recojet_subjet_jetindex              =ientry.recojet_subjet_jetindex
                recojet_subjet_Px                    =ientry.recojet_subjet_Px
                recojet_subjet_Py                    =ientry.recojet_subjet_Py
                recojet_subjet_Pz                    =ientry.recojet_subjet_Pz
                recojet_subjet_jetChargePt_kappa_0_30=ientry.recojet_subjet_jetChargePt_kappa_0_30


                if BTagMax==0:
                    continue;


                if performBTagCuts and BTagMax<fCut_sj_BTagMax_min :
                    continue;


                fill_jet_mass_1_subjet_histos=True
                fill_jet_mass_2_subjet_histos=True
                for i in range(len(recojet_subjet_E)):
                    #subjets filled in order, meaning if sj x is the first one with jetindex in question, sj x+1 is the other one of the pair
                    #check for exactly 4 subjets in collection happened earlier on
                    temp=TLorentzVector(0,0,0,0)
                    temp.SetPxPyPzE(recojet_subjet_Px[i],recojet_subjet_Py[i],recojet_subjet_Pz[i],recojet_subjet_E[i])
                    if(fill_jet_mass_1_subjet_histos and recojet_subjet_jetindex[i]==ind_jetM1):
                        temp2=TLorentzVector(0,0,0,0)
                        if recojet_subjet_jetindex[i]!=recojet_subjet_jetindex[i+1]:
                            print 'sj indices were supposed to be the same, mass1'
                        temp2.SetPxPyPzE(recojet_subjet_Px[i+1],recojet_subjet_Py[i+1],recojet_subjet_Pz[i],recojet_subjet_E[i+1])
                        if recojet_subjet_E[i+1]>recojet_subjet_E[i]:
                            hist_vec_reco_1D[40].Fill(recojet_subjet_E[i+1]/(recojet_subjet_E[i]+recojet_subjet_E[i+1]),weight) 
                        else: 
                            hist_vec_reco_1D[40].Fill(recojet_subjet_E[i]/(recojet_subjet_E[i]+recojet_subjet_E[i+1]),weight) 
                        hist_vec_reco_1D[43].Fill(degrees(temp.Angle(temp2.Vect())),weight)
                        fill_jet_mass_1_subjet_histos=False
                    if(fill_jet_mass_2_subjet_histos and recojet_subjet_jetindex[i]==ind_jetM2):
                        temp2=TLorentzVector(0,0,0,0)
                        temp2.SetPxPyPzE(recojet_subjet_Px[i+1],recojet_subjet_Py[i+1],recojet_subjet_Pz[i],recojet_subjet_E[i+1])
                        if recojet_subjet_jetindex[i]!=recojet_subjet_jetindex[i+1]:
                            print 'sj indices were supposed to be the same, mass2'
                        if recojet_subjet_E[i+1]>recojet_subjet_E[i]:
                            hist_vec_reco_1D[42].Fill(recojet_subjet_E[i+1]/(recojet_subjet_E[i]+recojet_subjet_E[i+1]),weight) 
                        else: 
                            hist_vec_reco_1D[42].Fill(recojet_subjet_E[i]/(recojet_subjet_E[i]+recojet_subjet_E[i+1]),weight) 
                        hist_vec_reco_1D[45].Fill(degrees(temp.Angle(temp2.Vect())),weight) 
                        fill_jet_mass_2_subjet_histos=False
               
                recojet_subjet_rfj_j_svtx_Charge   = ientry.recojet_subjet_rfj_j_svtx_Charge
                recojet_subjet_rfj_j_svtx_r        = ientry.recojet_subjet_rfj_j_svtx_r
                recojet_subjet_rfj_j_svtx_E        = ientry.recojet_subjet_rfj_j_svtx_E
                recojet_subjet_rfj_j_svtx_Mass     = ientry.recojet_subjet_rfj_j_svtx_Mass
                recojet_subjet_rfj_j_svtx_nTrack   = ientry.recojet_subjet_rfj_j_svtx_nTrack
                recojet_subjet_rfj_j_svtx_jetindex = ientry.recojet_subjet_rfj_j_svtx_jetindex




                if(recojet_subjet_jetindex[0]!=recojet_subjet_jetindex[1]) or (recojet_subjet_jetindex[2]!=recojet_subjet_jetindex[3]):
                    print "from first principle these subjet indices should be the same, but why aren't they "<<recojet_subjet_jetindex[0],recojet_subjet_jetindex[1],recojet_subjet_jetindex[2],recojet_subjet_jetindex[3]
              
                hist_vec_reco_1D[51].Fill(rj_m1.M()-rj_m2.M(),weight)
                hist_vec_reco_1D[0].Fill(rj_m1.M(),weight)
                hist_vec_reco_1D[1].Fill(rj_m2.M(),weight)
                hist_vec_reco_1D[2].Fill(rj_m1_rfj.M(),weight)
                hist_vec_reco_1D[3].Fill(rj_m2_rfj.M(),weight)

                hist_vec_reco_1D[65].Fill(ientry.recojet_beta1_N2[ind_jetM1],weight)
                hist_vec_reco_1D[66].Fill(ientry.recojet_beta2_N2[ind_jetM1],weight)
                hist_vec_reco_1D[67].Fill(ientry.recojet_beta0_5_N2[ind_jetM1],weight)

                hist_vec_reco_1D[68].Fill(ientry.recojet_beta1_N3[ind_jetM1],weight)
                hist_vec_reco_1D[69].Fill(ientry.recojet_beta2_N3[ind_jetM1],weight)
                hist_vec_reco_1D[70].Fill(ientry.recojet_beta0_5_N3[ind_jetM1],weight)

                hist_vec_reco_1D[71].Fill(ientry.recojet_beta1_C2[ind_jetM1],weight)
                hist_vec_reco_1D[72].Fill(ientry.recojet_beta2_C2[ind_jetM1],weight)
                hist_vec_reco_1D[73].Fill(ientry.recojet_beta0_5_C2[ind_jetM1],weight)

                hist_vec_reco_1D[74].Fill(ientry.recojet_beta1_C3[ind_jetM1],weight)
                hist_vec_reco_1D[75].Fill(ientry.recojet_beta2_C3[ind_jetM1],weight)
                hist_vec_reco_1D[76].Fill(ientry.recojet_beta0_5_C3[ind_jetM1],weight)

                hist_vec_reco_1D[77].Fill(ientry.recojet_beta1_D2[ind_jetM1],weight)
                hist_vec_reco_1D[78].Fill(ientry.recojet_beta2_D2[ind_jetM1],weight)
                hist_vec_reco_1D[79].Fill(ientry.recojet_beta0_5_D2[ind_jetM1],weight)

                hist_vec_reco_1D[80].Fill(ientry.recojet_nsubjettiness2[ind_jetM1]/ientry.recojet_nsubjettiness1[ind_jetM1],weight)
                if(ientry.recojet_nsubjettiness2[ind_jetM1]!=0):
                    hist_vec_reco_1D[81].Fill(ientry.recojet_nsubjettiness3[ind_jetM1]/ientry.recojet_nsubjettiness2[ind_jetM1],weight)

                hist_vec_reco_1D[128].Fill(ientry.recojet_beta1_N2_E_theta[ind_jetM1],weight)
                hist_vec_reco_1D[129].Fill(ientry.recojet_beta2_N2_E_theta[ind_jetM1],weight)
                hist_vec_reco_1D[130].Fill(ientry.recojet_beta0_5_N2_E_theta[ind_jetM1],weight)

                hist_vec_reco_1D[131].Fill(ientry.recojet_beta1_N3_E_theta[ind_jetM1],weight)
                hist_vec_reco_1D[132].Fill(ientry.recojet_beta2_N3_E_theta[ind_jetM1],weight)
                hist_vec_reco_1D[133].Fill(ientry.recojet_beta0_5_N3_E_theta[ind_jetM1],weight)

                hist_vec_reco_1D[134].Fill(ientry.recojet_beta1_C2_E_theta[ind_jetM1],weight)
                hist_vec_reco_1D[135].Fill(ientry.recojet_beta2_C2_E_theta[ind_jetM1],weight)
                hist_vec_reco_1D[136].Fill(ientry.recojet_beta0_5_C2_E_theta[ind_jetM1],weight)

                hist_vec_reco_1D[137].Fill(ientry.recojet_beta1_C3_E_theta[ind_jetM1],weight)
                hist_vec_reco_1D[138].Fill(ientry.recojet_beta2_C3_E_theta[ind_jetM1],weight)
                hist_vec_reco_1D[139].Fill(ientry.recojet_beta0_5_C3_E_theta[ind_jetM1],weight)

                hist_vec_reco_1D[140].Fill(ientry.recojet_beta1_D2_E_theta[ind_jetM1],weight)
                hist_vec_reco_1D[141].Fill(ientry.recojet_beta2_D2_E_theta[ind_jetM1],weight)
                hist_vec_reco_1D[142].Fill(ientry.recojet_beta0_5_D2_E_theta[ind_jetM1],weight)

                #print 'length in background hists',len(hist_vec_reco_1D),len(hist_vec_reco_2D),ind_jetM1


                #hist_vec_reco_1D[143].Fill(ientry.recojet_nsubjettiness2_lrz[ind_jetM1]/ientry.recojet_nsubjettiness1_lrz[ind_jetM1],weight)
                #if(ientry.recojet_nsubjettiness2_lrz[ind_jetM1]!=0):
                #    hist_vec_reco_1D[144].Fill(ientry.recojet_nsubjettiness3_lrz[ind_jetM1]/ientry.recojet_nsubjettiness2_lrz[ind_jetM1],weight)

                hist_vec_reco_1D[82].Fill(ientry.recojet_dij_21[ind_jetM1]/pow(recojet_E[ind_jetM1],2),weight)
                hist_vec_reco_1D[83].Fill(ientry.recojet_dij_32[ind_jetM1]/pow(recojet_E[ind_jetM1],2,),weight)
                hist_vec_reco_1D[84].Fill(ientry.recojet_dij_43[ind_jetM1]/pow(recojet_E[ind_jetM1],2,),weight)

                if(ientry.recojet_dij_21[ind_jetM1]!=0):
                    hist_vec_reco_1D[85].Fill(ientry.recojet_dij_32[ind_jetM1]/ientry.recojet_dij_21[ind_jetM1],weight)
                    hist_vec_reco_1D[86].Fill(ientry.recojet_dij_43[ind_jetM1]/ientry.recojet_dij_21[ind_jetM1],weight)
                if(ientry.recojet_dij_32[ind_jetM1]!=0):
                    hist_vec_reco_1D[87].Fill(ientry.recojet_dij_43[ind_jetM1]/ientry.recojet_dij_32[ind_jetM1],weight)


                hist_vec_reco_1D[88].Fill(ientry.recojet_beta1_N2[ind_jetM2],weight)
                hist_vec_reco_1D[89].Fill(ientry.recojet_beta2_N2[ind_jetM2],weight)
                hist_vec_reco_1D[90].Fill(ientry.recojet_beta0_5_N2[ind_jetM2],weight)

                hist_vec_reco_1D[91].Fill(ientry.recojet_beta1_N3[ind_jetM2],weight)
                hist_vec_reco_1D[92].Fill(ientry.recojet_beta2_N3[ind_jetM2],weight)
                hist_vec_reco_1D[93].Fill(ientry.recojet_beta0_5_N3[ind_jetM2],weight)

                hist_vec_reco_1D[94].Fill(ientry.recojet_beta1_C2[ind_jetM2],weight)
                hist_vec_reco_1D[95].Fill(ientry.recojet_beta2_C2[ind_jetM2],weight)
                hist_vec_reco_1D[96].Fill(ientry.recojet_beta0_5_C2[ind_jetM2],weight)

                hist_vec_reco_1D[97].Fill(ientry.recojet_beta1_C3[ind_jetM2],weight)
                hist_vec_reco_1D[98].Fill(ientry.recojet_beta2_C3[ind_jetM2],weight)
                hist_vec_reco_1D[99].Fill(ientry.recojet_beta0_5_C3[ind_jetM2],weight)

                hist_vec_reco_1D[100].Fill(ientry.recojet_beta1_D2[ind_jetM2],weight)
                hist_vec_reco_1D[101].Fill(ientry.recojet_beta2_D2[ind_jetM2],weight)
                hist_vec_reco_1D[102].Fill(ientry.recojet_beta0_5_D2[ind_jetM2],weight)

                hist_vec_reco_1D[103].Fill(ientry.recojet_nsubjettiness2[ind_jetM2]/ientry.recojet_nsubjettiness1[ind_jetM2],weight)
                if(ientry.recojet_nsubjettiness2[ind_jetM2]!=0):
                    hist_vec_reco_1D[104].Fill(ientry.recojet_nsubjettiness3[ind_jetM2]/ientry.recojet_nsubjettiness2[ind_jetM2],weight)

                hist_vec_reco_1D[145].Fill(ientry.recojet_beta1_N2_E_theta[ind_jetM2],weight)
                hist_vec_reco_1D[146].Fill(ientry.recojet_beta2_N2_E_theta[ind_jetM2],weight)
                hist_vec_reco_1D[147].Fill(ientry.recojet_beta0_5_N2_E_theta[ind_jetM2],weight)

                hist_vec_reco_1D[148].Fill(ientry.recojet_beta1_N3_E_theta[ind_jetM2],weight)
                hist_vec_reco_1D[149].Fill(ientry.recojet_beta2_N3_E_theta[ind_jetM2],weight)
                hist_vec_reco_1D[150].Fill(ientry.recojet_beta0_5_N3_E_theta[ind_jetM2],weight)

                hist_vec_reco_1D[151].Fill(ientry.recojet_beta1_C2_E_theta[ind_jetM2],weight)
                hist_vec_reco_1D[152].Fill(ientry.recojet_beta2_C2_E_theta[ind_jetM2],weight)
                hist_vec_reco_1D[153].Fill(ientry.recojet_beta0_5_C2_E_theta[ind_jetM2],weight)

                hist_vec_reco_1D[154].Fill(ientry.recojet_beta1_C3_E_theta[ind_jetM2],weight)
                hist_vec_reco_1D[155].Fill(ientry.recojet_beta2_C3_E_theta[ind_jetM2],weight)
                hist_vec_reco_1D[156].Fill(ientry.recojet_beta0_5_C3_E_theta[ind_jetM2],weight)

                hist_vec_reco_1D[157].Fill(ientry.recojet_beta1_D2_E_theta[ind_jetM2],weight)
                hist_vec_reco_1D[158].Fill(ientry.recojet_beta2_D2_E_theta[ind_jetM2],weight)
                hist_vec_reco_1D[159].Fill(ientry.recojet_beta0_5_D2_E_theta[ind_jetM2],weight)

                hist_vec_reco_1D[162].Fill(ientry.reco_y32,weight)

                #hist_vec_reco_1D[160].Fill(ientry.recojet_nsubjettiness2_lrz[ind_jetM2]/ientry.recojet_nsubjettiness1_lrz[ind_jetM2],weight)
                #if(ientry.recojet_nsubjettiness2_lrz[ind_jetM2]!=0):
                #    hist_vec_reco_1D[161].Fill(ientry.recojet_nsubjettiness3_lrz[ind_jetM2]/ientry.recojet_nsubjettiness2_lrz[ind_jetM2],weight)

                hist_vec_reco_1D[105].Fill(ientry.recojet_dij_21[ind_jetM2]/pow(recojet_E[ind_jetM2],2),weight)
                hist_vec_reco_1D[106].Fill(ientry.recojet_dij_32[ind_jetM2]/pow(recojet_E[ind_jetM2],2),weight)
                hist_vec_reco_1D[107].Fill(ientry.recojet_dij_43[ind_jetM2]/pow(recojet_E[ind_jetM2],2),weight)
                if(ientry.recojet_dij_21[ind_jetM2]!=0):
                    hist_vec_reco_1D[108].Fill(ientry.recojet_dij_32[ind_jetM2]/ientry.recojet_dij_21[ind_jetM2],weight)
                    hist_vec_reco_1D[109].Fill(ientry.recojet_dij_43[ind_jetM2]/ientry.recojet_dij_21[ind_jetM2],weight)
                if(ientry.recojet_dij_32[ind_jetM2]!=0):
                    hist_vec_reco_1D[110].Fill(ientry.recojet_dij_43[ind_jetM2]/ientry.recojet_dij_32[ind_jetM2],weight)

                hist_vec_reco_2D[0].Fill(rj_m1.M(),rj_m2.M(),weight)
                hist_vec_reco_2D[1].Fill(rj_m1_rfj.M(),rj_m2_rfj.M(),weight)
                hist_vec_reco_2D[2].Fill(degrees(rj_m1.Theta()),degrees(rj_m2.Theta()),weight)

                hist_vec_reco_2D[20].Fill(ientry.recojet_beta1_N2[ind_jetM1],ientry.recojet_beta1_N2[ind_jetM2],weight)

                hist_vec_reco_2D[21].Fill(ientry.recojet_beta2_N2[ind_jetM1],ientry.recojet_beta2_N2[ind_jetM2],weight)
                hist_vec_reco_2D[22].Fill(ientry.recojet_beta0_5_N2[ind_jetM1],ientry.recojet_beta0_5_N2[ind_jetM2],weight)

                hist_vec_reco_2D[23].Fill(ientry.recojet_beta1_N3[ind_jetM1],ientry.recojet_beta1_N3[ind_jetM2],weight)
                hist_vec_reco_2D[24].Fill(ientry.recojet_beta2_N3[ind_jetM1],ientry.recojet_beta2_N3[ind_jetM2],weight)
                hist_vec_reco_2D[25].Fill(ientry.recojet_beta0_5_N3[ind_jetM1],ientry.recojet_beta0_5_N3[ind_jetM2],weight)

                hist_vec_reco_2D[26].Fill(ientry.recojet_beta1_C2[ind_jetM1],ientry.recojet_beta1_C2[ind_jetM2],weight)
                hist_vec_reco_2D[27].Fill(ientry.recojet_beta2_C2[ind_jetM1],ientry.recojet_beta2_C2[ind_jetM2],weight)
                hist_vec_reco_2D[28].Fill(ientry.recojet_beta0_5_C2[ind_jetM1],ientry.recojet_beta0_5_C2[ind_jetM2],weight)

                hist_vec_reco_2D[29].Fill(ientry.recojet_beta1_C3[ind_jetM1],ientry.recojet_beta1_C3[ind_jetM2],weight)
                hist_vec_reco_2D[30].Fill(ientry.recojet_beta2_C3[ind_jetM1],ientry.recojet_beta2_C3[ind_jetM2],weight)
                hist_vec_reco_2D[31].Fill(ientry.recojet_beta0_5_C3[ind_jetM1],ientry.recojet_beta0_5_C3[ind_jetM2],weight)

                hist_vec_reco_2D[32].Fill(ientry.recojet_beta1_D2[ind_jetM1],ientry.recojet_beta1_D2[ind_jetM2],weight)
                hist_vec_reco_2D[33].Fill(ientry.recojet_beta2_D2[ind_jetM1],ientry.recojet_beta2_D2[ind_jetM2],weight)
                hist_vec_reco_2D[34].Fill(ientry.recojet_beta0_5_D2[ind_jetM1],ientry.recojet_beta0_5_D2[ind_jetM2],weight)

                hist_vec_reco_2D[35].Fill(ientry.recojet_nsubjettiness2[ind_jetM1]/ientry.recojet_nsubjettiness1[ind_jetM1],ientry.recojet_nsubjettiness2[ind_jetM2]/ientry.recojet_nsubjettiness1[ind_jetM2],weight)
                if(ientry.recojet_nsubjettiness2[ind_jetM2]!=0 and ientry.recojet_nsubjettiness2[ind_jetM1]!=0):
                    hist_vec_reco_2D[36].Fill(ientry.recojet_nsubjettiness3[ind_jetM1]/ientry.recojet_nsubjettiness2[ind_jetM1],ientry.recojet_nsubjettiness3[ind_jetM2]/ientry.recojet_nsubjettiness2[ind_jetM2],weight)
                    
                hist_vec_reco_2D[43].Fill(ientry.recojet_beta1_N2_E_theta[ind_jetM1],ientry.recojet_beta1_N2_E_theta[ind_jetM2],weight)

                hist_vec_reco_2D[44].Fill(ientry.recojet_beta2_N2_E_theta[ind_jetM1],ientry.recojet_beta2_N2_E_theta[ind_jetM2],weight)
                hist_vec_reco_2D[45].Fill(ientry.recojet_beta0_5_N2_E_theta[ind_jetM1],ientry.recojet_beta0_5_N2_E_theta[ind_jetM2],weight)

                hist_vec_reco_2D[46].Fill(ientry.recojet_beta1_N3_E_theta[ind_jetM1],ientry.recojet_beta1_N3_E_theta[ind_jetM2],weight)
                hist_vec_reco_2D[47].Fill(ientry.recojet_beta2_N3_E_theta[ind_jetM1],ientry.recojet_beta2_N3_E_theta[ind_jetM2],weight)
                hist_vec_reco_2D[48].Fill(ientry.recojet_beta0_5_N3_E_theta[ind_jetM1],ientry.recojet_beta0_5_N3_E_theta[ind_jetM2],weight)

                hist_vec_reco_2D[49].Fill(ientry.recojet_beta1_C2_E_theta[ind_jetM1],ientry.recojet_beta1_C2_E_theta[ind_jetM2],weight)
                hist_vec_reco_2D[50].Fill(ientry.recojet_beta2_C2_E_theta[ind_jetM1],ientry.recojet_beta2_C2_E_theta[ind_jetM2],weight)
                hist_vec_reco_2D[51].Fill(ientry.recojet_beta0_5_C2_E_theta[ind_jetM1],ientry.recojet_beta0_5_C2_E_theta[ind_jetM2],weight)

                hist_vec_reco_2D[52].Fill(ientry.recojet_beta1_C3_E_theta[ind_jetM1],ientry.recojet_beta1_C3_E_theta[ind_jetM2],weight)
                hist_vec_reco_2D[53].Fill(ientry.recojet_beta2_C3_E_theta[ind_jetM1],ientry.recojet_beta2_C3_E_theta[ind_jetM2],weight)
                hist_vec_reco_2D[54].Fill(ientry.recojet_beta0_5_C3_E_theta[ind_jetM1],ientry.recojet_beta0_5_C3_E_theta[ind_jetM2],weight)

                hist_vec_reco_2D[55].Fill(ientry.recojet_beta1_D2_E_theta[ind_jetM1],ientry.recojet_beta1_D2_E_theta[ind_jetM2],weight)
                hist_vec_reco_2D[56].Fill(ientry.recojet_beta2_D2_E_theta[ind_jetM1],ientry.recojet_beta2_D2_E_theta[ind_jetM2],weight)
                hist_vec_reco_2D[57].Fill(ientry.recojet_beta0_5_D2_E_theta[ind_jetM1],ientry.recojet_beta0_5_D2_E_theta[ind_jetM2],weight)

                #hist_vec_reco_2D[58].Fill(ientry.recojet_nsubjettiness2_lrz[ind_jetM1]/ientry.recojet_nsubjettiness1_lrz[ind_jetM1],ientry.recojet_nsubjettiness2_lrz[ind_jetM2]/ientry.recojet_nsubjettiness1_lrz[ind_jetM2],weight)
                #if(ientry.recojet_nsubjettiness2_lrz[ind_jetM2]!=0 and ientry.recojet_nsubjettiness2_lrz[ind_jetM1]!=0):
                #    hist_vec_reco_2D[59].Fill(ientry.recojet_nsubjettiness3_lrz[ind_jetM1]/ientry.recojet_nsubjettiness2_lrz[ind_jetM1],ientry.recojet_nsubjettiness3_lrz[ind_jetM2]/ientry.recojet_nsubjettiness2_lrz[ind_jetM2],weight)

                hist_vec_reco_2D[37].Fill(ientry.recojet_dij_21[ind_jetM1]/pow(recojet_E[ind_jetM1],2),ientry.recojet_dij_21[ind_jetM2]/pow(recojet_E[ind_jetM2],2),weight)
                hist_vec_reco_2D[38].Fill(ientry.recojet_dij_32[ind_jetM1]/pow(recojet_E[ind_jetM1],2),ientry.recojet_dij_32[ind_jetM2]/pow(recojet_E[ind_jetM2],2),weight)
                hist_vec_reco_2D[39].Fill(ientry.recojet_dij_43[ind_jetM1]/pow(recojet_E[ind_jetM1],2),ientry.recojet_dij_43[ind_jetM2]/pow(recojet_E[ind_jetM2],2),weight)
                if(ientry.recojet_dij_21[ind_jetM1]!=0 and ientry.recojet_dij_21[ind_jetM2]!=0):
                    hist_vec_reco_2D[40].Fill(ientry.recojet_dij_32[ind_jetM1]/ientry.recojet_dij_21[ind_jetM1],ientry.recojet_dij_32[ind_jetM2]/ientry.recojet_dij_21[ind_jetM2],weight)
                    hist_vec_reco_2D[41].Fill(ientry.recojet_dij_43[ind_jetM1]/ientry.recojet_dij_21[ind_jetM1],ientry.recojet_dij_43[ind_jetM2]/ientry.recojet_dij_21[ind_jetM2],weight)
                if(ientry.recojet_dij_32[ind_jetM1]!=0 and ientry.recojet_dij_32[ind_jetM2]!=0):
                    hist_vec_reco_2D[42].Fill(ientry.recojet_dij_43[ind_jetM1]/ientry.recojet_dij_32[ind_jetM1],ientry.recojet_dij_43[ind_jetM2]/ientry.recojet_dij_32[ind_jetM2],weight)


                #check if mass-ordered and max BTag subjet jet index are the same, use original mass
                if recojet_subjet_rfj_j_jetindex[ind_sj_BTagMax]==ind_jetM1:
                    hist_vec_reco_1D[4].Fill(rj_m1.M(),weight)
                    hist_vec_reco_1D[5].Fill(rj_m2.M(),weight)
                    hist_vec_reco_2D[3].Fill(rj_m1.M(),rj_m2.M(),weight)
                else:
                    hist_vec_reco_1D[4].Fill(rj_m2.M(),weight)
                    hist_vec_reco_1D[5].Fill(rj_m1.M(),weight)   
                    hist_vec_reco_2D[3].Fill(rj_m2.M(),rj_m1.M(),weight)
                #check if mass-ordered and max LTag subjet jet index are the same, use original mass
                if recojet_subjet_rfj_j_jetindex[ind_sj_LTagMax]==ind_jetM1:
                    hist_vec_reco_1D[6].Fill(rj_m1.M(),weight)
                    hist_vec_reco_1D[7].Fill(rj_m2.M(),weight)
                    #here the most light flavored subjet is supposed to come from Z, expect that to be worse
                    hist_vec_reco_2D[4].Fill(rj_m2.M(),rj_m1.M(),weight)
                else:
                    hist_vec_reco_1D[6].Fill(rj_m2.M(),weight)
                    hist_vec_reco_1D[7].Fill(rj_m1.M(),weight) 
                    hist_vec_reco_2D[4].Fill(rj_m1.M(),rj_m2.M(),weight)
                #check if mass-ordered and BTag sum jet index are the same, use original mass
                if ind_j_BTagSumMax==ind_jetM1:
                    hist_vec_reco_1D[8].Fill(rj_m1.M(),weight)
                    hist_vec_reco_1D[9].Fill(rj_m2.M(),weight)
                    hist_vec_reco_2D[5].Fill(rj_m1.M(),rj_m2.M(),weight)
                else:
                    hist_vec_reco_1D[8].Fill(rj_m2.M(),weight)
                    hist_vec_reco_1D[9].Fill(rj_m1.M(),weight)   
                    hist_vec_reco_2D[5].Fill(rj_m2.M(),rj_m1.M(),weight)
                #check if mass-ordered and max LTag subjet jet index are the same, use original mass
                if ind_j_LTagSumMax==ind_jetM1:
                    hist_vec_reco_1D[10].Fill(rj_m1.M(),weight)
                    hist_vec_reco_1D[11].Fill(rj_m2.M(),weight)
                    hist_vec_reco_2D[6].Fill(rj_m2.M(),rj_m1.M(),weight)
                else:
                    hist_vec_reco_1D[10].Fill(rj_m2.M(),weight)
                    hist_vec_reco_1D[11].Fill(rj_m1.M(),weight) 
                    hist_vec_reco_2D[6].Fill(rj_m1.M(),rj_m2.M(),weight)
                hist_vec_reco_1D[12].Fill(degrees(rj_m1.Theta()),weight)
                hist_vec_reco_1D[13].Fill(degrees(rj_m2.Theta()),weight)
                #up to here 13 1D filled, 7 2D
                hist_vec_reco_1D[14].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],weight)
                hist_vec_reco_1D[15].Fill(recojet_subjet_rfj_j_CTag[ind_sj_CTagMax],weight) 
                hist_vec_reco_1D[16].Fill(recojet_subjet_rfj_j_OTag[ind_sj_LTagMax],weight)
                hist_vec_reco_2D[7].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_OTag[ind_sj_LTagMax])  
                if ind_j_BTagSumMax==recojet_subjet_rfj_j_jetindex[0]:
                  hist_vec_reco_1D[17].Fill(recojet_subjet_rfj_j_BTag[0]+recojet_subjet_rfj_j_BTag[1],weight)
                  hist_vec_reco_1D[18].Fill(recojet_subjet_rfj_j_BTag[2]+recojet_subjet_rfj_j_BTag[3],weight)
                  hist_vec_reco_1D[19].Fill(recojet_subjet_rfj_j_OTag[2]+recojet_subjet_rfj_j_OTag[3],weight)
                  hist_vec_reco_1D[20].Fill(recojet_subjet_rfj_j_CTag[0]+recojet_subjet_rfj_j_CTag[1],weight)
                  hist_vec_reco_2D[8].Fill(recojet_subjet_rfj_j_BTag[0]+recojet_subjet_rfj_j_BTag[1],recojet_subjet_rfj_j_BTag[2]+recojet_subjet_rfj_j_BTag[3],weight)  
                  hist_vec_reco_2D[9].Fill(recojet_subjet_rfj_j_BTag[0]+recojet_subjet_rfj_j_BTag[1],recojet_subjet_rfj_j_OTag[2]+recojet_subjet_rfj_j_OTag[3],weight)  
                  hist_vec_reco_2D[10].Fill(recojet_subjet_rfj_j_BTag[0]+recojet_subjet_rfj_j_BTag[1],recojet_subjet_rfj_j_CTag[2]+recojet_subjet_rfj_j_CTag[3],weight)
                  #if(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[0]*recojet_subjet_rfj_j_jetChargePt_kappa_0_30[1]) <0:
                  #   if recojet_subjet_rfj_j_jetChargePt_kappa_0_30[0]>0:
                  #      hist_vec_reco_2D[11].Fill(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[0],recojet_subjet_rfj_j_jetChargePt_kappa_0_30[1],weight)
                  #   else:
                  #      hist_vec_reco_2D[11].Fill(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[1],recojet_subjet_rfj_j_jetChargePt_kappa_0_30[0],weight) 
                  #if(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[2]*recojet_subjet_rfj_j_jetChargePt_kappa_0_30[3]) <0:
                  #   if recojet_subjet_rfj_j_jetChargePt_kappa_0_30[2]>0:
                  #      hist_vec_reco_2D[12].Fill(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[2],recojet_subjet_rfj_j_jetChargePt_kappa_0_30[3],weight)
                  #   else:
                  #      hist_vec_reco_2D[12].Fill(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[3],recojet_subjet_rfj_j_jetChargePt_kappa_0_30[2],weight)   
                else: 
                  hist_vec_reco_1D[17].Fill(recojet_subjet_rfj_j_BTag[2]+recojet_subjet_rfj_j_BTag[3],weight)
                  hist_vec_reco_1D[18].Fill(recojet_subjet_rfj_j_BTag[0]+recojet_subjet_rfj_j_BTag[1],weight)
                  hist_vec_reco_1D[19].Fill(recojet_subjet_rfj_j_OTag[0]+recojet_subjet_rfj_j_OTag[1],weight)
                  hist_vec_reco_1D[20].Fill(recojet_subjet_rfj_j_CTag[0]+recojet_subjet_rfj_j_CTag[1],weight)
                  hist_vec_reco_2D[8].Fill(recojet_subjet_rfj_j_BTag[2]+recojet_subjet_rfj_j_BTag[3],recojet_subjet_rfj_j_BTag[0]+recojet_subjet_rfj_j_BTag[1],weight)
                  hist_vec_reco_2D[9].Fill(recojet_subjet_rfj_j_BTag[2]+recojet_subjet_rfj_j_BTag[3],recojet_subjet_rfj_j_OTag[0]+recojet_subjet_rfj_j_OTag[1],weight)
                  hist_vec_reco_2D[10].Fill(recojet_subjet_rfj_j_BTag[2]+recojet_subjet_rfj_j_BTag[3],recojet_subjet_rfj_j_CTag[0]+recojet_subjet_rfj_j_CTag[1],weight)
                  #if(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[2]*recojet_subjet_rfj_j_jetChargePt_kappa_0_30[3]) <0:
                  #   if recojet_subjet_rfj_j_jetChargePt_kappa_0_30[2]>0:
                  #      hist_vec_reco_2D[11].Fill(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[2],recojet_subjet_rfj_j_jetChargePt_kappa_0_30[3],weight)
                  #   else:
                  #      hist_vec_reco_2D[11].Fill(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[3],recojet_subjet_rfj_j_jetChargePt_kappa_0_30[2],weight)
                  #if(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[0]*recojet_subjet_rfj_j_jetChargePt_kappa_0_30[1]) <0:
                  #   if recojet_subjet_rfj_j_jetChargePt_kappa_0_30[0]>0:
                  #      hist_vec_reco_2D[12].Fill(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[0],recojet_subjet_rfj_j_jetChargePt_kappa_0_30[1],weight)
                  #   else:
                  #      hist_vec_reco_2D[12].Fill(recojet_subjet_rfj_j_jetChargePt_kappa_0_30[1],recojet_subjet_rfj_j_jetChargePt_kappa_0_30[0],weight) 
                #both quantities go from -180 to 180
                hist_vec_reco_1D[21].Fill(degrees(rj_m1.Theta()-rj_m2.Theta()),weight)
                hist_vec_reco_1D[22].Fill(degrees(rj_m1.DeltaPhi(rj_m2)),weight)
                hist_vec_reco_1D[127].Fill(degrees(rj_m1.Angle(rj_m2.Vect())),weight)
                #up to here BTag combinations have been checked, also max subjet indices and mass ordering according to them
                #next after checking the sum vs the sum of the opposite jet, now check BTag Sum subjet vs subjet tagging info from other fat jet
                hist_vec_reco_1D[23].Fill(recojet_subjet_rfj_j_CTag[ind_sj_BTagMax],weight)
                hist_vec_reco_1D[24].Fill(recojet_subjet_rfj_j_OTag[ind_sj_BTagMax],weight)
                hist_vec_reco_1D[25].Fill(recojet_subjet_rfj_j_BTag[ind_sj_LTagMax],weight)
                hist_vec_reco_1D[26].Fill(recojet_subjet_rfj_j_CTag[ind_sj_LTagMax],weight)
                ind_sj2_BTagMax=0
                if ind_sj_BTagMax==0:
                   ind_sj2_BTagMax=1
                elif ind_sj_BTagMax>1:
                   if ind_sj_BTagMax ==2: 
                      ind_sj2_BTagMax=3
                   else:
                      ind_sj2_BTagMax=2
                else:
                    if ind_sj_BTagMax!=1:
                      print "what is wrong, code should have been now 1 ",ind_sj_BTagMax
                ind_sj2_BTagMax=0
                if ind_sj_BTagMax==0:
                   ind_sj2_BTagMax=1
                elif ind_sj_BTagMax>1:
                   if ind_sj_BTagMax ==2: 
                      ind_sj2_BTagMax=3
                   else:
                      ind_sj2_BTagMax=2
                else:
                    if ind_sj_BTagMax!=1:
                      print "what is wrong, code should have been now 1",ind_sj_BTagMax
                #print "ind_sj2,ind_sj ",ind_sj2_BTagMax,ind_sj_BTagMax
                hist_vec_reco_1D[27].Fill(recojet_subjet_rfj_j_BTag[ind_sj2_BTagMax],weight)
                hist_vec_reco_1D[28].Fill(recojet_subjet_rfj_j_CTag[ind_sj2_BTagMax],weight)
                hist_vec_reco_1D[29].Fill(recojet_subjet_rfj_j_OTag[ind_sj2_BTagMax],weight)
                hist_vec_reco_2D[13].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_BTag[ind_sj2_BTagMax],weight)
                hist_vec_reco_2D[14].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_CTag[ind_sj2_BTagMax],weight)
                hist_vec_reco_2D[15].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_OTag[ind_sj2_BTagMax],weight)
                if ind_sj_BTagMax < 2:
                   if(recojet_subjet_rfj_j_BTag[2]>recojet_subjet_rfj_j_BTag[3]):
                      hist_vec_reco_1D[30].Fill(recojet_subjet_rfj_j_BTag[2],weight)
                      hist_vec_reco_2D[16].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_BTag[2],weight)
                      hist_vec_reco_2D[17].Fill(recojet_subjet_rfj_j_BTag[2],recojet_subjet_rfj_j_BTag[3],weight)
                   else:
                      hist_vec_reco_1D[30].Fill(recojet_subjet_rfj_j_BTag[3],weight)
                      hist_vec_reco_2D[16].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_BTag[3],weight)
                      hist_vec_reco_2D[17].Fill(recojet_subjet_rfj_j_BTag[3],recojet_subjet_rfj_j_BTag[2],weight)
                   if(recojet_subjet_rfj_j_CTag[2]>recojet_subjet_rfj_j_CTag[3]):
                      hist_vec_reco_1D[31].Fill(recojet_subjet_rfj_j_CTag[2],weight)
                      hist_vec_reco_2D[18].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_CTag[2],weight)
                   else:
                      hist_vec_reco_1D[31].Fill(recojet_subjet_rfj_j_CTag[3],weight)
                      hist_vec_reco_2D[18].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_CTag[3],weight)
                   if(recojet_subjet_rfj_j_OTag[2]>recojet_subjet_rfj_j_OTag[3]):
                      hist_vec_reco_1D[32].Fill(recojet_subjet_rfj_j_OTag[2],weight)
                      hist_vec_reco_2D[19].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_OTag[2],weight)
                   else:
                      hist_vec_reco_1D[32].Fill(recojet_subjet_rfj_j_OTag[3],weight)
                      hist_vec_reco_2D[19].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_OTag[3],weight)
                else:
                   if(recojet_subjet_rfj_j_BTag[0]>recojet_subjet_rfj_j_BTag[1]):
                      hist_vec_reco_1D[30].Fill(recojet_subjet_rfj_j_BTag[0],weight)
                      hist_vec_reco_2D[16].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_BTag[0],weight)
                      hist_vec_reco_2D[17].Fill(recojet_subjet_rfj_j_BTag[0],recojet_subjet_rfj_j_BTag[1],weight)
                   else:
                      hist_vec_reco_1D[30].Fill(recojet_subjet_rfj_j_BTag[1],weight)
                      hist_vec_reco_2D[16].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_BTag[1],weight)
                      hist_vec_reco_2D[17].Fill(recojet_subjet_rfj_j_BTag[1],recojet_subjet_rfj_j_BTag[0],weight)
                   if(recojet_subjet_rfj_j_CTag[0]>recojet_subjet_rfj_j_CTag[1]):
                      hist_vec_reco_1D[31].Fill(recojet_subjet_rfj_j_CTag[0],weight)
                      hist_vec_reco_2D[18].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_CTag[0],weight)
                   else:
                      hist_vec_reco_1D[31].Fill(recojet_subjet_rfj_j_CTag[1],weight)
                      hist_vec_reco_2D[18].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_CTag[1],weight)
                   if(recojet_subjet_rfj_j_OTag[0]>recojet_subjet_rfj_j_OTag[1]):
                      hist_vec_reco_1D[32].Fill(recojet_subjet_rfj_j_OTag[0],weight)
                      hist_vec_reco_2D[19].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_OTag[0],weight)
                   else:
                      hist_vec_reco_1D[32].Fill(recojet_subjet_rfj_j_OTag[1],weight)     
                      hist_vec_reco_2D[19].Fill(recojet_subjet_rfj_j_BTag[ind_sj_BTagMax],recojet_subjet_rfj_j_OTag[1],weight)

                veto_double_vertex=False;
                for s in range(len(recojet_subjet_rfj_j_svtx_Charge)):
                    #print "now in loop ",s
                    for s2 in range(len(recojet_subjet_rfj_j_svtx_Charge)):
                        if (s2!=s and abs(recojet_subjet_rfj_j_svtx_E[s]-recojet_subjet_rfj_j_svtx_E[s2])<1.e-6):
                            veto_double_vertex=True;
                            print "seems vertices are the same ",s2,s,recojet_subjet_rfj_j_svtx_E[s],recojet_subjet_rfj_j_svtx_E[s2],recojet_subjet_rfj_j_svtx_E[s]-recojet_subjet_rfj_j_svtx_E[s2], len(recojet_subjet_rfj_j_svtx_Charge),num_entry
                            break

                    if veto_double_vertex:
                        #print "we should break now ",s,num_entry
                        break
    print 'total events after all running BG histos ', hist_vec_reco_1D[0].Integral(0,hist_vec_reco_1D[0].GetNbinsX()+1),"masscuts/thetacuts?",performMassCuts,performThetaCuts,num_count
    return None


    

def process_event(i_final_histo_name_,i_input_file_name_,i_xsec_,i_lumi_, i_use_partonInfo_,i_fill_signal_histos_,i_bool_applyMassCuts_,i_bool_useRectMassCuts_,i_bool_applyThetaCuts_,i_performBTagCuts_,i_performC2Cuts_,i_performD2Cuts_):
    print "at start of process event"         
    input_file_=root.TFile.Open(i_input_file_name_)
    #input_file2_=root.TFile.Open(i_input_file_name2_)
    lumi=i_lumi_
    xsec_=i_xsec_
    use_partonInfo=i_use_partonInfo_
    bool_performMassCut=i_bool_applyMassCuts_
    bool_useRectangularMassCut=i_bool_useRectMassCuts_
    bool_performthetaCut=i_bool_applyThetaCuts_
    bool_performBTagCut=i_performBTagCuts_
    bool_performC2Cut=i_performC2Cuts_
    bool_performD2Cut=i_performD2Cuts_
    fill_signal_histos=i_fill_signal_histos_
    print 'process event, parton/signal/mass/massrect/theta/BTagCut/C2/D2 cut ',use_partonInfo,fill_signal_histos,bool_performMassCut,bool_useRectangularMassCut,bool_performthetaCut,bool_performBTagCut,bool_performC2Cut,bool_performD2Cut
          
                            
    file_histogram = root.TFile(i_final_histo_name_, "RECREATE")





    #limits are for the high energy selection
    n_bins_high=100
    lim_mass_low=0
    lim_mass_high=200

    lim_mass_delta_mass_low=0.
    lim_mass_delta_mass_high=73.5

    lim_Tagging_low=-0.05
    lim_Tagging_high=1.05

    n_bins_svtx_Charge=23
    lim_svtx_Charge_low=-11.5
    lim_svtx_Charge_high=11.5

    n_bins_svtx_nTrk=31
    lim_svtx_nTrk_low=-0.5
    lim_svtx_nTrk_high=30.5

    lim_subjet_Charge_low=-2.00
    lim_subjet_Charge_high=2.00

    lim_energy_low=-200;
    lim_energy_high=3500.;

    lim_dangle_low=-0.01
    lim_dangle_sj_high=75.
    lim_dangle_high=30.0
    lim_dangle_max=180.0

    lim_theta_min=0.0
    lim_theta_max=180.0

    lim_deltaE_low=-0.10
    lim_deltaE_high=0.10

    lim_vtx_mass_low=0
    lim_vtx_mass_high=10

    lim_BTag_low=0.
    lim_BTag_high=1.

    lim_energy_jet_low=400;
    lim_energy_jet_high=1750.

    lim_energy_rat_j2_low=0.
    lim_energy_rat_j1_low=0.5
    lim_energy_rat_high=1.

    lim_energy_jet_rel_low=-0.5
    lim_energy_jet_rel_high=0.5

    n_bins_low=50
    lim_energy_low_mediumSqrtS=1000.
    lim_energy_low_highSqrtS=2500.

    lim_energy_sqrtS_rel_low=-0.5
    lim_energy_sqrtS_rel_high=0.5

    lim_tag_low=0.0
    lim_tag_high=1.0


    lim_N2_low=0.0
    lim_b1_N2_high=0.50
    lim_b2_N2_high=0.40
    lim_b0_5_N2_high=0.55

    lim_N3_low=0.0
    lim_b1_N3_high=5.0
    lim_b2_N3_high=7.5
    lim_b0_5_N3_high=3.00

    lim_C2_low=0.0
    lim_b1_C2_high=0.55
    lim_b2_C2_high=0.75
    lim_b0_5_C2_high=0.50

    lim_C3_low=0.0
    lim_b1_C3_high=0.75
    lim_b2_C3_high=0.85
    lim_b0_5_C3_high=0.60

    lim_D2_low=0.0
    lim_b1_D2_high=10.0
    lim_b2_D2_high=55.0
    lim_b0_5_D2_high=4.0

    lim_tau21_low=0.0
    lim_tau21_high=1.2

    lim_tau32_low=0.0
    lim_tau32_high=1.2

    lim_dij_low=0.0
    lim_dij21_rat_high=0.06
    lim_dij32_rat_high=0.005
    lim_dij43_rat_high=0.001
    
    lim_dij32_over_dij21_rat_high=1.00
    lim_dij43_over_dij21_rat_high=0.30
    lim_dij43_over_dij32_rat_high=1.20

    if fill_signal_histos:
       h_HZqq_signal_mass_H_matched = TH1F( "h_HZqq_signal_mass_H_matched", "", n_bins_high, lim_mass_low, lim_mass_high )
       h_HZqq_signal_mass_Z_matched = TH1F( "h_HZqq_signal_mass_Z_matched", "", n_bins_high, lim_mass_low, lim_mass_high )

       h_HZqq_signal_Btag_b_H_matched = TH1F( "h_HZqq_signal_Btag_b_H_matched", "", n_bins_high, lim_Tagging_low, lim_Tagging_high )
       h_HZqq_signal_Btag_bbar_H_matched = TH1F( "h_HZqq_signal_Btag_bbar_matched", "", n_bins_high, lim_Tagging_low, lim_Tagging_high )

       h_HZqq_signal_svtx_Charge_H_matched = TH1F( "h_HZqq_signal_svtx_Charge_H_matched", "", n_bins_svtx_Charge, lim_svtx_Charge_low, lim_svtx_Charge_high )
       h_HZqq_signal_svtx_Charge_Z_matched = TH1F( "h_HZqq_signal_svtx_Charge_Z_matched", "", n_bins_svtx_Charge, lim_svtx_Charge_low, lim_svtx_Charge_high )
       h_HZqq_signal_svtx_nTrk_H_matched = TH1F( "h_HZqq_signal_svtx_nTrk_H_matched", "", n_bins_svtx_nTrk, lim_svtx_nTrk_low, lim_svtx_nTrk_high )
       h_HZqq_signal_svtx_nTrk_Z_matched = TH1F( "h_HZqq_signal_svtx_nTrk_Z_matched", "", n_bins_svtx_nTrk, lim_svtx_nTrk_low, lim_svtx_nTrk_high )

       h_HZqq_signal_Angle_sj_rf_sj_H_matched = TH1F( "h_HZqq_signal_Angle_sj_rf_sj_H_matched", "", n_bins_high, lim_dangle_low, lim_dangle_high )
       h_HZqq_signal_Angle_sj_rf_sj_Z_matched = TH1F( "h_HZqq_signal_Angle_sj_rf_sj_Z_matched", "", n_bins_high, lim_dangle_low, lim_dangle_high )
       h_HZqq_signal_deltaE_rel_sj_rf_sj_H_matched = TH1F( "h_HZqq_signal_deltaE_rel_sj_rf_sj_H_matched", "", n_bins_high, lim_deltaE_low, lim_deltaE_high )
       h_HZqq_signal_deltaE_rel_sj_rf_sj_Z_matched = TH1F( "h_HZqq_signal_deltaE_rel_sj_rf_sj_Z_matched", "", n_bins_high, lim_deltaE_low, lim_deltaE_high )

       h_HZqq_signal_Angle_rf_or_j_H_matched = TH1F( "h_HZqq_signal_Angle_rf_or_j_H_matched", "", 6*n_bins_high, lim_dangle_low, lim_dangle_max )
       h_HZqq_signal_Angle_rf_or_j_Z_matched = TH1F( "h_HZqq_signal_Angle_rf_or_j_Z_matched", "", 6*n_bins_high, lim_dangle_low, lim_dangle_max )
       h_HZqq_signal_deltaE_rel_rf_or_j_H_matched = TH1F( "h_HZqq_signal_deltaE_rel_rf_or_j_H_matched", "", n_bins_high, lim_deltaE_low, lim_deltaE_high )
       h_HZqq_signal_deltaE_rel_rf_or_j_Z_matched = TH1F( "h_HZqq_signal_deltaE_rel_rf_or_j_Z_matched", "", n_bins_high, lim_deltaE_low, lim_deltaE_high )


       h_HZqq_signal_BTag_rfj_H = TH1F( "h_HZqq_signal_BTag_rfj_H", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_CTag_rfj_H = TH1F( "h_HZqq_signal_CTag_rfj_H", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_OTag_rfj_H = TH1F( "h_HZqq_signal_OTag_rfj_H", "", n_bins_high, lim_BTag_low,lim_BTag_high)

       h_HZqq_signal_BTag_rfj_Z_b = TH1F( "h_HZqq_signal_BTag_rfj_Z_b", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_CTag_rfj_Z_b = TH1F( "h_HZqq_signal_CTag_rfj_Z_b", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_OTag_rfj_Z_b = TH1F( "h_HZqq_signal_OTag_rfj_Z_b", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_BTag_rfj_Z_c = TH1F( "h_HZqq_signal_BTag_rfj_Z_c", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_CTag_rfj_Z_c = TH1F( "h_HZqq_signal_CTag_rfj_Z_c", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_OTag_rfj_Z_c = TH1F( "h_HZqq_signal_OTag_rfj_Z_c", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_BTag_rfj_Z_l = TH1F( "h_HZqq_signal_BTag_rfj_Z_l", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_CTag_rfj_Z_l = TH1F( "h_HZqq_signal_CTag_rfj_Z_l", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_OTag_rfj_Z_l = TH1F( "h_HZqq_signal_OTag_rfj_Z_l", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)


       h_HZqq_signal_BTag_max_rfj_sj_H = TH1F( "h_HZqq_signal_BTag_max_rfj_sj_H", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_BTag_min_rfj_sj_H = TH1F( "h_HZqq_signal_BTag_min_rfj_sj_H", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_BTag_sum_rfj_sj_H = TH1F( "h_HZqq_signal_BTag_sum_rfj_sj_H", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)

       h_HZqq_signal_CTag_sum_rfj_sj_H = TH1F( "h_HZqq_signal_CTag_sum_rfj_sj_H", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)
       h_HZqq_signal_OTag_sum_rfj_sj_H = TH1F( "h_HZqq_signal_OTag_sum_rfj_sj_H", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)

       h_HZqq_signal_BTag_max_rfj_sj_Z_b = TH1F( "h_HZqq_signal_BTag_max_rfj_sj_Z_b", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_BTag_min_rfj_sj_Z_b = TH1F( "h_HZqq_signal_BTag_min_rfj_sj_Z_b", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_BTag_sum_rfj_sj_Z_b = TH1F( "h_HZqq_signal_BTag_sum_rfj_sj_Z_b", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)
       h_HZqq_signal_CTag_sum_rfj_sj_Z_b = TH1F( "h_HZqq_signal_CTag_sum_rfj_sj_Z_b", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)
       h_HZqq_signal_OTag_sum_rfj_sj_Z_b = TH1F( "h_HZqq_signal_OTag_sum_rfj_sj_Z_b", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)

       h_HZqq_signal_BTag_max_rfj_sj_Z_c = TH1F( "h_HZqq_signal_BTag_max_rfj_sj_Z_c", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_BTag_min_rfj_sj_Z_c = TH1F( "h_HZqq_signal_BTag_min_rfj_sj_Z_c", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_BTag_sum_rfj_sj_Z_c = TH1F( "h_HZqq_signal_BTag_sum_rfj_sj_Z_c", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)
       h_HZqq_signal_CTag_sum_rfj_sj_Z_c = TH1F( "h_HZqq_signal_CTag_sum_rfj_sj_Z_c", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)
       h_HZqq_signal_OTag_sum_rfj_sj_Z_c = TH1F( "h_HZqq_signal_OTag_sum_rfj_sj_Z_c", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)

       h_HZqq_signal_BTag_max_rfj_sj_Z_l = TH1F( "h_HZqq_signal_BTag_max_rfj_sj_Z_l", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_BTag_min_rfj_sj_Z_l = TH1F( "h_HZqq_signal_BTag_min_rfj_sj_Z_l", "", n_bins_high, lim_BTag_low,lim_BTag_high)
       h_HZqq_signal_BTag_sum_rfj_sj_Z_l = TH1F( "h_HZqq_signal_BTag_sum_rfj_sj_Z_l", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)
       h_HZqq_signal_CTag_sum_rfj_sj_Z_l = TH1F( "h_HZqq_signal_CTag_sum_rfj_sj_Z_l", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)
       h_HZqq_signal_OTag_sum_rfj_sj_Z_l = TH1F( "h_HZqq_signal_OTag_sum_rfj_sj_Z_l", "", n_bins_high, lim_BTag_low,2.*lim_BTag_high)

       h_HZqq_signal_Angle_rfj1_rfj2 = TH1F( "h_HZqq_signal_Angle_rfj1_rfj2", "", 60*n_bins_high, lim_dangle_low, lim_dangle_max )
       h_HZqq_signal_Angle_j1_j2 = TH1F( "h_HZqq_signal_Angle_j2_j2", "", 60*n_bins_high, lim_dangle_low, lim_dangle_max )

       h_HZqq_signal_o_mass_H_matched = TH1F( "h_HZqq_signal_o_mass_H_matched", "", n_bins_high, lim_mass_low, lim_mass_high )
       h_HZqq_signal_o_mass_Z_matched = TH1F( "h_HZqq_signal_o_mass_Z_matched", "", n_bins_high, lim_mass_low, lim_mass_high )

       h_HZqq_signal_rf_mass_H_matched = TH1F( "h_HZqq_signal_rf_mass_H_matched", "", n_bins_high, lim_mass_low, lim_mass_high )
       h_HZqq_signal_rf_mass_Z_matched = TH1F( "h_HZqq_signal_rf_mass_Z_matched", "", n_bins_high, lim_mass_low, lim_mass_high )


       h_HZqq_signal_o_E_H_matched = TH1F( "h_HZqq_signal_o_E_H_matched", "", n_bins_high, lim_energy_jet_low, lim_energy_jet_high )
       h_HZqq_signal_o_E_Z_matched = TH1F( "h_HZqq_signal_o_E_Z_matched", "", n_bins_high, lim_energy_jet_low, lim_energy_jet_high )

       h_HZqq_signal_corr_E_H_matched = TH1F( "h_HZqq_signal_corr_E_H_matched", "", n_bins_high, lim_energy_jet_low, lim_energy_jet_high )
       h_HZqq_signal_corr_E_Z_matched = TH1F( "h_HZqq_signal_corr_E_Z_matched", "", n_bins_high, lim_energy_jet_low, lim_energy_jet_high )

       h_HZqq_signal_o_deltaE_H_matched = TH1F( "h_HZqq_signal_o_deltaE_H_matched", "", n_bins_high, lim_energy_jet_rel_low, lim_energy_jet_rel_high )
       h_HZqq_signal_o_deltaE_Z_matched = TH1F( "h_HZqq_signal_o_deltaE_Z_matched", "", n_bins_high, lim_energy_jet_rel_low, lim_energy_jet_rel_high )

       h_HZqq_signal_corr_deltaE_H_matched = TH1F( "h_HZqq_signal_corr_deltaE_H_matched", "", n_bins_high, lim_energy_jet_rel_low, lim_energy_jet_rel_high )
       h_HZqq_signal_corr_deltaE_Z_matched = TH1F( "h_HZqq_signal_corr_deltaE_Z_matched", "", n_bins_high, lim_energy_jet_rel_low, lim_energy_jet_rel_high )


       h_HZqq_signal_Angle_jet1_H = TH1F( "h_HZqq_signal_Angle_jet1_H", "", 6*n_bins_high, lim_dangle_low, lim_dangle_max )
       h_HZqq_signal_Angle_jet2_Z = TH1F( "h_HZqq_signal_Angle_jet2_Z", "", 6*n_bins_high, lim_dangle_low, lim_dangle_max )

       h_HZqq_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj = TH1F("h_HZqq_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj","", n_bins_high, lim_subjet_Charge_low,lim_subjet_Charge_high);
       h_HZqq_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj.GetXaxis().SetTitle("subjet-charge_E (q-(Z)-matched)");
       h_HZqq_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj = TH1F("h_HZqq_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj","", n_bins_high, lim_subjet_Charge_low,lim_subjet_Charge_high);
       h_HZqq_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj.GetXaxis().SetTitle("subjet-charge_E (q-(Z)-matched)");
       h_HZqq_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj = TH1F("h_HZqq_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj","", n_bins_high, lim_subjet_Charge_low,lim_subjet_Charge_high);
       h_HZqq_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj.GetXaxis().SetTitle("subjet-charge_E (q-(Z)-matched)");
       h_HZqq_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj = TH1F("h_HZqq_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj","", n_bins_high, lim_subjet_Charge_low,lim_subjet_Charge_high);
       h_HZqq_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj.GetXaxis().SetTitle("subjet-charge_E (q-(Z)-matched)");
       h_HZqq_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj = TH1F("h_HZqq_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj","", n_bins_high, lim_subjet_Charge_low,lim_subjet_Charge_high);
       h_HZqq_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj.GetXaxis().SetTitle("subjet-charge_E (q-(Z)-matched)");

       h_HZqq_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj = TH1F("h_HZqq_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj","", n_bins_high, lim_subjet_Charge_low,lim_subjet_Charge_high);
       h_HZqq_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj.GetXaxis().SetTitle("subjet-charge_E (q-(Z)-matched)");
       h_HZqq_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj = TH1F("h_HZqq_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj","", n_bins_high, lim_subjet_Charge_low,lim_subjet_Charge_high);
       h_HZqq_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj.GetXaxis().SetTitle("subjet-charge_E (q-(Z)-matched)");
       h_HZqq_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj = TH1F("h_HZqq_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj","", n_bins_high, lim_subjet_Charge_low,lim_subjet_Charge_high);
       h_HZqq_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj.GetXaxis().SetTitle("subjet-charge_E (q-(Z)-matched)");
       h_HZqq_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj = TH1F("h_HZqq_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj","", n_bins_high, lim_subjet_Charge_low,lim_subjet_Charge_high);
       h_HZqq_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj.GetXaxis().SetTitle("subjet-charge_E (q-(Z)-matched)");
       h_HZqq_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj = TH1F("h_HZqq_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj","", n_bins_high, lim_subjet_Charge_low,lim_subjet_Charge_high);
       h_HZqq_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj.GetXaxis().SetTitle("subjet-charge_E (q-(Z)-matched)");

       h_HZqq_signal_1D_hist_list=[]
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_mass_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_mass_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_Btag_b_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_Btag_bbar_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_svtx_Charge_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_svtx_Charge_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_svtx_nTrk_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_svtx_nTrk_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_Angle_sj_rf_sj_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_Angle_sj_rf_sj_Z_matched)
       #[9] histos by then
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_deltaE_rel_sj_rf_sj_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_deltaE_rel_sj_rf_sj_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_Angle_rf_or_j_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_Angle_rf_or_j_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_deltaE_rel_rf_or_j_H_matched)
       #14 histos by then
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_deltaE_rel_rf_or_j_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_rfj_H)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_CTag_rfj_H)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_OTag_rfj_H)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_rfj_Z_b)
       #index [19] histos done by now
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_CTag_rfj_Z_b)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_OTag_rfj_Z_b)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_rfj_Z_c)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_CTag_rfj_Z_c)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_OTag_rfj_Z_c)
       #24 histos by now
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_rfj_Z_l)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_CTag_rfj_Z_l)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_OTag_rfj_Z_l)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_max_rfj_sj_H)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_min_rfj_sj_H)
       #29 histos by now
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_sum_rfj_sj_H)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_CTag_sum_rfj_sj_H)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_OTag_sum_rfj_sj_H)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_max_rfj_sj_Z_b)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_min_rfj_sj_Z_b)
       #34 histos by now
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_sum_rfj_sj_Z_b)    
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_CTag_sum_rfj_sj_Z_b)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_OTag_sum_rfj_sj_Z_b)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_max_rfj_sj_Z_c)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_min_rfj_sj_Z_c)
       #39 histos by now
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_sum_rfj_sj_Z_c)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_CTag_sum_rfj_sj_Z_c)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_OTag_sum_rfj_sj_Z_c)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_max_rfj_sj_Z_l)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_min_rfj_sj_Z_l)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_BTag_sum_rfj_sj_Z_l)
       #index [45] histos by now
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_CTag_sum_rfj_sj_Z_l)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_OTag_sum_rfj_sj_Z_l)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_Angle_rfj1_rfj2)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_Angle_j1_j2)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_o_mass_H_matched)
       #index [50] histos done by now
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_o_mass_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_rf_mass_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_rf_mass_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_o_E_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_o_E_Z_matched)
       #55 histos by now
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_corr_E_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_corr_E_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_o_deltaE_H_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_o_deltaE_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_corr_deltaE_H_matched)
       #60 histos by now
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_corr_deltaE_Z_matched)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_Angle_jet1_H)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_signal_Angle_jet2_Z)
       #63 hists by now
       h_HZqq_signal_1D_hist_list.append(h_HZqq_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_pos_gj)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_pos_gj)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_pos_gj)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_pos_gj)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_pos_gj)
       
       h_HZqq_signal_1D_hist_list.append(h_HZqq_subjet_charge_E_kappa_0_15_gen_sj_matched_Z_q_neg_gj)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_subjet_charge_E_kappa_0_20_gen_sj_matched_Z_q_neg_gj)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_subjet_charge_E_kappa_0_25_gen_sj_matched_Z_q_neg_gj)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_subjet_charge_E_kappa_0_30_gen_sj_matched_Z_q_neg_gj)
       h_HZqq_signal_1D_hist_list.append(h_HZqq_subjet_charge_E_kappa_0_50_gen_sj_matched_Z_q_neg_gj)
       #73 hists by now

       h_HZqq_signal_1D_sqrtS_part = TH1F( "h_HZqq_signal_1D_sqrtS_part", "", n_bins_high, lim_energy_low,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_genTot_isoPh = TH1F( "h_HZqq_signal_1D_sqrtS_genTot_isoPh", "", n_bins_high, lim_energy_low,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_genTot_isoPh_EMiss = TH1F( "h_HZqq_signal_1D_sqrtS_genTot_isoPh_EMiss", "", n_bins_high, lim_energy_low,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_genTot_isoPh_TrueEMiss = TH1F( "h_HZqq_signal_1D_sqrtS_genTot_isoPh_TrueEMiss", "", n_bins_high, lim_energy_low,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_genTot_isoPh_inv = TH1F( "h_HZqq_signal_1D_sqrtS_genTot_isoPh_inv", "", n_bins_high, lim_energy_low,lim_energy_high)

       h_HZqq_signal_1D_sqrtS_gen_j1_j2 = TH1F( "h_HZqq_signal_1D_sqrtS_gen_j1_j2", "", n_bins_high, lim_energy_low,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_gen_j1_j2_EMiss = TH1F( "h_HZqq_signal_1D_sqrtS_gen_j1_j2_EMiss", "", n_bins_high, lim_energy_low,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_gen_j1_j2_TrueEMiss = TH1F( "h_HZqq_signal_1D_sqrtS_gen_j1_j2_TrueEMiss", "", n_bins_high, lim_energy_low,lim_energy_high)

       h_HZqq_signal_1D_sqrtS_recoTot_isoPh = TH1F( "h_HZqq_signal_1D_sqrtS_recoTot_isoPh", "", n_bins_high, lim_energy_low,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss = TH1F( "h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss", "", n_bins_high, lim_energy_low,lim_energy_high)

       h_HZqq_signal_1D_sqrtS_reco_j1_j2 = TH1F( "h_HZqq_signal_1D_sqrtS_reco_j1_j2", "", n_bins_high, lim_energy_low,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss = TH1F( "h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss", "", n_bins_high, lim_energy_low,lim_energy_high)

       h_HZqq_signal_1D_sqrtS_recoTot_isoPh_high_sqrtS = TH1F( "h_HZqq_signal_1D_sqrtS_recoTot_isoPh_high_sqrtS", "", n_bins_low, lim_energy_low_highSqrtS,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss_high_sqrtS = TH1F( "h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss_high_sqrtS", "", n_bins_low, lim_energy_low_highSqrtS,lim_energy_high)

       h_HZqq_signal_1D_sqrtS_reco_j1_j2_high_sqrtS = TH1F( "h_HZqq_signal_1D_sqrtS_reco_j1_j2_high_sqrtS", "", n_bins_low, lim_energy_low_highSqrtS,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss_high_sqrtS = TH1F( "h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss_high_sqrtS", "", n_bins_low, lim_energy_low_highSqrtS,lim_energy_high)

       h_HZqq_signal_1D_sqrtS_recoTot_isoPh_high_sqrtS_part2500 = TH1F( "h_HZqq_signal_1D_sqrtS_recoTot_isoPh_high_sqrtS_part2500", "", n_bins_low, lim_energy_low_highSqrtS,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss_high_sqrtS_part2500 = TH1F( "h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss_high_sqrtS_part2500", "", n_bins_low, lim_energy_low_highSqrtS,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_reco_j1_j2_high_sqrtS_part2500 = TH1F( "h_HZqq_signal_1D_sqrtS_reco_j1_j2_high_sqrtS_part2500", "", n_bins_low, lim_energy_low_highSqrtS,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss_high_sqrtS_part2500 = TH1F( "h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss_high_sqrtS_part2500", "", n_bins_low, lim_energy_low_highSqrtS,lim_energy_high)

       h_HZqq_signal_1D_sqrtS_recoTot_isoPh_part2500 = TH1F( "h_HZqq_signal_1D_sqrtS_recoTot_isoPh_part2500", "", n_bins_high, lim_energy_low,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss_part2500 = TH1F( "h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss_part2500", "", n_bins_high, lim_energy_low,lim_energy_high)

       h_HZqq_signal_1D_sqrtS_reco_j1_j2_part2500 = TH1F( "h_HZqq_signal_1D_sqrtS_reco_j1_j2_part2500", "", n_bins_high, lim_energy_low,lim_energy_high)
       h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss_part2500 = TH1F( "h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss_part2500", "", n_bins_high, lim_energy_low,lim_energy_high)

       h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_part2500 = TH1F( "h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_part2500", "", n_bins_high, lim_energy_sqrtS_rel_low,lim_energy_sqrtS_rel_high)
       h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_EMiss_part2500 = TH1F( "h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_EMiss_part2500", "", n_bins_high, lim_energy_sqrtS_rel_low,lim_energy_sqrtS_rel_high)

       h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500 = TH1F( "h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500", "", n_bins_high, lim_energy_sqrtS_rel_low,lim_energy_sqrtS_rel_high)
       h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500 = TH1F( "h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500", "", n_bins_high, lim_energy_sqrtS_rel_low,lim_energy_sqrtS_rel_high)

       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco=[]
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_part)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_genTot_isoPh)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_genTot_isoPh_EMiss)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_genTot_isoPh_TrueEMiss)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_genTot_isoPh_inv)
       #index [4] histo up to here
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_gen_j1_j2)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_gen_j1_j2_EMiss)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_gen_j1_j2_TrueEMiss)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_recoTot_isoPh)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss)
       #9 histos up to here
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_reco_j1_j2)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_recoTot_isoPh_high_sqrtS)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss_high_sqrtS)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_reco_j1_j2_high_sqrtS)
       #14 histos up to here
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss_high_sqrtS)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_recoTot_isoPh_high_sqrtS_part2500)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss_high_sqrtS_part2500)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_reco_j1_j2_high_sqrtS_part2500)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss_high_sqrtS_part2500)
       #19 histos up to here
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_recoTot_isoPh_part2500)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_recoTot_isoPh_EMiss_part2500)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_reco_j1_j2_part2500)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_sqrtS_reco_j1_j2_EMiss_part2500)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_part2500)
       #24 histos up to here
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_delta_sqrtS_rel_recoTot_isoPh_EMiss_part2500)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_part2500)
       h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco.append(h_HZqq_signal_1D_delta_sqrtS_rel_reco_j1_j2_EMiss_part2500)


       for hist in h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco:
          hist.Sumw2()

       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )
       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh_EMiss = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh_EMiss", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )
       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh_TrueEMiss = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh_TrueEMiss", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )
       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh_inv = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh_inv", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )

       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_gen_j1_j2 = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_gen_j1_j2", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )
       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_gen_j1_j2_EMiss = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_gen_j1_j2_EMiss", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )
       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_gen_j1_j2_TrueEMiss = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_gen_j1_j2_TrueEMiss", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )

       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_recoTot_isoPh = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_recoTot_isoPh", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )
       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_recoTot_isoPh_EMiss = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_recoTot_isoPh_EMiss", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )

       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_reco_j1_j2 = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_reco_j1_j2", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )
       h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_reco_j1_j2_EMiss = TH2F( "h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_reco_j1_j2_EMiss", "", n_bins_high, lim_energy_low,lim_energy_high, n_bins_high, lim_energy_low,lim_energy_high )


       h_HZqq_signal_2D_hist_list_part_vs_gen_reco=[]
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh)
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh_EMiss)
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh_TrueEMiss)
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_genTot_isoPh_inv)
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_gen_j1_j2)
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_gen_j1_j2_EMiss)
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_gen_j1_j2_TrueEMiss)
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_recoTot_isoPh)
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_recoTot_isoPh_EMiss)
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_reco_j1_j2)
       h_HZqq_signal_2D_hist_list_part_vs_gen_reco.append(h_HZqq_signal_2D_sqrtS_part_vs_sqrtS_reco_j1_j2_EMiss)

       for hist in h_HZqq_signal_2D_hist_list_part_vs_gen_reco:
          hist.Sumw2()

       for hist in h_HZqq_signal_1D_hist_list:
          hist.Sumw2()
          hist.SetLineWidth(2)
          hist.SetLineColor(root.kBlack)
          hist.GetYaxis().SetTitle("Events")

       fill_HZ_histograms(input_file_,xsec_,use_partonInfo,h_HZqq_signal_1D_hist_list,h_HZqq_signal_1D_hist_list_sqrtS_part_gen_reco,h_HZqq_signal_2D_hist_list_part_vs_gen_reco,lumi,bool_performMassCut,bool_useRectangularMassCut,bool_performthetaCut,bool_performBTagCut,bool_performC2Cut,bool_performD2Cut)

    h_mass_jet1 = TH1F( "h_mass_jet1", "", n_bins_high, lim_mass_low, lim_mass_high )
    h_mass_jet2 = TH1F( "h_mass_jet2", "", n_bins_high, lim_mass_low, lim_mass_high )
    h_mass_rfjet1 = TH1F( "h_mass_rfjet1", "", n_bins_high, lim_mass_low, lim_mass_high )
    h_mass_rfjet2 = TH1F( "h_mass_rfjet2", "", n_bins_high, lim_mass_low, lim_mass_high )
    h_mass_jet_sj_BTagMax = TH1F( "h_mass_jet_sj_BTagMax", "", n_bins_high, lim_mass_low, lim_mass_high )
    h_mass_jet_not_sj_BTagMax = TH1F( "h_mass_jet_not_sj_BTagMax", "", n_bins_high, lim_mass_low, lim_mass_high )
    #index 5 histos done
    h_mass_jet_sj_LTagMax = TH1F( "h_mass_jet_sj_LTagMax", "", n_bins_high, lim_mass_low, lim_mass_high )
    h_mass_jet_not_sj_LTagMax = TH1F( "h_mass_jet_not_sj_LTagMax", "", n_bins_high, lim_mass_low, lim_mass_high )
    h_mass_jet_BTagSumMax = TH1F( "h_mass_jet_BTagSumMax", "", n_bins_high, lim_mass_low, lim_mass_high )
    h_mass_jet_not_BTagSumMax = TH1F( "h_mass_jet_not_BTagSumMax", "", n_bins_high, lim_mass_low, lim_mass_high )
    h_mass_jet_LTagSumMax = TH1F( "h_mass_jet_LTagSumMax", "", n_bins_high, lim_mass_low, lim_mass_high )
    #index 10 histos done
    h_mass_jet_not_LTagSumMax = TH1F( "h_mass_jet_not_LTagSumMax", "", n_bins_high, lim_mass_low, lim_mass_high )
    h_theta_jet1 = TH1F( "h_theta_jet1", "", n_bins_high, lim_theta_min, lim_theta_max)
    h_theta_jet2 = TH1F( "h_theta_jet2", "", n_bins_high, lim_theta_min, lim_theta_max)
    h_BTag_sj_BTagMax = TH1F( "h_BTag_sj_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    h_CTag_sj_CTagMax = TH1F( "h_CTag_sj_CTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    #index 15 histos done
    h_LTag_sj_LTagMax = TH1F( "h_LTag_sj_LTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    h_BTagSum_jet_BTagSumMax = TH1F( "h_BTagSum_jet_BTagSumMax", "", n_bins_high, lim_tag_low, 2.*lim_tag_high )
    h_BTagSum_jet_BTagSumMin = TH1F( "h_BTagSum_jet_BTagSumMin", "", n_bins_high, lim_tag_low, 2.*lim_tag_high )
    h_CTagSum_jet_BTagSumMin = TH1F( "h_CTagSum_jet_BTagSumMin", "", n_bins_high, lim_tag_low, 2.*lim_tag_high )
    h_OTagSum_jet_BTagSumMin = TH1F( "h_LTagSum_jet_BTagSumMin", "", n_bins_high, lim_tag_low, 2.*lim_tag_high )
    #index 20 histos done
    h_theta_jet_m1_min_jet_m2 = TH1F( "h_theta_jet_m1_min_jet_m2", "", n_bins_high, -lim_dangle_max,lim_dangle_max)
    h_dphi_jet_m1_jet_m2 = TH1F( "h_dphi_jet_m1_jet_m2", "", n_bins_high, -lim_dangle_max,lim_dangle_max)
    h_CTag_sj_BTagMax = TH1F( "h_CTag_sj_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    h_LTag_sj_BTagMax = TH1F( "h_LTag_sj_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    h_BTag_sj_LTagMax = TH1F( "h_BTag_sj_LTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    #index 25 histos done
    h_CTag_sj_LTagMax = TH1F( "h_CTag_sj_LTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    h_BTag_sj2_BTagMax = TH1F( "h_BTag_sj2_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    h_CTag_sj2_BTagMax = TH1F( "h_CTag_sj2_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    h_LTag_sj2_BTagMax = TH1F( "h_LTag_sj2_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    h_BTag_sj1_j2_sj_BTagMax = TH1F( "h_BTag_sj1_BTag_j2_sj_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    #index 30 histos done
    h_CTag_sj1_BTag_j2_sj_BTagMax = TH1F( "h_CTag_sj1_BTag_j2_sj_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    h_LTag_sj1_BTag_j2_sj_BTagMax = TH1F( "h_LTag_sj1_BTag_j2_sj_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high )
    h_BTagMax_j1_mass_noBTagCut = TH1F( "h_BTagMax_j1_mass_noBTagCut", "", n_bins_high, lim_tag_low, lim_tag_high )   
    h_BTag_rfj1_E_j2_mass_noBTagCut = TH1F( "h_BTag_rfj1_E_j2_mass_noBTagCut", "", n_bins_high, lim_tag_low, lim_tag_high )   
    h_CTag_rfj1_E_j2_mass_noBTagCut = TH1F( "h_CTag_rfj1_E_j2_mass_noBTagCut", "", n_bins_high, lim_tag_low, lim_tag_high )    
    #index 35 done
    h_LTag_rfj1_E_j2_mass_noBTagCut = TH1F( "h_LTag_rfj1_E_j2_mass_noBTagCut", "", n_bins_high,  lim_tag_low, lim_tag_high )  
    h_BTagMax_j2_mass_noBTagCut = TH1F( "h_BTagMax_j2_mass_noBTagCut", "", n_bins_high, lim_tag_low, lim_tag_high )   
    h_E_rfj_BTagMax_over_E_rfj_comb_j1_mass = TH1F( "h_E_rfj_BTagMax_over_E_rfj_comb_j1_mass", "", n_bins_high, lim_energy_rat_j2_low, lim_energy_rat_high)
    h_E_rfsj1_over_E_rfj_j1_mass = TH1F( "h_E_rfsj1_over_E_rfj_j1_mass", "", n_bins_high, lim_energy_rat_j1_low, lim_energy_rat_high)   
    h_E_sj1_over_E_j1_mass = TH1F( "h_E_sj1_over_E_j1_mass", "", n_bins_high, lim_energy_rat_j1_low, lim_energy_rat_high)    
    #index 40 done
    h_E_rfsj1_over_E_rfj_j2_mass = TH1F( "h_E_rfsj1_over_E_rfj_j2_mass", "", n_bins_high, lim_energy_rat_j1_low, lim_energy_rat_high)   
    h_E_sj1_over_E_j2_mass = TH1F( "h_E_sj1_over_E_j2_mass", "", n_bins_high, lim_energy_rat_j1_low, lim_energy_rat_high)   
    h_Angle_sj1_sj2_j1_mass = TH1F( "h_Angle_sj1_sj2_j1_mass", "", n_bins_high, lim_dangle_low,  lim_dangle_sj_high)
    h_Angle_srfj1_srfj2_j1_mass = TH1F( "h_Angle_srfj1_srfj2_j1_mass", "", n_bins_high, lim_dangle_low,  lim_dangle_sj_high)   
    h_Angle_sj1_sj2_j2_mass = TH1F( "h_Angle_sj1_sj2_j2_mass", "", n_bins_high, lim_dangle_low,  lim_dangle_sj_high)   
    #index 45 done
    h_Angle_srfj1_srfj2_j2_mass = TH1F( "h_Angle_srfj1_srfj2_j2_mass", "", n_bins_high, lim_dangle_low,  lim_dangle_sj_high)

    h_BTag_rfj1_E_j2_mass_BTagCut = TH1F( "h_BTag_rfj1_E_j2_mass_BTagCut", "", n_bins_high, lim_tag_low, lim_tag_high )   
    h_CTag_rfj1_E_j2_mass_BTagCut = TH1F( "h_CTag_rfj1_E_j2_mass_BTagCut", "", n_bins_high, lim_tag_low, lim_tag_high )    
    h_LTag_rfj1_E_j2_mass_BTagCut = TH1F( "h_LTag_rfj1_E_j2_mass_BTagCut", "", n_bins_high,  lim_tag_low, lim_tag_high ) 
    h_BTagMax_j2_mass_BTagCut = TH1F( "h_BTagMax_j2_mass_BTagCut", "", n_bins_high, lim_tag_low, lim_tag_high )  
    #index 50 done by here
    h_m1_min_m2 = TH1F( "h_m1_min_m2", "", n_bins_high, lim_mass_delta_mass_low,lim_mass_delta_mass_high )  
    h_sqrtS_j1_j2_EMiss_high_sqrtS = TH1F( "h_sqrtS_j1_j2_EMiss_high_sqrtS", "", n_bins_high, lim_energy_low_mediumSqrtS,lim_energy_high)
    h_sqrtS_j1_j2_orig_high_sqrtS = TH1F( "h_sqrtS_j1_j2_orig_high_sqrtS", "", n_bins_high, lim_energy_low_mediumSqrtS,lim_energy_high)
    h_Esum_j1_j2_EMiss_sqrtS_2500 = TH1F( "h_Esum_j1_j2_EMiss_sqrtS_2500", "", n_bins_high, lim_energy_low_mediumSqrtS,lim_energy_high)
    h_Esum_j1_j2_EMiss_sqrtS_2500_orig = TH1F( "h_Esum_j1_j2_EMiss_sqrtS_2500_orig", "", n_bins_high, lim_energy_low_mediumSqrtS,lim_energy_high)
    #index 55 done now
    h_Esum_j1_j2_orig_sqrtS_2500_orig = TH1F( "h_Esum_j1_j2_orig_sqrtS_2500_orig", "", n_bins_high, lim_energy_low_mediumSqrtS,lim_energy_high)
    h_MET_over_Esum_j1_j2_orig_sqrtS_2500_orig = TH1F( "h_MET_over_Esum_j1_j2_orig_sqrtS_2500_orig", "", n_bins_high, lim_energy_rat_j2_low, lim_energy_rat_j1_low)
    h_MET_over_sqrtS_j1_j2_orig_sqrtS_2500_orig = TH1F( "h_MET_over_sqrtS_j1_j2_orig_sqrtS_2500_orig", "", n_bins_high, lim_energy_rat_j2_low, lim_energy_rat_j1_low)
    h_MET_over_Esum_j1_j2_orig_sqrtS_2500_EMiss = TH1F( "h_MET_over_Esum_j1_j2_orig_sqrtS_2500_EMiss", "", n_bins_high, lim_energy_rat_j2_low, lim_energy_rat_j1_low)
    h_MET_over_sqrtS_j1_j2_orig_sqrtS_2500_EMiss = TH1F( "h_MET_over_sqrtS_j1_j2_orig_sqrtS_2500_EMiss", "", n_bins_high, lim_energy_rat_j2_low, lim_energy_rat_j1_low)
    #index 60 done now
    h_MET_over_E_totPFO_sqrtS_2500_orig = TH1F( "h_MET_over_E_totPFO_sqrtS_2500_orig", "", n_bins_high, lim_energy_rat_j2_low, lim_energy_rat_j1_low)
    h_MET_over_E_totPFO_sqrtS_2500_EMiss = TH1F( "h_MET_over_E_totPFO_sqrtS_2500_EMiss", "", n_bins_high, lim_energy_rat_j2_low, lim_energy_rat_j1_low)
    h_Esum_j1_j2_orig_sqrtS_2500_EMiss = TH1F( "h_Esum_j1_j2_orig_EMiss_sqrtS_2500_EMiss", "", n_bins_high, lim_energy_low_mediumSqrtS,lim_energy_high)
    h_sqrtS_j1_j2_orig_sqrtS_2500_EMiss = TH1F( "h_sqrtS_j1_j2_EMiss_sqrtS_2500_EMiss", "", n_bins_high, lim_energy_low_mediumSqrtS,lim_energy_high)
    h_jet1_beta1_N2 = TH1F( "h_jet1_beta1_N2", "", n_bins_high, lim_N2_low,lim_b1_N2_high)
    #index 65 done now
    h_jet1_beta2_N2 = TH1F( "h_jet1_beta2_N2", "", n_bins_high, lim_N2_low,lim_b2_N2_high)
    h_jet1_beta0_5_N2 = TH1F( "h_jet1_beta0_5_N2", "", n_bins_high, lim_N2_low,lim_b0_5_N2_high)
    h_jet1_beta1_N3 = TH1F( "h_jet1_beta1_N3", "", n_bins_high, lim_N3_low,lim_b1_N3_high)
    h_jet1_beta2_N3 = TH1F( "h_jet1_beta2_N3", "", n_bins_high, lim_N3_low,lim_b2_N3_high)
    h_jet1_beta0_5_N3 = TH1F( "h_jet1_beta0_5_N3", "", n_bins_high, lim_N3_low,lim_b0_5_N3_high)
    #index 70 done now
    h_jet1_beta1_C2 = TH1F( "h_jet1_beta1_C2", "", n_bins_high, lim_C2_low,lim_b1_C2_high)
    h_jet1_beta2_C2 = TH1F( "h_jet1_beta2_C2", "", n_bins_high, lim_C2_low,lim_b2_C2_high)
    h_jet1_beta0_5_C2 = TH1F( "h_jet1_beta0_5_C2", "", n_bins_high, lim_C2_low,lim_b0_5_C2_high)
    h_jet1_beta1_C3 = TH1F( "h_jet1_beta1_C3", "", n_bins_high, lim_C3_low,lim_b1_C3_high)
    h_jet1_beta2_C3 = TH1F( "h_jet1_beta2_C3", "", n_bins_high, lim_C3_low,lim_b2_C3_high)
    #index 75 done now
    h_jet1_beta0_5_C3 = TH1F( "h_jet1_beta0_5_C3", "", n_bins_high, lim_C3_low,lim_b0_5_C3_high)
    h_jet1_beta1_D2 = TH1F( "h_jet1_beta1_D2", "", n_bins_high, lim_D2_low,lim_b1_D2_high)
    h_jet1_beta2_D2 = TH1F( "h_jet1_beta2_D2", "", n_bins_high, lim_D2_low,lim_b2_D2_high)
    h_jet1_beta0_5_D2 = TH1F( "h_jet1_beta0_5_D2", "", n_bins_high, lim_D2_low,lim_b0_5_D2_high)
    h_jet1_tau21 = TH1F( "h_jet1_tau21", "", n_bins_high, lim_tau21_low,lim_tau21_high)
    #index 80 done now
    h_jet1_tau32 = TH1F( "h_jet1_tau32", "", n_bins_high, lim_tau32_low,lim_tau32_high)
    h_jet1_dij21_over_E = TH1F( "h_jet1_dij21_over_E", "", n_bins_high, lim_dij_low,lim_dij21_rat_high)
    h_jet1_dij32_over_E = TH1F( "h_jet1_dij32_over_E", "", n_bins_high, lim_dij_low,lim_dij32_rat_high)
    h_jet1_dij43_over_E = TH1F( "h_jet1_dij43_over_E", "", n_bins_high, lim_dij_low,lim_dij43_rat_high)
    h_jet1_dij32_over_dij21 = TH1F( "h_jet1_dij32_over_dij21", "", n_bins_high, lim_dij_low,lim_dij32_over_dij21_rat_high)
    #index 85 done now
    h_jet1_dij43_over_dij21 = TH1F( "h_jet1_dij43_over_dij21", "", n_bins_high, lim_dij_low,lim_dij43_over_dij21_rat_high)
    h_jet1_dij43_over_dij32 = TH1F( "h_jet1_dij43_over_dij32", "", n_bins_high, lim_dij_low,lim_dij43_over_dij32_rat_high)
    h_jet2_beta1_N2 = TH1F( "h_jet2_beta1_N2", "", n_bins_high, lim_N2_low,lim_b1_N2_high)
    h_jet2_beta2_N2 = TH1F( "h_jet2_beta2_N2", "", n_bins_high, lim_N2_low,lim_b2_N2_high)
    h_jet2_beta0_5_N2 = TH1F( "h_jet2_beta0_5_N2", "", n_bins_high, lim_N2_low,lim_b0_5_N2_high)
    #index 90 done now
    h_jet2_beta1_N3 = TH1F( "h_jet2_beta1_N3", "", n_bins_high, lim_N3_low,lim_b1_N3_high)
    h_jet2_beta2_N3 = TH1F( "h_jet2_beta2_N3", "", n_bins_high, lim_N3_low,lim_b2_N3_high)
    h_jet2_beta0_5_N3 = TH1F( "h_jet2_beta0_5_N3", "", n_bins_high, lim_N3_low,lim_b0_5_N3_high)
    h_jet2_beta1_C2 = TH1F( "h_jet2_beta1_C2", "", n_bins_high, lim_C2_low,lim_b1_C2_high)
    h_jet2_beta2_C2 = TH1F( "h_jet2_beta2_C2", "", n_bins_high, lim_C2_low,lim_b2_C2_high)
    #index 95 done now
    h_jet2_beta0_5_C2 = TH1F( "h_jet2_beta0_5_C2", "", n_bins_high, lim_C2_low,lim_b0_5_C2_high)
    h_jet2_beta1_C3 = TH1F( "h_jet2_beta1_C3", "", n_bins_high, lim_C3_low,lim_b1_C3_high)
    h_jet2_beta2_C3 = TH1F( "h_jet2_beta2_C3", "", n_bins_high, lim_C3_low,lim_b2_C3_high)
    h_jet2_beta0_5_C3 = TH1F( "h_jet2_beta0_5_C3", "", n_bins_high, lim_C3_low,lim_b0_5_C3_high)
    h_jet2_beta1_D2 = TH1F( "h_jet2_beta1_D2", "", n_bins_high, lim_D2_low,lim_b1_D2_high)
    #index 100 done now
    h_jet2_beta2_D2 = TH1F( "h_jet2_beta2_D2", "", n_bins_high, lim_D2_low,lim_b2_D2_high)
    h_jet2_beta0_5_D2 = TH1F( "h_jet2_beta0_5_D2", "", n_bins_high, lim_D2_low,lim_b0_5_D2_high)
    h_jet2_tau21 = TH1F( "h_jet2_tau21", "", n_bins_high, lim_tau21_low,lim_tau21_high)
    h_jet2_tau32 = TH1F( "h_jet2_tau32", "", n_bins_high, lim_tau32_low,lim_tau32_high)
    h_jet2_dij21_over_E = TH1F( "h_jet2_dij21_over_E", "", n_bins_high, lim_dij_low,lim_dij21_rat_high)
    #index 105 done now
    h_jet2_dij32_over_E = TH1F( "h_jet2_dij32_over_E", "", n_bins_high, lim_dij_low,lim_dij32_rat_high)
    h_jet2_dij43_over_E = TH1F( "h_jet2_dij43_over_E", "", n_bins_high, lim_dij_low,lim_dij43_rat_high)
    h_jet2_dij32_over_dij21 = TH1F( "h_jet2_dij32_over_dij21", "", n_bins_high, lim_dij_low,lim_dij32_over_dij21_rat_high)
    h_jet2_dij43_over_dij21 = TH1F( "h_jet2_dij43_over_dij21", "", n_bins_high, lim_dij_low,lim_dij43_over_dij21_rat_high)
    h_jet2_dij43_over_dij32 = TH1F( "h_jet2_dij43_over_dij32", "", n_bins_high, lim_dij_low,lim_dij43_over_dij32_rat_high)
    #index 110 done now
    #comment these out later, just for checking the diquark, fourquark and 6 quark sqrtS for reconstructed sqrtS of 2500 with and without EMiss correction
    h_diquark_sqrtS_reco_sqrtS_2500_orig = TH1F( "h_diquark_sqrtS_reco_sqrtS_2500_orig", "", n_bins_high,  0,lim_energy_high)
    h_diquark_sqrtS_reco_sqrtS_2500_EMiss = TH1F( "h_diquark_sqrtS_reco_sqrtS_2500_EMiss", "", n_bins_high, 0, lim_energy_high)
    h_fourquark_sqrtS_reco_sqrtS_2500_orig = TH1F( "h_fourquark_sqrtS_reco_sqrtS_2500_orig", "", n_bins_high, 0,lim_energy_high)
    h_fourquark_sqrtS_reco_sqrtS_2500_EMiss = TH1F( "h_fourquark_sqrtS_reco_sqrtS_2500_EMiss", "", n_bins_high, 0, lim_energy_high)
    h_sixquark_sqrtS_reco_sqrtS_2500_orig = TH1F( "h_sixquark_sqrtS_reco_sqrtS_2500_orig", "", n_bins_high,  0,lim_energy_high)
    #index 115 done now
    h_sixquark_sqrtS_reco_sqrtS_2500_EMiss = TH1F( "h_sixquark_sqrtS_reco_sqrtS_2500_EMiss", "", n_bins_high, 0, lim_energy_high)
    h_ee_out_sqrtS_reco_sqrtS_2500_orig = TH1F( "h_ee_out_sqrtS_reco_sqrtS_2500_orig", "", n_bins_high,  0,lim_energy_high)
    h_ee_out_sqrtS_reco_sqrtS_2500_EMiss = TH1F( "h_ee_out_sqrtS_reco_sqrtS_2500_EMiss", "", n_bins_high, 0, lim_energy_high)
    h_diquark_E_sum_reco_sqrtS_2500_orig = TH1F( "h_diquark_E_sum_reco_sqrtS_2500_orig", "", n_bins_high,  0,lim_energy_high)
    h_diquark_E_sum_reco_sqrtS_2500_EMiss = TH1F( "h_diquark_E_sum_reco_sqrtS_2500_EMiss", "", n_bins_high, 0, lim_energy_high)
    #index 120 done now
    h_fourquark_E_sum_reco_sqrtS_2500_orig = TH1F( "h_fourquark_E_sum_reco_sqrtS_2500_orig", "", n_bins_high, 0,lim_energy_high)
    h_fourquark_E_sum_reco_sqrtS_2500_EMiss = TH1F( "h_fourquark_E_sum_reco_sqrtS_2500_EMiss", "", n_bins_high, 0, lim_energy_high)
    h_sixquark_E_sum_reco_sqrtS_2500_orig = TH1F( "h_sixquark_E_sum_reco_sqrtS_2500_orig", "", n_bins_high,  0,lim_energy_high)
    h_sixquark_E_sum_reco_sqrtS_2500_EMiss = TH1F( "h_sixquark_E_sum_reco_sqrtS_2500_EMiss", "", n_bins_high, 0, lim_energy_high)
    h_ee_out_E_sum_reco_sqrtS_2500_orig = TH1F( "h_ee_out_E_sum_reco_sqrtS_2500_orig", "", n_bins_high,  0,lim_energy_high)
    #index 125 done now
    h_ee_out_E_sum_reco_sqrtS_2500_EMiss = TH1F( "h_ee_out_E_sum_reco_sqrtS_2500_EMiss", "", n_bins_high, 0, lim_energy_high)
    h_angle_jet_m1_jet_m2 = TH1F( "h_angle_jet_m1_jet_m2", "", n_bins_high, -lim_dangle_max,lim_dangle_max)

    h_jet1_beta1_N2_E_theta = TH1F( "h_jet1_beta1_N2_E_theta", "", n_bins_high, lim_N2_low,lim_b1_N2_high)
    h_jet1_beta2_N2_E_theta = TH1F( "h_jet1_beta2_N2_E_theta", "", n_bins_high, lim_N2_low,lim_b2_N2_high)
    h_jet1_beta0_5_N2_E_theta = TH1F( "h_jet1_beta0_5_N2_E_theta", "", n_bins_high, lim_N2_low,lim_b0_5_N2_high)
    #index 130 done now
    h_jet1_beta1_N3_E_theta = TH1F( "h_jet1_beta1_N3_E_theta", "", n_bins_high, lim_N3_low,lim_b1_N3_high)
    h_jet1_beta2_N3_E_theta = TH1F( "h_jet1_beta2_N3_E_theta", "", n_bins_high, lim_N3_low,lim_b2_N3_high)
    h_jet1_beta0_5_N3_E_theta = TH1F( "h_jet1_beta0_5_N3_E_theta", "", n_bins_high, lim_N3_low,lim_b0_5_N3_high)
    h_jet1_beta1_C2_E_theta = TH1F( "h_jet1_beta1_C2_E_theta", "", n_bins_high, lim_C2_low,lim_b1_C2_high)
    h_jet1_beta2_C2_E_theta = TH1F( "h_jet1_beta2_C2_E_theta", "", n_bins_high, lim_C2_low,lim_b2_C2_high)
    #index 135 done now
    h_jet1_beta0_5_C2_E_theta = TH1F( "h_jet1_beta0_5_C2_E_theta", "", n_bins_high, lim_C2_low,lim_b0_5_C2_high)
    h_jet1_beta1_C3_E_theta = TH1F( "h_jet1_beta1_C3_E_theta", "", n_bins_high, lim_C3_low,lim_b1_C3_high)
    h_jet1_beta2_C3_E_theta = TH1F( "h_jet1_beta2_C3_E_theta", "", n_bins_high, lim_C3_low,lim_b2_C3_high)
    h_jet1_beta0_5_C3_E_theta = TH1F( "h_jet1_beta0_5_C3_E_theta", "", n_bins_high, lim_C3_low,lim_b0_5_C3_high)
    h_jet1_beta1_D2_E_theta = TH1F( "h_jet1_beta1_D2_E_theta", "", n_bins_high, lim_D2_low,lim_b1_D2_high)
    #index 140 done now
    h_jet1_beta2_D2_E_theta = TH1F( "h_jet1_beta2_D2_E_theta", "", n_bins_high, lim_D2_low,lim_b2_D2_high)
    h_jet1_beta0_5_D2_E_theta = TH1F( "h_jet1_beta0_5_D2_E_theta", "", n_bins_high, lim_D2_low,lim_b0_5_D2_high)
    h_jet1_tau21_lrz = TH1F( "h_jet1_tau21_lrz", "", n_bins_high, lim_tau21_low,lim_tau21_high)
    h_jet1_tau32_lrz = TH1F( "h_jet1_tau32_lrz", "", n_bins_high, lim_tau32_low,lim_tau32_high)
    h_jet2_beta1_N2_E_theta = TH1F( "h_jet2_beta1_N2_E_theta", "", n_bins_high, lim_N2_low,lim_b1_N2_high)
    #index 145 done now
    h_jet2_beta2_N2_E_theta = TH1F( "h_jet2_beta2_N2_E_theta", "", n_bins_high, lim_N2_low,lim_b2_N2_high)
    h_jet2_beta0_5_N2_E_theta = TH1F( "h_jet2_beta0_5_N2_E_theta", "", n_bins_high, lim_N2_low,lim_b0_5_N2_high)
    h_jet2_beta1_N3_E_theta = TH1F( "h_jet2_beta1_N3_E_theta", "", n_bins_high, lim_N3_low,lim_b1_N3_high)
    h_jet2_beta2_N3_E_theta = TH1F( "h_jet2_beta2_N3_E_theta", "", n_bins_high, lim_N3_low,lim_b2_N3_high)
    h_jet2_beta0_5_N3_E_theta = TH1F( "h_jet2_beta0_5_N3_E_theta", "", n_bins_high, lim_N3_low,lim_b0_5_N3_high)
    #index 150 done now
    h_jet2_beta1_C2_E_theta = TH1F( "h_jet2_beta1_C2_E_theta", "", n_bins_high, lim_C2_low,lim_b1_C2_high)
    h_jet2_beta2_C2_E_theta = TH1F( "h_jet2_beta2_C2_E_theta", "", n_bins_high, lim_C2_low,lim_b2_C2_high)
    h_jet2_beta0_5_C2_E_theta = TH1F( "h_jet2_beta0_5_C2_E_theta", "", n_bins_high, lim_C2_low,lim_b0_5_C2_high)
    h_jet2_beta1_C3_E_theta = TH1F( "h_jet2_beta1_C3_E_theta", "", n_bins_high, lim_C3_low,lim_b1_C3_high)
    h_jet2_beta2_C3_E_theta = TH1F( "h_jet2_beta2_C3_E_theta", "", n_bins_high, lim_C3_low,lim_b2_C3_high)
    #index 155 done now
    h_jet2_beta0_5_C3_E_theta = TH1F( "h_jet2_beta0_5_C3_E_theta", "", n_bins_high, lim_C3_low,lim_b0_5_C3_high)
    h_jet2_beta1_D2_E_theta = TH1F( "h_jet2_beta1_D2_E_theta", "", n_bins_high, lim_D2_low,lim_b1_D2_high)
    h_jet2_beta2_D2_E_theta = TH1F( "h_jet2_beta2_D2_E_theta", "", n_bins_high, lim_D2_low,lim_b2_D2_high)
    h_jet2_beta0_5_D2_E_theta = TH1F( "h_jet2_beta0_5_D2_E_theta", "", n_bins_high, lim_D2_low,lim_b0_5_D2_high)
    h_jet2_tau21_lrz = TH1F( "h_jet2_tau21_lrz", "", n_bins_high, lim_tau21_low,lim_tau21_high)
    #index 160 done now
    h_jet2_tau32_lrz = TH1F( "h_jet2_tau32_lrz", "", n_bins_high, lim_tau32_low,lim_tau32_high)
    lim_y32_low=0.
    lim_y32_high=0.006
    h_reco_y32 = TH1F( "h_reco_y32", "", n_bins_high, lim_y32_low,lim_y32_high)

 
    h_signal_background_1D_hist_list=[]
    h_signal_background_1D_hist_list.append(h_mass_jet1)
    h_signal_background_1D_hist_list.append(h_mass_jet2)
    h_signal_background_1D_hist_list.append(h_mass_rfjet1)
    h_signal_background_1D_hist_list.append(h_mass_rfjet2)
    h_signal_background_1D_hist_list.append(h_mass_jet_sj_BTagMax)
    h_signal_background_1D_hist_list.append(h_mass_jet_not_sj_BTagMax)
    #index 5 histos done
    h_signal_background_1D_hist_list.append(h_mass_jet_sj_LTagMax)
    h_signal_background_1D_hist_list.append(h_mass_jet_not_sj_LTagMax)
    h_signal_background_1D_hist_list.append(h_mass_jet_BTagSumMax)
    h_signal_background_1D_hist_list.append(h_mass_jet_not_BTagSumMax)
    h_signal_background_1D_hist_list.append(h_mass_jet_LTagSumMax)
    #index 10 histos done
    h_signal_background_1D_hist_list.append(h_mass_jet_not_LTagSumMax)
    h_signal_background_1D_hist_list.append(h_theta_jet1)
    h_signal_background_1D_hist_list.append(h_theta_jet2)
    h_signal_background_1D_hist_list.append(h_BTag_sj_BTagMax)
    h_signal_background_1D_hist_list.append(h_CTag_sj_CTagMax)
    #index 15 histos done
    h_signal_background_1D_hist_list.append(h_LTag_sj_LTagMax)
    h_signal_background_1D_hist_list.append(h_BTagSum_jet_BTagSumMax)
    h_signal_background_1D_hist_list.append(h_BTagSum_jet_BTagSumMin)
    h_signal_background_1D_hist_list.append(h_CTagSum_jet_BTagSumMin)
    h_signal_background_1D_hist_list.append(h_OTagSum_jet_BTagSumMin)
    #index 20 histos done
    h_signal_background_1D_hist_list.append(h_theta_jet_m1_min_jet_m2)
    h_signal_background_1D_hist_list.append(h_dphi_jet_m1_jet_m2)
    h_signal_background_1D_hist_list.append(h_CTag_sj_BTagMax)
    h_signal_background_1D_hist_list.append(h_LTag_sj_BTagMax)
    h_signal_background_1D_hist_list.append(h_BTag_sj_LTagMax)
    #index 25 histos done
    h_signal_background_1D_hist_list.append(h_CTag_sj_LTagMax)
    h_signal_background_1D_hist_list.append(h_BTag_sj2_BTagMax)
    h_signal_background_1D_hist_list.append(h_CTag_sj2_BTagMax)
    h_signal_background_1D_hist_list.append(h_LTag_sj2_BTagMax)
    h_signal_background_1D_hist_list.append(h_BTag_sj1_j2_sj_BTagMax)
    #index 30 histos done
    h_signal_background_1D_hist_list.append(h_CTag_sj1_BTag_j2_sj_BTagMax)
    h_signal_background_1D_hist_list.append(h_LTag_sj1_BTag_j2_sj_BTagMax)
    h_signal_background_1D_hist_list.append(h_BTagMax_j1_mass_noBTagCut)
    h_signal_background_1D_hist_list.append(h_BTag_rfj1_E_j2_mass_noBTagCut)
    h_signal_background_1D_hist_list.append(h_CTag_rfj1_E_j2_mass_noBTagCut)
    #ubdex 35 done
    h_signal_background_1D_hist_list.append(h_LTag_rfj1_E_j2_mass_noBTagCut)
    h_signal_background_1D_hist_list.append(h_BTagMax_j2_mass_noBTagCut)
    h_signal_background_1D_hist_list.append(h_E_rfj_BTagMax_over_E_rfj_comb_j1_mass)
    h_signal_background_1D_hist_list.append(h_E_rfsj1_over_E_rfj_j1_mass)
    h_signal_background_1D_hist_list.append(h_E_sj1_over_E_j1_mass)
    #index 40 done
    h_signal_background_1D_hist_list.append(h_E_rfsj1_over_E_rfj_j2_mass)
    h_signal_background_1D_hist_list.append(h_E_sj1_over_E_j2_mass)
    h_signal_background_1D_hist_list.append(h_Angle_sj1_sj2_j1_mass)
    h_signal_background_1D_hist_list.append(h_Angle_srfj1_srfj2_j1_mass)
    h_signal_background_1D_hist_list.append(h_Angle_sj1_sj2_j2_mass)
    #index 45 done
    h_signal_background_1D_hist_list.append(h_Angle_srfj1_srfj2_j2_mass)
    h_signal_background_1D_hist_list.append(h_BTag_rfj1_E_j2_mass_BTagCut)
    h_signal_background_1D_hist_list.append(h_CTag_rfj1_E_j2_mass_BTagCut)
    h_signal_background_1D_hist_list.append(h_LTag_rfj1_E_j2_mass_BTagCut)
    h_signal_background_1D_hist_list.append(h_BTagMax_j2_mass_BTagCut)
    #index 50 done by here
    h_signal_background_1D_hist_list.append(h_m1_min_m2)
    h_signal_background_1D_hist_list.append(h_sqrtS_j1_j2_EMiss_high_sqrtS)
    h_signal_background_1D_hist_list.append(h_sqrtS_j1_j2_orig_high_sqrtS)
    h_signal_background_1D_hist_list.append(h_Esum_j1_j2_EMiss_sqrtS_2500)
    h_signal_background_1D_hist_list.append(h_Esum_j1_j2_EMiss_sqrtS_2500_orig)
    #index 55 done now
    h_signal_background_1D_hist_list.append(h_Esum_j1_j2_orig_sqrtS_2500_orig)
    h_signal_background_1D_hist_list.append(h_MET_over_Esum_j1_j2_orig_sqrtS_2500_orig)
    h_signal_background_1D_hist_list.append(h_MET_over_sqrtS_j1_j2_orig_sqrtS_2500_orig)
    h_signal_background_1D_hist_list.append(h_MET_over_Esum_j1_j2_orig_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_MET_over_sqrtS_j1_j2_orig_sqrtS_2500_EMiss)
    #index 60 done now
    h_signal_background_1D_hist_list.append(h_MET_over_E_totPFO_sqrtS_2500_orig)
    h_signal_background_1D_hist_list.append(h_MET_over_E_totPFO_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_Esum_j1_j2_orig_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_sqrtS_j1_j2_orig_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_jet1_beta1_N2)
    #index 65 done now
    h_signal_background_1D_hist_list.append(h_jet1_beta2_N2)
    h_signal_background_1D_hist_list.append(h_jet1_beta0_5_N2)
    h_signal_background_1D_hist_list.append(h_jet1_beta1_N3)
    h_signal_background_1D_hist_list.append(h_jet1_beta2_N3)
    h_signal_background_1D_hist_list.append(h_jet1_beta0_5_N3)
    #index 70 done now
    h_signal_background_1D_hist_list.append(h_jet1_beta1_C2)
    h_signal_background_1D_hist_list.append(h_jet1_beta2_C2)
    h_signal_background_1D_hist_list.append(h_jet1_beta0_5_C2)
    h_signal_background_1D_hist_list.append(h_jet1_beta1_C3)
    h_signal_background_1D_hist_list.append(h_jet1_beta2_C3)
    #index 75 done now
    h_signal_background_1D_hist_list.append(h_jet1_beta0_5_C3)
    h_signal_background_1D_hist_list.append(h_jet1_beta1_D2)
    h_signal_background_1D_hist_list.append(h_jet1_beta2_D2)
    h_signal_background_1D_hist_list.append(h_jet1_beta0_5_D2)
    h_signal_background_1D_hist_list.append(h_jet1_tau21)
    #index 80 done now
    h_signal_background_1D_hist_list.append(h_jet1_tau32)
    h_signal_background_1D_hist_list.append(h_jet1_dij21_over_E)
    h_signal_background_1D_hist_list.append(h_jet1_dij32_over_E)
    h_signal_background_1D_hist_list.append(h_jet1_dij43_over_E)
    h_signal_background_1D_hist_list.append(h_jet1_dij32_over_dij21)
    #index 85 done now
    h_signal_background_1D_hist_list.append(h_jet1_dij43_over_dij21)
    h_signal_background_1D_hist_list.append(h_jet1_dij43_over_dij32)
    h_signal_background_1D_hist_list.append(h_jet2_beta1_N2)
    h_signal_background_1D_hist_list.append(h_jet2_beta2_N2)
    h_signal_background_1D_hist_list.append(h_jet2_beta0_5_N2)
    #index 90 done now
    h_signal_background_1D_hist_list.append(h_jet2_beta1_N3)
    h_signal_background_1D_hist_list.append(h_jet2_beta2_N3)
    h_signal_background_1D_hist_list.append(h_jet2_beta0_5_N3)
    h_signal_background_1D_hist_list.append(h_jet2_beta1_C2)
    h_signal_background_1D_hist_list.append(h_jet2_beta2_C2)
    #index 95 done now
    h_signal_background_1D_hist_list.append(h_jet2_beta0_5_C2)
    h_signal_background_1D_hist_list.append(h_jet2_beta1_C3)
    h_signal_background_1D_hist_list.append(h_jet2_beta2_C3)
    h_signal_background_1D_hist_list.append(h_jet2_beta0_5_C3)
    h_signal_background_1D_hist_list.append(h_jet2_beta1_D2)
    #index 100 done now
    h_signal_background_1D_hist_list.append(h_jet2_beta2_D2)
    h_signal_background_1D_hist_list.append(h_jet2_beta0_5_D2)
    h_signal_background_1D_hist_list.append(h_jet2_tau21)
    h_signal_background_1D_hist_list.append(h_jet2_tau32)
    h_signal_background_1D_hist_list.append(h_jet2_dij21_over_E)
    #index 105 done now
    h_signal_background_1D_hist_list.append(h_jet2_dij32_over_E)
    h_signal_background_1D_hist_list.append(h_jet2_dij43_over_E)
    h_signal_background_1D_hist_list.append(h_jet2_dij32_over_dij21)
    h_signal_background_1D_hist_list.append(h_jet2_dij43_over_dij21)
    h_signal_background_1D_hist_list.append(h_jet2_dij43_over_dij32)
    #index 110 done now
    h_signal_background_1D_hist_list.append(h_diquark_sqrtS_reco_sqrtS_2500_orig)
    h_signal_background_1D_hist_list.append(h_diquark_sqrtS_reco_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_fourquark_sqrtS_reco_sqrtS_2500_orig)
    h_signal_background_1D_hist_list.append(h_fourquark_sqrtS_reco_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_sixquark_sqrtS_reco_sqrtS_2500_orig)
    #index 115 done now
    h_signal_background_1D_hist_list.append(h_sixquark_sqrtS_reco_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_ee_out_sqrtS_reco_sqrtS_2500_orig)
    h_signal_background_1D_hist_list.append(h_ee_out_sqrtS_reco_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_diquark_E_sum_reco_sqrtS_2500_orig)
    h_signal_background_1D_hist_list.append(h_diquark_E_sum_reco_sqrtS_2500_EMiss)
    #index 120 done now
    h_signal_background_1D_hist_list.append(h_fourquark_E_sum_reco_sqrtS_2500_orig)
    h_signal_background_1D_hist_list.append(h_fourquark_E_sum_reco_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_sixquark_E_sum_reco_sqrtS_2500_orig)
    h_signal_background_1D_hist_list.append(h_sixquark_E_sum_reco_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_ee_out_E_sum_reco_sqrtS_2500_orig)
    #index 125 done now
    h_signal_background_1D_hist_list.append(h_ee_out_E_sum_reco_sqrtS_2500_EMiss)
    h_signal_background_1D_hist_list.append(h_angle_jet_m1_jet_m2)
    h_signal_background_1D_hist_list.append(h_jet1_beta1_N2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta2_N2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta0_5_N2_E_theta)
    #index 130 done now
    h_signal_background_1D_hist_list.append(h_jet1_beta1_N3_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta2_N3_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta0_5_N3_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta1_C2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta2_C2_E_theta)
    #index 135 done now     
    h_signal_background_1D_hist_list.append(h_jet1_beta0_5_C2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta1_C3_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta2_C3_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta0_5_C3_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta1_D2_E_theta)
    #index 140 done now
    h_signal_background_1D_hist_list.append(h_jet1_beta2_D2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_beta0_5_D2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet1_tau21_lrz)
    h_signal_background_1D_hist_list.append(h_jet1_tau32_lrz)
    h_signal_background_1D_hist_list.append(h_jet2_beta1_N2_E_theta)
    #index 145 done now
    h_signal_background_1D_hist_list.append(h_jet2_beta2_N2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta0_5_N2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta1_N3_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta2_N3_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta0_5_N3_E_theta)
    #index 150 done now
    h_signal_background_1D_hist_list.append(h_jet2_beta1_C2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta2_C2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta0_5_C2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta1_C3_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta2_C3_E_theta)
    #index 155 done now
    h_signal_background_1D_hist_list.append(h_jet2_beta0_5_C3_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta1_D2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta2_D2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_beta0_5_D2_E_theta)
    h_signal_background_1D_hist_list.append(h_jet2_tau21_lrz)
    #index 160 done now
    h_signal_background_1D_hist_list.append(h_jet2_tau32_lrz)
    h_signal_background_1D_hist_list.append(h_reco_y32)



    for hist in h_signal_background_1D_hist_list:
        hist.Sumw2()
        
    h_2D_mass_jet1_vs_jet2 = TH2F( "h_2D_mass_jet1_vs_jet2", "", n_bins_high, lim_mass_low, lim_mass_high, n_bins_high, lim_mass_low, lim_mass_high )
    h_2D_mass_rfjet1_vs_rfjet2 = TH2F( "h_2D_mass_rfjet1_vs_rfjet2", "", n_bins_high, lim_mass_low, lim_mass_high, n_bins_high, lim_mass_low, lim_mass_high )
    h_2D_theta_jet1_vs_jet2 = TH2F( "h_2D_theta_jet1_vs_jet2", "", n_bins_high, lim_theta_min, lim_theta_max, n_bins_high, lim_theta_min, lim_theta_max )
    h_2D_mass_jet_sj_BTagMax_vs_jet_sj_not_BTagMax = TH2F( "h_2D_mass_jet_sj_BTagMax_vs_jet_sj_not_BTagMax", "", n_bins_high, lim_mass_low, lim_mass_high, n_bins_high, lim_mass_low, lim_mass_high )
    h_2D_mass_jet_sj_LTagMax_vs_jet_sj_not_LTagMax = TH2F( "h_2D_mass_jet_sj_LTagMax_vs_jet_sj_not_LTagMax", "", n_bins_high, lim_mass_low, lim_mass_high, n_bins_high, lim_mass_low, lim_mass_high )
    h_2D_mass_jet_BTagSumMax_vs_jet_not_BTagSumMax = TH2F( "h_2D_mass_jet_BTagSumMax_vs_jet_not_BTagSumMax", "", n_bins_high, lim_mass_low, lim_mass_high, n_bins_high, lim_mass_low, lim_mass_high )
    #index 5 histos done
    h_2D_mass_jet_LTagSumMax_vs_jet_not_LTagSumMax = TH2F( "h_2D_mass_jet_LTagSumMax_vs_jet_not_LTagSumMax", "", n_bins_high, lim_mass_low, lim_mass_high, n_bins_high, lim_mass_low, lim_mass_high )
    h_2D_BTag_sj_BTagMax_vs_LTag_sj_LTagMax = TH2F( "h_2D_BTag_sj_BTagMax_vs_LTag_sj_LTagMax", "", n_bins_high, lim_tag_low, lim_tag_high, n_bins_high, lim_tag_low, lim_tag_high )
    h_2D_BTagSum_jet_BTagSumMax_vs_BTagSum_jet_BTagSumMin = TH2F( "h_2D_BTagSum_jet_BTagSumMax_vs_BTagSum_jet_BTagSumMin", "", n_bins_high, lim_tag_low, 2.*lim_tag_high, n_bins_high, lim_tag_low, 2.*lim_tag_high )    
    h_2D_BTagSum_jet_BTagSumMax_vs_LTagSum_jet_BTagSumMin = TH2F( "h_2D_BTagSum_jet_BTagSumMax_vs_LTagSum_jet_BTagSumMin", "", n_bins_high, lim_tag_low, 2.*lim_tag_high, n_bins_high, lim_tag_low, 2.*lim_tag_high )  
    h_2D_BTagSum_jet_BTagSumMax_vs_CTagSum_jet_BTagSumMin = TH2F( "h_2D_BTagSum_jet_BTagSumMax_vs_CTagSum_jet_BTagSumMin", "", n_bins_high, lim_tag_low, 2.*lim_tag_high, n_bins_high, lim_tag_low, 2.*lim_tag_high )  
    #index 10 histos done
    h_2D_sjcharge_neg_vs_sjcharge_pos_jet_BTagSumMax = TH2F( "h_2D_sjcharge_neg_vs_sjcharge_pos_jet_BTagSumMax", "", n_bins_high, lim_subjet_Charge_low, lim_subjet_Charge_high, n_bins_high, lim_subjet_Charge_low, lim_subjet_Charge_high) 
    h_2D_sjcharge_neg_vs_sjcharge_pos_jet_BTagSumMin = TH2F( "h_2D_sjcharge_neg_vs_sjcharge_pos_jet_BTagSumMin", "", n_bins_high, lim_subjet_Charge_low, lim_subjet_Charge_high, n_bins_high, lim_subjet_Charge_low, lim_subjet_Charge_high) 
    h_2D_BTag_sj_BTagMax_vs_BTag_sj2_BTagMax = TH2F( "h_2D_BTag_sj_BTagMax_vs_BTag_sj2_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high, n_bins_high, lim_tag_low, lim_tag_high )
    h_2D_BTag_sj_BTagMax_vs_CTag_sj2_BTagMax = TH2F( "h_2D_BTag_sj_BTagMax_vs_CTag_sj2_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high, n_bins_high, lim_tag_low, lim_tag_high )    
    h_2D_BTag_sj_BTagMax_vs_LTag_sj2_BTagMax = TH2F( "h_2D_BTag_sj_BTagMax_vs_LTag_sj2_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high, n_bins_high, lim_tag_low, lim_tag_high )
    #index 15 histos done by now
    h_2D_BTag_sj_BTagMax_vs_BTag_sj1_j2_sj_BTagMax = TH2F( "h_2D_BTag_sj_BTagMax_vs_BTag_sj1_j2_sj_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high, n_bins_high, lim_tag_low, lim_tag_high )
    h_2D_BTag_sj1_vs_BTag_sj2_j2_sj_BTagMax = TH2F( "h_2D_BTag_sj1_vs_BTag_sj2_BTag_j2_sj_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high, n_bins_high, lim_tag_low, lim_tag_high )
    h_2D_BTag_sj_BTagMax_vs_CTag_sj1_j2_sj_BTagMax = TH2F( "h_2D_BTag_sj_BTagMax_vs_CTag_sj1_j2_sj_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high, n_bins_high, lim_tag_low, lim_tag_high )
    h_2D_BTag_sj_BTagMax_vs_LTag_sj1_j2_sj_BTagMax = TH2F( "h_2D_BTag_sj_BTagMax_vs_LTag_sj1_j2_sj_BTagMax", "", n_bins_high, lim_tag_low, lim_tag_high, n_bins_high, lim_tag_low, lim_tag_high )

    h_2D_beta1_N2_jet1_vs_jet2 = TH2F( "h_2D_beta1_N2_jet1_vs_jet2", "", n_bins_high, lim_N2_low,lim_b1_N2_high, n_bins_high, lim_N2_low,lim_b1_N2_high)
    #index 20 done now
    h_2D_beta2_N2_jet1_vs_jet2 = TH2F( "h_2D_beta2_N2_jet1_vs_jet2", "", n_bins_high, lim_N2_low,lim_b2_N2_high, n_bins_high, lim_N2_low,lim_b2_N2_high)
    h_2D_beta0_5_N2_jet1_vs_jet2 = TH2F( "h_2D_beta0_5_N2_jet1_vs_jet2", "", n_bins_high, lim_N2_low,lim_b0_5_N2_high, n_bins_high, lim_N2_low,lim_b0_5_N2_high)
    h_2D_beta1_N3_jet1_vs_jet2 = TH2F( "h_2D_beta1_N3_jet1_vs_jet2", "", n_bins_high, lim_N3_low,lim_b1_N3_high, n_bins_high, lim_N3_low,lim_b1_N3_high)
    h_2D_beta2_N3_jet1_vs_jet2 = TH2F( "h_2D_beta2_N3_jet1_vs_jet2", "", n_bins_high, lim_N3_low,lim_b2_N3_high, n_bins_high, lim_N3_low,lim_b2_N3_high)
    h_2D_beta0_5_N3_jet1_vs_jet2 = TH2F( "h_2D_beta0_5_N3_jet1_vs_jet2", "", n_bins_high, lim_N3_low,lim_b0_5_N3_high, n_bins_high, lim_N3_low,lim_b0_5_N3_high)
    #index 25 done now
    h_2D_beta1_C2_jet1_vs_jet2 = TH2F( "h_2D_beta1_C2_jet1_vs_jet2", "", n_bins_high, lim_C2_low,lim_b1_C2_high, n_bins_high, lim_C2_low,lim_b1_C2_high)
    h_2D_beta2_C2_jet1_vs_jet2 = TH2F( "h_2D_beta2_C2_jet1_vs_jet2", "", n_bins_high, lim_C2_low,lim_b2_C2_high, n_bins_high, lim_C2_low,lim_b2_C2_high)
    h_2D_beta0_5_C2_jet1_vs_jet2 = TH2F( "h_2D_beta0_5_C2_jet1_vs_jet2", "", n_bins_high, lim_C2_low,lim_b0_5_C2_high, n_bins_high, lim_C2_low,lim_b0_5_C2_high)
    h_2D_beta1_C3_jet1_vs_jet2 = TH2F( "h_2D_beta1_C3_jet1_vs_jet2", "", n_bins_high, lim_C3_low,lim_b1_C3_high, n_bins_high, lim_C3_low,lim_b1_C3_high)
    h_2D_beta2_C3_jet1_vs_jet2 = TH2F( "h_2D_beta2_C3_jet1_vs_jet2", "", n_bins_high, lim_C3_low,lim_b2_C3_high, n_bins_high, lim_C3_low,lim_b2_C3_high)
    #index 30 done now
    h_2D_beta0_5_C3_jet1_vs_jet2 = TH2F( "h_2D_beta0_5_C3_jet1_vs_jet2", "", n_bins_high, lim_C3_low,lim_b0_5_C3_high, n_bins_high, lim_C3_low,lim_b0_5_C3_high)
    h_2D_beta1_D2_jet1_vs_jet2 = TH2F( "h_2D_beta1_D2_jet1_vs_jet2", "", n_bins_high, lim_D2_low,lim_b1_D2_high, n_bins_high, lim_D2_low,lim_b1_D2_high)
    h_2D_beta2_D2_jet1_vs_jet2 = TH2F( "h_2D_beta2_D2_jet1_vs_jet2", "", n_bins_high, lim_D2_low,lim_b2_D2_high, n_bins_high, lim_D2_low,lim_b2_D2_high)
    h_2D_beta0_5_D2_jet1_vs_jet2 = TH2F( "h_2D_beta0_5_D2_jet1_vs_jet2", "", n_bins_high, lim_D2_low,lim_b0_5_D2_high, n_bins_high, lim_D2_low,lim_b0_5_D2_high)
    h_2D_tau21_jet1_vs_jet2 = TH2F( "h_2D_tau21_jet1_vs_jet2", "", n_bins_high, lim_tau21_low,lim_tau21_high, n_bins_high, lim_tau21_low,lim_tau21_high)
    #index 35 done now
    h_2D_tau32_jet1_vs_jet2 = TH2F( "h_2D_tau32_jet1_vs_jet2", "", n_bins_high, lim_tau32_low,lim_tau32_high, n_bins_high, lim_tau32_low,lim_tau32_high)
    h_2D_dij21_over_E_jet1_vs_jet2 = TH2F( "h_2D_dij21_over_E_jet1_vs_jet2", "", n_bins_high, lim_dij_low,lim_dij21_rat_high, n_bins_high, lim_dij_low,lim_dij21_rat_high)
    h_2D_dij32_over_E_jet1_vs_jet2 = TH2F( "h_2D_dij32_over_E_jet1_vs_jet2", "", n_bins_high, lim_dij_low,lim_dij32_rat_high, n_bins_high, lim_dij_low,lim_dij32_rat_high)
    h_2D_dij43_over_E_jet1_vs_jet2 = TH2F( "h_2D_dij43_over_E_jet1_vs_jet2", "", n_bins_high, lim_dij_low,lim_dij43_rat_high, n_bins_high, lim_dij_low,lim_dij43_rat_high)
    h_2D_dij32_over_dij21_jet1_vs_jet2 = TH2F( "h_2D_dij32_over_dij21_jet1_vs_jet2", "", n_bins_high, lim_dij_low,lim_dij32_over_dij21_rat_high, n_bins_high, lim_dij_low,lim_dij32_over_dij21_rat_high)
    #index 40 done now
    h_2D_dij43_over_dij21_jet1_vs_jet2 = TH2F( "h_2D_dij43_over_dij21_jet1_vs_jet2", "", n_bins_high, lim_dij_low,lim_dij43_over_dij21_rat_high, n_bins_high, lim_dij_low,lim_dij43_over_dij21_rat_high)
    h_2D_dij43_over_dij32_jet1_vs_jet2 = TH2F( "h_2D_dij43_over_dij32_jet1_vs_jet2", "", n_bins_high, lim_dij_low,lim_dij43_over_dij32_rat_high, n_bins_high, lim_dij_low,lim_dij43_over_dij32_rat_high)
    h_2D_beta1_N2_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta1_N2_E_theta_jet1_vs_jet2", "", n_bins_high, lim_N2_low,lim_b1_N2_high, n_bins_high, lim_N2_low,lim_b1_N2_high)
    h_2D_beta2_N2_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta2_N2_E_theta_jet1_vs_jet2", "", n_bins_high, lim_N2_low,lim_b2_N2_high, n_bins_high, lim_N2_low,lim_b2_N2_high)
    h_2D_beta0_5_N2_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta0_5_N2_E_theta_jet1_vs_jet2", "", n_bins_high, lim_N2_low,lim_b0_5_N2_high, n_bins_high, lim_N2_low,lim_b0_5_N2_high)
    #index 45 done now
    h_2D_beta1_N3_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta1_N3_E_theta_jet1_vs_jet2", "", n_bins_high, lim_N3_low,lim_b1_N3_high, n_bins_high, lim_N3_low,lim_b1_N3_high)
    h_2D_beta2_N3_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta2_N3_E_theta_jet1_vs_jet2", "", n_bins_high, lim_N3_low,lim_b2_N3_high, n_bins_high, lim_N3_low,lim_b2_N3_high)
    h_2D_beta0_5_N3_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta0_5_N3_E_theta_jet1_vs_jet2", "", n_bins_high, lim_N3_low,lim_b0_5_N3_high, n_bins_high, lim_N3_low,lim_b0_5_N3_high)
    h_2D_beta1_C2_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta1_C2_E_theta_jet1_vs_jet2", "", n_bins_high, lim_C2_low,lim_b1_C2_high, n_bins_high, lim_C2_low,lim_b1_C2_high)
    h_2D_beta2_C2_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta2_C2_E_theta_jet1_vs_jet2", "", n_bins_high, lim_C2_low,lim_b2_C2_high, n_bins_high, lim_C2_low,lim_b2_C2_high)
    #index 50 done now
    h_2D_beta0_5_C2_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta0_5_C2_E_theta_jet1_vs_jet2", "", n_bins_high, lim_C2_low,lim_b0_5_C2_high, n_bins_high, lim_C2_low,lim_b0_5_C2_high)
    h_2D_beta1_C3_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta1_C3_E_theta_jet1_vs_jet2", "", n_bins_high, lim_C3_low,lim_b1_C3_high, n_bins_high, lim_C3_low,lim_b1_C3_high)
    h_2D_beta2_C3_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta2_C3_E_theta_jet1_vs_jet2", "", n_bins_high, lim_C3_low,lim_b2_C3_high, n_bins_high, lim_C3_low,lim_b2_C3_high)
    h_2D_beta0_5_C3_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta0_5_C3_E_theta_jet1_vs_jet2", "", n_bins_high, lim_C3_low,lim_b0_5_C3_high, n_bins_high, lim_C3_low,lim_b0_5_C3_high)
    h_2D_beta1_D2_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta1_D2_E_theta_jet1_vs_jet2", "", n_bins_high, lim_D2_low,lim_b1_D2_high, n_bins_high, lim_D2_low,lim_b1_D2_high)
    #index 55 done now
    h_2D_beta2_D2_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta2_D2_E_theta_jet1_vs_jet2", "", n_bins_high, lim_D2_low,lim_b2_D2_high, n_bins_high, lim_D2_low,lim_b2_D2_high)
    h_2D_beta0_5_D2_E_theta_jet1_vs_jet2 = TH2F( "h_2D_beta0_5_D2_E_theta_jet1_vs_jet2", "", n_bins_high, lim_D2_low,lim_b0_5_D2_high, n_bins_high, lim_D2_low,lim_b0_5_D2_high)
    h_2D_tau21_lrz_jet1_vs_jet2 = TH2F( "h_2D_tau21_lrz_jet1_vs_jet2", "", n_bins_high, lim_tau21_low,lim_tau21_high, n_bins_high, lim_tau21_low,lim_tau21_high)
    h_2D_tau32_lrz_jet1_vs_jet2 = TH2F( "h_2D_tau32_lrz_jet1_vs_jet2", "", n_bins_high, lim_tau32_low,lim_tau32_high, n_bins_high, lim_tau32_low,lim_tau32_high)





    h_signal_background_2D_hist_list=[]
    h_signal_background_2D_hist_list.append(h_2D_mass_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_mass_rfjet1_vs_rfjet2)
    h_signal_background_2D_hist_list.append(h_2D_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_mass_jet_sj_BTagMax_vs_jet_sj_not_BTagMax)
    h_signal_background_2D_hist_list.append(h_2D_mass_jet_sj_LTagMax_vs_jet_sj_not_LTagMax)
    h_signal_background_2D_hist_list.append(h_2D_mass_jet_BTagSumMax_vs_jet_not_BTagSumMax)
    #index 5 histos done
    h_signal_background_2D_hist_list.append(h_2D_mass_jet_LTagSumMax_vs_jet_not_LTagSumMax)
    h_signal_background_2D_hist_list.append(h_2D_BTag_sj_BTagMax_vs_LTag_sj_LTagMax)
    h_signal_background_2D_hist_list.append(h_2D_BTagSum_jet_BTagSumMax_vs_BTagSum_jet_BTagSumMin)
    h_signal_background_2D_hist_list.append(h_2D_BTagSum_jet_BTagSumMax_vs_LTagSum_jet_BTagSumMin)
    h_signal_background_2D_hist_list.append(h_2D_BTagSum_jet_BTagSumMax_vs_CTagSum_jet_BTagSumMin)
    #index 10 histos done
    h_signal_background_2D_hist_list.append(h_2D_sjcharge_neg_vs_sjcharge_pos_jet_BTagSumMax)
    h_signal_background_2D_hist_list.append(h_2D_sjcharge_neg_vs_sjcharge_pos_jet_BTagSumMin)
    h_signal_background_2D_hist_list.append(h_2D_BTag_sj_BTagMax_vs_BTag_sj2_BTagMax)
    h_signal_background_2D_hist_list.append(h_2D_BTag_sj_BTagMax_vs_CTag_sj2_BTagMax)
    h_signal_background_2D_hist_list.append(h_2D_BTag_sj_BTagMax_vs_LTag_sj2_BTagMax)
    #index 15 histos done by now
    h_signal_background_2D_hist_list.append(h_2D_BTag_sj_BTagMax_vs_BTag_sj1_j2_sj_BTagMax)
    h_signal_background_2D_hist_list.append(h_2D_BTag_sj1_vs_BTag_sj2_j2_sj_BTagMax)
    h_signal_background_2D_hist_list.append(h_2D_BTag_sj_BTagMax_vs_CTag_sj1_j2_sj_BTagMax)
    h_signal_background_2D_hist_list.append(h_2D_BTag_sj_BTagMax_vs_LTag_sj1_j2_sj_BTagMax)
    h_signal_background_2D_hist_list.append(h_2D_beta1_N2_jet1_vs_jet2)
    #index 20 done now
    h_signal_background_2D_hist_list.append(h_2D_beta2_N2_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta0_5_N2_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta1_N3_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta2_N3_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta0_5_N3_jet1_vs_jet2)
    #index 25 done now
    h_signal_background_2D_hist_list.append(h_2D_beta1_C2_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta2_C2_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta0_5_C2_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta1_C3_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta2_C3_jet1_vs_jet2)
    #index 30 done now
    h_signal_background_2D_hist_list.append(h_2D_beta0_5_C3_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta1_D2_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta2_D2_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta0_5_D2_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_tau21_jet1_vs_jet2)
    #index 35 done now
    h_signal_background_2D_hist_list.append(h_2D_tau32_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_dij21_over_E_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_dij32_over_E_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_dij43_over_E_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_dij32_over_dij21_jet1_vs_jet2)
    #index 40 done now
    h_signal_background_2D_hist_list.append(h_2D_dij43_over_dij21_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_dij43_over_dij32_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta1_N2_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta2_N2_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta0_5_N2_E_theta_jet1_vs_jet2)
    #index 45 done now
    h_signal_background_2D_hist_list.append(h_2D_beta1_N3_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta2_N3_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta0_5_N3_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta1_C2_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta2_C2_E_theta_jet1_vs_jet2)
    #index 50 done now
    h_signal_background_2D_hist_list.append(h_2D_beta0_5_C2_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta1_C3_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta2_C3_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta0_5_C3_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta1_D2_E_theta_jet1_vs_jet2)
    #index 55 done now
    h_signal_background_2D_hist_list.append(h_2D_beta2_D2_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_beta0_5_D2_E_theta_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_tau21_lrz_jet1_vs_jet2)
    h_signal_background_2D_hist_list.append(h_2D_tau32_lrz_jet1_vs_jet2)




    for hist in h_signal_background_2D_hist_list:
       hist.Sumw2()

    print 'length of hist list/hist2D list ',len(h_signal_background_1D_hist_list),len(h_signal_background_2D_hist_list)
    fill_background_histograms(input_file_,xsec_,use_partonInfo,h_signal_background_1D_hist_list,h_signal_background_2D_hist_list,lumi,bool_performMassCut,bool_useRectangularMassCut,bool_performthetaCut,bool_performBTagCut,bool_performC2Cut,bool_performD2Cut)
    
    CLICdpStyle()


    file_histogram.Write()
    file_histogram.Close()
  
    return None

def process_files():

    lumi_=4000.
    use_partonInfo_=True
    fill_signal_histos_=False

    #performMassCuts_=True
    #useMassRectCuts_=False
    #performThetaCuts_=True
    #performBTagCuts_=True

    useMassRectCuts_=False
    performMassCuts_=True
    #change to ellipse cuts for next run
    performThetaCuts_=False
    performBTagCuts_=False
    #C2 and D2 are very correlated, prefer D2 cuts for now
    performC2Cuts_=False
    performD2Cuts_=False

    print 'start processing of files'

    cross_section_= 3.83
    fill_signal_histos_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_hzqq_13391_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_hzqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar_withSignalHistos.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    fill_signal_histos_=False
    use_partonInfo_=False

    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_13394_to_13397_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qqqq_13394_to_13397_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qqqq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qqqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 902.
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_mqqqq_2_TeV_13696_to_13699_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_mqq_1TeV_13425_to_13428_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qqqq_mqqqq_2TeV_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=  369.8 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_13399_to_13402_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_13399_to_13402_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 1269.
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_


    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_mqq_1_TeV_13425_to_13428_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_mqq_1TeV_13425_to_13428_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_mqq_1TeV_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 170.8
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbcbbc_13094_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_bbcbbc_13094_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_bbcbbc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 9.2271753e-3 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbubbu_13095_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_bbubbu_13095_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_bbubbu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=  9.1731760e-3 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ddcyyc_13096_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ddcyyc_13096_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_ddcyyc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=  1.3757137
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_dduyyu_13097_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_dduyyu_13097_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_dduyyu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=  14.498909
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscbbc_13098_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_sscbbc_13098_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_sscbbc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=  12.499614
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscssc_13099_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_sscssc_13099_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_sscssc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=  1.1651315
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssussu_13123_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ssussu_13123_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_ssussu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=  1.2615661e-2
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssubbu_13292_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ssubbu_13292_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_ssubbu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=   5.4145233e-2 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycbbu_13318_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycbbu_13318_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yycbbu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=  13.394883
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycddu_13326_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycddu_13326_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yycddu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=   	2.0054737
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycssu_13323_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycssu_13323_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yycssu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 2.0248353
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyubbc_13320_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyubbc_13320_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yyubbc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=13.330064
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyuddc_13328_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyuddc_13328_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yyuddc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=2.0034170
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyussc_13325_polm80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyussc_13325_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/test_ee_yyussc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=2.0189010
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_


    #NOW WE GO TO THE POSITIVE POLARIZATION CASES

    lumi_=1000.

    use_partonInfo_=True
    fillGenInfo_=True
    fill_signal_histos_=True
    #use_partonInfo_=False
    cross_section_= 2.67
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_hzqq_13392_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/hzqq_13392_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_hzqq_noMassCut_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_OnlyPartonSelection.root"  
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_hzqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar_withSignalHistos.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    fill_signal_histos_=False
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_13393_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qqqq_13393_RunEventStatisticsHistogram.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qqqq_noSqrtS.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qqqq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root"
    cross_section_= 120.
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qqqq_mqqqq_2_TeV_13700_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7ee_qq_mqq_1_TeV_13429_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qqqq_mqqqq_2TeV_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root"
    cross_section_=   49.2
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_13398_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_qq_13398_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root"
    cross_section_= 786.
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_qq_mqq_1_TeV_13429_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7ee_qq_mqq_1_TeV_13429_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_mqq_1TeV_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root"
    cross_section_=  73.5 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_



    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbcbbc_13071_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_bbcbbc_13071_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_bbcbbc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 2.9986901e-3 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_bbubbu_13072_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_bbubbu_13072_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_bbubbu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 2.9825397e-3 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ddcyyc_13073_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ddcyyc_13073_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_ddcyyc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=1.7824610e-1
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_dduyyu_13074_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_dduyyu_13074_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_dduyyu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 5.0109474
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscbbc_13075_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_sscbbc_13075_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_sscbbc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=4.8938333
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_sscssc_13076_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_sscssc_13076_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_sscssc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 1.3677677e-1 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssussu_13077_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ssussu_13077_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_ssussu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 3.3776171e-3 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_ssubbu_13293_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_ssubbu_13293_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_ssubbu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 2.3216638e-2 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycbbu_13322_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycbbu_13322_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yycbbu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 5.2101109
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycddu_13319_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycddu_13319_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yycddu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 4.0984879e-1 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yycssu_13327_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yycssu_13327_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yycssu_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=  4.1853929e-1 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyubbc_13324_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyubbc_13324_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yyubbc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 5.2070149
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyuddc_13321_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyuddc_13321_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yyuddc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_= 4.1203686e-1 
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_

    use_partonInfo_=False
    fillGenInfo_=True
    input_file_name_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/HZStudy_ee_yyussc_13329_polp80_3TeV_wO_CLIC_o3_v14_DR7.root"
    #final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_qq_m_j1_80_160_theta1_20_160_m_j2_50_135.root"
    #input_file_name2_="/eos/user/w/weberma2/data/HZAnalyzerFiles/190709_gcc62_CT/VtxVLC7RFJVLC7/ee_yyussc_13329_RunEventStatisticsHistogram.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/test_ee_yyussc_AnalysisBaselineHistos_ellipse_m1_126_dm_35_m2_92_5_dm_35_EThetaVar.root" 
    cross_section_=4.2245034e-1
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_,use_partonInfo_,fill_signal_histos_,performMassCuts_,useMassRectCuts_,performThetaCuts_,performBTagCuts_,performC2Cuts_,performD2Cuts_)
    print 'finished file', final_histo_name_





    return None

process_files()



