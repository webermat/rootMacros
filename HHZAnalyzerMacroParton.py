from ROOT import gROOT, TCanvas, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor
import ROOT as root
from math import cos, sin, pi, degrees, radians, pow, sqrt,acos

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

def setUpperCanvas(canvas_name) :
    c1= TCanvas(canvas_name,canvas_name,10,50,600,500)
    c1.cd()
    gPad.SetTopMargin(0.06)
    return c1

def squared_mass(vl) :
    m_H=125.1
    m_Z=91.188
    if len(vl)!=3:
        print 'expect 3 vectors exactly for mass calculation, something could be wrong'
        return -1
    else:
        vector_squared_sum=(vl[0].M()-m_H)**2+(vl[1].M()-m_H)**2+(vl[2].M()-m_Z)**2
        vector_squared_sum=min(vector_squared_sum,(vl[0].M()-m_Z)**2+(vl[1].M()-m_H)**2+(vl[2].M()-m_H)**2)
        vector_squared_sum=min(vector_squared_sum,(vl[0].M()-m_H)**2+(vl[1].M()-m_Z)**2+(vl[2].M()-m_H)**2)
        return vector_squared_sum
 
 
    
def orderedThreeVector(vl) :
    result_vector=[]
    m_H=125.1
    m_Z=91.188
    if len(vl)<3:
        print 'expect at least 3 vectors, something needs to be fixed'
        return None;
    if len(vl)==3:
        result_vector = vl
    else:
        if len(vl)==4:
            minrange=-1
            for vect in (range(len(vl))):
                #print 'in index',vect,'of',len(vl)
                for vect2 in range(vect+1,len(vl),1): 
                    temp_list=vl[:]
                    #print 'vect,vect2 at start',vect,vect2,len(vl),len(temp_list)
                    temp_list[vect]=temp_list[vect]+temp_list[vect2]
                    temp_list.pop(vect2)      
                    #print 'vect,vect2',vect,vect2,len(temp_list)
                #for vect in len(range(vl-1)):
                    if(minrange<0):
                        result_vector =  temp_list
                        #print 'get in here once',vect,vect2,squared_mass(temp_list),minrange
                        minrange=squared_mass(temp_list)
                        #print 'get in here once 2',vect,vect2,squared_mass(temp_list),minrange
                    elif squared_mass(temp_list)<minrange:
                        #print 'get in here maybe',vect,vect2,squared_mass(temp_list),minrange
                        minrange=squared_mass(temp_list)
                        result_vector =  temp_list
        elif len(vl)==5:
            minrange=-1
            #print 'begin stuff here'
            for vect in (range(len(vl))):
                #print 'in index',vect,'of',len(vl)
                for vect2 in range(vect+1,len(vl),1): 
                    #print 'zero check',vect,vect+1,vect2
                    #print 'vect,vect2 at start',vect,vect2,len(vl)
                    for vect3 in range(vect+1,len(vl),1): 
                        #print 'first draft',vect,vect2,vect3
                        if vect3!=vect2 :
                            #print 'first check',vect,vect2,vect3
                            for vect4 in range(vect3+1,len(vl)):
                                #print 'second draft',vect,vect2,vect3,vect4
                                if vect4!=vect2:
                                    #print 'vect,vect2 at start',vect,vect2,len(vl),len(temp_list)
                                    temp_list=vl[:]
                                    #print 'second check',vect,vect2,vect3,vect4,len(vl),len(temp_list)
                                    temp_list[vect]=temp_list[vect]+temp_list[vect2]
                                    temp_list[vect3]=temp_list[vect3]+temp_list[vect4]
                                    if(vect2<vect4):
                                        #print 'pop1',vect,vect2,vect3,vect4,len(temp_list)
                                        temp_list.pop(vect4)      
                                        #print 'pop1 vect4 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                        temp_list.pop(vect2) 
                                        #print 'pop1 vect2 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                    else:
                                        #print 'pop2',vect,vect2,vect3,vect4,len(temp_list)
                                        temp_list.pop(vect2)      
                                        #print 'pop2 vect2 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                        temp_list.pop(vect4)  
                                        #print 'pop2 vect4 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                    if(minrange<0):
                                        result_vector =  temp_list
                                        #print 'get in here once',vect,vect2,squared_mass(temp_list),minrange
                                        minrange=squared_mass(temp_list)
                                        #print 'get in here once 2',vect,vect2,squared_mass(temp_list),minrange
                                    elif squared_mass(temp_list)<minrange:
                                        #print 'get in here maybe',vect,vect2,squared_mass(temp_list),minrange
                                        minrange=squared_mass(temp_list)
                                        result_vector =  temp_list
        elif len(vl)==6:
            minrange=-1
            #print 'begin stuff here'
            for vect in (range(len(vl))):
                #print 'in index',vect,'of',len(vl)
                for vect2 in range(vect+1,len(vl),1): 
                    #print 'zero check',vect,vect+1,vect2
                    #print 'vect,vect2 at start',vect,vect2,len(vl)
                    for vect3 in range(vect+1,len(vl),1): 
                        #print 'first draft',vect,vect2,vect3
                        if vect3!=vect2 :
                            #print 'first check',vect,vect2,vect3
                            for vect4 in range(vect3+1,len(vl)):
                                #print 'second draft',vect,vect2,vect3,vect4
                                if vect4!=vect2:
                                    #print 'vect,vect2 at start',vect,vect2,len(vl),len(temp_list)
                                    for vect5 in range(vect3+1,len(vl),1): 
                                        #print 'third draft',vect,vect2,vect3,vect4,vect5
                                        if vect5!=vect and vect5!=vect2 and vect5!=vect3 and vect5!= vect4:
                                            #print 'third check',vect,vect2,vect3,vect4,vect5
                                            for vect6 in range(vect5+1,len(vl)):
                                                #print 'fourth draft',vect,vect2,vect3,vect4,vect5,vect6
                                                if vect6!=vect2 and vect6!=vect4:
                                                    temp_list=vl[:]
                                                    #print 'vect,vect2,vect3,vect4,vect5,vect6 at at end',vect,vect2,vect3,vect4,vect5,vect6,len(vl),len(temp_list)
                                                    temp_list=vl[:]
                                                    #print 'second check',vect,vect2,vect3,vect4,len(vl),len(temp_list)
                                                    temp_list[vect]=temp_list[vect]+temp_list[vect2]
                                                    temp_list[vect3]=temp_list[vect3]+temp_list[vect4]
                                                    temp_list[vect5]=temp_list[vect5]+temp_list[vect6]
                                                    array_second=[vect2,vect4,vect6]
                                                    array_second.sort()
                                                    #print 'sorted array', array_second
                                                    temp_list.pop(array_second[2]) 
                                                    temp_list.pop(array_second[1])   
                                                    temp_list.pop(array_second[0])   
                                                    #print 'pop1 vect4 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                                    #temp_list.pop(vect2) 
                                                    #print 'pop1 vect2 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                                    #else:
                                                    #print 'pop2',vect,vect2,vect3,vect4,len(temp_list)
                                                    #temp_list.pop(vect2)      
                                                    #print 'pop2 vect2 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                                    #temp_list.pop(vect4)  
                                                    #print 'pop2 vect4 gone ',vect,vect2,vect3,vect4,len(temp_list)
                                                    if(minrange<0):
                                                        result_vector =  temp_list
                                                        #print 'get in here once',vect,vect2,squared_mass(temp_list),minrange
                                                        minrange=squared_mass(temp_list)
                                                        #print 'get in here once 2',vect,vect2,squared_mass(temp_list),minrange
                                                    elif squared_mass(temp_list)<minrange:
                                                        #print 'get in here maybe',vect,vect2,squared_mass(temp_list),minrange
                                                        minrange=squared_mass(temp_list)
                                                        result_vector =  temp_list                                                        
    result_vector.sort(key=lambda x: x.M(), reverse=True)
    return result_vector
 
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

def fill_HHZ_histograms(file,xsec,hist_vec_reco_1D,lumi):
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
    print "tree-entries ",tree.GetEntries(), " weight ",weight, "xsec",xsec,"lumi",lumi,"total event",weight*tree.GetEntries(), 'to', xsec*lumi, 'size of vector',len(hist_vec_reco_1D)

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

        #quarks ordered by energy, no cut on flavour 
        tempH1_q1=TLorentzVector(0,0,0,0)
        tempH1_q2=TLorentzVector(0,0,0,0)

        tempH2_q1=TLorentzVector(0,0,0,0)
        tempH2_q2=TLorentzVector(0,0,0,0)

        tempH1_V1_q1=TLorentzVector(0,0,0,0)
        tempH1_V1_q2=TLorentzVector(0,0,0,0)

        tempH1_V2_q1=TLorentzVector(0,0,0,0)
        tempH1_V2_q2=TLorentzVector(0,0,0,0)

        tempH2_V1_q1=TLorentzVector(0,0,0,0)
        tempH2_V1_q2=TLorentzVector(0,0,0,0)

        tempH2_V2_q1=TLorentzVector(0,0,0,0)
        tempH2_V2_q2=TLorentzVector(0,0,0,0)

        temp_true_inv=TLorentzVector(0,0,0,0)

        true_inv_E=ientry.true_inv_E
        true_inv_Px=ientry.true_inv_Px
        true_inv_Py=ientry.true_inv_Py
        true_inv_Pz=ientry.true_inv_Pz

        temp_true_inv.SetPxPyPzE(true_inv_Px,true_inv_Py,true_inv_Pz,true_inv_E)

        H1H2_decays_bbar = True
        H1H2_decays_qqbar = True
        Z_decays_qqbar = True

        trueME_E=ientry.trueME_E
        trueME_Px=ientry.trueME_Px
        trueME_Py=ientry.trueME_Py
        trueME_Pz=ientry.trueME_Pz
        trueME_PDGID=ientry.trueME_PDGID
        trueME_ParentID=ientry.trueME_ParentID
        trueME_PDGID_D1=ientry.trueME_PDGID_D1
        
        index_firstH=-1
        index_secondH=-1
        index_firstZ=-1
        
        #indices 0-1 outgoing electrons after ISR and beam strahlung (i.e. the "real collision")
        #index 2,3 ISR photons of incoming electrons
        #index 4 is H1, index 5 is H2, 6 is Z
        #index 7 and 8 are daughters from H1
        #index 9 and 10 are decay products of H2
        #index 11 and 12 are decay products of Z 
        for i in range(len(trueME_E)):
            if trueME_PDGID[4]!=25 or trueME_PDGID[5]!=25 or trueME_PDGID[6]!=23:
                num_total_exception+=1
            if i<2:
                temp=TLorentzVector(0,0,0,0)
                temp.SetPxPyPzE(trueME_Px[i],trueME_Py[i],trueME_Pz[i],trueME_E[i])
                tempTotEventP4+=temp
            if trueME_PDGID[i]==25 and index_firstH==-1:
                index_firstH=i
            if trueME_PDGID[i]==25 and i!=index_firstH and index_secondH==-1 :
                index_secondH=i
            if trueME_PDGID[i]==23 and index_firstZ==-1:
                index_firstZ=i
            if index_firstZ!=-1 and index_firstH!=-1 and index_secondH!=-1 :
                break
        if (abs(trueME_PDGID[index_firstH+3])!=5 or abs(trueME_PDGID[index_firstH+4])!=5) or (abs(trueME_PDGID[index_secondH+4])!=5 or abs(trueME_PDGID[index_secondH+4])!=5):
            H1H2_decays_bbar=False

        if (abs(trueME_PDGID[index_firstH+3])>6 or abs(trueME_PDGID[index_firstH+4])>6) or (abs(trueME_PDGID[index_secondH+4])>6 or abs(trueME_PDGID[index_secondH+4])>6):
            H1H2_decays_qqbar=False

        #H1 decays into quarks or gluons, aka jets
        if abs(trueME_PDGID[index_firstH+3])<7 or trueME_PDGID[index_firstH+3]==21 :
            if trueME_ParentID[index_firstH+3]!=1 or trueME_ParentID[index_firstH+4]!=1:
                print 'should have been daughters from H1',trueME_ParentID[index_firstH+3],trueME_ParentID[index_firstH+4]
            if trueME_E[index_firstH+3]> trueME_E[index_firstH+4]:
                tempH1_q1.SetPxPyPzE(trueME_Px[index_firstH+3],trueME_Py[index_firstH+3],trueME_Pz[index_firstH+3],trueME_E[index_firstH+3])
                tempH1_q2.SetPxPyPzE(trueME_Px[index_firstH+4],trueME_Py[index_firstH+4],trueME_Pz[index_firstH+4],trueME_E[index_firstH+4])
            else :
                tempH1_q2.SetPxPyPzE(trueME_Px[index_firstH+3],trueME_Py[index_firstH+3],trueME_Pz[index_firstH+3],trueME_E[index_firstH+3])
                tempH1_q1.SetPxPyPzE(trueME_Px[index_firstH+4],trueME_Py[index_firstH+4],trueME_Pz[index_firstH+4],trueME_E[index_firstH+4])
            if trueME_PDGID[index_firstH+3]==5 :
                tempH1_b.SetPxPyPzE(trueME_Px[index_firstH+3],trueME_Py[index_firstH+3],trueME_Pz[index_firstH+3],trueME_E[index_firstH+3])
                tempH1_bbar.SetPxPyPzE(trueME_Px[index_firstH+4],trueME_Py[index_firstH+4],trueME_Pz[index_firstH+4],trueME_E[index_firstH+4])
                if abs(trueME_PDGID[index_firstH+3]) != abs(trueME_PDGID[index_firstH+4]):
                    print 'pdg id should for H1 be the same in abs value, what is wrong ',trueME_PDGID[index_firstH+3],-trueME_PDGID[index_firstH+4]
            elif trueME_PDGID[index_firstH+3]==-5 :
                tempH1_bbar.SetPxPyPzE(trueME_Px[index_firstH+3],trueME_Py[index_firstH+3],trueME_Pz[index_firstH+3],trueME_E[index_firstH+3])
                tempH1_b.SetPxPyPzE(trueME_Px[index_firstH+4],trueME_Py[index_firstH+4],trueME_Pz[index_firstH+4],trueME_E[index_firstH+4])             
            if trueME_PDGID[index_firstH+3] != -trueME_PDGID[index_firstH+4]:
                if H1H2_decays_bbar:
                    print 'pdg id should for H1 be the same, what is wrong bbbar',trueME_PDGID[index_firstH+3],-trueME_PDGID[index_firstH+4]
        #H2 decays into quarks or gluons, aka jets (after index,+1 -->Z, +2 -->H1 d1, +3 --> H1 d2
        if abs(trueME_PDGID[index_secondH+4])<7 or trueME_PDGID[index_secondH+5]==21 :
            if trueME_ParentID[index_secondH+4]!=2 or trueME_ParentID[index_secondH+5]!=2:
                print 'should have been daughters from H2',trueME_ParentID[index_firstH+3],trueME_ParentID[index_firstH+4]
            if trueME_E[index_secondH+4]> trueME_E[index_secondH+5]:
                tempH2_q1.SetPxPyPzE(trueME_Px[index_secondH+4],trueME_Py[index_secondH+4],trueME_Pz[index_secondH+4],trueME_E[index_secondH+4])
                tempH2_q2.SetPxPyPzE(trueME_Px[index_secondH+5],trueME_Py[index_secondH+5],trueME_Pz[index_secondH+5],trueME_E[index_secondH+5])
            else :
                tempH2_q2.SetPxPyPzE(trueME_Px[index_secondH+4],trueME_Py[index_secondH+4],trueME_Pz[index_secondH+4],trueME_E[index_secondH+4])
                tempH2_q1.SetPxPyPzE(trueME_Px[index_secondH+5],trueME_Py[index_secondH+5],trueME_Pz[index_secondH+5],trueME_E[index_secondH+5])
            if trueME_PDGID[index_secondH+4]==5 :
                tempH2_b.SetPxPyPzE(trueME_Px[index_secondH+4],trueME_Py[index_secondH+4],trueME_Pz[index_secondH+4],trueME_E[index_secondH+4])
                tempH2_bbar.SetPxPyPzE(trueME_Px[index_secondH+5],trueME_Py[index_secondH+5],trueME_Pz[index_secondH+5],trueME_E[index_secondH+5])
                if abs(trueME_PDGID[index_secondH+4]) != abs(trueME_PDGID[index_secondH+5]):
                    print 'pdg id should for H2 be the same, what is wrong ',trueME_PDGID[index_secondH+4],-trueME_PDGID[index_secondH+5]
            elif trueME_PDGID[index_secondH+4]==-5 :
                tempH2_bbar.SetPxPyPzE(trueME_Px[index_secondH+4],trueME_Py[index_secondH+4],trueME_Pz[index_secondH+4],trueME_E[index_secondH+4])
                tempH2_b.SetPxPyPzE(trueME_Px[index_secondH+5],trueME_Py[index_secondH+5],trueME_Pz[index_secondH+5],trueME_E[index_secondH+5])             
            if trueME_PDGID[index_secondH+4] != -trueME_PDGID[index_secondH+5]:
                if H1H2_decays_bbar:
                    print 'pdg id should for H2 be the same, what is wrong bbbar ',trueME_PDGID[index_secondH+4],-trueME_PDGID[index_secondH+5]

        if (trueME_PDGID[index_firstZ+5]>6 or  trueME_PDGID[index_firstZ+6]>6):
            Z_decays_qqbar=False
        
        #positive charge, aka up type quarks or down bar type quarks
        if(trueME_PDGID[index_firstZ+5]==-1 or trueME_PDGID[index_firstZ+5]==2 or trueME_PDGID[index_firstZ+5]==-3 or trueME_PDGID[index_firstZ+5]==4 or trueME_PDGID[index_firstZ+5]==-5 or trueME_PDGID[index_firstZ+5]==6):
            if trueME_PDGID[index_firstZ+5] != -trueME_PDGID[index_firstZ+6]:
                print 'pdg id should be the same for Z, what is wrong ',trueME_PDGID[index_firstZ+5],-trueME_PDGID[index_firstZ+6]
            if trueME_ParentID[index_firstZ+5]!=3 or trueME_ParentID[index_firstZ+6]!=3:
                print 'should have been daughters from Z',trueME_ParentID[index_firstZ+5],trueME_ParentID[index_firstZ+6]
            tempZ_q_pos.SetPxPyPzE(trueME_Px[index_firstZ+5],trueME_Py[index_firstZ+5],trueME_Pz[index_firstZ+5],trueME_E[index_firstZ+5])
            tempZ_q_neg.SetPxPyPzE(trueME_Px[index_firstZ+6],trueME_Py[index_firstZ+6],trueME_Pz[index_firstZ+6],trueME_E[index_firstZ+6])
        else:
            #quark index 6 is negatively charged
            if trueME_PDGID[index_firstZ+5] != -trueME_PDGID[index_firstZ+6]:
                print 'pdg id should be the same for Z, what is wrong ',trueME_PDGID[index_firstZ+5],-trueME_PDGID[index_firstZ+6]
            if trueME_ParentID[index_firstZ+5]!=3 or trueME_ParentID[index_firstZ+6]!=3:
                print 'should have been daughters from Z',trueME_ParentID[index_firstZ+5],trueME_ParentID[index_firstZ+6]
            tempZ_q_neg.SetPxPyPzE(trueME_Px[index_firstZ+5],trueME_Py[index_firstZ+5],trueME_Pz[index_firstZ+5],trueME_E[index_firstZ+5])
            tempZ_q_pos.SetPxPyPzE(trueME_Px[index_firstZ+6],trueME_Py[index_firstZ+6],trueME_Pz[index_firstZ+6],trueME_E[index_firstZ+6])
            
        tempZ_first=TLorentzVector(0,0,0,0)
        tempZ_first.SetPxPyPzE(trueME_Px[index_firstZ],trueME_Py[index_firstZ],trueME_Pz[index_firstZ],trueME_E[index_firstZ])

        tempH1P4.SetPxPyPzE(trueME_Px[index_firstH],trueME_Py[index_firstH],trueME_Pz[index_firstH],trueME_E[index_firstH])
        tempH2P4.SetPxPyPzE(trueME_Px[index_secondH],trueME_Py[index_secondH],trueME_Pz[index_secondH],trueME_E[index_secondH])
        tempHHP4=tempH1P4+tempH2P4
        tempZP4=tempZ_q_pos+tempZ_q_neg
        tempTotEventP4HHZ=tempHHP4+tempZ_first
        hist_vec_reco_1D[0].Fill(tempTotEventP4.M(),weight)
        hist_vec_reco_1D[42].Fill(tempTotEventP4HHZ.M(),weight)
        hist_vec_reco_1D[1].Fill(tempH1P4.E()+tempH2P4.E(),weight)

        if(tempH1P4.E()>tempH2P4.E()):
            hist_vec_reco_1D[2].Fill(tempH1P4.E(),weight)
            hist_vec_reco_1D[3].Fill(tempH2P4.E(),weight)
            hist_vec_reco_1D[118].Fill(tempH1P4.P(),weight)
            hist_vec_reco_1D[119].Fill(tempH2P4.P(),weight)
            hist_vec_reco_1D[62].Fill(degrees(tempH1P4.Theta()),weight)
            hist_vec_reco_1D[63].Fill(degrees(tempH2P4.Theta()),weight)
            hist_vec_reco_1D[121].Fill(1-cos(tempH1P4.Angle(tempZP4.Vect())),weight)
            hist_vec_reco_1D[122].Fill(1-cos(tempH2P4.Angle(tempZP4.Vect())),weight)
        else:
            hist_vec_reco_1D[2].Fill(tempH2P4.E(),weight)
            hist_vec_reco_1D[3].Fill(tempH1P4.E(),weight)
            hist_vec_reco_1D[118].Fill(tempH2P4.P(),weight)
            hist_vec_reco_1D[119].Fill(tempH1P4.P(),weight)
            hist_vec_reco_1D[62].Fill(degrees(tempH2P4.Theta()),weight)
            hist_vec_reco_1D[63].Fill(degrees(tempH1P4.Theta()),weight)
            hist_vec_reco_1D[121].Fill(1-cos(tempH2P4.Angle(tempZP4.Vect())),weight)
            hist_vec_reco_1D[122].Fill(1-cos(tempH1P4.Angle(tempZP4.Vect())),weight)
        hist_vec_reco_1D[4].Fill(tempZ_first.E(),weight)
        hist_vec_reco_1D[120].Fill(tempZ_first.P(),weight)

        hist_vec_reco_1D[5].Fill((tempH1P4+tempH2P4).M(),weight)
        hist_vec_reco_1D[6].Fill(degrees((tempH1P4+tempH2P4).Angle(tempZ_first.Vect())),weight)
        hist_vec_reco_1D[7].Fill(degrees(DeltaPhi((tempH1P4+tempH2P4).Phi(),tempZ_first.Phi())),weight)
        hist_vec_reco_1D[8].Fill(degrees(abs((tempH1P4+tempH2P4).Theta()-tempZ_first.Theta())),weight)
        hist_vec_reco_1D[9].Fill(degrees(tempH1P4.Angle(tempH2P4.Vect())),weight)
        hist_vec_reco_1D[10].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempH2P4.Phi())),weight)
        hist_vec_reco_1D[11].Fill(degrees(abs(tempH1P4.Theta()-tempH2P4.Theta())),weight)
        hist_vec_reco_1D[39].Fill(1.-cos(tempH1P4.Angle(tempH2P4.Vect())),weight)
        hist_vec_reco_1D[64].Fill(degrees(tempZP4.Theta()),weight)

        #for histogram filling: H1 energy> H2 energy
        if H1H2_decays_bbar:
            if(tempH1P4.E()>tempH2P4.E()):
                hist_vec_reco_1D[12].Fill(degrees(tempH1_bbar.Angle(tempH1_b.Vect())),weight)
                hist_vec_reco_1D[13].Fill(degrees(DeltaPhi((tempH1_bbar).Phi(),tempH1_b.Phi())),weight)
                hist_vec_reco_1D[14].Fill(degrees(abs((tempH1_bbar).Theta()-tempH1_b.Theta())),weight)
                hist_vec_reco_1D[60].Fill(1.-cos(tempH1_bbar.Angle(tempH1_b.Vect())),weight)

                hist_vec_reco_1D[15].Fill(degrees(tempH2_bbar.Angle(tempH2_b.Vect())),weight)
                hist_vec_reco_1D[16].Fill(degrees(DeltaPhi((tempH2_bbar).Phi(),tempH2_b.Phi())),weight)
                hist_vec_reco_1D[17].Fill(degrees(abs((tempH2_bbar).Theta()-tempH2_b.Theta())),weight)
                hist_vec_reco_1D[61].Fill(1.-cos(tempH2_bbar.Angle(tempH2_b.Vect())),weight)
            else:
                hist_vec_reco_1D[12].Fill(degrees(tempH2_bbar.Angle(tempH2_b.Vect())),weight)
                hist_vec_reco_1D[13].Fill(degrees(DeltaPhi((tempH2_bbar).Phi(),tempH2_b.Phi())),weight)
                hist_vec_reco_1D[14].Fill(degrees(abs(tempH2_bbar.Theta()-tempH2_b.Theta())),weight)
                hist_vec_reco_1D[60].Fill(1.-cos(tempH2_bbar.Angle(tempH2_b.Vect())),weight)

                hist_vec_reco_1D[15].Fill(degrees(tempH1_bbar.Angle(tempH1_b.Vect())),weight)
                hist_vec_reco_1D[16].Fill(degrees(DeltaPhi((tempH1_bbar).Phi(),tempH1_b.Phi())),weight)
                hist_vec_reco_1D[17].Fill(degrees(abs((tempH1_bbar).Theta()-tempH1_b.Theta())),weight)
                hist_vec_reco_1D[61].Fill(1.-cos(tempH1_bbar.Angle(tempH1_b.Vect())),weight)


        if H1H2_decays_qqbar:
            if(tempH1P4.E()>tempH2P4.E()):
                hist_vec_reco_1D[18].Fill(degrees(tempH1_q2.Angle(tempH1_q1.Vect())),weight)
                hist_vec_reco_1D[19].Fill(degrees(DeltaPhi((tempH1_q2).Phi(),tempH1_q1.Phi())),weight)
                hist_vec_reco_1D[20].Fill(degrees(abs((tempH1_q2).Theta()-tempH1_q1.Theta())),weight)
                hist_vec_reco_1D[40].Fill(1.-cos(tempH1_q2.Angle(tempH1_q1.Vect())),weight)

                hist_vec_reco_1D[21].Fill(degrees(tempH2_q2.Angle(tempH2_q1.Vect())),weight)
                hist_vec_reco_1D[22].Fill(degrees(DeltaPhi((tempH2_q2).Phi(),tempH2_q1.Phi())),weight)
                hist_vec_reco_1D[23].Fill(degrees(abs((tempH2_q2).Theta()-tempH2_q1.Theta())),weight)
                hist_vec_reco_1D[41].Fill(1.-cos(tempH2_q2.Angle(tempH2_q1.Vect())),weight)
            else:
                hist_vec_reco_1D[18].Fill(degrees(tempH2_q2.Angle(tempH2_q1.Vect())),weight)
                hist_vec_reco_1D[19].Fill(degrees(DeltaPhi((tempH2_q2).Phi(),tempH2_q1.Phi())),weight)
                hist_vec_reco_1D[20].Fill(degrees(abs(tempH2_q2.Theta()-tempH2_q1.Theta())),weight)
                hist_vec_reco_1D[40].Fill(1.-cos(tempH2_q2.Angle(tempH2_q1.Vect())),weight)

                hist_vec_reco_1D[21].Fill(degrees(tempH1_q2.Angle(tempH1_q1.Vect())),weight)
                hist_vec_reco_1D[22].Fill(degrees(DeltaPhi((tempH1_q2).Phi(),tempH1_q1.Phi())),weight)
                hist_vec_reco_1D[23].Fill(degrees(abs((tempH1_q2).Theta()-tempH1_q1.Theta())),weight)
                hist_vec_reco_1D[41].Fill(1.-cos(tempH1_q2.Angle(tempH1_q1.Vect())),weight)

        if Z_decays_qqbar:
                hist_vec_reco_1D[24].Fill(degrees(tempZ_q_neg.Angle(tempZ_q_pos.Vect())),weight)
                hist_vec_reco_1D[25].Fill(degrees(DeltaPhi((tempZ_q_neg).Phi(),tempZ_q_pos.Phi())),weight)
                hist_vec_reco_1D[26].Fill(degrees(abs((tempZ_q_neg).Theta()-tempZ_q_pos.Theta())),weight)
                hist_vec_reco_1D[59].Fill(1.-cos(tempZ_q_neg.Angle(tempZ_q_pos.Vect())),weight)


        tempHClosest=TLorentzVector(0,0,0,0)
        tempHClosest=tempH1P4
        tempHFar=TLorentzVector(0,0,0,0)
        tempHFar=tempH2P4
        if tempH2P4.Angle(tempZ_first.Vect())<tempH1P4.Angle(tempZ_first.Vect()):
            hist_vec_reco_1D[43].Fill(degrees(tempZ_first.Angle(tempH2P4.Vect())),weight)
            hist_vec_reco_1D[44].Fill(degrees(DeltaPhi(tempZ_first.Phi(),tempH2P4.Phi())),weight)
            hist_vec_reco_1D[45].Fill(degrees(abs(tempZ_first.Theta()-tempH2P4.Theta())),weight)
            hist_vec_reco_1D[46].Fill(degrees(tempZ_first.Angle(tempH1P4.Vect())),weight)
            hist_vec_reco_1D[47].Fill(degrees(DeltaPhi(tempZ_first.Phi(),tempH1P4.Phi())),weight)
            hist_vec_reco_1D[48].Fill(degrees(abs(tempZ_first.Theta()-tempH1P4.Theta())),weight)
            hist_vec_reco_1D[49].Fill((tempZ_first+tempH2P4).E(),weight)
            hist_vec_reco_1D[50].Fill((tempZ_first+tempH2P4).E()-tempH1P4.E(),weight)

            tempHClosest=tempH2P4
            tempHFar=tempH1P4
            #checked that H2 closer to Z than H1 to Z
            if tempH2P4.Angle(tempH1P4.Vect())< tempH2P4.Angle(tempZ_first.Vect()):
                #H1 closer to H2 than H2 to Z
                hist_vec_reco_1D[51].Fill(degrees(tempH1P4.Angle(tempH2P4.Vect())),weight)
                hist_vec_reco_1D[52].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempH2P4.Phi())),weight)
                hist_vec_reco_1D[53].Fill(degrees(abs(tempH1P4.Theta()-tempH2P4.Theta())),weight)
                hist_vec_reco_1D[57].Fill(1-cos(tempH1P4.Angle(tempH2P4.Vect())),weight)
            else:
                #H2 to Z closest angle
                hist_vec_reco_1D[51].Fill(degrees(tempH2P4.Angle(tempZ_first.Vect())),weight)
                hist_vec_reco_1D[52].Fill(degrees(DeltaPhi(tempH2P4.Phi(),tempZ_first.Phi())),weight)
                hist_vec_reco_1D[53].Fill(degrees(abs(tempH2P4.Theta()-tempZ_first.Theta())),weight)
                hist_vec_reco_1D[57].Fill(1-cos(tempZ_first.Angle(tempH2P4.Vect())),weight)
                #H1 and H2 further appart than H1 to Z
                if tempH2P4.Angle(tempH1P4.Vect())> tempH1P4.Angle(tempZ_first.Vect()):
                    #H1 and H2 further appart than H1 to Z
                    hist_vec_reco_1D[54].Fill(degrees(tempH1P4.Angle(tempH2P4.Vect())),weight)
                    hist_vec_reco_1D[55].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempH2P4.Phi())),weight)
                    hist_vec_reco_1D[56].Fill(degrees(abs(tempH1P4.Theta()-tempH2P4.Theta())),weight)
                    hist_vec_reco_1D[58].Fill(1-cos(tempH1P4.Angle(tempH2P4.Vect())),weight)
                else:
                    #then H1 to Z > H1,H2 > H2,Z
                    hist_vec_reco_1D[54].Fill(degrees(tempH1P4.Angle(tempZ_first.Vect())))
                    hist_vec_reco_1D[55].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempZ_first.Phi())),weight)
                    hist_vec_reco_1D[56].Fill(degrees(abs(tempH1P4.Theta()-tempZ_first.Theta())),weight)
                    hist_vec_reco_1D[58].Fill(1-cos(tempH1P4.Angle(tempZ_first.Vect())),weight)
        else:
            #clear H1 closer to Z than H2 to Z
            hist_vec_reco_1D[43].Fill(degrees(tempZ_first.Angle(tempH1P4.Vect())),weight)
            hist_vec_reco_1D[44].Fill(degrees(DeltaPhi(tempZ_first.Phi(),tempH1P4.Phi())),weight)
            hist_vec_reco_1D[45].Fill(degrees(abs(tempZ_first.Theta()-tempH1P4.Theta())),weight)
            hist_vec_reco_1D[46].Fill(degrees(tempZ_first.Angle(tempH2P4.Vect())),weight)
            hist_vec_reco_1D[47].Fill(degrees(DeltaPhi(tempZ_first.Phi(),tempH2P4.Phi())),weight)
            hist_vec_reco_1D[48].Fill(degrees(abs(tempZ_first.Theta()-tempH2P4.Theta())),weight)
            hist_vec_reco_1D[49].Fill((tempZ_first+tempH1P4).E(),weight)
            hist_vec_reco_1D[50].Fill((tempZ_first+tempH1P4).E()-tempH2P4.E(),weight)
            if tempH2P4.Angle(tempH1P4.Vect())< tempH1P4.Angle(tempZ_first.Vect()):
                #H1 closer to H2 than H1 to Z
                hist_vec_reco_1D[51].Fill(degrees(tempH1P4.Angle(tempH2P4.Vect())))
                hist_vec_reco_1D[52].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempH2P4.Phi())),weight)
                hist_vec_reco_1D[53].Fill(degrees(abs(tempH1P4.Theta()-tempH2P4.Theta())),weight)
                hist_vec_reco_1D[57].Fill(1-cos(tempH1P4.Angle(tempH2P4.Vect())),weight)
            else:
                #H1 to Z closest angle
                hist_vec_reco_1D[51].Fill(degrees(tempH1P4.Angle(tempZ_first.Vect())))
                hist_vec_reco_1D[52].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempZ_first.Phi())),weight)
                hist_vec_reco_1D[53].Fill(degrees(abs(tempH1P4.Theta()-tempZ_first.Theta())),weight)
                hist_vec_reco_1D[57].Fill(1-cos(tempH1P4.Angle(tempZ_first.Vect())),weight)
                #H1 and H2 further appart than H2 to Z
                if tempH2P4.Angle(tempH1P4.Vect())> tempH2P4.Angle(tempZ_first.Vect()):
                    #H1 and H2 further appart than H2 to Z, and H2,Z > H1,Z
                    hist_vec_reco_1D[54].Fill(degrees(tempH1P4.Angle(tempH2P4.Vect())))
                    hist_vec_reco_1D[55].Fill(degrees(DeltaPhi(tempH1P4.Phi(),tempH2P4.Phi())),weight)
                    hist_vec_reco_1D[56].Fill(degrees(abs(tempH1P4.Theta()-tempH2P4.Theta())),weight)
                    hist_vec_reco_1D[58].Fill(1-cos(tempH1P4.Angle(tempH2P4.Vect())),weight)
                else:
                    #then H2 to Z > H1,H2 > H1,Z
                    hist_vec_reco_1D[54].Fill(degrees(tempH2P4.Angle(tempZ_first.Vect())))
                    hist_vec_reco_1D[55].Fill(degrees(DeltaPhi(tempH2P4.Phi(),tempZ_first.Phi())),weight)
                    hist_vec_reco_1D[56].Fill(degrees(abs(tempH2P4.Theta()-tempZ_first.Theta())),weight)
                    hist_vec_reco_1D[58].Fill(1-cos(tempH2P4.Angle(tempZ_first.Vect())),weight)

        genjet_E=ientry.genjet_E
        genjet_Px=ientry.genjet_Px
        genjet_Py=ientry.genjet_Py
        genjet_Pz=ientry.genjet_Pz

        temp_gj_sum=TLorentzVector(0,0,0,0)

        for ind in range(len(genjet_E)):
            temp_gj=TLorentzVector(0,0,0,0)
            temp_gj.SetPxPyPzE(genjet_Px[ind],genjet_Py[ind],genjet_Pz[ind],genjet_E[ind])
            temp_gj_sum+=temp_gj
            #print' inex/E/sumE ',ind,temp_gj.E(),temp_gj_sum.E()
        

        hist_vec_reco_1D[113].Fill(temp_gj_sum.M(),weight)
        hist_vec_reco_1D[115].Fill((temp_true_inv+temp_gj_sum).M(),weight)

        #fill now bbar only histos
        if(H1H2_decays_bbar==False):
            continue;
        if(Z_decays_qqbar==False):
            continue;
        num_count+=1

        hist_vec_reco_1D[117].Fill(tempTotEventP4.M(),weight)

        recojet_E=ientry.recojet_E
        recojet_Px=ientry.recojet_Px
        recojet_Py=ientry.recojet_Py
        recojet_Pz=ientry.recojet_Pz

        if len(recojet_E)<3:
            print 'too small recojet length',len(recojetE)<3

        recojet_vector=[]

        for ind in range(len(recojet_E)):
            temp_=TLorentzVector(0,0,0,0)
            temp_.SetPxPyPzE(recojet_Px[ind],recojet_Py[ind],recojet_Pz[ind],recojet_E[ind])
            recojet_vector.append(temp_)
        #for vect in range(len(recojet_vector)):
        #    print 'mass vect before ',vect,recojet_vector[vect].M()
        
        recojet_vector.sort(key=lambda x: x.M(), reverse=True)
        #for vect in range(len(recojet_vector)):
        #    print 'mass vect after ',vect,recojet_vector[vect].M()
        recojet_vector_combto3=orderedThreeVector(recojet_vector)
        #for vect in recojet_vector_combto3:
        #   print num_entry,'reco mass vect after combination to 3 ',vect.M(),vect.E()
        if len(genjet_E)<3:
            print 'too small genjet length',len(genjetE)<3
        genjet_vector=[]
        for ind in range(len(genjet_E)):
            temp_=TLorentzVector(0,0,0,0)
            temp_.SetPxPyPzE(genjet_Px[ind],genjet_Py[ind],genjet_Pz[ind],genjet_E[ind])
            genjet_vector.append(temp_)
        #for vect in range(len(genjet_vector)):
        #    print 'gen mass vect before ',vect,genjet_vector[vect].M()
        
        genjet_vector.sort(key=lambda x: x.M(), reverse=True)
        #for vect in range(len(genjet_vector)):
        #    print num_entry,'gen mass vect after, index ',vect,genjet_vector[vect].M(),genjet_vector[vect].E()
        genjet_vector_combto3=orderedThreeVector(genjet_vector)
        #for vect in genjet_vector_combto3:
        #   print num_entry,'gen mass vect after combination to 3 ',vect.M(),vect.E()

        if(tempH1P4.E()>tempH2P4.E()):
            hist_vec_reco_1D[77].Fill(degrees(recojet_vector[0].Angle(tempH1P4.Vect())),weight)
            hist_vec_reco_1D[78].Fill(degrees(recojet_vector[0].Angle(tempH2P4.Vect())),weight)

            hist_vec_reco_1D[80].Fill(degrees(recojet_vector[1].Angle(tempH1P4.Vect())),weight)
            hist_vec_reco_1D[81].Fill(degrees(recojet_vector[1].Angle(tempH2P4.Vect())),weight)

            hist_vec_reco_1D[83].Fill(degrees(recojet_vector[2].Angle(tempH1P4.Vect())),weight)
            hist_vec_reco_1D[84].Fill(degrees(recojet_vector[2].Angle(tempH2P4.Vect())),weight)
        else:
            hist_vec_reco_1D[77].Fill(degrees(recojet_vector[0].Angle(tempH2P4.Vect())),weight)
            hist_vec_reco_1D[78].Fill(degrees(recojet_vector[0].Angle(tempH1P4.Vect())),weight)

            hist_vec_reco_1D[80].Fill(degrees(recojet_vector[1].Angle(tempH2P4.Vect())),weight)
            hist_vec_reco_1D[81].Fill(degrees(recojet_vector[1].Angle(tempH1P4.Vect())),weight)

            hist_vec_reco_1D[83].Fill(degrees(recojet_vector[2].Angle(tempH2P4.Vect())),weight)
            hist_vec_reco_1D[84].Fill(degrees(recojet_vector[2].Angle(tempH1P4.Vect())),weight)

        hist_vec_reco_1D[79].Fill(degrees(recojet_vector[1].Angle(tempZP4.Vect())),weight)
        hist_vec_reco_1D[82].Fill(degrees(recojet_vector[1].Angle(tempZP4.Vect())),weight)
        hist_vec_reco_1D[85].Fill(degrees(recojet_vector[2].Angle(tempZP4.Vect())),weight)

        hist_vec_reco_1D[86].Fill(recojet_vector[0].M(),weight)
        hist_vec_reco_1D[87].Fill(recojet_vector[1].M(),weight)
        hist_vec_reco_1D[88].Fill(recojet_vector[2].M(),weight)
        hist_vec_reco_1D[133].Fill(recojet_vector[0].E(),weight)
        hist_vec_reco_1D[134].Fill(recojet_vector[1].E(),weight)
        hist_vec_reco_1D[135].Fill(recojet_vector[2].E(),weight)
        if len(recojet_vector)>3:
            hist_vec_reco_1D[127].Fill(degrees(min(recojet_vector[3].Angle(tempH1P4.Vect()),recojet_vector[3].Angle(tempH2P4.Vect()))),weight)
            hist_vec_reco_1D[128].Fill(degrees(recojet_vector[3].Angle(tempZP4.Vect())),weight)
            hist_vec_reco_1D[124].Fill(recojet_vector[3].M(),weight)
            hist_vec_reco_1D[136].Fill(recojet_vector[2].E(),weight)
            if len(recojet_vector)>4:
                hist_vec_reco_1D[139].Fill(degrees(min(recojet_vector[4].Angle(tempH1P4.Vect()),recojet_vector[4].Angle(tempH2P4.Vect()))),weight)
                hist_vec_reco_1D[140].Fill(degrees(recojet_vector[4].Angle(tempZP4.Vect())),weight)
                if len(recojet_vector)>5:
                    hist_vec_reco_1D[143].Fill(degrees(min(recojet_vector[5].Angle(tempH1P4.Vect()),recojet_vector[5].Angle(tempH2P4.Vect()))),weight)
                    hist_vec_reco_1D[144].Fill(degrees(recojet_vector[5].Angle(tempZP4.Vect())),weight)


        for rj in range(len(recojet_vector)):
            min_angle = 200.
            for gj in genjet_vector:
                if( degrees(recojet_vector[rj].Angle(gj.Vect()))<min_angle) :
                    min_angle=degrees(recojet_vector[rj].Angle(gj.Vect()))
            if rj==0:
                hist_vec_reco_1D[145].Fill(min_angle,weight)
            elif rj==1:
               hist_vec_reco_1D[146].Fill(min_angle,weight)
            elif rj==2:
               hist_vec_reco_1D[147].Fill(min_angle,weight)
            elif rj==3:
               hist_vec_reco_1D[148].Fill(min_angle,weight)
            elif rj==4:
               hist_vec_reco_1D[149].Fill(min_angle,weight)
            elif rj==5:
               hist_vec_reco_1D[150].Fill(min_angle,weight)

        hist_vec_reco_1D[114].Fill(temp_gj_sum.M(),weight)
        hist_vec_reco_1D[116].Fill((temp_true_inv+temp_gj_sum).M(),weight)
        hist_vec_reco_1D[154].Fill(recojet_vector_combto3[0].M(),weight)
        hist_vec_reco_1D[155].Fill(recojet_vector_combto3[1].M(),weight)
        hist_vec_reco_1D[156].Fill(recojet_vector_combto3[2].M(),weight)

        for rj in range(len(recojet_vector_combto3)):
            min_angle = 200.
            for gj in genjet_vector_combto3:
                if( degrees(recojet_vector_combto3[rj].Angle(gj.Vect()))<min_angle) :
                    min_angle=degrees(recojet_vector_combto3[rj].Angle(gj.Vect()))
            if rj==0:
                hist_vec_reco_1D[157].Fill(min_angle,weight)
            elif rj==1:
               hist_vec_reco_1D[158].Fill(min_angle,weight)
            elif rj==2:
               hist_vec_reco_1D[159].Fill(min_angle,weight)

        hist_vec_reco_1D[160].Fill(degrees(min(recojet_vector_combto3[0].Angle(tempH1P4.Vect()),recojet_vector_combto3[0].Angle(tempH2P4.Vect()))),weight)
        hist_vec_reco_1D[161].Fill(degrees(recojet_vector_combto3[0].Angle(tempZP4.Vect())),weight)
        hist_vec_reco_1D[162].Fill(degrees(min(recojet_vector_combto3[1].Angle(tempH1P4.Vect()),recojet_vector_combto3[1].Angle(tempH2P4.Vect()))),weight)
        hist_vec_reco_1D[163].Fill(degrees(recojet_vector_combto3[1].Angle(tempZP4.Vect())),weight)
        hist_vec_reco_1D[164].Fill(degrees(min(recojet_vector_combto3[2].Angle(tempH1P4.Vect()),recojet_vector_combto3[2].Angle(tempH2P4.Vect()))),weight)
        hist_vec_reco_1D[165].Fill(degrees(recojet_vector_combto3[2].Angle(tempZP4.Vect())),weight)

        if(tempH1P4.E()>tempH2P4.E()):
            hist_vec_reco_1D[101].Fill(degrees(genjet_vector[0].Angle(tempH1P4.Vect())),weight)
            hist_vec_reco_1D[102].Fill(degrees(genjet_vector[0].Angle(tempH2P4.Vect())),weight)

            hist_vec_reco_1D[107].Fill(degrees(genjet_vector[1].Angle(tempH1P4.Vect())),weight)
            hist_vec_reco_1D[108].Fill(degrees(genjet_vector[1].Angle(tempH2P4.Vect())),weight)

            hist_vec_reco_1D[83].Fill(degrees(genjet_vector[2].Angle(tempH1P4.Vect())),weight)
            hist_vec_reco_1D[84].Fill(degrees(genjet_vector[2].Angle(tempH2P4.Vect())),weight)
        else:
            hist_vec_reco_1D[101].Fill(degrees(genjet_vector[0].Angle(tempH2P4.Vect())),weight)
            hist_vec_reco_1D[102].Fill(degrees(genjet_vector[0].Angle(tempH1P4.Vect())),weight)

            hist_vec_reco_1D[104].Fill(degrees(genjet_vector[1].Angle(tempH2P4.Vect())),weight)
            hist_vec_reco_1D[105].Fill(degrees(genjet_vector[1].Angle(tempH1P4.Vect())),weight)

            hist_vec_reco_1D[107].Fill(degrees(genjet_vector[2].Angle(tempH2P4.Vect())),weight)
            hist_vec_reco_1D[108].Fill(degrees(genjet_vector[2].Angle(tempH1P4.Vect())),weight)

        hist_vec_reco_1D[103].Fill(degrees(genjet_vector[1].Angle(tempZP4.Vect())),weight)
        hist_vec_reco_1D[106].Fill(degrees(genjet_vector[1].Angle(tempZP4.Vect())),weight)
        hist_vec_reco_1D[109].Fill(degrees(genjet_vector[2].Angle(tempZP4.Vect())),weight)

        hist_vec_reco_1D[110].Fill(genjet_vector[0].M(),weight)
        hist_vec_reco_1D[111].Fill(genjet_vector[1].M(),weight)
        hist_vec_reco_1D[112].Fill(genjet_vector[2].M(),weight)

        hist_vec_reco_1D[129].Fill(genjet_vector[0].E(),weight)
        hist_vec_reco_1D[130].Fill(genjet_vector[1].E(),weight)
        hist_vec_reco_1D[131].Fill(genjet_vector[2].E(),weight)

        hist_vec_reco_1D[151].Fill(genjet_vector_combto3[0].M(),weight)
        hist_vec_reco_1D[152].Fill(genjet_vector_combto3[1].M(),weight)
        hist_vec_reco_1D[153].Fill(genjet_vector_combto3[2].M(),weight)

        if len(genjet_vector)>3:
            hist_vec_reco_1D[123].Fill(genjet_vector[3].M(),weight)
            hist_vec_reco_1D[125].Fill(degrees(min(genjet_vector[3].Angle(tempH1P4.Vect()),genjet_vector[3].Angle(tempH2P4.Vect()))),weight)
            hist_vec_reco_1D[126].Fill(degrees(genjet_vector[3].Angle(tempZP4.Vect())),weight)
            hist_vec_reco_1D[132].Fill(genjet_vector[3].E(),weight)
            if len(genjet_vector)>4:
                hist_vec_reco_1D[137].Fill(degrees(min(genjet_vector[4].Angle(tempH1P4.Vect()),genjet_vector[4].Angle(tempH2P4.Vect()))),weight)
                hist_vec_reco_1D[138].Fill(degrees(genjet_vector[4].Angle(tempZP4.Vect())),weight)
                if len(genjet_vector)>5:
                    hist_vec_reco_1D[141].Fill(degrees(min(genjet_vector[5].Angle(tempH1P4.Vect()),genjet_vector[5].Angle(tempH2P4.Vect()))),weight)
                    hist_vec_reco_1D[142].Fill(degrees(genjet_vector[5].Angle(tempZP4.Vect())),weight)

    print 'total events after all running signal histos', hist_vec_reco_1D[0].Integral(0,hist_vec_reco_1D[0].GetNbinsX()+1)," ",hist_vec_reco_1D[110].Integral(0,hist_vec_reco_1D[110].GetNbinsX()+1),num_entry,hist_vec_reco_1D[0].GetEntries(),num_count, num_total_exception
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
    
    h_sqrtS_e1_e2_effective = TH1F("h_sqrtS_e1_e2_effective","", n_bins_high, lim_energy_low,lim_energy_high);
    h_sqrtS_e1_e2_effective_HHZ_bbbbqq = TH1F("h_sqrtS_e1_e2_effective_HHZ_bbbbqq","", n_bins_high, lim_energy_low,lim_energy_high);
    h_sqrtS_HHZ = TH1F("h_sqrtS_HHZ","", n_bins_high, lim_energy_low,lim_energy_high);
    h_H1_H2_E_sum = TH1F("h_H1_H2_E_sum","", n_bins_high, lim_energy_low,lim_energy_high);
    h_H1_E = TH1F("h_H1_E","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_H2_E = TH1F("h_H2_E","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_Z_E = TH1F("h_Z_E","", n_bins_high, lim_energy_low,0.5*lim_energy_high);

    h_H1_P = TH1F("h_H1_P","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_H2_P = TH1F("h_H2_P","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_Z_P = TH1F("h_Z_P","", n_bins_high, lim_energy_low,0.5*lim_energy_high);

    h_sqrtS_gj_sum_all = TH1F("h_sqrtS_gj_sum_all","", n_bins_high, lim_energy_low,lim_energy_high);
    h_sqrtS_gj_sum_HHZ_bbbbqq = TH1F("h_sqrtS_gj_sum_HHZ_bbbbqq ","", n_bins_high, lim_energy_low,lim_energy_high);
    h_sqrtS_gj_sum_plus_InvVec_all = TH1F("h_sqrtS_gj_sum_plus_InvVec_all","", n_bins_high, lim_energy_low,lim_energy_high);
    h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq = TH1F("h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq ","", n_bins_high, lim_energy_low,lim_energy_high);

    lim_mass_low=0;
    lim_mass_high=500;
    
    h_H1_H2_mass = TH1F("h_H1_H2_mass","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    
    lim_dalpha_low=0.;
    lim_dalpha_high=180.;

    h_theta_H1 = TH1F("h_theta_H1","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_theta_H2 = TH1F("h_theta_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_theta_Z = TH1F("h_theta_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    
    h_dalpha_H1_H2_comb_vs_Z = TH1F("h_dalpha_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dphi_H1_H2_comb_vs_Z = TH1F("h_dphi_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dtheta_H1_H2_comb_vs_Z = TH1F("h_dtheta_H1_H2_comb_vs_Z","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    
    h_dalpha_H1_H2 = TH1F("h_dalpha_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dphi_H1_H2 = TH1F("h_dphi_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dtheta_H1_H2 = TH1F("h_dtheta_H1_H2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
   
    h_dalpha_H_Z_min = TH1F("h_dalpha_H_Z_min","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dphi_H_Z_min = TH1F("h_dphi_H_Z_min","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dtheta_H_Z_min = TH1F("h_dtheta_H_Z_min","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dalpha_H_Z_max = TH1F("h_dalpha_H_Z_max","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dphi_H_Z_max = TH1F("h_dphi_H_Z_max","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dtheta_H_Z_max = TH1F("h_dtheta_H_Z_max","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_Esum_Z_H_minAngle = TH1F("Esum_Z_H_minAngle","", n_bins_high, lim_energy_low,lim_energy_high);
    h_delta_Esum_Z_H_minAngle_min_E_H_Z_maxAngle = TH1F("h_delta_Esum_Z_H_minAngle_min_E_H_Z_maxAngle","", n_bins_high, -0.5*lim_energy_high,0.5*lim_energy_high);

    lim_dalpha_qqbar_low=0.;
    lim_dalpha_qqbar_high=45.;
    
    h_dalpha_H1_bbbar = TH1F("h_dalpha_H1_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dphi_H1_bbbar = TH1F("h_dphi_H1_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dtheta_H1_bbbar = TH1F("h_dtheta_H1_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    
    h_dalpha_H2_bbbar = TH1F("h_dalpha_H2_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dphi_H2_bbbar = TH1F("h_dphi_H2_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dtheta_H2_bbbar = TH1F("h_dtheta_H2_bbbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    
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
    
    h_dalpha_max_allH_bbbar = TH1F("h_dalpha_max_allH_bbbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
    h_dphi_max_allH_bbbar = TH1F("h_dphi_max_allH_bbbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
    h_dtheta_max_allH_bbbar = TH1F("h_dtheta_max_allH_bbbar","", n_bins_high, lim_dalpha_qqbar_allH_low,lim_dalpha_qqbar_allH_high);
    
    h_dalpha_allH_diboson_qqbar = TH1F("h_dalpha_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dphi_allH_diboson_qqbar = TH1F("h_dphi_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    h_dtheta_allH_diboson_qqbar = TH1F("h_dtheta_allH_diboson_qqbar","", n_bins_high, lim_dalpha_qqbar_low,lim_dalpha_qqbar_high);
    
    lim_H1_E_over_allH_E_low=0.5;
    lim_H1_E_over_allH_E_high=1.0;
    
    h_E_H1_over_E_allH = TH1F("h_E_H1_over_E_allH","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
    h_E_H1_over_E_allH_bbbar = TH1F("h_E_H1_over_E_allH_bbbar","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
    h_E_H1_over_E_allH_qqbar = TH1F("h_E_H1_over_E_allH_qqbar","", n_bins_high, lim_H1_E_over_allH_E_low,lim_H1_E_over_allH_E_high);
    
    lim_d_ij_qqbar_low=0.;
    lim_d_ij_qqbar_high=5.0;

    h_d_ij_H1_bbbar = TH1F("h_d_ij_H1_bbbar","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
    h_d_ij_H2_bbbar = TH1F("h_d_ij_H2_bbbar","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);

    h_d_ij_H1_H2 = TH1F("h_d_ij_H1_H2","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
    h_d_ij_H1_Z = TH1F("h_d_ij_H1_Z","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
    h_d_ij_H2_Z = TH1F("h_d_ij_H2_Z","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
    h_d_ij_H1_qqbar = TH1F("h_d_ij_H1_qqbar","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
    h_d_ij_H2_qqbar = TH1F("h_d_ij_H2_qqbar","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
    h_d_ij_Z_qqbar = TH1F("h_d_ij_Z_qqbar","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);

    h_dalpha_min_V1_V2 = TH1F("h_dalpha_min_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dphi_min_V1_V2 = TH1F("h_dphi_min_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dtheta_min_V1_V2 = TH1F("h_dtheta_min_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

    h_dalpha_max_V1_V2 = TH1F("h_dalpha_max_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dphi_max_V1_V2 = TH1F("h_dphi_max_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);
    h_dtheta_max_V1_V2 = TH1F("h_dtheta_max_V1_V2","", n_bins_high, lim_dalpha_low,lim_dalpha_high);

    h_d_ij_min_V1_V2 = TH1F("h_d_ij_min_V1_V2","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);
    h_d_ij_max_V1_V2 = TH1F("h_d_ij_max_V1_V2","", n_bins_high, lim_d_ij_qqbar_low,lim_d_ij_qqbar_high);

    lim_angle_low=0
    lim_angle_high=180

    h_dalpha_H1_rj_HHZ_bbbbqq = TH1F("h_dalpha_H1_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dphi_H1_rj_HHZ_bbbbqq = TH1F("h_dphi_H1_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dtheta_H1_rj_HHZ_bbbbqq = TH1F("h_dtheta_H1_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_H2_rj_HHZ_bbbbqq = TH1F("h_dalpha_H2_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dphi_H2_rj_HHZ_bbbbqq = TH1F("h_dphi_H2_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dtheta_H2_rj_HHZ_bbbbqq = TH1F("h_dtheta_H2_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_Z_rj_HHZ_bbbbqq = TH1F("h_dalpha_Z_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dphi_Z_rj_HHZ_bbbbqq = TH1F("h_dphi_Z_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dtheta_Z_rj_HHZ_bbbbqq = TH1F("h_dtheta_Z_rj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_H1_gj_HHZ_bbbbqq = TH1F("h_dalpha_H1_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dphi_H1_gj_HHZ_bbbbqq = TH1F("h_dphi_H1_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dtheta_H1_gj_HHZ_bbbbqq = TH1F("h_dtheta_H1_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_H2_gj_HHZ_bbbbqq = TH1F("h_dalpha_H2_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dphi_H2_gj_HHZ_bbbbqq = TH1F("h_dphi_H2_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dtheta_H2_gj_HHZ_bbbbqq = TH1F("h_dtheta_H2_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_Z_gj_HHZ_bbbbqq = TH1F("h_dalpha_Z_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dphi_Z_gj_HHZ_bbbbqq = TH1F("h_dphi_Z_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dtheta_Z_gj_HHZ_bbbbqq = TH1F("h_dtheta_Z_gj_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    lim_mass_low=0
    lim_mass_high=250
    h_mass_Z_rj_HHZ_bbbbqq = TH1F("h_mass_Z_rj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_H1_rj_HHZ_bbbbqq = TH1F("h_mass_H1_rj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_H2_rj_HHZ_bbbbqq = TH1F("h_mass_H2_rj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);

    h_dalpha_rj1_H1_HHZ_bbbbqq = TH1F("h_dalpha_rj1_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj1_H2_HHZ_bbbbqq = TH1F("h_dalpha_rj1_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj1_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj1_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_rj2_H1_HHZ_bbbbqq = TH1F("h_dalpha_rj2_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj2_H2_HHZ_bbbbqq = TH1F("h_dalpha_rj2_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj2_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj2_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_rj3_H1_HHZ_bbbbqq = TH1F("h_dalpha_rj3_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj3_H2_HHZ_bbbbqq = TH1F("h_dalpha_rj3_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj3_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj3_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_mass_rj1_HHZ_bbbbqq = TH1F("h_mass_rj1_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_rj2_HHZ_bbbbqq = TH1F("h_mass_rj2_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_rj3_HHZ_bbbbqq = TH1F("h_mass_rj3_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);

    h_mass_Z_gj_HHZ_bbbbqq = TH1F("h_mass_Z_gj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_H1_gj_HHZ_bbbbqq = TH1F("h_mass_H1_gj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_H2_gj_HHZ_bbbbqq = TH1F("h_mass_H2_gj_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);

    h_dalpha_gj1_H1_HHZ_bbbbqq = TH1F("h_dalpha_gj1_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_gj1_H2_HHZ_bbbbqq = TH1F("h_dalpha_gj1_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_gj1_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj1_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_gj2_H1_HHZ_bbbbqq = TH1F("h_dalpha_gj2_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_gj2_H2_HHZ_bbbbqq = TH1F("h_dalpha_gj2_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_gj2_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj2_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_gj3_H1_HHZ_bbbbqq = TH1F("h_dalpha_gj3_H1_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_gj3_H2_HHZ_bbbbqq = TH1F("h_dalpha_gj3_H2_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_gj3_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj3_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_mass_gj1_HHZ_bbbbqq = TH1F("h_mass_gj1_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_gj2_HHZ_bbbbqq = TH1F("h_mass_gj2_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_gj3_HHZ_bbbbqq = TH1F("h_mass_gj3_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);

    h_mass_gj4_HHZ_bbbbqq = TH1F("h_mass_gj4_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_rj4_HHZ_bbbbqq = TH1F("h_mass_rj4_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);

    h_dalpha_gj4_Hs_HHZ_bbbbqq = TH1F("h_dalpha_gj4_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_gj4_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj4_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_rj4_Hs_HHZ_bbbbqq = TH1F("h_dalpha_rj4_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj4_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj4_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_gj1_E_bbbbqq = TH1F("h_gj1_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_gj2_E_bbbbqq = TH1F("h_gj2_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_gj3_E_bbbbqq = TH1F("h_gj3_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_gj4_E_bbbbqq = TH1F("h_gj4_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);

    h_rj1_E_bbbbqq = TH1F("h_rj1_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_rj2_E_bbbbqq = TH1F("h_rj2_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_rj3_E_bbbbqq = TH1F("h_rj3_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);
    h_rj4_E_bbbbqq = TH1F("h_rj4_E_bbbbqq","", n_bins_high, lim_energy_low,0.5*lim_energy_high);

    h_dalpha_gj5_Hs_HHZ_bbbbqq = TH1F("h_dalpha_gj5_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_gj5_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj5_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_rj5_Hs_HHZ_bbbbqq = TH1F("h_dalpha_rj5_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj5_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj5_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_gj6_Hs_HHZ_bbbbqq = TH1F("h_dalpha_gj6_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_gj6_Z_HHZ_bbbbqq = TH1F("h_dalpha_gj6_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_rj6_Hs_HHZ_bbbbqq = TH1F("h_dalpha_rj6_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj6_Z_HHZ_bbbbqq = TH1F("h_dalpha_rj6_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_rj1_gjs_bbbbqq = TH1F("h_dalpha_rj1_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj2_gjs_bbbbqq = TH1F("h_dalpha_rj2_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj3_gjs_bbbbqq = TH1F("h_dalpha_rj3_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj4_gjs_bbbbqq = TH1F("h_dalpha_rj4_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj5_gjs_bbbbqq = TH1F("h_dalpha_rj5_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_rj6_gjs_bbbbqq = TH1F("h_dalpha_rj6_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_mass_comb_gj1_HHZ_bbbbqq = TH1F("h_mass_comb_gj1_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_comb_gj2_HHZ_bbbbqq = TH1F("h_mass_comb_gj2_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_comb_gj3_HHZ_bbbbqq = TH1F("h_mass_comb_gj3_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);

    h_mass_comb_rj1_HHZ_bbbbqq = TH1F("h_mass_comb_rj1_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_comb_rj2_HHZ_bbbbqq = TH1F("h_mass_comb_rj2_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);
    h_mass_comb_rj3_HHZ_bbbbqq = TH1F("h_mass_comb_rj3_HHZ_bbbbqq","", n_bins_high, lim_mass_low,lim_mass_high);

    h_dalpha_comb_rj1_gjs_bbbbqq = TH1F("h_dalpha_comb_rj1_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_comb_rj2_gjs_bbbbqq = TH1F("h_dalpha_comb_rj2_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_comb_rj3_gjs_bbbbqq = TH1F("h_dalpha_comb_rj3_gjs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    h_dalpha_comb_rj1_Hs_bbbbqq = TH1F("h_dalpha_comb_rj1_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_comb_rj1_Z_bbbbqq = TH1F("h_dalpha_comb_rj1_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_comb_rj2_Hs_bbbbqq = TH1F("h_dalpha_comb_rj2_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_comb_rj2_Z_bbbbqq = TH1F("h_dalpha_comb_rj2_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_comb_rj3_Hs_bbbbqq = TH1F("h_dalpha_comb_rj3_Hs_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);
    h_dalpha_comb_rj3_Z_bbbbqq = TH1F("h_dalpha_comb_rj3_Z_HHZ_bbbbqq","", n_bins_high, lim_angle_low,lim_angle_high);

    hist_vec_HHZ_parton_list=[]
    
    hist_vec_HHZ_parton_list.append(h_sqrtS_e1_e2_effective)
    hist_vec_HHZ_parton_list.append(h_H1_H2_E_sum)
    hist_vec_HHZ_parton_list.append(h_H1_E)
    hist_vec_HHZ_parton_list.append(h_H2_E)
    hist_vec_HHZ_parton_list.append(h_Z_E)
    #4 done here
    hist_vec_HHZ_parton_list.append(h_H1_H2_mass)
    hist_vec_HHZ_parton_list.append(h_dalpha_H1_H2_comb_vs_Z)
    hist_vec_HHZ_parton_list.append(h_dphi_H1_H2_comb_vs_Z)
    hist_vec_HHZ_parton_list.append(h_dtheta_H1_H2_comb_vs_Z)
    hist_vec_HHZ_parton_list.append(h_dalpha_H1_H2)
    #9 done here
    hist_vec_HHZ_parton_list.append(h_dphi_H1_H2)
    hist_vec_HHZ_parton_list.append(h_dtheta_H1_H2)
    hist_vec_HHZ_parton_list.append(h_dalpha_H1_bbbar)
    hist_vec_HHZ_parton_list.append(h_dphi_H1_bbbar)
    hist_vec_HHZ_parton_list.append(h_dtheta_H1_bbbar)
    #14 done here
    hist_vec_HHZ_parton_list.append(h_dalpha_H2_bbbar)
    hist_vec_HHZ_parton_list.append(h_dphi_H2_bbbar)
    hist_vec_HHZ_parton_list.append(h_dtheta_H2_bbbar)
    hist_vec_HHZ_parton_list.append(h_dalpha_H1_qqbar)
    hist_vec_HHZ_parton_list.append(h_dphi_H1_qqbar)
    #19 done here
    hist_vec_HHZ_parton_list.append(h_dtheta_H1_qqbar)
    hist_vec_HHZ_parton_list.append(h_dalpha_H2_qqbar)
    hist_vec_HHZ_parton_list.append(h_dphi_H2_qqbar)
    hist_vec_HHZ_parton_list.append(h_dtheta_H2_qqbar)
    hist_vec_HHZ_parton_list.append(h_dalpha_Z_qqbar)
    #24 done here
    hist_vec_HHZ_parton_list.append(h_dphi_Z_qqbar)
    hist_vec_HHZ_parton_list.append(h_dtheta_Z_qqbar)
    hist_vec_HHZ_parton_list.append(h_dalpha_max_allH_qqbar)
    hist_vec_HHZ_parton_list.append(h_dphi_max_allH_qqbar)
    hist_vec_HHZ_parton_list.append(h_dtheta_max_allH_qqbar)
    #29 done here
    hist_vec_HHZ_parton_list.append(h_dalpha_max_allH_bbbar)
    hist_vec_HHZ_parton_list.append(h_dphi_max_allH_bbbar)
    hist_vec_HHZ_parton_list.append(h_dtheta_max_allH_bbbar)
    hist_vec_HHZ_parton_list.append(h_dalpha_allH_diboson_qqbar)
    hist_vec_HHZ_parton_list.append(h_dphi_allH_diboson_qqbar)
    #34 done here
    hist_vec_HHZ_parton_list.append(h_dtheta_allH_diboson_qqbar)
    hist_vec_HHZ_parton_list.append(h_E_H1_over_E_allH)
    hist_vec_HHZ_parton_list.append(h_E_H1_over_E_allH_bbbar)
    hist_vec_HHZ_parton_list.append(h_E_H1_over_E_allH_qqbar)
    hist_vec_HHZ_parton_list.append(h_d_ij_H1_H2)
    #39 done up to here
    hist_vec_HHZ_parton_list.append(h_d_ij_H1_qqbar)
    hist_vec_HHZ_parton_list.append(h_d_ij_H2_qqbar)
    hist_vec_HHZ_parton_list.append(h_sqrtS_HHZ)
    hist_vec_HHZ_parton_list.append(h_dalpha_H_Z_min)
    hist_vec_HHZ_parton_list.append(h_dphi_H_Z_min)
    #44 done up to here
    hist_vec_HHZ_parton_list.append(h_dtheta_H_Z_min)
    hist_vec_HHZ_parton_list.append(h_dalpha_H_Z_max)
    hist_vec_HHZ_parton_list.append(h_dphi_H_Z_max)
    hist_vec_HHZ_parton_list.append(h_dtheta_H_Z_max)
    hist_vec_HHZ_parton_list.append(h_Esum_Z_H_minAngle)
    #49 done up to here --> dphi and dtheta for maximally angularly separated, but still check if it is more on phi or theta
    hist_vec_HHZ_parton_list.append(h_delta_Esum_Z_H_minAngle_min_E_H_Z_maxAngle)
    hist_vec_HHZ_parton_list.append(h_dalpha_min_V1_V2)
    hist_vec_HHZ_parton_list.append(h_dphi_min_V1_V2)
    hist_vec_HHZ_parton_list.append(h_dtheta_min_V1_V2)
    hist_vec_HHZ_parton_list.append(h_dalpha_max_V1_V2)
    #54 done up to now
    hist_vec_HHZ_parton_list.append(h_dphi_max_V1_V2)
    hist_vec_HHZ_parton_list.append(h_dtheta_max_V1_V2)
    hist_vec_HHZ_parton_list.append(h_d_ij_min_V1_V2)
    hist_vec_HHZ_parton_list.append(h_d_ij_max_V1_V2)
    hist_vec_HHZ_parton_list.append(h_d_ij_Z_qqbar)
    #59 done up to now
    hist_vec_HHZ_parton_list.append(h_d_ij_H1_bbbar)
    hist_vec_HHZ_parton_list.append(h_d_ij_H2_bbbar)
    hist_vec_HHZ_parton_list.append(h_theta_H1)
    hist_vec_HHZ_parton_list.append(h_theta_H2)
    hist_vec_HHZ_parton_list.append(h_theta_Z)
    #64 done up to now

    hist_vec_HHZ_parton_list.append(h_dalpha_H1_rj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dphi_H1_rj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dtheta_H1_rj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_H2_rj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dphi_H2_rj_HHZ_bbbbqq)
    #69 done up to now
    hist_vec_HHZ_parton_list.append(h_dtheta_H2_rj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_Z_rj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dphi_Z_rj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dtheta_Z_rj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_Z_rj_HHZ_bbbbqq)
    #74 done up to now
    hist_vec_HHZ_parton_list.append(h_mass_H1_rj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_H2_rj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj1_H1_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj1_H2_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj1_Z_HHZ_bbbbqq) 
    #79 done up to now
    hist_vec_HHZ_parton_list.append(h_dalpha_rj2_H1_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj2_H2_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj2_Z_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj3_H1_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj3_H2_HHZ_bbbbqq)
    #84 done up to now
    hist_vec_HHZ_parton_list.append(h_dalpha_rj3_Z_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_rj1_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_rj2_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_rj3_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_H1_gj_HHZ_bbbbqq)
    #89 up to here
    hist_vec_HHZ_parton_list.append(h_dphi_H1_gj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dtheta_H1_gj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_H2_gj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dphi_H2_gj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dtheta_H2_gj_HHZ_bbbbqq)
    #94 done up to now    
    hist_vec_HHZ_parton_list.append(h_dalpha_Z_gj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dphi_Z_gj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dtheta_Z_gj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_Z_gj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_H1_gj_HHZ_bbbbqq)
    #99 done up to now   
    hist_vec_HHZ_parton_list.append(h_mass_H2_gj_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj1_H1_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj1_H2_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj1_Z_HHZ_bbbbqq) 
    hist_vec_HHZ_parton_list.append(h_dalpha_gj2_H1_HHZ_bbbbqq)
    #104 done up to now
    hist_vec_HHZ_parton_list.append(h_dalpha_gj2_H2_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj2_Z_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj3_H1_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj3_H2_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj3_Z_HHZ_bbbbqq)
    #109 done up to now
    hist_vec_HHZ_parton_list.append(h_mass_gj1_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_gj2_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_gj3_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_sqrtS_gj_sum_all)
    hist_vec_HHZ_parton_list.append(h_sqrtS_gj_sum_HHZ_bbbbqq)
    #114 done until now
    hist_vec_HHZ_parton_list.append(h_sqrtS_gj_sum_plus_InvVec_all)
    hist_vec_HHZ_parton_list.append(h_sqrtS_gj_sum_plus_InvVec_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_sqrtS_e1_e2_effective_HHZ_bbbbqq)

    hist_vec_HHZ_parton_list.append(h_H1_P)
    hist_vec_HHZ_parton_list.append(h_H2_P)
    #119 done until now
    hist_vec_HHZ_parton_list.append(h_Z_P)
    hist_vec_HHZ_parton_list.append(h_d_ij_H1_Z)
    hist_vec_HHZ_parton_list.append(h_d_ij_H2_Z)

    hist_vec_HHZ_parton_list.append(h_mass_gj4_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_rj4_HHZ_bbbbqq)
    #124 done until now
    hist_vec_HHZ_parton_list.append(h_dalpha_gj4_Hs_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj4_Z_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj4_Hs_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj4_Z_HHZ_bbbbqq)

    hist_vec_HHZ_parton_list.append(h_gj1_E_bbbbqq)
    #129 done until now
    hist_vec_HHZ_parton_list.append(h_gj2_E_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_gj3_E_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_gj4_E_bbbbqq)

    hist_vec_HHZ_parton_list.append(h_rj1_E_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_rj2_E_bbbbqq)
    #134 done until now
    hist_vec_HHZ_parton_list.append(h_rj3_E_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_rj4_E_bbbbqq)

    hist_vec_HHZ_parton_list.append(h_dalpha_gj5_Hs_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj5_Z_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj5_Hs_HHZ_bbbbqq)
    #139 done
    hist_vec_HHZ_parton_list.append(h_dalpha_rj5_Z_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj6_Hs_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_gj6_Z_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj6_Hs_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj6_Z_HHZ_bbbbqq)
    #144 done 
    hist_vec_HHZ_parton_list.append(h_dalpha_rj1_gjs_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj2_gjs_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj3_gjs_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj4_gjs_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_rj5_gjs_bbbbqq)
    #149 done
    hist_vec_HHZ_parton_list.append(h_dalpha_rj6_gjs_bbbbqq)

    hist_vec_HHZ_parton_list.append(h_mass_comb_gj1_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_comb_gj2_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_comb_gj3_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_comb_rj1_HHZ_bbbbqq)
    #154 DONE
    hist_vec_HHZ_parton_list.append(h_mass_comb_rj2_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_mass_comb_rj3_HHZ_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj1_gjs_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj2_gjs_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj3_gjs_bbbbqq)
    #159 done
    hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj1_Hs_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj1_Z_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj2_Hs_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj2_Z_bbbbqq)
    hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj3_Hs_bbbbqq)
    #165 done
    hist_vec_HHZ_parton_list.append(h_dalpha_comb_rj3_Z_bbbbqq)



    for hist in hist_vec_HHZ_parton_list:
        hist.Sumw2()


    print 'length of hist list ',len(hist_vec_HHZ_parton_list)
    fill_HHZ_histograms(input_file_,xsec_,hist_vec_HHZ_parton_list,lumi)
    
    CLICdpStyle()


    file_histogram.Write()
    file_histogram.Close()
  
    return None

def process_files():

    lumi_=4000.
    print 'start processing of files'

#9,8,5
    cross_section_= 6.06e-02
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC14_NJet6_with_yij/polm80/HHZStudy_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14.root"
    #input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC10_NJet3_noIsoPhLep/polm80/HHZStudy_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14.root"
    #eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC15_NJet3/polm80/HHZStudy_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_yij/polm80/test_hhz_14343_polm80_3TeV_wO_CLIC_o3_v14_partonlevelOnly_noIsoPh_TrueMCJets.root"  
    process_event(final_histo_name_,input_file_name_,cross_section_,lumi_)
    print 'finished file', final_histo_name_

    lumi_=1000.
    cross_section_= 4.23e-02
    input_file_name_="/eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC10_NJet3_noIsoPhLep/polp80/HHZStudy_hhz_14344_polp80_3TeV_wO_CLIC_o3_v14.root"
    #eos/user/w/weberma2/data/HHZAnalyzerFiles/VLC15_NJet3/polp80/HHZStudy_hhz_14344_polp80_3TeV_wO_CLIC_o3_v14.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC10_NJets3/polp80/test_hhz_14344_polp80_3TeV_wO_CLIC_o3_v14_partonlevelOnly_noIsoPh.root"  
    #process_event(final_histo_name_,input_file_name_,cross_section_,lumi_)
    print 'finished file', final_histo_name_

    return None

process_files()



