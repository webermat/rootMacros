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


    

def process_event(i_final_histo_name_,i_input_file_name_,i_weight_file):
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

    s_jet1_C2_beta1 = array('f',[0])
    reader.AddVariable('jet1_C2_beta1',s_jet1_C2_beta1)
    s_jet2_C2_beta1 = array('f',[0])
    reader.AddVariable('jet2_C2_beta1',s_jet2_C2_beta1)
    s_jet1_tau21 = array('f',[0])
    reader.AddVariable('jet1_tau21',s_jet1_tau21)
    s_jet2_tau21 = array('f',[0])
    reader.AddVariable('jet2_tau21',s_jet2_tau21)

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
    tree.SetBranchAddress('jet1_C2_beta1',s_jet1_C2_beta1)
    tree.SetBranchAddress('jet2_C2_beta1',s_jet2_C2_beta1)
    tree.SetBranchAddress('jet1_tau21',s_jet1_tau21)
    tree.SetBranchAddress('jet2_tau21',s_jet2_tau21)
    v_eventWeight = array('f',[0])
    tree.SetBranchAddress('eventWeight',v_eventWeight)
    reader.BookMVA('BDT',TString(i_weight_file))


    n_bins_high=100
    lim_BDT_low=-1.0
    lim_BDT_high=1.0
    lim_mass_low=0
    lim_mass_high=200
    lim_theta_low=0
    lim_theta_high=180

    lim_theta_low=0
    lim_theta_high=180

    lim_dtheta_low=-180
    lim_dtheta_high=180

    lim_BTag_low=0;
    lim_BTag_high=0;

    lim_D2_low=0;
    lim_D2_high=10;

    print 'do i get here maybe'

    BDT_cuts=[]
    #print 'do i get here maybe 1'
    BDT_cuts.append(-0.20)
    #print 'do i get here maybe 2'
    BDT_cuts.append(-0.15)
    BDT_cuts.append(-0.10)
    BDT_cuts.append(-0.05)
    BDT_cuts.append(0.00)
    BDT_cuts.append(0.05)
    BDT_cuts.append(0.10)
    BDT_cuts.append(0.15)
    BDT_cuts.append(0.20)
    BDT_cuts.append(0.25)
    BDT_cuts.append(0.30)
    BDT_cuts.append(0.35)
    BDT_cuts.append(0.40)
    BDT_cuts.append(0.45)
    BDT_cuts.append(0.50)
    BDT_cuts.append(0.55)
    BDT_cuts.append(0.60)
    BDT_cuts.append(0.65)
    BDT_cuts.append(0.70)
    BDT_cuts.append(0.75)
    BDT_cuts.append(0.80)
    BDT_cuts.append(0.85)
    BDT_cuts.append(0.90)
    BDT_cuts.append(0.95)
    BDT_cuts.append(1.00)
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
        h_signal_background_1D_hist_list.append(h_BDT_output_Eff)
        for hist in h_signal_background_1D_hist_list:
            hist.Sumw2()
            
            
        num_entry=-1
        BDTCut=0.50
        for ientry in range(tree.GetEntries()):
            num_entry+=1;
            if num_entry%(int(tree.GetEntries()/5.)) == 0:
                print "sig BG in entry ",num_entry,bdt_value
                    
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

        print 'at end of file processing of ',i_final_histo_name_,BDTCut,h_jet1_mass.Integral()
                
        fileout.Write()
    fileout.Close()

 
    return None

def process_files():


    #files_weights_='/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/dataset/filesPolp80GiniIndexNoNormSkipNormalization/dataset/weights/TMVAClassification_BDT.weights.xml'
    files_weights_='/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage1_25NCuts20_qqqq2TeV_allVar/dataset/weights/TMVAClassification_BDT.weights.xml'


    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"


    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/MVTrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage1_25NCuts20.root"  

 

    process_event(final_histo_name_,input_file_,files_weights_)

    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/MVTrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage1_25NCuts20.root" 
    process_event(final_histo_name_,input_file_,files_weights_)

 
    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/test_ee_qqqq_mqqqq_2TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/MVTrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage1_25NCuts20.root" 
    process_event(final_histo_name_,input_file_,files_weights_)

    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/test_ee_qqqqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/MVTrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage1_25NCuts20.root" 
    process_event(final_histo_name_,input_file_,files_weights_)


    files_weights_='/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage1_25NCuts20_qqqq2TeV_allVar/dataset/weights/TMVAClassification_BDT.weights.xml'


    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/MVTrainingReader_hzqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage1_25NCuts20.root"  
    process_event(final_histo_name_,input_file_,files_weights_)

    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/MVTrainingReader_ee_qq_mqq_1TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage1_25NCuts20.root" 
    process_event(final_histo_name_,input_file_,files_weights_)

 
    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/test_ee_qqqq_mqqqq_2TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/MVTrainingReader_ee_qqqq_mqqqq_2TeV_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage1_25NCuts20.root" 
    process_event(final_histo_name_,input_file_,files_weights_)

    input_file_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/test_ee_qqqqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/MVTrainingReader_ee_qqqqqq_histofiles_BDT_GiniIndexNormNumEventsMaxDepth3NTrees400Shrinkage1_25NCuts20.root" 
    process_event(final_histo_name_,input_file_,files_weights_)






    return None

process_files()



