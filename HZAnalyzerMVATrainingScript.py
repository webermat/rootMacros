from ROOT import gROOT, TCanvas, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor,TMVA,TCut
import ROOT as root
from math import cos, sin, pi, degrees, radians, pow, sqrt,acos
from array import array
import os



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


    

def process_event(i_final_histo_name_,i_input_signal_file_name_,i_input_bkg_file_name_):
    print "at start of process event"     
    root.TMVA.Tools.Instance()
    #root.TMVA.PyMethodBase.PyInitialize()

    _output_training_name_=i_final_histo_name_
    _input_signal_name_=i_input_signal_file_name_
    _input_bkg_ee_qq_name_=i_input_bkg_file_name_[0]
    _input_bkg_ee_qqqq_name_=i_input_bkg_file_name_[1]
    _input_bkg_ee_qqqqqq_name_=i_input_bkg_file_name_[2]


    fileout = root.TFile(_output_training_name_,"RECREATE")

    input_file_signal_=root.TFile.Open(_input_signal_name_)
    input_file_bkg_ee_qq_=root.TFile.Open(_input_bkg_ee_qq_name_)
    input_file_bkg_ee_qqqq_=root.TFile.Open(_input_bkg_ee_qqqq_name_)
    input_file_bkg_ee_qqqqqq_=root.TFile.Open(_input_bkg_ee_qqqqqq_name_)

    tree_sig = input_file_signal_.Get("MVATrainingVariables")
    tree_bkg_ee_qq = input_file_bkg_ee_qq_.Get("MVATrainingVariables")
    tree_bkg_ee_qqqq = input_file_bkg_ee_qqqq_.Get("MVATrainingVariables")
    tree_bkg_ee_qqqqqq = input_file_bkg_ee_qqqqqq_.Get("MVATrainingVariables")

    #flags V is if verbose set or not (typically false) , if it should be silent or not (silent== no output from mva, default false) ,transformations can be  
    #I;D;P;U;G,D,for identity, decorrelation,PCA,Uniform   and Gaussianisation  followed  by  decorrelation transformations
    #analysistype can be Classification,Regression,Multiclass or Auto (default Auto)
    factory = TMVA.Factory("TMVAClassification", fileout,
                            ":".join([
                                "!V",
                                "!Silent",
                                "Color",
                                "DrawProgressBar",
                                "Transformations=I;D;P;G,D",
                                "AnalysisType=Classification"]
                                     ))


    dataloader = TMVA.DataLoader('dataset')
    dataloader.SetWeightExpression("eventWeight")
    dataloader.SetSignalWeightExpression("eventWeight")
    dataloader.SetBackgroundWeightExpression("eventWeight")
    dataloader.AddVariable("jet1_mass", 'F')
    dataloader.AddVariable("jet2_mass", 'F')
    dataloader.AddVariable("jet1_theta", 'F')
    dataloader.AddVariable("jet2_theta", 'F')
    dataloader.AddVariable("jet1_D2_beta1", 'F')
    dataloader.AddVariable("jet2_D2_beta1", 'F')
    dataloader.AddVariable("jet1_BTag_rfj_BTagMax", 'F')
    dataloader.AddSpectator("jet1_CTag_rfj_CTagMax", 'F')
    #dataloader.AddVariable("jet1_BTag_rfj_BTagMax", 'F')
    dataloader.AddVariable("jet1_C2_beta1", 'F')
    dataloader.AddVariable("jet2_C2_beta1", 'F')
    dataloader.AddSpectator("jet1_tau21", 'F')
    dataloader.AddSpectator("jet2_tau21", 'F')
    dataloader.AddSpectator("costheta1_for_Atheta1", 'F')
    dataloader.AddSpectator("costheta2_for_Atheta1theta2",'F')
    dataloader.AddSpectator("phi_for_Aphis",'F')
    dataloader.AddVariable("reco_y32",'F')
    dataloader.AddSpectator("reco_y43",'F')
    dataloader.AddSpectator("dphi_j1j2",'F')
    dataloader.AddSpectator("jet1_d21",'F')
    dataloader.AddSpectator("jet1_d32",'F')
    dataloader.AddSpectator("jet1_d43",'F')
    dataloader.AddSpectator("jet2_d21",'F')
    dataloader.AddSpectator("jet2_d32",'F')
    dataloader.AddSpectator("jet2_d43",'F')
 
    dataloader.AddVariable("jet1_C3_beta1", 'F')
    dataloader.AddVariable("jet2_C3_beta1", 'F')
    dataloader.AddSpectator("jet1_N2_beta1", 'F')
    dataloader.AddSpectator("jet2_N2_beta1", 'F')
    dataloader.AddSpectator("jet1_N3_beta1", 'F')
    dataloader.AddSpectator("jet2_N3_beta1", 'F')

    def_weight=1

    dataloader.AddSignalTree(tree_sig, def_weight,"Training and Testing")
    dataloader.AddBackgroundTree(tree_bkg_ee_qq,def_weight,"Training and Testing")
    dataloader.AddBackgroundTree(tree_bkg_ee_qqqq, def_weight,"Training and Testing")
    dataloader.AddBackgroundTree(tree_bkg_ee_qqqqqq, def_weight,"Training and Testing")
 

    #method = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT",
    #               ":".join([
    #                   "!H",
    #                   "!V",
    #                   "NTrees=850",-->default 800
    #                   "nEventsMin=150", -->default is 0
    #                   "MaxDepth=3",-->default
    #                   "BoostType=AdaBoost",-->default
    #                   "AdaBoostBeta=0.5", -->default
    #                   "SeparationType=GiniIndex",-->default, e.g. SDivSqrtSPlusB"
    #                   "nCuts=20",-->default
    #                   "PruneMethod=NoPruning",-->default
    #                   ]))
 
    #factory.TrainAllMethods()
    #factory.TestAllMethods()
    #factory.EvaluateAllMethods()


    cut_S = TCut("")
    cut_B = TCut("")
    dataloader.PrepareTrainingAndTestTree( cut_S, cut_B,
                                        ":".join([ "V",
                                                   "nTrain_Signal=0",
                                                   "nTrain_Background=0",
                                                   "SplitMode=Random",
                                                   "NormMode=NumEvents",
                                                   #"NormMode=None",
                                                   ])
                                        ) 

    factory.BookMethod(dataloader, root.TMVA.Types.kBDT, "BDT",
                       ":".join(["!H",
                                 "V",
                                 "NTrees=250",
                                 "MaxDepth=3",
                                 #"BoostType=Grad",
                                 #"Shrinkage=1.00",
                                 "BoostType=AdaBoost",
                                 "AdaBoostBeta=0.20",
                                 #"SeparationType=GiniIndexWithLaplace",
                                 "SeparationType=GiniIndex",
                                 "nCuts=-1",
                                 "PruneMethod=NoPruning",
                                 #"SkipNormalization",
                                 ]))

    #factory.WriteDataInformation()

    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()
    fileout.Close()


  
    return None

def process_files():


    input_file_signal_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_withPartonHistos_AllEvents.root"
    #/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_June24/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root"

    input_files_bkg=[] 
    input_file_ee_qq_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta.root"
    input_file_ee_qqqq_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/test_ee_qqqq_mqqqq_2TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta.root"
    input_file_ee_qqqqqq_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/test_ee_qqqqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta.root"
    input_files_bkg.append(input_file_ee_qq_)
    input_files_bkg.append(input_file_ee_qqqq_)
    input_files_bkg.append(input_file_ee_qqqqqq_)


    os.chdir('/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/')
    #os.rmdir('dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees250AdaBoostBeta020NCuts_m1')
    os.mkdir('dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees250AdaBoostBeta020NCuts_hzqq_AllEvents_m1_qqqq2TeV_y32_jet1_D2_C2_C3_jet2_D2_C2_C3')
    os.chdir('dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees250AdaBoostBeta020NCuts_hzqq_AllEvents_m1_qqqq2TeV_y32_jet1_D2_C2_C3_jet2_D2_C2_C3')

                       
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35_Aug7/MVATrainingWeights_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees250NCuts_m1__hzqq_AllEvents__ee_qq_mqq_1TeV__ee_qqqq_mqqqq_2TeV__ee_qqqqqq_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_y32_jet1_D2_C2_C3_jet2_D2_C2_C3.root"  
    process_event(final_histo_name_,input_file_signal_,input_files_bkg)

    input_file_signal_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta_withPartonHistos_AllEvents.root"
    #/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta.root"

    input_files_bkg=[]
    input_file_ee_qq_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta.root"
    input_file_ee_qqqq_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_ee_qqqq_mqqqq_2TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta.root"
    input_file_ee_qqqqqq_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/test_ee_qqqqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree_RecoSelectionOnly_JetSub_with_E_theta.root"
    input_files_bkg.append(input_file_ee_qq_)
    input_files_bkg.append(input_file_ee_qqqq_)
    input_files_bkg.append(input_file_ee_qqqqqq_)
     
    final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/MVATrainingWeights_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees250NCuts_m1__hzqq_AllEvents__ee_qq_mqq_1TeV__ee_qqqq_mqqqq_2TeV__ee_qqqqqq_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_y32_jet1_D2_C2_C3_jet2_D2_C2_C3.root"  

    os.chdir('/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190709Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35_Aug7/')
    #os.rmdir('dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees250AdaBoostBeta020NCuts_m1')
    os.mkdir('dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees250AdaBoostBeta020NCuts_hzqq_AllEvents_m1_qqqq2TeV_y32_jet1_D2_C2_C3_jet2_D2_C2_C3')
    os.chdir('dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees250AdaBoostBeta020NCuts_hzqq_AllEvents_m1_qqqq2TeV_y32_jet1_D2_C2_C3_jet2_D2_C2_C3')

    process_event(final_histo_name_,input_file_signal_,input_files_bkg)


    return None

process_files()



