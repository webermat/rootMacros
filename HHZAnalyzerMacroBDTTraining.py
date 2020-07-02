from ROOT import gROOT, TCanvas, TF1, TH1F, TH1, TH2,  TH2F, TGraph, TCanvas, TLegend, TTree, TLorentzVector, TVector3, TStyle, gPad,gStyle,TColor,TMVA,TCut
import ROOT as root
from math import cos, sin, pi, degrees, radians, pow, sqrt,acos
from array import array
import sys
from array import array
import os

def process_event(i_final_histo_name_,i_input_signal_file_name_,i_input_bkg_file_name_):
    print "at start of process event"    
    root.TMVA.Tools.Instance()

    _output_training_name_=i_final_histo_name_
    _input_signal_name_=i_input_signal_file_name_
    _input_bkg_hzqq_name_=i_input_bkg_file_name_[0]
    _input_bkg_ee_qq_name_=i_input_bkg_file_name_[1]
    _input_bkg_ee_qqqq_name_=i_input_bkg_file_name_[2]
    _input_bkg_ee_qqqqqq_name_=i_input_bkg_file_name_[3]
    _input_bkg_WWH_qqqqH_name_=i_input_bkg_file_name_[4]
    _input_bkg_ZZH_qqqqH_name_=i_input_bkg_file_name_[5]

    fileout = root.TFile(_output_training_name_,"RECREATE")

    input_file_signal_=root.TFile.Open(_input_signal_name_)
    input_file_bkg_ee_qqqq_=root.TFile.Open(_input_bkg_ee_qqqq_name_)
    input_file_bkg_hzqq_=root.TFile.Open(_input_bkg_hzqq_name_)
    input_file_bkg_ee_qq_=root.TFile.Open(_input_bkg_ee_qq_name_)
    input_file_bkg_ee_qqqqqq_=root.TFile.Open(_input_bkg_ee_qqqqqq_name_)
    input_file_bkg_WWH_qqqqH_=root.TFile.Open(_input_bkg_WWH_qqqqH_name_)
    input_file_bkg_ZZH_qqqqH_=root.TFile.Open(_input_bkg_ZZH_qqqqH_name_)

    tree_sig = input_file_signal_.Get("MVATrainingVariables")
    tree_bkg_hzqq = input_file_bkg_hzqq_.Get("MVATrainingVariables")
    tree_bkg_ee_qq = input_file_bkg_ee_qq_.Get("MVATrainingVariables")
    tree_bkg_ee_qqqq = input_file_bkg_ee_qqqq_.Get("MVATrainingVariables")
    tree_bkg_ee_qqqqqq = input_file_bkg_ee_qqqqqq_.Get("MVATrainingVariables")
    tree_bkg_WWH_qqqqH = input_file_bkg_WWH_qqqqH_.Get("MVATrainingVariables")
    tree_bkg_ZZH_qqqqH = input_file_bkg_ZZH_qqqqH_.Get("MVATrainingVariables")

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
    dataloader.SetWeightExpression("weight")
    dataloader.SetSignalWeightExpression("weight")
    dataloader.SetBackgroundWeightExpression("weight")
    dataloader.AddVariable("comb_jet1_mass", 'F')
    dataloader.AddVariable("comb_jet2_mass", 'F')
    dataloader.AddVariable("comb_jet3_mass", 'F')
    dataloader.AddVariable("comb_jet1_dalpha", 'F')
    dataloader.AddVariable("comb_jet2_dalpha", 'F')
    dataloader.AddSpectator("comb_jet3_dalpha", 'F')
    dataloader.AddSpectator("comb_jet1_E1_over_Etot", 'F')
    dataloader.AddSpectator("comb_jet2_E1_over_Etot", 'F')
    dataloader.AddSpectator("comb_jet3_E1_over_Etot", 'F')
    dataloader.AddVariable("comb_jet1_BTagMax", 'F')
    dataloader.AddVariable("comb_jet2_BTagMax", 'F')
    dataloader.AddVariable("comb_jet3_BTagMax", 'F')
    dataloader.AddSpectator("comb_jet1_cosThetaHel_absmin", 'F')
    dataloader.AddSpectator("comb_jet2_cosThetaHel_absmin", 'F')
    dataloader.AddSpectator("comb_jet3_cosThetaHel_absmin", 'F')
    dataloader.AddSpectator("comb_jet1_E", 'F')
    dataloader.AddSpectator("comb_jet2_E", 'F')
    dataloader.AddSpectator("comb_jet3_E", 'F')
    dataloader.AddSpectator("comb_jet1_theta", 'F')
    dataloader.AddSpectator("comb_jet2_theta", 'F')
    dataloader.AddSpectator("comb_jet3_theta", 'F')
    dataloader.AddVariable("BTag_sum_max2", 'F')
    dataloader.AddVariable("BTag_sum_max3", 'F')
    dataloader.AddSpectator("BTag_sum_max4", 'F')
    dataloader.AddSpectator("BTag_sum_all", 'F')
    #CTag and LTag variables, only present for tight mass selection
    dataloader.AddSpectator("CTag_sum_max2", 'F')
    dataloader.AddSpectator("CTag_sum_max3", 'F')
    dataloader.AddVariable("LTag_sum_all", 'F')
    dataloader.AddVariable("comb_jet1_CTagMax", 'F')
    dataloader.AddVariable("comb_jet2_CTagMax", 'F')
    dataloader.AddSpectator("comb_jet3_CTagMax", 'F')
    dataloader.AddVariable("comb_jet1_LTagMax", 'F')
    dataloader.AddVariable("comb_jet2_LTagMax", 'F')
    dataloader.AddVariable("comb_jet3_LTagMax", 'F')
    dataloader.AddSpectator("comb_jet1_LTagMin", 'F')
    dataloader.AddSpectator("comb_jet2_LTagMin", 'F')
    dataloader.AddSpectator("comb_jet3_LTagMin", 'F')
    dataloader.AddVariable("jet1_theta", 'F')
    dataloader.AddVariable("jet2_theta", 'F')
    dataloader.AddVariable("jet3_theta", 'F')
    dataloader.AddVariable("jet4_theta", 'F')
    dataloader.AddVariable("jet5_theta", 'F')
    dataloader.AddVariable("jet6_theta", 'F')
    dataloader.AddVariable("jet1_E", 'F')
    dataloader.AddVariable("jet2_E", 'F')
    dataloader.AddVariable("jet3_E", 'F')
    dataloader.AddVariable("jet4_E", 'F')
    dataloader.AddVariable("jet5_E", 'F')
    dataloader.AddVariable("jet6_E", 'F')
    dataloader.AddSpectator("jet1_BTag", 'F')
    dataloader.AddSpectator("jet2_BTag", 'F')
    dataloader.AddSpectator("jet3_BTag", 'F')
    dataloader.AddSpectator("jet4_BTag", 'F')
    dataloader.AddSpectator("jet5_BTag", 'F')
    dataloader.AddSpectator("jet6_BTag", 'F')
    dataloader.AddVariable("y12", 'F')
    dataloader.AddVariable("y23", 'F')
    dataloader.AddVariable("y34", 'F')
    dataloader.AddVariable("y45", 'F')
    dataloader.AddVariable("y56", 'F')

    def_weight=1

    dataloader.AddSignalTree(tree_sig, def_weight,"Training and Testing")
    dataloader.AddBackgroundTree(tree_bkg_ee_qqqq, def_weight,"Training and Testing")
    dataloader.AddBackgroundTree(tree_bkg_ee_qq,def_weight,"Training and Testing")
    dataloader.AddBackgroundTree(tree_bkg_ee_qqqqqq, def_weight,"Training and Testing")
    dataloader.AddBackgroundTree(tree_bkg_hzqq,def_weight,"Training and Testing")
    dataloader.AddBackgroundTree(tree_bkg_WWH_qqqqH, def_weight,"Training and Testing")
    dataloader.AddBackgroundTree(tree_bkg_ZZH_qqqqH, def_weight,"Training and Testing")

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
                                 "NTrees=300",
                                 "MaxDepth=3",
                                 #"BoostType=Grad",
                                 #"Shrinkage=1.50",
                                 "BoostType=AdaBoost",
                                 "AdaBoostBeta=0.20",
                                 #"SeparationType=GiniIndexWithLaplace",
                                 "SeparationType=GiniIndex",
                                 "nCuts=20",
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
    #noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins-->correct name now
    #noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTag_noLTagMins_ --> has also no noj3CMax
    #noEratios_noCosThetaHel 
    #hhqq_vs_hzqq__ee_qq__ee_qqqq__ee_qqqqqq

    print 'start VLC11 polm80'
    input_file_signal_VLC11_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hhqq_14364.root"

    input_files_bkg_VLC11_polm80_=[] 
    input_file_hzqq_VLC11_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hzqq_13391.root"
    input_file_ee_qq_VLC11_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qq_13399_to_13402.root"
    input_file_ee_qqqq_VLC11_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqq_13394_to_13397.root"
    input_file_ee_qqqqqq_VLC11_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqqqq.root"
    input_file_WWH_qqqqH_VLC11_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_WWH_qqqqH_14734.root"
    input_file_ZZH_qqqqH_VLC11_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ZZH_qqqqH_14726.root"
    input_files_bkg_VLC11_polm80_.append(input_file_hzqq_VLC11_polm80_)
    input_files_bkg_VLC11_polm80_.append(input_file_ee_qq_VLC11_polm80_)
    input_files_bkg_VLC11_polm80_.append(input_file_ee_qqqq_VLC11_polm80_)
    input_files_bkg_VLC11_polm80_.append(input_file_ee_qqqqqq_VLC11_polm80_)
    input_files_bkg_VLC11_polm80_.append(input_file_WWH_qqqqH_VLC11_polm80_)
    input_files_bkg_VLC11_polm80_.append(input_file_ZZH_qqqqH_VLC11_polm80_)

    os.chdir('/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/')
    os.mkdir('dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_BTag3_1_00_BTag4_2_20___hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ')
    os.chdir('dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_BTag3_1_00_BTag4_2_20___hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ')

    final_histo_name_VLC11_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polm80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/MVATrainingWeights_ee_qqqqqq_only_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ.root"  
    process_event(final_histo_name_VLC11_polm80_,input_file_signal_VLC11_polm80_,input_files_bkg_VLC11_polm80_)

    print 'start VLC11 polp80'
    input_file_signal_VLC11_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hhqq_14365.root"

    input_files_bkg_VLC11_polp80_=[] 
    input_file_hzqq_VLC11_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hzqq_13392.root"
    input_file_ee_qq_VLC11_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qq_13398.root"
    input_file_ee_qqqq_VLC11_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqq_13393.root"
    input_file_ee_qqqqqq_VLC11_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqqqq.root"
    input_file_WWH_qqqqH_VLC11_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_WWH_qqqqH_14735.root"
    input_file_ZZH_qqqqH_VLC11_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ZZH_qqqqH_14727.root"
    input_files_bkg_VLC11_polp80_.append(input_file_hzqq_VLC11_polp80_)
    input_files_bkg_VLC11_polp80_.append(input_file_ee_qq_VLC11_polp80_)
    input_files_bkg_VLC11_polp80_.append(input_file_ee_qqqq_VLC11_polp80_)
    input_files_bkg_VLC11_polp80_.append(input_file_ee_qqqqqq_VLC11_polp80_)
    input_files_bkg_VLC11_polp80_.append(input_file_WWH_qqqqH_VLC11_polp80_)
    input_files_bkg_VLC11_polp80_.append(input_file_ZZH_qqqqH_VLC11_polp80_)

    os.chdir('/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/')
    os.mkdir('dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_BTag3_1_00_BTag4_2_20___hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ')
    os.chdir('dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins_BTag3_1_00_BTag4_2_20___hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ')

    final_histo_name_VLC11_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC11_NJets6_finalAnalysis/polp80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/MVATrainingWeights_ee_qqqqqq_only_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_hzqq_qq_qqqq_qqqqqq_WWH_WWZ.root"  
    process_event(final_histo_name_VLC11_polp80_,input_file_signal_VLC11_polp80_,input_files_bkg_VLC11_polp80_)
    ''''
    print 'start VLC14 polm80'
    input_file_signal_VLC14_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hhqq_14364.root"

    input_files_bkg_VLC14_polm80_=[] 
    input_file_hzqq_VLC14_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hzqq_13391.root"
    input_file_ee_qq_VLC14_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qq_13399_to_13402.root"
    input_file_ee_qqqq_VLC14_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqq_13394_to_13397.root"
    input_file_ee_qqqqqq_VLC14_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqqqq.root"
    input_files_bkg_VLC14_polm80_.append(input_file_hzqq_VLC14_polm80_)
    input_files_bkg_VLC14_polm80_.append(input_file_ee_qq_VLC14_polm80_)
    input_files_bkg_VLC14_polm80_.append(input_file_ee_qqqq_VLC14_polm80_)
    input_files_bkg_VLC14_polm80_.append(input_file_ee_qqqqqq_VLC14_polm80_)
 


    os.chdir('/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/')
    os.mkdir('dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq')
    os.chdir('dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq')

    final_histo_name_VLC14_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polm80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/MVATrainingWeights_ee_qqqqqq_only_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq__ee_qqqqqq.root"  
    process_event(final_histo_name_VLC14_polm80_,input_file_signal_VLC14_polm80_,input_files_bkg_VLC14_polm80_)

    print 'start VLC14 polp80'
    input_file_signal_VLC14_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hhqq_14365.root"

    input_files_bkg_VLC14_polp80_=[] 
    input_file_hzqq_VLC14_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hzqq_13392.root"
    input_file_ee_qq_VLC14_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qq_13398.root"
    input_file_ee_qqqq_VLC14_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqq_13393.root"
    input_file_ee_qqqqqq_VLC14_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqqqq.root"
    input_files_bkg_VLC14_polp80_.append(input_file_hzqq_VLC14_polp80_)
    input_files_bkg_VLC14_polp80_.append(input_file_ee_qq_VLC14_polp80_)
    input_files_bkg_VLC14_polp80_.append(input_file_ee_qqqq_VLC14_polp80_)
    input_files_bkg_VLC14_polp80_.append(input_file_ee_qqqqqq_VLC14_polp80_)


    os.chdir('/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/')
    os.mkdir('dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq')
    os.chdir('dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq')

    final_histo_name_VLC14_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC14_NJets6_finalAnalysis/polp80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/MVATrainingWeights_ee_qqqqqq_only_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq__vs_hzqq__ee_qq__ee_qqqq.root"  
    process_event(final_histo_name_VLC14_polp80_,input_file_signal_VLC14_polp80_,input_files_bkg_VLC14_polp80_)

    print 'start VLC7 polm80'
    input_file_signal_VLC7_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hhqq_14364.root"

    input_files_bkg_VLC7_polm80_=[] 
    input_file_hzqq_VLC7_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hzqq_13391.root"
    input_file_ee_qq_VLC7_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qq_13399_to_13402.root"
    input_file_ee_qqqq_VLC7_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqq_13394_to_13397.root"
    input_file_ee_qqqqqq_VLC7_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqqqq.root"
    input_files_bkg_VLC7_polm80_.append(input_file_hzqq_VLC7_polm80_)
    input_files_bkg_VLC7_polm80_.append(input_file_ee_qq_VLC7_polm80_)
    input_files_bkg_VLC7_polm80_.append(input_file_ee_qqqq_VLC7_polm80_)
    input_files_bkg_VLC7_polm80_.append(input_file_ee_qqqqqq_VLC7_polm80_)


    os.chdir('/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/')
    os.mkdir('dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq')
    os.chdir('dataset/filesPolm80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq')

    final_histo_name_VLC7_polm80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polm80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/MVATrainingWeights_ee_qqqqqq_only_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq__ee_qqqqqq.root"  
    process_event(final_histo_name_VLC7_polm80_,input_file_signal_VLC7_polm80_,input_files_bkg_VLC7_polm80_)

    print 'start VLC7 polp80'
    input_file_signal_VLC7_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hhqq_14365.root"

    input_files_bkg_VLC7_polp80_=[] 
    input_file_hzqq_VLC7_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_hzqq_13392.root"
    input_file_ee_qq_VLC7_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qq_13398.root"
    input_file_ee_qqqq_VLC7_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqq_13393.root"
    input_file_ee_qqqqqq_VLC7_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/NTuplesAnalysis_E_theta_BTag3_1_00_BTag4_2_20_tight_Mass_Cuts/ntuple_HHZ_ECut_j1_150_j2_100_j4_50_dTheta12_80_BTag3_1_00_BTag4_2_20_M2_75_M3_50_150_ee_qqqqqq.root"
    input_files_bkg_VLC7_polp80_.append(input_file_hzqq_VLC7_polp80_)
    input_files_bkg_VLC7_polp80_.append(input_file_ee_qq_VLC7_polp80_)
    input_files_bkg_VLC7_polp80_.append(input_file_ee_qqqq_VLC7_polp80_)
    input_files_bkg_VLC7_polp80_.append(input_file_ee_qqqqqq_VLC7_polp80_)

    os.chdir('/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/')
    os.mkdir('dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq')
    os.chdir('dataset/filesPolp80GiniIndexNormNumEventsMaxDepth3NTrees300AdaBoostBeta020NCuts_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq')

    final_histo_name_VLC7_polp80_="/eos/user/w/weberma2/HistoFiles/HHZAnalyzer/190904Prod/VLC7_NJets6_finalAnalysis/polp80/MVATrainingTrees_Jan2020_E_theta_BTag_tight_Mass_Cuts/MVATrainingWeights_ee_qqqqqq_only_GiniIndexNormNumEventsMaxDepth3AdaBoostBeta020NTrees300NCuts20_noEratios_noCosThetaHel_noCTagSum_noj3dalphaCTagCTagMax_noLTagMins___hhqq_vs_ee_qqqqqq__ee_qqqqqq.root"  
    process_event(final_histo_name_VLC7_polp80_,input_file_signal_VLC7_polp80_,input_files_bkg_VLC7_polp80_)
    '''

    return None

process_files()



