#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TFile.h"


void process_event(TString i_final_histo_name_,TString i_input_signal_file_name_,std::vector <TString > i_input_bkg_file_name_){
  TMVA::Tools::Instance();
  
  TFile* input_file_signal_=TFile::Open(i_input_signal_file_name_);
  TFile* input_file_bkg_ee_qq_=TFile::Open(i_input_bkg_file_name_[0]);
  TFile* input_file_bkg_ee_qqqq_=TFile::Open(i_input_bkg_file_name_[1]);
  TFile* input_file_bkg_ee_qqqqqq_=TFile::Open(i_input_bkg_file_name_[2]);
  
  TTree* tree_sig = (TTree*)input_file_signal_->Get("MVATrainingVariables");
  TTree* tree_bkg_ee_qq = (TTree*)input_file_bkg_ee_qq_->Get("MVATrainingVariables");
  TTree* tree_bkg_ee_qqqq = (TTree*)input_file_bkg_ee_qqqq_->Get("MVATrainingVariables");
  TTree* tree_bkg_ee_qqqqqq = (TTree*)input_file_bkg_ee_qqqqqq_->Get("MVATrainingVariables");
  
  
  TFile* fileout = new TFile(i_final_histo_name_,"RECREATE");
  
  //flags V is if verbose set or not (typically false) , if it should be silent or not (silent== no output from mva, default false) ,transformations can be  
  //I;D;P;U;G,D,for identity, decorrelation,PCA,Uniform   and Gaussianisation  followed  by  decorrelation transformations
  //analysistype can be Classification,Regression,Multiclass or Auto (default Auto)
  TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", fileout,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
  
  TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");
  //factory->SetWeightExpression("eventWeight");
  dataloader->SetSignalWeightExpression("eventWeight");
  dataloader->SetBackgroundWeightExpression("eventWeight");
  dataloader->AddVariable("jet1_mass", 'F');
  dataloader->AddVariable("jet2_mass", 'F');
  dataloader->AddVariable("jet1_theta", 'F');
  dataloader->AddVariable("deltatheta:=jet1_theta-jet2_theta", 'F');
  dataloader->AddVariable("jet1_D2_beta1", 'F');
  dataloader->AddVariable("jet2_D2_beta1", 'F');
  dataloader->AddVariable("jet1_BTag_rfj_BTagMax", 'F');
  dataloader->AddSpectator("jet1_C2_beta1", 'F');
  dataloader->AddSpectator("jet2_C2_beta1", 'F');
  dataloader->AddSpectator("jet1_tau21", 'F');
  dataloader->AddSpectator("jet2_tau21", 'F');
  
  float def_weight=1;
  
  dataloader->AddSignalTree(tree_sig, def_weight);
  dataloader->AddBackgroundTree(tree_bkg_ee_qq,def_weight);
  dataloader->AddBackgroundTree(tree_bkg_ee_qqqq, def_weight);
  dataloader->AddBackgroundTree(tree_bkg_ee_qqqqqq, def_weight);
  
  /*
    #method = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT",
    #               ":".join([
    #                   "!H",
    #                   "!V",
    #                   "NTrees=850",-->default 800
    #                   "nEventsMin=150", -->default is 0
    #                   "MaxDepth=3",-->default
    #                   "BoostType=AdaBoost",-->default
    #                   "AdaBoostBeta=0.5", -->default
    #                   "SeparationType=GiniIndex",-->default
    #                   "nCuts=20",-->default
    #                   "PruneMethod=NoPruning",-->default
    #                   ]))
    
    #factory.TrainAllMethods()
    #factory.TestAllMethods()
    #factory.EvaluateAllMethods()
  */
  
  TCut cut_S = "";
  TCut cut_B = "";
  dataloader->PrepareTrainingAndTestTree( cut_S, cut_B,
					  "!V:nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents" );
  
  factory->BookMethod( dataloader,TMVA::Types::kBDT, "BDT","!H:!V");
  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  fileout->Close();
}

  
void process_files(){


  TString input_file_signal_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root";
  std::vector< TString >input_files_bkg;
  TString input_file_ee_qq_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root";
  TString input_file_ee_qqqq_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/test_ee_qqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root";
  TString input_file_ee_qqqqqq_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/test_ee_qqqqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root";
  input_files_bkg.push_back(input_file_ee_qq_);
  input_files_bkg.push_back(input_file_ee_qqqq_);
  input_files_bkg.push_back(input_file_ee_qqqqqq_);
    
  TString final_histo_name_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polm80/MVATrainingTrees_dm35/MVTrainingWeights__hzqq__ee_qq_mqq_1TeV__ee_qqqq__ee_qqqqqq_m1_126_dm_35_m2_92_5_dm_35_noThetaCutC++.root";  
  process_event(final_histo_name_,input_file_signal_,input_files_bkg);
 
  TString input_file_signal_polp_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/test_hzqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root";
  std::vector< TString >input_files_bkg_polp_;
  TString input_file_ee_qq_polp_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/test_ee_qq_mqq_1TeV_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root";
  TString input_file_ee_qqqq_polp_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/test_ee_qqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root";
  TString input_file_ee_qqqqqq_polp_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/test_ee_qqqqqq_ellipse_m1_126_dm_35_m2_92_5_dm_35_noThetaCut_MVATrainingTree.root";
  input_files_bkg_polp_.push_back(input_file_ee_qq_polp_);
  input_files_bkg_polp_.push_back(input_file_ee_qqqq_polp_);
  input_files_bkg_polp_.push_back(input_file_ee_qqqqqq_polp_);
    
  TString final_histo_name_polp_="/eos/user/w/weberma2/HistoFiles/HZAnalyzer/190417Prod/VLC7VtxRFJVLC7/polp80/MVATrainingTrees_dm35/MVTrainingWeights__hzqq__ee_qq_mqq_1TeV__ee_qqqq__ee_qqqqqq_m1_126_dm_35_m2_92_5_dm_35_noThetaCutC++.root";  
  process_event(final_histo_name_polp_,input_file_signal_polp_,input_files_bkg_polp_);

}




