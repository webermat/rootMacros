# rootMacros
root macros processing

submissionScript: contains python scripts to submit jobs to the distributed computing system of ILCDirac


xml-files to run the ExtractionSteps (and first transformation of features) from CLICPerformance/Calorimetry

C++ macros for an early version of the analysis of physics data for Higgs+Z boson all-hadronic analysis

several python scripts for HZAnalysis:

HZAnalyzerMVATrainingTree.py:
creates out of the first tree a transformed tree after modifications on the Higgs jet

HZAnalyzerMVATrainingScript.py:
training and testing of the boosted decision tree (BDT) to classify HZ events

HZAnalyzerMVAReadingScript.py
final selection based on the output of the BDT

HZAnalyzerMVAPlottingScript.py
final histograms and plots, plots after preselection

several script to evaluate the impact of systematic uncertainties
