##Example to submit Marlin job: MarlinExample.py
from DIRAC.Core.Base import Script
Script.parseCommandLine()
from ILCDIRAC.Interfaces.API.DiracILC import  DiracILC
from ILCDIRAC.Interfaces.API.NewInterface.UserJob import *
from ILCDIRAC.Interfaces.API.NewInterface.Applications import *
jobGroup = "ph10_sim_CLIC_o3_v05_2016_07_04"
dirac = DiracILC(True,jobGroup+".rep")
job = UserJob()
#all input files to be send here
job.setInputSandbox( ["LFN:/ilc/user/w/webermat/160620/photonReReco/lib.tar.gz" , "LFN:/ilc/user/w/webermat/160502/photonReReco/PandoraSettingsFast.xml", "LFN:/ilc/user/w/webermat/160502/photonReReco/PandoraLikelihoodData9EBin_CLIC_ILD.xml" ] )
job.setOutputSandbox ( [ "*.log","*.out","hnunu_rec_6265_photonrereco_3_0_TeV_ntuple160620_233_25.root"] )
job.setOutputData(["hnunu_rec_6265_photonrereco_3_0_TeV_ntuple160620_233_25.root"], "160620_"+jobGroup, "CERN-SRM")
job.setInputData([
#job.setParametricInputData([
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_233.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_234.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_235.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_236.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_237.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_238.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_239.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_24.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_240.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_241.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_242.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_243.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_244.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_245.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_246.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_247.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_248.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_249.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_25.slcio"
])
job.setName("h_nunu_rec_6265_3TeV")
#job.setDestination("LCG.CERN.ch")
#job.setCPUTime(300000)
ma = Marlin()
ma.setVersion('ILCSoft-01-17-09_gcc48')
#ILCSoft-01-17-08_gcc48
ma.setSteeringFile("test_fastjet_new_reprocess_drop.xml")
ma.setGearFile("/afs/cern.ch/user/w/weberma2/public/clic_ild_cdr.gear")
ma.setInputFile([
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_233.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_234.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_235.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_236.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_237.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_238.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_239.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_24.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_240.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_241.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_242.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_243.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_244.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_245.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_246.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_247.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_248.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_249.slcio",
"/ilc/prod/clic/3tev/h_nunu/ILD/REC/00006265/000/h_nunu_rec_6265_25.slcio"
])
ma.setOutputFile( "hnunu_rec_6265_photonrereco_3_0_TeV_ntuple160620_233_25.root")
ma.setProcessorsToExclude(["libDDMarlinPandora.so"])
res = job.append(ma)
if not res['OK']:
  print res['Message']
  exit()
#dirac.submit(job)
print job.submit(dirac)
