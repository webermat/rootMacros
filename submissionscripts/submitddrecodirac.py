import sys, os, subprocess
from string import Template
from DIRAC.Core.Base import Script
Script.parseCommandLine()
from ILCDIRAC.Interfaces.API.DiracILC import  DiracILC
from ILCDIRAC.Interfaces.API.NewInterface.UserJob import *
from ILCDIRAC.Interfaces.API.NewInterface.Applications import *
  
#beware end index now one lower  than previously,  pay attention also start from 0 and not 1
for ind in range(0,200,1):
   jobGroup = "pim1500_CLIC_o3_v14_2018-11-01_CLIC_gcc62"
   dirac = DiracILC(True,jobGroup+".rep")
   job = UserJob()
   job.setJobGroup(jobGroup)
   job.setOutputSandbox ( [ "*.log","*.out","*.py"] )
   job.setBannedSites(['LCG.IN2P3-CC.fr','OSG.UConn.us','LCG.Cracow.pl','OSG.MIT.us','LCG.Glasgow.uk','OSG.CIT.us','OSG.BNL.us','LCG.Brunel.uk','LCG.QMUL.uk'])
   job = UserJob()
   job.setJobGroup(jobGroup)
   job.setOutputSandbox ( [ "*.log","*.out","*.py"] )
   job.setInputSandbox( ["LFN:/ilc/user/w/webermat/181011/PandoraSettings/PandoraSettingsDefault.xml","LFN:/ilc/user/w/webermat/181011/PandoraSettings/PandoraLikelihoodData12EBin.xml"] ) 
   #job.setInputSandbox( ["LFN:/ilc/user/w/webermat/171221/PandoraSettingsDefault_171221_SWC_CLIC_n_K0L_12bins_Ph_gcc62.xml","LFN:/ilc/user/w/webermat/171221/PandoraLikelihoodData12EBin_Zuds500_CLIC_o3_v14_171221_gcc62.xml"] )   
   job.setBannedSites(['LCG.INP3-CC.fr','OSG.UConn.us','LCG.Cracow.pl','OSG.MIT.us','LCG.Glasgow.uk','OSG.CIT.us','OSG.BNL.us','LCG.Brunel.uk','LCG.QMUL.uk','LCG.Oxford.uk'])
   ma = Marlin()
   ma.setVersion('ILCSoft-2018-11-01_gcc62') 
   ma.setInputFile("LFN:/ilc/user/w/webermat/ddsim/ILC18-10-11_gcc62_CLIC_o3_v14/pim1500_split100evts/ddsim_ILC18-10-11_gcc62_CLIC_o3_v14_pim1500_0%s.slcio"%((str)(ind)))

   #ma.setSteeringFile("/afs/cern.ch/work/w/weberma2/steeringFiles/testCalibPerformance_171221_CLIC_o3_v14_PhLikelihood12Ebin_TauFinder_SWC_CLIC_n_K0L.xml")
   ma.setSteeringFile("/afs/cern.ch/work/w/weberma2/steeringFiles/clicReconstruction.xml")
   ma.setDetectorModel("CLIC_o3_v14")
   #taurootfilename2="TauSignalStatistics_pim1500__CT_CLIC_o3_v13_0%s.root"%((str)(ind))
   #ma.setExtraCLIArguments("--MyTauFinder.FileName_Signal={taurootfilename}".format(taurootfilename=taurootfilename2))
   lcoutputreco="pim1500_ILC18-11-01_gcc62_CLIC_o3_v14__0%s.slcio"%((str)(ind))
   ma.setOutputFile(lcoutputreco)
   res=job.append(ma)
   if not res['OK']:
      print res['Message']
      exit()
   #job.setOutputData([lcoutputreco],"ddsimrecovalidation/ILC18-10-11_gcc62_CT_CLIC_o3_v14_SWC_CLIC_n_K0L_newTune/pim1500_500_pim15002EBin","CERN-DST-EOS")
   job.setOutputData([lcoutputreco],"ddsimrecovalidation/ILC18-11-01_gcc62_CLIC_o3_v14/pim1500","CERN-DST-EOS")
   print lcoutputreco
   job.setName(lcoutputreco)
   job.dontPromptMe()
   job.submit(dirac)


