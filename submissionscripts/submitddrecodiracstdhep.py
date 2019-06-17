import sys, os, subprocess
from string import Template

from DIRAC.Core.Base import Script
Script.parseCommandLine()
from ILCDIRAC.Interfaces.API.DiracILC import  DiracILC
from ILCDIRAC.Interfaces.API.NewInterface.UserJob import *
from ILCDIRAC.Interfaces.API.NewInterface.Applications import *

#0 to 10 for all, but 500, there it is 0 to 41
#for input_ind in range(9495,9496,1):
for input_ind in range(0,10,1):
   for ind in range(0,20,1):#for low stat 91-200 1 job, for 380,500,750,1000,2000 10 jobs, for 1500 and 3000 do 20 jobs per file, 1-401 for bbar
      #for 91-200 split 1000, for 500 slit 250, for 380,750,1000,2000 split 100,  same for 1500 and 3000 split125
      #jobGroup = "bb_180518_gcc62_SWC_CLIC_DR200"
      jobGroup = "Zuds3000_180518_gcc62_SWC_CLIC_DR200"
      dirac = DiracILC(True,jobGroup+".rep")
      job = UserJob()
      job.setJobGroup(jobGroup)
      job.setOutputSandbox ( [ "*.log","*.out","*.py"] )
      job.setBannedSites(['LCG.IN2P3-CC.fr','OSG.UConn.us','LCG.Cracow.pl','OSG.MIT.us','LCG.Glasgow.uk','OSG.CIT.us','OSG.BNL.us','LCG.Brunel.uk','LCG.QMUL.uk'])
      #pay attention that the Zuds 500 here is NOT changed
      job.setInputSandbox( ["LFN:/ilc/user/w/webermat/171109/PandoraSettingsDefault_171109_CLIC_o3_v13_12EBin.xml","LFN:/ilc/user/w/webermat/171109/PandoraLikelihoodData12EBin_Zuds500_CLIC_o3_v13_171109.xml"] ) 
      job.setBannedSites(['LCG.INP3-CC.fr','OSG.UConn.us','LCG.Cracow.pl','OSG.MIT.us','LCG.Glasgow.uk','OSG.CIT.us','OSG.BNL.us','LCG.Brunel.uk','LCG.QMUL.uk','LCG.Oxford.uk'])
      ma = Marlin()
      ma.setVersion('ILCSoft-2018-05-18_gcc62')
      #ma.setInputFile("LFN:/ilc/prod/clic/3tev/bb/CLIC_o3_v14/SIM/00009495/000/bb_sim_%s_%s.slcio"%(str(input_ind),(str)(ind)))
      ma.setInputFile("LFN:/ilc/user/w/webermat/ddsimstdhep/ILC18-03-02_CLIC_o3_v14/Zuds3000_split125/ddsim_ILC18-03-02_Zuds3000_CLIC_o3_v14_0%s_%s.slcio"%(str(input_ind),(str)(ind)))
      ma.setSteeringFile("/afs/cern.ch/work/w/weberma2/steeringFiles/testCalibPerformance_180514_CLIC_o3_v14_PhLikelihood12Ebin_TauFinder_SWC_CLIC.xml")
      ma.setDetectorModel("CLIC_o3_v14")
      #taurootfilename2="TauSignalStatistics_Zuds3000_CT_CLIC_o3_v14_0%s_%s.root"%(str(input_ind),(str)(ind))
      #ma.setExtraCLIArguments("--MyTauFinder.FileName_Signal={taurootfilename}".format(taurootfilename=taurootfilename2))
      #lcoutputreco ="ddsmp_ILC180518_gcc62_bb_9495_CLIC_o3_v14_SWC_CLIC_DR200_0%s_%s.slcio"%(str(input_ind),(str)(ind))
      lcoutputreco ="ddsmp_ILC180518_gcc62_Zuds3000_CLIC_o3_v14_SWC_CLIC_DR200_0%s_%s.slcio"%(str(input_ind),(str)(ind))
      ma.setOutputFile(lcoutputreco)
      res=job.append(ma)
      if not res['OK']:
         print res['Message']
         exit()
      #job.setOutputData([lcoutputreco],"ddsimrecostdhep/ILC17-05-18_gcc62_CT_CLIC_o3_v14_SWC/bb_SIM_9495_D0_Z0_DR200","CERN-DST-EOS")
      job.setOutputData([lcoutputreco],"ddsimrecostdhep/ILC18-05-18_gcc62_CT_CLIC_o3_v14_SWC/Zuds3000_D0_Z0_DR200","CERN-DST-EOS")
      print lcoutputreco
      job.setName(lcoutputreco)
      job.dontPromptMe()
      job.submit(dirac)

