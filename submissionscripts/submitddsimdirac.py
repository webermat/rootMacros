import sys, os, subprocess
from string import Template
from DIRAC.Core.Base import Script
Script.parseCommandLine() 
from ILCDIRAC.Interfaces.API.DiracILC import  DiracILC
from ILCDIRAC.Interfaces.API.NewInterface.UserJob import *
from ILCDIRAC.Interfaces.API.NewInterface.Applications import *

jobGroup = "sim_K0L2_190220"
dirac = DiracILC(True,jobGroup+".rep")
job = UserJob()
job.setJobGroup(jobGroup)
job.setOutputSandbox ( [ "*.log","*.out","*.py"] )
job.setBannedSites(['LCG.IK0L2P3-CC.fr','OSG.UConn.us','LCG.Cracow.pl','OSG.MIT.us','LCG.Glasgow.uk','OSG.CIT.us','OSG.BNL.us','LCG.Brunel.uk','LCG.QMUL.uk'])
#pay attention that the Zuds200 here is NOT changed
job.setBannedSites(['LCG.INP3-CC.fr','OSG.UConn.us','LCG.Cracow.pl','OSG.MIT.us','LCG.Glasgow.uk','OSG.CIT.us','OSG.BNL.us','LCG.Brunel.uk','LCG.QMUL.uk','LCG.Oxford.uk'])
job.setSplitFilesAcrossJobs("LFN:/ilc/user/w/webermat/slcioFiles/mcparticlesK0L2_uniformCosTheta_70000.slcio",eventsPerFile=70000,eventsPerJob=500)
ddsim = DDSim()
ddsim.setVersion("ILCSoft-2019-02-20_gcc62")
ddsim.setDetectorModel("CLIC_o3_v14")
#filelist=[]
ddsimoutputname="ddsim_ILC19-02-20_gcc62_CLIC_o3_v14_K0L2.slcio"
ddsim.setOutputFile(ddsimoutputname)
res=job.append(ddsim)
if not res['OK']:
   print res['Message']
   exit()
job.setOutputData([ddsimoutputname],"ddsim/ILC19-02-20_gcc62_CT/CLIC_o3_v14/K0L2_70000evts","CERN-DST-EOS")
print ddsimoutputname
job.setName(ddsimoutputname)
job.dontPromptMe()
job.submit(dirac)
