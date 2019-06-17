import sys, os, subprocess
from string import Template
from DIRAC.Core.Base import Script
Script.parseCommandLine() 
from ILCDIRAC.Interfaces.API.DiracILC import  DiracILC
from ILCDIRAC.Interfaces.API.NewInterface.UserJob import *
from ILCDIRAC.Interfaces.API.NewInterface.Applications import *

from DIRAC.Resources.Catalog.FileCatalogClient import FileCatalogClient

jobGroup = "Zuds3000_12957_181101_allrangetuning_SWC_allSamples_nPFO_1"

#12959 1500
#12954,12955,12956,12957 are 100,200,500,3000

meta = {}
fc = FileCatalogClient()

meta = {}
meta['ProdID']='12957' 
meta['Datatype']='SIM'

res = fc.findFilesByMetadata(meta)
if not res['OK']:
   print res['Message']



lfns = res['Value']
print "Found %s files" % len(lfns)
filelist=[]
for lfn in lfns:
   filelist.append("LFN:"+lfn)
#print filelist



dirac = DiracILC(True,jobGroup+".rep")
job = UserJob()
job.setJobGroup(jobGroup)
job.setOutputSandbox ( [ "*.log","*.out","*.py","*.xml"] )
job.setBannedSites(['OSG.UCon.us','LCG.Cracow.pl','OSG.MIT.us','LCG.Glasgow.uk','OSG.CIT.us','OSG.BNL.us','LCG.Brunel.uk','LCG.QMUL.uk'])
#default settings and default likelihood file
job.setInputSandbox( ["LFN:/ilc/user/w/webermat/190222/PandoraSettings/PandoraLikelihoodData12EBin.xml","LFN:/ilc/user/w/webermat/190222/PandoraSettings/PandoraSettingsDefault.xml","LFN:/ilc/user/w/webermat/190222/libClicPerformance.tar.gz"] ) 
job.setBannedSites(['LCG.INP3-CC.fr','OSG.UCon.us','LCG.Cracow.pl','OSG.MIT.us','LCG.Glasgow.uk','OSG.CIT.us','OSG.BNL.us','LCG.Brunel.uk','LCG.QMUL.uk','LCG.Oxford.uk'])
job.setSplitInputData(filelist, numberOfFilesPerJob=30)
ma = Marlin()
ma.setVersion('ILCSoft-2018-11-01_gcc62')
#ma.setSteeringFile("LFN:/ilc/user/w/webermat/190222/clicReconstruction20181101JetAnalyzer.xml")
ma.setSteeringFile("/afs/cern.ch/user/w/weberma2/work/steeringFiles/clicReconstruction20181101JetAnalyzer.xml")
#ma.setSteeringFile("/home/weberma2/clicReconstruction20181101JetAnalyzer.xml")
ma.setDetectorModel("CLIC_o3_v14")
JetRootfilename="JetStudy_Zuds3000_12957_noOverlay_CLIC_o3_v14_CT_PandoraPFOs.root"
ma.setExtraCLIArguments("--MyJetAnalyzer.OutputRootFileName={Jetrootfilename2}".format(Jetrootfilename2=JetRootfilename))
res=job.append(ma)
if not res['OK']:
   print res['Message']
   exit()
job.setOutputData([JetRootfilename],"ddrecoZuds/CLIC_o3_v14/ILC18-11-01_gcc62_CT/SWC_allRange_allSamples_nPFO_1/Zuds3000_V2/12957","CERN-DST-EOS")
#job.setOutputData([JetRootfilename],"ddrecoZuds/ILC18-11-01_gcc62_CT/CLIC_o3_v14/SWCEndcap_0_85_0_92_redSamples_nPFO_1/Zuds3000/12957","CERN-DST-EOS")
#job.setOutputData([JetRootfilename],"ddrecoZuds/ILC18-11-01_gcc62_CT/CLIC_o3_v14/SWCEndcap_0_85_0_92_allSamples_Erat_0_99/Zuds3000/12957","CERN-DST-EOS")
#job.setOutputData([JetRootfilename],"ddrecoZuds/ILC18-11-01_gcc62_CT/CLIC_o3_v14/SWCEndcap_0_85_0_92_redSamples_Erat_0_99/Zuds3000/12957","CERN-DST-EOS")
print jobGroup
job.setName(jobGroup)
job.dontPromptMe()
job.submit(dirac)
