import sys, os, subprocess
from string import Template

from DIRAC.Core.Base import Script
Script.parseCommandLine()
from ILCDIRAC.Interfaces.API.DiracILC import  DiracILC
from ILCDIRAC.Interfaces.API.NewInterface.UserJob import *
from ILCDIRAC.Interfaces.API.NewInterface.Applications import *

from DIRAC.Resources.Catalog.FileCatalogClient import FileCatalogClient

fc = FileCatalogClient()

#14343 hhz polm80
#14344 hhz polp80


#13392 hzqq polp80
#13391 hzqq polm80

#13390 hzll polp80
#13389 hzll polm80

#ee qq samples 13398 polp80, ee_qq_mqq_1TeV 13429 with 1 TeV m(qq) cut
#ee qq samples 13399,13400,13401,13402 polm80  ee_qq_mqq_1TeV 13425,13426,13427,13428 with 1 TeV m(qq) cut

#ee_qqqq samples 13393 polp80, ee_qqqq_mqqqq_2TeV 13700 with 2 TeV m(qqqq) cut
#ee_qqqq samples 13394,13395,13396,13397 polm80 ee_qqqq_mqqqq_2TeV 13696,13697,13698,13699 with 2 TeV m(qqqq) cut

#ee_bbcbbc 13094 polm  #ee_bbcbbc 13071 polp       #ee_bbubbu 13095 polm #ee_bbubbu 13072 polp      #ee_ddcyyc 13096 polm  #ee_ddcyyc 13073 polp        #ee_dduyyu 13097 polm #ee_dduyyu 13074 polp
#ee_sscbbc 13098 polm #ee_sscbbc 13075 polp   #ee_sscssc 13099 polm #ee_sscssc 13076 polp      #ee_ssussu 13123 polm #ee_ssussu 13077 polp    #ee_ssubbu 13292 polm #ee_ssubbu 13293 polp      #ee_yycbbu 13318 polm #ee_yycbbu 13322 polp
#ee_yycddu 13326 polm #ee_yycddu 13319 polp     #ee_yycssu 13323 polm #ee_yycssu 13327 polp   #ee_yyubbc 13320 polm #ee_yyubbc 13324 polp     #ee_yyuddc 13328 polm #ee_yyuddc 13321 polp     #ee_yyussc 13325 polm #ee_yyussc 13329 polm

file_processor_list=[]
#process those hzqq files
file_processor_list.append(['13392','hzqq','polp80'])
file_processor_list.append(['13391','hzqq','polm80'])
#process those ee_qq files
file_processor_list.append(['13398','ee_qq','polp80'])
file_processor_list.append(['13399','ee_qq','polm80'])
file_processor_list.append(['13400','ee_qq','polm80'])
file_processor_list.append(['13401','ee_qq','polm80'])
file_processor_list.append(['13402','ee_qq','polm80'])
#process those ee_qqqq files
file_processor_list.append(['13393','ee_qqqq','polp80'])
file_processor_list.append(['13394','ee_qqqq','polm80'])
file_processor_list.append(['13395','ee_qqqq','polm80'])
file_processor_list.append(['13396','ee_qqqq','polm80'])
file_processor_list.append(['13397','ee_qqqq','polm80'])
#process those ee_qqqqqq files
file_processor_list.append(['13094','ee_bbcbbc','polm80'])
file_processor_list.append(['13071','ee_bbcbbc','polp80'])
file_processor_list.append(['13095','ee_bbubbu','polm80'])
file_processor_list.append(['13072','ee_bbubbu','polp80'])
file_processor_list.append(['13096','ee_ddcyyc','polm80'])
file_processor_list.append(['13073','ee_ddcyyc','polp80'])
file_processor_list.append(['13097','ee_dduyyu','polm80'])
file_processor_list.append(['13074','ee_dduyyu','polp80'])
file_processor_list.append(['13098','ee_sscbbc','polm80'])
file_processor_list.append(['13075','ee_sscbbc','polp80'])
file_processor_list.append(['13099','ee_sscssc','polm80'])
file_processor_list.append(['13076','ee_sscssc','polp80'])
file_processor_list.append(['13123','ee_ssussu','polm80'])
file_processor_list.append(['13077','ee_ssussu','polp80'])
file_processor_list.append(['13292','ee_ssubbu','polm80'])
file_processor_list.append(['13293','ee_ssubbu','polp80'])
file_processor_list.append(['13318','ee_yycbbu','polm80'])
file_processor_list.append(['13322','ee_yycbbu','polp80'])
file_processor_list.append(['13326','ee_yycddu','polm80'])
file_processor_list.append(['13319','ee_yycddu','polp80'])
file_processor_list.append(['13323','ee_yycssu','polm80'])
file_processor_list.append(['13327','ee_yycssu','polp80'])
file_processor_list.append(['13320','ee_yyubbc','polm80'])
file_processor_list.append(['13324','ee_yyubbc','polp80'])
file_processor_list.append(['13328','ee_yyuddc','polm80'])
file_processor_list.append(['13321','ee_yyuddc','polp80'])
file_processor_list.append(['13325','ee_yyussc','polm80'])
file_processor_list.append(['13329','ee_yyussc','polp80'])

for file in file_processor_list:
   meta = {} 
   meta['ProdID']=file[0] 
   meta['Datatype']='DST'
 
   res = fc.findFilesByMetadata(meta)
   if not res['OK']:
      print res['Message']

   lfns = res['Value']
   #print "Found %s files" % len(lfns)
   filelist=[]
   for lfn in lfns:
      filelist.append(lfn)
   #print filelist
   #filelist2=filelist[0]
   #print file[1],file[0],file[2]
   #check if big radius in fastjetanalyzer indeed 0.7 or 1.0
   jobGroup = "HHZAnalyzer_19-09-04_%s_%s_VLC10PFOs_NJet3_%s"%(file[1],file[0],file[2])
   dirac = DiracILC(True,jobGroup+".rep")
   job = UserJob()
   #job.setExecutionEnv({'ROOT_INCLUDE_PATH':'./'})
   job.setJobGroup(jobGroup)
   job.setOutputSandbox ( [ "*.log","*.out","*.py"] )
   job.setBannedSites(['LCG.IN2P3-CC.fr','OSG.UConn.us','LCG.Cracow.pl','OSG.MIT.us','LCG.Glasgow.uk','OSG.CIT.us','OSG.BNL.us','LCG.Brunel.uk','LCG.QMUL.uk'])
   #pay attention that the Zuds200 here is NOT changed
   #job.setInputSandbox( ["LFN:/ilc/user/w/webermat/190904/libClicPerformanceNJet2.tar.gz", "LFN:/ilc/user/w/webermat/190412/vtxprob.tar.gz","LFN:/ilc/user/w/webermat/190412/flavourTagging04-01_ct_90deg/lcfiweights.tar.gz"] ) 
   job.setInputSandbox( ["LFN:/ilc/user/w/webermat/190904/libClicPerformance.tar.gz", "LFN:/ilc/user/w/webermat/190412/vtxprob.tar.gz","LFN:/ilc/user/w/webermat/190412/flavourTagging04-01_ct_90deg/lcfiweights.tar.gz"] ) 
   job.setBannedSites(['LCG.INP3-CC.fr','OSG.UConn.us','LCG.Cracow.pl','OSG.MIT.us','LCG.Glasgow.uk','OSG.CIT.us','OSG.BNL.us','LCG.Brunel.uk','LCG.QMUL.uk','LCG.Oxford.uk'])
   job.setSplitInputData(filelist, numberOfFilesPerJob=50)
   ma = Marlin()
   ma.setVersion('ILCSoft-2019-09-04_gcc62')
   #ma.setInputFile("LFN:/ilc/user/w/webermat/ddsimstdheptautau/ILC18-10-11_gcc62_CLIC_o3_v14/tautau200/ddsim_ILC181011_gcc62_tautau_200_CLIC_o3_v14_0%s.slcio"%(str(input_ind)))
   ma.setSteeringFile("/eos/user/w/weberma2/steeringFiles/testHHZAnalyzer.xml")
   ma.setDetectorModel("CLIC_o3_v14")
   HHZrootfilename2="HHZStudy_%s_%s_%s_3TeV_wO_CLIC_o3_v14.root"%(file[1],file[0],file[2])
   #print HHZrootfilename2
   #TrackPtMin default -1, saveMEInfo=true and Rmax=0.7
   ma.setExtraCLIArguments("--MyHHZAnalyzer.saveMEInfo=false --MyHHZAnalyzer.OutputRootFileName={HHZrootfilename}".format(HHZrootfilename=HHZrootfilename2))
   res=job.append(ma)
   if not res['OK']:
      print res['Message']
      exit()
   outputdirectory="HHZAnalyzer/ILC19-09-04_gcc62_CT/CLIC_o3_v14/VtxRFJVLC10_NJet3/%s/%s/%s"%(file[2],file[1],file[0])
   #print outputdirectory
   job.setOutputData([HHZrootfilename2],outputdirectory,"CERN-DST-EOS")
   print jobGroup,meta['ProdID'],file[2],file[1]
   job.setName(jobGroup)
   job.dontPromptMe()
   job.submit(dirac)

