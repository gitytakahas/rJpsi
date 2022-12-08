from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config
#from PhysicsTools.NanoAODTools.crab.BcToJpsiMuNu_new_nanoAOD_filesPath import BcToJpsiMuNu_files
#from BcToJpsiMuNu_new_nanoAOD_filesPath import BcToJpsiMuNu_files

config = Configuration()

config.section_("General")
config.General.requestName = 'NanoPost2'
config.General.transferLogs=True
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['crab_script.py','../scripts/haddnano.py'] #hadd nano will not be needed once nano tools are in cmssw
config.JobType.sendPythonFolder	 = True
config.section_("Data")
#config.Data.inputDataset = '/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v4/MINIAODSIM'
config.Data.userInputFiles = ['root://cms-xrd-global.cern.ch//store/user/friti/crab_job_2020Jul28/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/crab_BcToJpsiMuNu_new/200728_084005/0000/RJPsi_mc_2020Jul28_1.root',
'root://cms-xrd-global.cern.ch//store/user/friti/crab_job_2020Jul28/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/crab_BcToJpsiMuNu_new/200728_084005/0000/RJPsi_mc_2020Jul28_2.root',
'root://cms-xrd-global.cern.ch//store/user/friti/crab_job_2020Jul28/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/crab_BcToJpsiMuNu_new/200728_084005/0000/RJPsi_mc_2020Jul28_3.root',
'root://cms-xrd-global.cern.ch//store/user/friti/crab_job_2020Jul28/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/crab_BcToJpsiMuNu_new/200728_084005/0000/RJPsi_mc_2020Jul28_4.root',
'root://cms-xrd-global.cern.ch//store/user/friti/crab_job_2020Jul28/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/crab_BcToJpsiMuNu_new/200728_084005/0000/RJPsi_mc_2020Jul28_5.root']
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 2

config.Data.outLFNDirBase = '/store/user/friti/NanoPileUp/'
config.Data.publication = False
config.Data.outputDatasetTag = 'BcToJpsiTauNu'
config.section_("Site")
config.Site.storageSite = "T3_CH_PSI"

#config.Site.storageSite = "T2_CH_CERN"
#config.section_("User")
#config.User.voGroup = 'dcms'

