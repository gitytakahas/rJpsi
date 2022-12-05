import os
import sys
import optparse
import datetime

from nanoAOD_filesPaths.BcToJpsiMuNu_files_path import BcToJpsiMuNu_files
from nanoAOD_filesPaths.BcToJpsiTauNu_files_path import BcToJpsiTauNu_files
from nanoAOD_filesPaths.OniaX_files_path import OniaX_files
from nanoAOD_filesPaths.data_files_path import data_files

def cfg_writer(sampleLabel, sampleDatasetPrivate, outdir, add):
    f = open("crab_cfg.py", "w")
    f.write("from WMCore.Configuration import Configuration\n")
    #f.write("from CRABClient.UserUtilities import config, getUsernameFromSiteDB\n")
    f.write("\nconfig = Configuration()\n")
    f.write("config.section_('General')\n")
    f.write("config.General.workArea = 'RJPsiNanoTools_%s' %('"+add+"')\n")
    f.write("config.General.requestName = '"+sampleLabel+"_nanotools'\n")
    f.write("config.General.transferLogs=True\n")
    f.write("config.section_('JobType')\n")
    f.write("config.JobType.pluginName = 'Analysis'\n")
    f.write("config.JobType.psetName = 'PSet.py'\n")
    f.write("config.JobType.scriptExe = 'crab_script.sh'\n")
    f.write("config.JobType.inputFiles = ['crab_script.py','../scripts/haddnano.py']\n")
    f.write("config.JobType.sendPythonFolder = True\n")
    f.write("config.section_('Data')\n")
    f.write("config.Data.userInputFiles = "+str(sampleDatasetPrivate)+"\n")
    f.write("config.Data.allowNonValidInputDataset = True\n")
    #f.write("config.Data.inputDBS = 'phys03'")
    #    f.write("config.Data.inputDBS = 'global'\n")
    f.write("config.Data.splitting = 'FileBased'\n")
    f.write("config.Data.unitsPerJob = 10\n")
    #config.Data.runRange = ''
    #f.write("config.Data.splitting = 'EventAwareLumiBased'")
    #f.write("config.Data.totalUnits = 10\n")
    f.write("config.Data.outLFNDirBase = '/store/user/friti/NanoPU_%s/%s' % ( '"+add+"','" +outdir+ "')\n")
    #f.write("config.Data.outLFNDirBase = '/store/user/friti/NanoPU/%s' % ('" +outdir+ "')\n")
    f.write("config.Data.publication = False\n")
    f.write("config.Data.outputDatasetTag = '"+sampleLabel+"'\n")
    f.write("config.section_('Site')\n")
    f.write("config.Site.storageSite = 'T3_CH_PSI'\n")
    #f.write("config.Site.storageSite = "T2_CH_CERN"
    #f.write("config.section_("User")
    #f.write("config.User.voGroup = 'dcms'
    f.close()

production_tag = datetime.date.today().strftime('%Y%b%d')

datasets = [BcToJpsiMuNu_files, BcToJpsiTauNu_files, OniaX_files]
names = ['BcToJpsiMuNu', 'BcToJpsiTauNu', 'OniaX']
for name,dataset in zip(names,datasets):
    print "Producing crab cfg file"
    cfg_writer(name, dataset, name, production_tag)

    print "Submitting crab jobs"
    os.system("crab submit crab_cfg.py")
    
