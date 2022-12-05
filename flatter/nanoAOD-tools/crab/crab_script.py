#!/usr/bin/env python
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

from  PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puAutoWeight_2018
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import puWeightProducer

'''
from nanoAOD_filesPaths.BcToJpsiMuNu_new_nanoAOD_filesPath import BcToJpsiMuNu_files
from nanoAOD_filesPaths.BcToJpsiTauNu_nanoAOD_filesPath import BcToJpsiTauNu_files
from nanoAOD_filesPaths.OniaX_nanoAOD_filesPath import OniaX_files
from nanoAOD_filesPaths.data_nanoAOD_filesPath import data_files
'''
#p=PostProcessor(".",['root://t3dcachedb03.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/friti/crab_job_2020Jul28/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/crab_BcToJpsiMuNu_new/200728_084005/0000/RJPsi_mc_2020Jul28_1.root'],"",modules=[puAutoWeight_2018()],provenance=True,fwkJobReport=True,maxEntries=10000)
'''
datasets = [BcToJpsiMuNu_files, BcToJpsiTauNu_files, OniaX_files, data_files]

for dataset,name in zip(datasets,['BcMu','BcTau','Onia','data']):
    files_path = dataset
    p=PostProcessor(name,files_path,"",modules=[puAutoWeight_2018()],provenance=True,fwkJobReport=True,haddFileName=name+'_pu.root')

    p.run()

print "DONE"
'''
p=PostProcessor('.',inputFiles(),"",modules=[puAutoWeight_2018()],provenance=True,haddFileName='tree.root', fwkJobReport=True)
#p=PostProcessor('.',['root://cms-xrd-global.cern.ch//store/user/friti/crab_job_2020Jul28/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/crab_BcToJpsiMuNu_new/200728_084005/0000/RJPsi_mc_2020Jul28_1.root'],"",modules=[puAutoWeight_2018()],provenance=True,haddFileName='tree.root',maxEntries=100)

p.run()
