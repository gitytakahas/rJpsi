#! /usr/bin/env python
import uproot
#import uproot_methods
import pandas as pd
import numpy as np
import os
#import time
#from helper import *
import xgboost as xgb
#from root_pandas.root_pandas.readwrite import to_root
import root_pandas
#from root_numpy import fill_hist, array2root, array2tree
#fromt 
#import root_numpy
from ROOT import TFile, TTree
import argparse
parser = argparse.ArgumentParser(description="A simple ttree plotter")
parser.add_argument("-f", "--file", dest="file", default="/work/ytakahas/work/BsTauTau/CMSSW_10_2_10/src/job/BcJpsiTauNu_truthmc/Myroot_new.root", help="Signal file", required=True)
parser.add_argument("-p", "--prefix", dest="prefix", default="None", help="Prefix", required=True)
parser.add_argument("-m", "--model", dest="model", default="normal", help="model", required=True)

parser.add_argument("-o", "--outdir", dest="outdir", default="None", help="output dir.", required=True)
#parser.add_argument("-n", "--name", dest="name", default="comb", help="output name", required=False)
args = parser.parse_args()

pd.set_option('display.max_colwidth', -1)



def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

fdir = args.outdir
ensureDir(fdir)

# normal
#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'b_iso', 'b_iso_ntracks', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_sumofdnn', 'ncand', 'estar']

#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'b_iso_0p7', 'b_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_sumofdnn', 'tau_sumofdnn_1prong', 'tau_sumofdnn_otherB', 'tau_sumofdnn_pu', 'perEVT_dnn', 'ncand', 'estar']
#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'b_iso_0p7', 'b_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_sumofdnn', 'tau_sumofdnn_1prong', 'tau_sumofdnn_otherB', 'tau_sumofdnn_pu', 'ncand', 'estar']

#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7', 'tau_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi','pi1_pt', 'ncand', 'estar']

#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7',  'dr_b_pv', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi', 'pi1_dnn',  'pi2_dnn', 'pi3_dnn',  'estar']

#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7', 'tau_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi', 'tau_sumofdnn', 'ncand', 'estar']

#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7', 'tau_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi', 'tau_sumofdnn', 'ncand', 'estar']

#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7', 'tau_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi','pi1_pt', 'pi1_dnn', 'pi2_pt', 'pi2_dnn', 'pi3_pt', 'pi3_dnn', 'ncand', 'estar']
#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7', 'tau_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi','pi1_pt', 'pi1_dnn', 'pi2_pt', 'pi2_dnn', 'pi3_pt', 'pi3_dnn', 'ncand', 'estar', 'tau_mass', 'b_mass', 'tau_rhomass1', 'tau_rhomass2']


#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7',  'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi',  'tau_sumofdnn', 'ncand', 'estar', 'perEVT_mc']
features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7',  'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi',  'tau_sumofdnn','tau_sumofdnn_1prong', 'tau_sumofdnn_otherB', 'tau_sumofdnn_pu', 'ncand', 'estar']


#features = ['b_pt', 'b_eta', 'b_alpha', 'b_vprob', 'tau_iso_0p7', 'tau_iso_ntracks_0p7', 'b_lips', 'b_pvips', 'b_mindoca', 'dr_b_pv', 'dr_jpsi_tau', 'tau_fls3d', 'tau_vprob', 'tau_fls3d_wjpsi', 'tau_sumofdnn', 'tau_sumofdnn_1prong', 'tau_sumofdnn_otherB', 'tau_sumofdnn_pu', 'ncand', 'estar', 'mm2', 'mu1_isLoose', 'mu2_isLoose']

training_branches = sorted(features)
ntree_limit = 800
_model = xgb.Booster({'nthread': 6})

#mpath='../mva/model_' + args.model + '/xgb_fulldata_None.model'
mpath='/work/ytakahas/work/analysis/CMSSW_10_2_10/src/rJpsi/mva/model_' + args.model + '/xgb_fulldata_None.model'
os.system('ls -lart ' + mpath)
print('model path = ', mpath)
#_model.load_model('/work/ytakahas/work/analysis/CMSSW_10_2_10/src/BcJpsiTauNu/mva/model_' + args.model + '/xgb_fulldata_None.model')
_model.load_model(mpath)

#_model.load_model('/work/ytakahas/work/mva/BcJpsiTau/model/xgb_fulldata_None.model')

print(args.file)
os.system('ls -lart ' + args.file)

events = uproot.open(args.file)['tree']

#print 'setup'
#ofile = TFile('Myroot.root', 'recreate')
#otree = TTree('tree', 'tree')

#mass = np.zeros(1, dtype=float)
#xgbs = np.zeros(1, dtype=float)
#otree.Branch('mass', mass, 'mass/F')
#otree.Branch('xgbs', xgbs, 'xgbs/F')

for i, params in enumerate(events.iterate(outputtype=pd.DataFrame, entrysteps=1000000)):
    print 
    print(i, 'making ... ' + fdir + '/Myroot_' + args.prefix + '_' + str(i) + '.root')
    print 

    #  print 'chheck1: ', i
    _branches = params.copy()
#    print(_branches)
    #  print 'chheck2: ', i
    xgbs = _model.predict(xgb.DMatrix(_branches[training_branches].sort_index(axis=1)), ntree_limit=ntree_limit)
    #  print(xgb)
    _branches['xgbs'] = xgbs.tolist()
    #  print(type(xgb))
    #  print(type(pd.DataFrame(xgb)))
#    print(_branches)
    
    
#    import pdb; pdb.set_trace()

#    print('after filtering ...')
#    print(_branches)
#    print(branches_filtered)

#  print(type(_branches))

#  for index, row in _branches.iterrows():
##    print(row['b_mass'], row['xgb'])
#
#    mass[0] = row['b_mass']
#    xgb[0] = row['xgb']
#
#    otree.Fill()
#     print(row['b_mass']), print(row['xgb'])
#    print(index, row['b_mass'])
#  _branches['xgb']

#    xgbsgt6 =  _branches['xgbs'] > 4.52
    xgbsgt =  _branches['xgbs'] > -100.
#    xgbsgt =  _branches['xgbs'] > 6.
#    xgbsgt =  _branches['tau_index'] ==0

#    import pdb; pdb.set_trace()

    branches_filtered = _branches[xgbsgt]

    print(branches_filtered)
    # by putting "first", it only make the first one to be false, namely the highest in pT!
#    final = branches_filtered[branches_filtered.duplicated(subset=['evt','run','lumi','jpsi_pt','jpsi_eta', 'jpsi_phi'], keep='first')==False]
#    final = branches_filtered
#    print(final)
#    print(branches_filtered)

#    import pdb; pdb.set_trace()

#    final.to_root(fdir + '/Myroot_' + args.prefix + '_' + str(i) + '.root', key='tree', store_index=True)
    branches_filtered.to_root(fdir + '/Myroot_' + args.prefix + '_' + str(i) + '.root', key='tree', store_index=True)

#ofile.Write()
#ofile.Close()
