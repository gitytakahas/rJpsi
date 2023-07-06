import copy, os,  math, sys, shutil
from numpy import array
from common.officialStyle import officialStyle
from array import array
import common.MultiDraw
import numpy as np
from varConfig import vardir
from common.DataMCPlot import *
from common.DisplayManager_postfit import DisplayManager
from common.DisplayManager_compare import DisplayManager_compare
from common.helper import *
from common.H2TauStyle import *
import ConfigParser
from optparse import OptionParser, OptionValueError
import json 
from ROOT import TChain

#gROOT.SetBatch(False)
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)


gROOT.ProcessLine(".L ~/tool//MultiDraw.cc+");
gROOT.Macro('common/functionmacro.C+')



#gStyle.SetPadRightMargin (0.14)
def comparisonPlots(hist, pname='sync.pdf', isLog=False, isRatio=True, clabel=''):

    display = DisplayManager(pname, isRatio, isLog, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)

def comparisonPlots_alt(hists, titles, norm=False, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he', prefix=None):

    display = DisplayManager_compare(pname, isLog, isRatio, 0.6, 0.7, doption)
    display.draw_legend = isLegend
    
    display.Draw(hists, titles, norm, prefix)


#pnfs_prefix='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/results_simultaneous/'
 
def add_lumi(year):
    lowX=0.7
    lowY=0.842
    lumi  = TPaveText(lowX, lowY+0.06, lowX+0.65, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.055)
    lumi.SetTextFont (   42 )
#    lumi.AddText("2016, " + str(luminumber) + " fb^{-1} (13TeV)")
    lumi.AddText(year)
    return lumi


def add_CMS():
    lowX=0.15
    lowY=0.87
    lumi  = TPaveText(lowX, lowY, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.05)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS")
    return lumi


def add_channel(cl):
    lowX=0.30
    lowY=0.81
#    lumi  = ROOT.TLatex(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
#    lumi  = ROOT.TLatex(lowX, lowY,  cl)
    lumi  = TLatex(lowX, lowY,  cl)
    lumi.SetNDC(True)
#    lumi.SetTextFont(61)
    lumi.SetTextSize(0.08)
#    lumi.SetBorderSize(   0 )
#    lumi.SetFillStyle(    0 )
#    lumi.SetTextAlign(   12 )
#    lumi.SetTextColor(    1 )
#    lumi.AddText(cl)
    return lumi


def add_Preliminary():
    lowX=0.26
    lowY=0.832
    lumi  = TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
#    lumi.SetTextFont(52)
    lumi.SetTextSize(0.05)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
#    lumi.AddText("#it{Simulation} Preliminary")
    lumi.AddText("#it{Preliminary}")
    return lumi

prefix ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/'


filedir={
#    'normal':{'filename':'root://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas/BcToJPsiMuMu_MuonPhys_2018_20230623/BcToJPsiMuMu_inclusive_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/230623_170305/0000/flatTuple_2*.root'},
#    'CF':{'filename':'root://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas/BcToJPsiMuMu_MuonPhys_2018_20230703_CF/BcToJPsiMuMu_inclusive_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/230703_140057/0000/flatTuple_2*.root'},
#    'CFv2':{'filename':'root://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas/BcToJPsiMuMu_MuonPhys_2018_20230704_CFv2/BcToJPsiMuMu_inclusive_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/230704_150336/0000/flatTuple_2*.root'}
    
    'normal':{'filename':prefix + 'job_pt_2018_MuonPhys/BcJpsiTau_inclusive/sig.root'},
    'CF':{'filename':prefix + 'job_pt_2018_MuonPhys_CF/BcJpsiTau_inclusive/Myroot.root'},
    'CFv2':{'filename':prefix + 'job_pt_2018_MuonPhys_CFv2/BcJpsiTau_inclusive/Myroot.root'}
}

finaldiscriminant = ['q2_simple', 'b_mass', 'b_vprob', 'tau_rhomass1', 'tau_rhomass2', 'b_eta', 'b_fls3d', 'b_alpha_zoom', 'b_pt', 'b_phi', 'b_lip', 'b_lips', 'b_pvip', 'b_pvips', 'b_fl3d', 'q2', 'mm2', 'ptmiss', 'estar', 'B_ptback']


for vkey, ivar in vardir.items():
    if vkey not in finaldiscriminant: 
        vardir.pop(vkey)




for varkey, vardir in vardir.items():

    
    hists = []
    titles = []
    
    idx = 0

    for key, var in filedir.items():
        
        print key
        
#        chain = TChain('tree', 'tree')
#        chain.Add(var['filename'])
        file = TFile(var['filename'])
        tree = file.Get('tree')
        
        
        #    file = TFile.Open(var['filename'])
        #    tree = file.Get('ntuplizer/tree')
        
        hist_ = TH1F(key + '_' + varkey, key + '_' + varkey, vardir['nbin'], vardir['xmin'], vardir['xmax'])
        hist_.GetXaxis().SetTitle(vardir['xtitle'])
        hist_.GetYaxis().SetTitle('a.u.')
        
        if 'var' in vardir:        
            tree.Draw(vardir['var'] + ' >> ' + hist_.GetName(), 'tau_isRight_3prong==1')
        else:
            tree.Draw(varkey + ' >> ' + hist_.GetName(), 'tau_isRight_3prong==1')

        
        hist_.SetLineColor(idx+1)
        hist_.SetMarkerColor(idx+1)
        idx += 1
        
        hists.append(copy.deepcopy(hist_))
        titles.append(key)
        
    comparisonPlots_alt(hists, titles, True, False, 'sergey_check/check_' + varkey + '.pdf', True, True, 'pE')

