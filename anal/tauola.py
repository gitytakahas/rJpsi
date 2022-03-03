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


lumi=59.6

#gROOT.SetBatch(False)
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
#gStyle.SetOptStat(0)



def comparisonPlots(hist, pname='sync.pdf', isLog=False, isRatio=True, clabel=''):

    display = DisplayManager(pname, isRatio, isLog, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)


def comparisonPlots_alt(hists, titles, norm=False, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he', prefix=None):

    display = DisplayManager_compare(pname, isLog, isRatio, 0.6, 0.7, doption)
    display.draw_legend = isLegend
    
    display.Draw(hists, titles, norm, prefix)

    

prefix='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Tauola/TauolaVariation/'

#ref = ['h3', 'h6']
#ref = ['h3']

maps = {
    'h1':{'title':'EvtGen2 + set1', 'color':2},
    'h2':{'title':'EvtGen1.6 + set0', 'color':3},
    'h3':{'title':'EvtGen1.6 + set1', 'color':1},
    'h04':{'title':'EvtGen2 + set0', 'color':4},
    'h4':{'title':'EvtGen2 + set1', 'color':2},
    'h5':{'title':'EvtGen1.6 + set0', 'color':3},
    'h6':{'title':'EvtGen1.6 + set1', 'color':1},
    'h07':{'title':'EvtGen2 + set0', 'color':4}
}

for ref, filename in zip(['h3', 'h6'], ['ntuple_test', 'ntuple_test_2']):
    
    file_prefix = ''
    file_down = 'h2'

    if ref=='h6':
        file_prefix = '_coarse'        
        file_down = 'h5'

    file = TFile(prefix + filename + '.root')

#    histnames = [key.GetName() for key in gDirectory.GetListOfKeys()]
    histnames = []

    if filename.find('2')==-1:
        histnames = ['h3', 'h2', 'h1', 'h04']
    else:
        histnames = ['h6', 'h5', 'h4', 'h07']

    hists = []
    titles = []

    for histname in histnames:

        print('check', filename, histname)
        
        hist = file.Get(histname)
        hist.SetLineColor(maps[histname]['color'])
        hist.SetMarkerColor(maps[histname]['color'])

        hists.append(copy.deepcopy(hist))
        titles.append(maps[histname]['title'])


    ensureDir('plots/tauola')
    comparisonPlots_alt(hists, titles, False, False, 'plots/tauola/compare' + file_prefix + '.pdf', True, True, 'hpE')


    # ratio
    
    ratios = []
    ratios_dir = {}
    titles = []
    
    for hist in hists:
        

        if hist.GetName() in ref: continue
        print hist.GetName()

        ratio = copy.deepcopy(hist)

        ratio.Divide(hists[0])
        ratio.GetXaxis().SetLabelSize(0.04)
        ratio.GetYaxis().SetRangeUser(0,10)
        ratio.SetLineColor(maps[hist.GetName()]['color'])
        ratio.SetMarkerColor(maps[hist.GetName()]['color'])
        
        for ibin in range(1, ratio.GetXaxis().GetNbins()+1):
            ratio.SetBinError(ibin, 0)


        ratios.append(copy.deepcopy(ratio))
        titles.append(maps[hist.GetName()]['title'])
        ratios_dir[hist.GetName()] = copy.deepcopy(ratio)

    comparisonPlots_alt(ratios, titles, False, False, 'plots/tauola/ratio' + file_prefix + '.pdf', False, True, 'hpE')


    envelope_up = TH1F('envelope_up', 'envelope_up', ratios[-1].GetXaxis().GetNbins(), ratios[-1].GetXaxis().GetXmin(), ratios[-1].GetXaxis().GetXmax())
    envelope_down = TH1F('envelope_down', 'envelope_down', ratios[-1].GetXaxis().GetNbins(), ratios[-1].GetXaxis().GetXmin(), ratios[-1].GetXaxis().GetXmax())

    os_envelope_up = TH1F('os_envelope_up', 'os_envelope_up', ratios[-1].GetXaxis().GetNbins(), ratios[-1].GetXaxis().GetXmin(), ratios[-1].GetXaxis().GetXmax())
    os_envelope_down = TH1F('os_envelope_down', 'os_envelope_down', ratios[-1].GetXaxis().GetNbins(), ratios[-1].GetXaxis().GetXmin(), ratios[-1].GetXaxis().GetXmax())

    for ibin in range(1, ratios[-1].GetXaxis().GetNbins()+1):
#        yvals = [ratio.GetBinContent(ibin) for ratio in ratios]


        val = ratios_dir[file_down].GetBinContent(ibin)

        envelope_down.SetBinContent(ibin, val)
        os_envelope_down.SetBinContent(ibin, val)

        os_envelope_up.SetBinContent(ibin, 1)

        val_up = -1
        
        if val >= 1:
            val_up = 1 - (val - 1)
        else:
            val_up = 1 + (1 - val)


        print 'check', ibin, val, val_up

#        envelope_up.SetBinContent(ibin, ratios_dir[file_down].GetBinContent(ibin))
        envelope_up.SetBinContent(ibin, val_up)

#        envelope_up.SetBinContent(ibin, 1)





#        print yvals, min(yvals), max(yvals)
#
#        envelope.SetBinContent(ibin, min(yvals))



    ensureDir('datacard/tauola')
    datacards = TFile('datacard/tauola/correction' + file_prefix + '.root','recreate')

    for ratio in ratios:
        ratio.Write()

    envelope_up.Write()
    envelope_down.Write()
    os_envelope_up.Write()
    os_envelope_down.Write()

    

