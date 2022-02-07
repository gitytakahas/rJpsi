import copy, os,  math, sys
from numpy import array
from officialStyle import officialStyle
from array import array
import common.MultiDraw
import numpy as np
from varConfig import vardir
from DataMCPlot import *
from DisplayManager_postfit import DisplayManager
from common.DisplayManager_compare import DisplayManager_compare
from helper import *
from H2TauStyle import *


lumi=59.6

#gROOT.SetBatch(False)
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
#gStyle.SetOptStat(0)


gROOT.ProcessLine(".L ~/tool//MultiDraw.cc+");


systs = []

for hammer in range(0, 9):
    systs.append('hammer_ebe_e' + str(hammer) + '_up')
    systs.append('hammer_ebe_e' + str(hammer) + '_down')


def applyHists(hists):

    
    colors = [1, 2, 4]

    for idx, hist in enumerate(hists):
        hist.SetFillStyle(0)
        hist.SetFillColor(0)
        hist.SetLineColor(colors[idx])
        hist.SetMarkerColor(colors[idx])
        hist.SetMarkerSize(0.5)
        hist.SetLineWidth(idx+1)
        hist.SetLineStyle(idx+1)


def comparisonPlots_alt(hists, titles, norm=False, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he'):

    display = DisplayManager_compare(pname, isLog, isRatio, 0.4, 0.7, doption)
    display.draw_legend = isLegend
    
    display.Draw(hists, titles, norm)


def getHist(vkey, channel, target, sys=None):

    filename = 'datacard/' + channel + '/' + vkey + '.root'
        
    if sys!=None:
        filename = 'datacard/' + channel + '/' + vkey + '_' + sys + '.root'

    file = TFile(filename)

    _hist = copy.deepcopy(file.Get(channel + '/' + target))

#    file.Close()
    
    return _hist


def setNameTitle(hist, name):
    hist.SetTitle(name)
    hist.SetName(name)

def draw(vkey, channels, target, sys=None, subtract=False, saveFig=False):

    hists = []
    titles = []
    
    for ii, channel in enumerate(channels):

        _hist = getHist(vkey, channel, target, sys)

#        filename = 'datacard/' + channel + '/' + vkey + '.root'
#        
#        if sys!=None:
#            filename = 'datacard/' + channel + '/' + vkey + '_' + sys + '.root'
#
#        file = TFile(filename)
#
#        _hist = copy.deepcopy(file.Get(channel + '/' + target))

        if subtract:

#            ['bg_bc', 'sig_others', 'sig_3p', 'bg_ul', 'data_obs']

            for proc in ['bg_bc', 'sig_others', 'sig_3p']:
            
#                if channel=='sr' and target=='data_obs' and proc=='sig_3p':
#                    print('THIS IS NOT CONSIDERED !!!!')
#                    continue

                _hist2 = getHist(vkey, channel, proc, sys) #copy.deepcopy(file.Get(channel + '/' + proc))
                _hist.Add(_hist2, -1)
        
        _hist.SetFillStyle(0)
        _hist.SetMarkerColor(ii+1)
        _hist.SetLineColor(ii+1)
        setNameTitle(_hist, target + '_' + channel)
        hists.append(_hist)
        titles.append(channel)

#        file.Close()



    hists2return = copy.deepcopy(hists)    


    if saveFig:
        _dirname = 'plots/compare/' + '_'.join(channels) + '_' + target
        ensureDir(_dirname)

        comparisonPlots_alt(hists, titles, True, False, _dirname + '/' + vkey +  '.pdf', True, True, 'hpe')


    _dirname = 'datacard/compare/' + '_'.join(channels) + '_' + target
    ensureDir(_dirname)

    file_output = TFile(_dirname + '/' + vkey +  '.root', 'recreate')
    for _hist in hists2return:
        _hist.Write()
    file_output.Write()
    file_output.Close()

        

    return hists2return


#finaldiscriminant = ['tau_rhomass_unrolled', 'tau_rhomass_unrolled_coarse']
finaldiscriminant = ['tau_rhomass_unrolled_coarse']
#finaldiscriminant = ['q2_simple']


for vkey, ivar in vardir.items():
    if vkey not in finaldiscriminant: 
        vardir.pop(vkey)

print '-'*80

for vkey, ivar in vardir.items():
    print '->', vkey.ljust(20), ivar

print '-'*80


ratio = 0.136

fitCat = 'sr'
ensureDir('syscompare/' + fitCat)


disc = 'tau_rhomass_unrolled_coarse'

for cr in ['sb_sr_bg_ul', 'lp_sb_data_obs']:
    file = TFile('datacard/compare/' + cr + '/' + disc + '.root')

    listofprocs = [key.GetName() for key in gDirectory.GetListOfKeys()]
    
    hists = []
    
    for proc in listofprocs:
        
        print(proc)

        hist = copy.deepcopy(file.Get(proc))
        hist.Scale(1./hist.GetSumOfWeights())

        hists.append(hist)


    ratio = copy.deepcopy(hists[1])
    ratio.Divide(hists[0])
    ratio.GetYaxis().SetTitle('ratio')
    ratio.GetYaxis().SetRangeUser(0,2)
    ratio.GetYaxis().SetNdivisions(506)
#    ratio.Fit('pol5')

    line = TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
    line.SetLineStyle(2)

    ensureDir('plots/closure/')

    canvas = TCanvas('closure')

    ratio.Draw()
    line.Draw('same')
    
    canvas.SaveAs('plots/closure/' + cr + '_' + disc + '.gif')
    



