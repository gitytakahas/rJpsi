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
        hists.append(_hist)
        titles.append(channel)

#        file.Close()



    hists2return = copy.deepcopy(hists)    

    if saveFig:
        _dirname = 'plots/compare/' + '_'.join(channels) + '_' + target
        ensureDir(_dirname)
        comparisonPlots_alt(hists, titles, True, False, _dirname + '/' + vkey +  '.pdf', True, True, 'hpe')

    return hists2return


finaldiscriminant = ['tau_rhomass_unrolled']
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



for vkey, ivar in vardir.items():

    # first compare the mc in sr and sb ... 
#    ratio1 = draw(vkey, ['hp_sb', 'hp_sr'], 'bg_ul')
#
#    ratio2 = draw(vkey, ['lp_sb', 'lp_sr'], 'data_obs', True)
#
#    ratio3 = draw(vkey, ['slp_sb', 'slp_sr'], 'data_obs', True)
#
#    ratio4 = draw(vkey, ['cr_hp_sb', 'cr_hp_sr'], 'data_obs', True)
#
#    ratio5 = draw(vkey, ['cr_lp_sb', 'cr_lp_sr'], 'data_obs', True)
#
#    ratio6 = draw(vkey, ['cr_slp_sb', 'cr_slp_sr'], 'data_obs', True)
#
#    ratio7 = draw(vkey, ['lp_sb', 'slp_sb', 'hp_sb'], 'data_obs', True)
#
#    ratio8 = draw(vkey, ['cr_lp_sb', 'cr_slp_sb'], 'data_obs', True)
#
#    ratio9 = draw(vkey, ['cr_lp_sr', 'cr_slp_sr'], 'data_obs', True)
#
#    ratio10 = draw(vkey, ['lp_sb', 'cr_lp_sb'], 'data_obs', True)
#
#    ratio11 = draw(vkey, ['lp_sr', 'cr_lp_sr'], 'data_obs', True)
#
#    ratio12 = draw(vkey, ['hp_sb'], 'data_obs', True, ratio2)

#    hists4ratio = draw(vkey, ['sr', 'sb'], 'bg_ul')
    
#    draw(vkey, ['cr_lp',  'cr_sb'], 'data_obs', True)

#    draw(vkey, ['cr_sb',  'cr_sr'], 'data_obs', True)






    draw(vkey, ['sb', 'sr'], 'bg_ul', None, False, True)

    draw(vkey, ['lp', 'sb'], 'data_obs', None, True, True)
    

#    if vkey=='tau_rhomass_unrolled':
    if vkey in finaldiscriminant:

        hists2write = []


        file = TFile('datacard/' + fitCat + '/' + vkey + '.root')
        file.cd(fitCat)
        
        listofprocs = [key.GetName() for key in gDirectory.GetListOfKeys()]
    
        for proc in listofprocs:
            _tmp = file.Get(fitCat + '/' + proc)
            hists2write.append(copy.deepcopy(_tmp))

        hists4ddbkg = draw(vkey, ['sr', 'sb'], 'data_obs', None, True, False)
        bkgHist = hists4ddbkg[1]        
        bkgHist.Scale(ratio)
        setNameTitle(bkgHist, 'dd_bkg')
        hists2write.append(bkgHist)

        ### reference
#        sig_3p_cent = getHist(vkey, fitCat, 'sig_3p')
#        sig_others_cent = getHist(vkey, fitCat, 'sig_others')
#        dd_bkg_cent = copy.deepcopy(bkgHist)


        #### shape variation

        for sys in systs:

            name_sys = sys.replace('_up','Up').replace('_down', 'Down')

            hists_sys = draw(vkey, ['sr', 'sb'], 'data_obs', sys, True, False)

            bkgHist_sys = hists_sys[1]
            bkgHist_sys.Scale(ratio)

            setNameTitle(bkgHist_sys, 'dd_bkg_' + name_sys)
            hists2write.append(bkgHist_sys)

            sig_3p_sys = getHist(vkey, fitCat, 'sig_3p', sys)
            sig_others_sys = getHist(vkey, fitCat, 'sig_others', sys)

            setNameTitle(sig_3p_sys, 'sig_3p_' + name_sys)
            setNameTitle(sig_others_sys, 'sig_others_' + name_sys)
            
            hists2write.append(sig_3p_sys)
            hists2write.append(sig_others_sys)

            
            
#            comparisonPlots(hists, titles, cat + ', ' + syst.replace('_up',''), 'syscompare/' + cat + '/' + syst + '.gif')
            

        ### CREATE shape comparisons 

        filename_new = 'datacard/'+ fitCat + '/' + vkey + '_new.root'
        file_new = TFile(filename_new, 'recreate')
        file_new.mkdir(fitCat)
        file_new.cd(fitCat)
        
        for hist2write in hists2write:
            hist2write.Write()

        file_new.Write()
        file_new.Close()



        ### shape comparison!
        file = TFile(filename_new)
        file.cd(fitCat)
        
        listofprocs = [key.GetName() for key in gDirectory.GetListOfKeys()]
    
        for proc in listofprocs:

            if proc.find('Up')==-1: continue

            strs = proc.split('_')

            procname = '_'.join(strs[:2])
            sysname = '_'.join(strs[2:]).replace('Up', '')

            print(procname, sysname)

            _cent = fitCat + '/' + procname
            _down = fitCat + '/' + proc.replace('Up','Down')
            _up = fitCat + '/' + proc

            print('comparing,', _cent, _down, _up)

            cent = copy.deepcopy(file.Get(_cent))
            down = copy.deepcopy(file.Get(_down))
            up = copy.deepcopy(file.Get(_up))
            
            hists = [cent, down, up]
            titles = [procname, proc, proc.replace('Up','Down')]

            applyHists(hists)

            _dirname = 'plots/syscompare'
            ensureDir(_dirname)
            comparisonPlots_alt(hists, titles, False, False, _dirname + '/' + vkey + '_' + procname + '_' + sysname + '.pdf', True, True, 'hpe')


        




