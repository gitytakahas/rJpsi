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
import json 

init = ConfigParser.SafeConfigParser()
init.read("./settings.ini")
binning_rhomass1 = json.loads(init.get('common','binning_rhomass1'))
binning_rhomass2 = json.loads(init.get('common','binning_rhomass2'))

print 'len(binning_rhomass1)', len(binning_rhomass1)
print 'len(binning_rhomass2)', len(binning_rhomass2)

lumi=59.6

def comparisonPlots(hist, pname='sync.pdf', isLog=False, isRatio=True, clabel=''):

    display = DisplayManager(pname, isRatio, isLog, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)


def comparisonPlots_alt(hists, titles, norm=False, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he', prefix=None):

    display = DisplayManager_compare(pname, isLog, isRatio, 0.6, 0.7, doption)
    display.draw_legend = isLegend
    
    display.Draw(hists, titles, norm, prefix)




#gROOT.SetBatch(False)
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
#gStyle.SetOptStat(0)

prefix = '/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Tauola/TauolaVariation/'
hdict = {
    'evt2_set0':{'file':'BPH-RunIIFall18GS-00271_set0_EvtGen2_new.root', 'hist':None, 'hist2d':None},
    'evt2_set1':{'file':'BPH-RunIIFall18GS-00271_set1_EvtGen2_new.root', 'hist':None, 'hist2d':None},
    'evt1_set0':{'file':'BPH-RunIIFall18GS-00271_set0_new.root', 'hist':None, 'hist2d':None},
    'evt1_set1':{'file':'BPH-RunIIFall18GS-00271_set1_new.root', 'hist':None, 'hist2d':None},
}


for key, var in hdict.items():
    file = TFile(prefix + var['file'])
    
    tree = file.Get('ZmuonAnalyzer/treemc')
    
    hname_2d = key + '_2d'
    
    hist_2d = TH2D(hname_2d, hname_2d, len(binning_rhomass2)-1, array('d',binning_rhomass2), len(binning_rhomass1)-1, array('d', binning_rhomass1))
    hist_2d.Sumw2()

    tree.Draw('truth_tau_dipion1_mass:truth_tau_dipion2_mass >> ' + hist_2d.GetName())

    # make it 1d ... 
    nbins_var = hist_2d.GetXaxis().GetNbins()*hist_2d.GetYaxis().GetNbins()
#    nbins_var2 = int((len(binning_rhomass1)-1)*(len(binning_rhomass2)-1))

    print '# of bins =', nbins_var #nbins_var2
    
    hname_1d = key + '_1d'

    hist_unrolled = TH1D(hname_1d, hname_1d, nbins_var, 0, nbins_var)

    idx_var = 1

    for iy in range(1, hist_2d.GetYaxis().GetNbins()+1):
        for ix in range(1, hist_2d.GetXaxis().GetNbins()+1):
            hist_unrolled.SetBinContent(idx_var, hist_2d.GetBinContent(ix, iy))
            hist_unrolled.SetBinError(idx_var, hist_2d.GetBinError(ix, iy))

#            print ix, iy, hist_2d.GetBinContent(ix,iy), hist_2d.GetBinError(ix,iy)
            
            idx_var += 1

#            multihists[hname_1d] = copy.deepcopy(hist_unrolled)


    hdict[key]['hist'] = copy.deepcopy(hist_unrolled)
    hdict[key]['hist2d'] = copy.deepcopy(hist_2d)


    

#prefix='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/Tauola/TauolaVariation/'

#ref = ['h3', 'h6']
#ref = ['h3']

#maps = {
#    'h1':{'title':'EvtGen2 + set1', 'color':2},
#    'h2':{'title':'EvtGen1.6 + set0', 'color':3},
#    'h3':{'title':'EvtGen1.6 + set1', 'color':1},
#    'h04':{'title':'EvtGen2 + set0', 'color':4},
#    'h4':{'title':'EvtGen2 + set1', 'color':2},
#    'h5':{'title':'EvtGen1.6 + set0', 'color':3},
#    'h6':{'title':'EvtGen1.6 + set1', 'color':1},
#    'h07':{'title':'EvtGen2 + set0', 'color':4}
#}

##for ref, filename in zip(['h3', 'h6'], ['ntuple_test', 'ntuple_test_2']):
##    
##    file_prefix = ''
##    file_down = 'h2'
##
##    if ref=='h6':
##        file_prefix = '_coarse'        
##        file_down = 'h5'
##
##    file = TFile(prefix + filename + '.root')
##
###    histnames = [key.GetName() for key in gDirectory.GetListOfKeys()]
##    histnames = []
##
##    if filename.find('2')==-1:
##        histnames = ['h3', 'h2', 'h1', 'h04']
##    else:
##        histnames = ['h6', 'h5', 'h4', 'h07']
##
##    hists = []
##    titles = []
##
##    for histname in histnames:
##
##        print('check', filename, histname)
##        
##        hist = file.Get(histname)
##        hist.SetLineColor(maps[histname]['color'])
##        hist.SetMarkerColor(maps[histname]['color'])
##
##        hists.append(copy.deepcopy(hist))
##        titles.append(maps[histname]['title'])
##
##
##    ensureDir('plots/tauola')
##    comparisonPlots_alt(hists, titles, False, False, 'plots/tauola/compare' + file_prefix + '.pdf', True, True, 'hpE')


    # ratio
    
##    ratios = []
##    ratios_dir = {}
##    titles = []
##    
##    for hist in hists:
##        
##
##        if hist.GetName() in ref: continue
##        print hist.GetName()
##
##        ratio = copy.deepcopy(hist)
##
##        ratio.Divide(hists[0])
##        ratio.GetXaxis().SetLabelSize(0.04)
##        ratio.GetYaxis().SetRangeUser(0,10)
##        ratio.SetLineColor(maps[hist.GetName()]['color'])
##        ratio.SetMarkerColor(maps[hist.GetName()]['color'])
##        
##        for ibin in range(1, ratio.GetXaxis().GetNbins()+1):
##            ratio.SetBinError(ibin, 0)
##
##
##        ratios.append(copy.deepcopy(ratio))
##        titles.append(maps[hist.GetName()]['title'])
##        ratios_dir[hist.GetName()] = copy.deepcopy(ratio)
##
##    comparisonPlots_alt(ratios, titles, False, False, 'plots/tauola/ratio' + file_prefix + '.pdf', False, True, 'hpE')


ratio = copy.deepcopy(hdict['evt1_set0']['hist'])
ratio.Divide(copy.deepcopy(hdict['evt1_set1']['hist'])) ### This is a reference i.e. default option

envelope_up = TH1F('envelope_up', 'envelope_up', ratio.GetXaxis().GetNbins(), ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())
envelope_down = TH1F('envelope_down', 'envelope_down', ratio.GetXaxis().GetNbins(), ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())


# one-sided
os_envelope_up = TH1F('os_envelope_up', 'os_envelope_up', ratio.GetXaxis().GetNbins(), ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())
os_envelope_down = TH1F('os_envelope_down', 'os_envelope_down', ratio.GetXaxis().GetNbins(), ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())

for ibin in range(1, ratio.GetXaxis().GetNbins()+1):

    val = ratio.GetBinContent(ibin)
    
    envelope_down.SetBinContent(ibin, val)
    os_envelope_down.SetBinContent(ibin, val)

    os_envelope_up.SetBinContent(ibin, 1)

    val_up = -1
        
    if val >= 1:
        val_up = 1 - (val - 1)
    else:
        val_up = 1 + (1 - val)
            

#    print 'check', ibin, val, val_up

    envelope_up.SetBinContent(ibin, val_up)



ensureDir('datacard/tauola')
datacards = TFile('datacard/tauola/correction.root','recreate')

for key, hist_ in hdict.items():
    hist_['hist'].Write()
    hist_['hist2d'].Write()

ratio.Write()

envelope_up.Write()
envelope_down.Write()
os_envelope_up.Write()
os_envelope_down.Write()

    

