import copy, os, math, sys
from numpy import array
#from CMGTools.H2TauTau.proto.plotter.categories_TauMu import cat_Inc
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, Double, TCanvas, kRed, TH1D, TRandom3, kBird, THStack
from common.DisplayManager import DisplayManager
from officialStyle import officialStyle
from array import array
import common.MultiDraw
import numpy as np
#from varConfig import vardir
from DataMCPlot import *

lumi=59.6
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def comparisonPlots(hists, titles, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he', prefix=None):

    display = DisplayManager(pname, isLog, isRatio, 0.2, 0.7, doption)
    display.draw_legend = isLegend
#    display.isEff = isEff
    
    display.Draw(hists, titles, prefix)




def calc(nsig, ndata):


    sig = 0
    sig_err = 0
    frac = 0
    frac_err = 0


    if ndata!=0:

        sig = nsig/math.sqrt(ndata)
        sig_err = math.sqrt(nsig/ndata + nsig*nsig/(4*ndata*ndata))       

        frac = nsig/ndata
        frac_err = math.sqrt(nsig/(ndata*ndata) + nsig*nsig/(ndata*ndata*ndata))

    return sig, sig_err, frac, frac_err


#processes = ['signal', 'bg_ul', 'bg_bc', 'data_obs']
processes = ['signal', 'data_obs']


datacards = {
    'mass':{'file':'datacard_mass/xgbs_fine_sr.root', 'leg':'mass', 'col':1},
    'nomass':{'file':'datacard_nomass/xgbs_fine_sr.root', 'leg':'no mass', 'col':2},
    'multiple':{'file':'datacard/xgbs_fine_sr.root', 'leg':'multiple', 'col':3},
}

ensureDir('plots/sb')


h_purity = []
h_sig = []
titles = []

dir = 'sr'

for dkey, var in datacards.items():
    
    file = TFile(var['file'])

    hists = {}

    for proc in processes:
        
        hist = file.Get(dir + '/' + proc)
        hist.SetName(proc)
        hist.SetTitle(proc)
        hists[proc + '_' + dir] = copy.deepcopy(hist)


    h_sig_truth = copy.deepcopy(hists['signal_' + dir])      
    h_sig_truth.SetMarkerColor(var['col'])
    h_sig_truth.SetLineColor(var['col'])

    h_significance_truth = copy.deepcopy(h_sig_truth)
    h_significance_truth.SetMarkerColor(var['col'])
    h_significance_truth.SetLineColor(var['col'])

    for ibin in range(1, hists['data_obs_' + dir].GetXaxis().GetNbins()+1):
        ndata = hists['data_obs_' + dir].Integral(ibin, hists['data_obs_' + dir].GetXaxis().GetNbins()+1)
        ntruth = h_sig_truth.Integral(ibin, hists['data_obs_' + dir].GetXaxis().GetNbins()+1)

        sig, sig_err, frac, frac_err = calc(ntruth, ndata)

        h_significance_truth.SetBinContent(ibin, sig)
        h_significance_truth.SetBinError(ibin, sig_err)
        
        h_sig_truth.SetBinContent(ibin, frac)
        h_sig_truth.SetBinError(ibin, frac_err)
        
        h_sig_truth.GetXaxis().SetLabelSize(0.05)
        h_significance_truth.GetXaxis().SetLabelSize(0.05)
        h_sig_truth.GetYaxis().SetTitle('S/(S+B)')
        h_significance_truth.GetYaxis().SetTitle('S/sqrt(S+B)')
            
#        h_sig_truth.Scale(2./h_sig_truth.GetMaximum())
#        h_significance_truth.Scale(2./h_significance_truth.GetMaximum())
            

    h_purity.append(h_sig_truth)
    h_sig.append(h_significance_truth)
    titles.append(dkey)
    


######################
### for multiple option
######################


h_sig_truth = copy.deepcopy(h_purity[-1])
h_sig_truth.SetMarkerColor(4)
h_sig_truth.SetLineColor(4)
h_sig_truth.SetLineStyle(2)

h_significance_truth = copy.deepcopy(h_sig[-1])
h_significance_truth.SetMarkerColor(4)
h_significance_truth.SetLineColor(4)
h_significance_truth.SetLineStyle(2)

h_sig_truth.GetXaxis().SetLabelSize(0.05)
h_sig_truth.GetYaxis().SetTitle('S/(S+B)')

h_significance_truth.GetXaxis().SetLabelSize(0.05)
h_significance_truth.GetYaxis().SetTitle('S/sqrt(S+B)')


ibin = 1

for xgbs in range(0, 31):
#for xgbs in range(1, 2):

    xgbs = 6 + xgbs*0.2
    xgbs = str(xgbs).replace('.','p')

    print('multiple reading file = ', 'datacard_' + str(xgbs) + '/xgbs_fine_sr.root')
    file = TFile('datacard_' + str(xgbs) + '/xgbs_fine_sr.root')

    hists = {}

    dir = 'sr'

    for proc in processes:
        
        hist = file.Get(dir + '/' + proc)
        hist.SetName(proc)
        hist.SetTitle(proc)
        hists[proc + '_' + dir] = copy.deepcopy(hist)

    ndata = hists['data_obs_' + dir].Integral()
    ntruth = hists['signal_' + dir].Integral()

    sig, sig_err, frac, frac_err = calc(ntruth, ndata)

    h_significance_truth.SetBinContent(ibin, sig)
    h_significance_truth.SetBinError(ibin, sig_err)
    
    h_sig_truth.SetBinContent(ibin, frac)
    h_sig_truth.SetBinError(ibin, frac_err)

    ibin += 1
        
h_purity.append(h_sig_truth)
h_sig.append(h_significance_truth)
titles.append("multiple highest")
    




print h_purity
print h_sig
print titles
comparisonPlots(h_purity, titles, False, 'plots/sb/purity.pdf', False, True, 'pzle')

comparisonPlots(h_sig, titles, False, 'plots/sb/significance.pdf', False, True, 'pzle')
