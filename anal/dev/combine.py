import copy, os, math, sys
from numpy import array
#from CMGTools.H2TauTau.proto.plotter.categories_TauMu import cat_Inc
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, Double, TCanvas, kRed, TH1D, TRandom3, kBird, THStack
#from common.DisplayManager_stack import DisplayManager
from officialStyle import officialStyle
from array import array
import common.MultiDraw
import numpy as np
from varConfig import vardir
from DataMCPlot import *
from DisplayManager_postfit import DisplayManager

lumi=59.6
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def comparisonPlots(hist, pname='sync.pdf', clabel='', isRatio=True):

    display = DisplayManager(pname, isRatio, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)


processes = ['signal', 'bg_ul', 'bg_bc', 'data_obs']

finaldiscriminant = ['q2', 'q2_simple', 'b_mass', 'b_mcorr', 'xgbs_zoom', 'estar', 'estar_simple', 'mm2', 'mm2_simple', 'tau_mass_zoom']

ensureDir('plots/combine')
ensureDir('datacard')



for discriminant in finaldiscriminant:

    hists = {}

    for dir in ['sr', 'sb', 'cr1', 'cr2']:
    
        file = TFile('datacard/' + discriminant + '_' + dir + '.root')

        for proc in processes:

            hist = file.Get(dir + '/' + proc)
            hist.SetName(proc)
            hist.SetTitle(proc)
            hists[proc + '_' + dir] = copy.deepcopy(hist)


    hist_cr2 = copy.deepcopy(hists['data_obs_cr2'])
    hist_cr2.Add(copy.deepcopy(hists['signal_cr2']), -1)
    hist_cr2.Add(copy.deepcopy(hists['bg_bc_cr2']), -1)
    
    hist_cr1 = copy.deepcopy(hists['data_obs_cr1'])
    hist_cr1.Add(copy.deepcopy(hists['signal_cr1']), -1)
    hist_cr1.Add(copy.deepcopy(hists['bg_bc_cr1']), -1)
    
    sf = hist_cr1.GetSumOfWeights()/hist_cr2.GetSumOfWeights()
    
    print 'sf = ', sf

    # prepare now the sideband. 


    hists4datacard = {}

    for dir in ['sr', 'sb']:

        if dir == 'sr':
            hist_sb = copy.deepcopy(hists['data_obs_sb'])
            hist_sb.Add(copy.deepcopy(hists['signal_sb']), -1)
            hist_sb.Add(copy.deepcopy(hists['bg_bc_sb']), -1)
            hist_sb.Scale(sf)
            hist_sb.SetName('bg_ul_v2')
            hist_sb.SetTitle('bg_ul_v2')

        Histo = DataMCPlot(discriminant + '_' + dir)

        hists4datacard_ = []

        for proc in processes:

            hist2add = copy.deepcopy(hists[proc + '_' + dir])

            hists4datacard_.append(hist2add)
        
            if proc=='bg_ul':
                hist2add = copy.deepcopy(hist_sb)
            

            hist2add.SetName(proc)
            hist2add.GetXaxis().SetLabelColor(1)
            hist2add.GetXaxis().SetLabelSize(0.0)

            if proc=='data_obs':
                hist2add.SetFillStyle(0)
                hist2add.Sumw2(False)
                hist2add.SetBinErrorOption(1)

            Histo.AddHistogram(proc, hist2add, 1)

            if proc=='data_obs':
                Histo.Hist(proc).stack = False


                
        # add BG histograms ... 
        for ibin in range(1, hist2add.GetXaxis().GetNbins()+1):

            dummy = TH1F('bg_bin' + str(ibin), 'bg_bin' + str(ibin),
                         hist2add.GetXaxis().GetNbins(), 
                         hist2add.GetXaxis().GetXmin(),
                         hist2add.GetXaxis().GetXmax())


            dummy.SetBinContent(ibin, 1)
        
            dummy.SetBinError(ibin, 0)
            
            hists4datacard_.append(copy.deepcopy(dummy))
            del dummy


        hists4datacard[dir] = hists4datacard_



        Histo._ApplyPrefs()
        print Histo

        comparisonPlots(Histo, 'plots/combine/' + discriminant  + '_' + dir + '.gif')
        

#        Histo.WriteDataCard('datacard/comb_' + discriminant + '.root', True, 'RECREATE', 'sr')

    ofile = TFile('datacard/comb_' + discriminant + '.root', 'recreate')

    for cat, lhist in hists4datacard.items():

        ofile.mkdir(cat)
        ofile.cd(cat)        
        
        for ihist in lhist:
            ihist.Write()
        


    ofile.Write()
    ofile.Close()
