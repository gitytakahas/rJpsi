import copy, os,  math, sys, shutil 
from numpy import array
from common.officialStyle import officialStyle
from array import array
import numpy as np
from varConfig import vardir
from common.DataMCPlot import *
from common.DisplayManager_postfit import DisplayManager
from common.DisplayManager_compare import DisplayManager_compare
from common.helper import *
from common.H2TauStyle import *
import ConfigParser
from optparse import OptionParser, OptionValueError

#gROOT.SetBatch(False)
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
#gStyle.SetOptStat(0)


def comparisonPlots(hist, lumi, pname='sync.pdf', isLog=False, isRatio=True, clabel=''):

    display = DisplayManager(pname, isRatio, isLog, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)

def comparisonPlots_alt(hists, titles, norm=False, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he'):

    display = DisplayManager_compare(pname, isLog, isRatio, 0.55, 0.75, doption)
    display.draw_legend = isLegend
    
    display.Draw(hists, titles, norm)

def add_lumi(luminumber):
    lowX=0.64
    lowY=0.842
    lumi  = TPaveText(lowX, lowY+0.06, lowX+0.65, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.055)
    lumi.SetTextFont (   42 )
#    lumi.AddText("2016, " + str(luminumber) + " fb^{-1} (13TeV)")
    lumi.AddText(str(luminumber) + " fb^{-1} (13TeV)")
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




hists = []

file = TFile('/work/ytakahas/work/analysis/CMSSW_10_2_10/src/rJpsi/anal/combine_sb3p5_sr4_simultaneous/tau_rhomass_unrolled_var.root')



#KEY: TH1Dbc_others;1bc_others
#KEY: TH1Djpsi_hc;1jpsi_hc
#KEY: TH1Djpsi_tau;1jpsi_tau
#KEY: TH1Dbg_ul;1bg_ul
#KEY: TH1Ddata_obs;1


#for proc in ['bc_others', 'bc_jpsi_dst', 'bc_jpsi_tau', 'bg_ul']:
#for proc in ['bc_others', 'bc_jpsi_dst', 'bc_jpsi_tau', 'bg_ul', 'dd_bkg', 'data_obs']:
for region in ['sr', 'sb']:
    for proc in ['bc_others', 'jpsi_hc', 'jpsi_tau', 'bg_ul', 'data_obs']:
        hists = []
        titles = []
    
        hist = file.Get('tauhad_' + region + '_2018/' + proc)

    
        hist_ = TH1F(hist.GetName(), hist.GetName(), hist.GetXaxis().GetNbins(), hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
        for ibin in range(1, hist.GetXaxis().GetNbins()+1):
        #        hist_.SetBinContent(ibin, float(hist.GetBinContent(ibin)/hist.GetSumOfWeights()))
            #        hist_.SetBinError(ibin, float(hist.GetBinError(ibin)/hist.GetSumOfWeights()))
            hist_.SetBinContent(ibin, float(hist.GetBinContent(ibin)))
            hist_.SetBinError(ibin, float(hist.GetBinError(ibin)))
        

        hist_.GetXaxis().SetTitle('Bin ID of the 2-dim. #rho_{1} vs #rho_{2}')
        hist_.SetFillColor(0)
        hist_.SetLineColor(hist.GetMarkerColor())
        hist_.SetMarkerColor(hist.GetMarkerColor())
        #    hist_.Scale(1./hist_.GetSumOfWeights())
        hists.append(copy.deepcopy(hist_))
        titles.append(proc)

        comparisonPlots_alt(hists, titles, False, False, 'display/' + proc + '_' + region + '.pdf', False, False, 'hpe')
