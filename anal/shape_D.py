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
from ROOT import TPaveText, gPad

#gROOT.SetBatch(False)
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
#gStyle.SetOptStat(0)
lumi = 59.5

def add_lumi(year):
    lowX=0.8
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



def normalize(hist, ii):

    hist_ = TH1F(hist.GetName(), hist.GetName(), hist.GetXaxis().GetNbins(), hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())

    for ibin in range(1, hist.GetXaxis().GetNbins()+1):
        hist_.SetBinContent(ibin, float(hist.GetBinContent(ibin)/hist.GetSumOfWeights()))
        hist_.SetBinError(ibin, float(hist.GetBinError(ibin)/hist.GetSumOfWeights()))

    hist_.GetXaxis().SetTitle('Bin ID of the 2-dim. #rho_{1} vs #rho_{2}')
    hist_.SetFillColor(0)
    hist_.SetLineColor(ii+1)
    hist_.SetMarkerColor(ii+1)

    return copy.deepcopy(hist_)


def modify(hist, ii):

    hist_ = TH1F(hist.GetName(), hist.GetName(), hist.GetXaxis().GetNbins(), hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())

    for ibin in range(1, hist.GetXaxis().GetNbins()+1):
        hist_.SetBinContent(ibin, float(hist.GetBinContent(ibin)))
        hist_.SetBinError(ibin, float(hist.GetBinError(ibin)))

    hist_.GetXaxis().SetTitle('Bin ID of the 2-dim. #rho_{1} vs #rho_{2}')
    hist_.SetFillColor(0)
    hist_.SetLineColor(ii+1)
    hist_.SetMarkerColor(ii+1)
    hist_.GetXaxis().SetLabelSize(0.0)

    return copy.deepcopy(hist_)

def comparisonPlots(hist, lumi, pname='sync.pdf', isLog=False, isRatio=True, clabel=''):

    display = DisplayManager(pname, isRatio, isLog, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)

def comparisonPlots_alt(hists, titles, norm=False, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he'):

    display = DisplayManager_compare(pname, isLog, isRatio, 0.5, 0.7, doption)
    display.draw_legend = isLegend
    
    display.Draw(hists, titles, norm)

##def add_lumi(luminumber):
##    lowX=0.64
##    lowY=0.842
##    lumi  = TPaveText(lowX, lowY+0.06, lowX+0.65, lowY+0.16, "NDC")
##    lumi.SetBorderSize(   0 )
##    lumi.SetFillStyle(    0 )
##    lumi.SetTextAlign(   12 )
##    lumi.SetTextColor(    1 )
##    lumi.SetTextSize(0.055)
##    lumi.SetTextFont (   42 )
###    lumi.AddText("2016, " + str(luminumber) + " fb^{-1} (13TeV)")
##    lumi.AddText(str(luminumber) + " fb^{-1} (13TeV)")
##    return lumi
##
##
##def add_CMS():
##    lowX=0.15
##    lowY=0.87
##    lumi  = TPaveText(lowX, lowY, lowX+0.15, lowY+0.16, "NDC")
##    lumi.SetTextFont(61)
##    lumi.SetTextSize(0.05)
##    lumi.SetBorderSize(   0 )
##    lumi.SetFillStyle(    0 )
##    lumi.SetTextAlign(   12 )
##    lumi.SetTextColor(    1 )
##    lumi.AddText("CMS")
##    return lumi
##
##
##def add_channel(cl):
##    lowX=0.30
##    lowY=0.81
###    lumi  = ROOT.TLatex(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
###    lumi  = ROOT.TLatex(lowX, lowY,  cl)
##    lumi  = TLatex(lowX, lowY,  cl)
##    lumi.SetNDC(True)
###    lumi.SetTextFont(61)
##    lumi.SetTextSize(0.08)
###    lumi.SetBorderSize(   0 )
###    lumi.SetFillStyle(    0 )
###    lumi.SetTextAlign(   12 )
###    lumi.SetTextColor(    1 )
###    lumi.AddText(cl)
##    return lumi
##
##
##def add_Preliminary():
##    lowX=0.26
##    lowY=0.832
##    lumi  = TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
###    lumi.SetTextFont(52)
##    lumi.SetTextSize(0.05)
##    lumi.SetBorderSize(   0 )
##    lumi.SetFillStyle(    0 )
##    lumi.SetTextAlign(   12 )
##    lumi.SetTextColor(    1 )
###    lumi.AddText("#it{Simulation} Preliminary")
##    lumi.AddText("#it{Preliminary}")
##    return lumi



ensureDir('plots')

ratios = []
envelopes = []

#for year in ['2016', '2017', '2018']:
for year in ['2018']:


    filenames = {'sr':year + '_sr_None/datacard/tau_rhomass_unrolled_var.root',
                 'sb':year + '_sb_None/datacard/tau_rhomass_unrolled_var.root'}

    hists = []
    titles = []

    
    file = TFile(filenames['sb'])
    hist = copy.deepcopy(file.Get('sb/data_obs'))

    for sub in ['bc_jpsi_tau', 'bc_jpsi_dst', 'bc_others']:
        hist2subtract = file.Get('sb/' + sub)
        if sub=='bc_jpsi_tau':
            hist.Add(hist2subtract,-0.25)
        else:
            hist.Add(hist2subtract,-1.)

    hist_sb = modify(hist, 0)
    hist_sb.SetName('from_sb')
    hist_sb.SetTitle('from_sb')

    # calculate the extrapolation ratio

    hist2ratios = {}

    for key, filename in filenames.items():
        
        file = TFile(filename)
        
        _tmp = copy.deepcopy(file.Get(key + '/bg_ul'))
        
        hist2ratios[key] = _tmp


    ratio = float(hist2ratios['sr'].GetSumOfWeights()/hist2ratios['sb'].GetSumOfWeights())


    hist_sb.Scale(ratio)
    
    for ibin in range(1, hist_sb.GetXaxis().GetNbins()+1):
        
        err = hist_sb.GetBinError(ibin)
        
        tot = math.pow(err,2) + math.pow(0.3*hist_sb.GetBinContent(ibin), 2)

#        print ibin, err, tot

        hist_sb.SetBinError(ibin, math.sqrt(tot))

    print('ratio=', ratio)

    hists2write = []

    file = TFile(year + '_sr_None/datacard/tau_rhomass_unrolled_var.root')

    for sub in ['data_obs', 'bc_jpsi_tau', 'bc_jpsi_dst', 'bc_others']:
        hist_ = file.Get('sr/' + sub)

        
        if sub=='bc_jpsi_dst':
            for ibin in range(1, hist_.GetXaxis().GetNbins()+1):
                
                err = hist_.GetBinError(ibin)
                
                tot = math.pow(err,2) + math.pow(0.38*hist_.GetBinContent(ibin), 2)
                
                print ibin, err, tot
                
                hist_.SetBinError(ibin, math.sqrt(tot))



        hists2write.append(copy.deepcopy(hist_))


    hists2write.append(hist_sb)
    Histo = DataMCPlot('Dbkg')
    
    for _hist in hists2write:

        _name = _hist.GetName()

        if _name.find('data')!=-1:
            _hist.SetFillStyle(0)
            _hist.Sumw2(False)
            _hist.SetBinErrorOption(1)
            

        _hist.SetName(_name)
        _hist.SetTitle(_name)
        _hist.GetXaxis().SetLabelColor(1)
        _hist.GetXaxis().SetLabelSize(0)

        Histo.AddHistogram(_hist.GetName(), _hist)
        print ( " Histo = DataMCPlot(vkey) with _hist ", _hist.GetName() ," bins ", _hist.GetNbinsX() )
        if _name.find('data')!=-1:
            Histo.Hist(_hist.GetName()).stack = False



    Histo._ApplyPrefs()
    print(Histo)
    comparisonPlots(Histo, lumi, 'plots/Dcheck_' + year +'.gif')



#    comparisonPlots_alt(hists, titles, False, False, 'plots/compare_data_obs_' + year + '.pdf', True, True, 'hpe')

#    comparisonPlots_alt(hists2norm, titles2norm, False, False, 'plots/compare_data_obs_' + year + '_yield.pdf', True, True, 'hpe')


