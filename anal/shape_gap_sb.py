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



#ensureDir('plots')
ensureDir('bkgcorr')

ratios = []
envelopes = []

for year in ['2016', '2017', '2018']:


    filenames = {'sb':'/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/results_simultaneous/' + year + '_sb_None/datacard/tau_rhomass_unrolled_var.root',
                 'gap':'/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/results_simultaneous/' + year + '_gap_None/datacard/tau_rhomass_unrolled_var.root'}

    hists = []
    titles = []

    hists2norm = []
    titles2norm = []
        
    
    ii = 0
    for key, filename in filenames.items():
    
        file = TFile(filename)
        hist = copy.deepcopy(file.Get(key + '/data_obs'))

        for sub in ['bc_jpsi_tau', 'bc_jpsi_dst', 'bc_others']:
            hist2subtract = file.Get(key + '/' + sub)
            if sub=='bc_jpsi_tau':
                hist.Add(hist2subtract,-0.25)
            else:
                hist.Add(hist2subtract,-1.)


        hist_ = normalize(hist, ii)

        hists.append(copy.deepcopy(hist_))
        titles.append(key + '(' + year + ')')


        hist2_ = modify(hist, ii)
#        print 'test=', hist.GetSumOfWeights(), hist2_.GetSumOfWeights()

        if key == 'sb':

            hist2ratios = {}

            for key2, filename2 in filenames.items():
    
                file2 = TFile(filename2)
                
                _tmp = copy.deepcopy(file2.Get(key2 + '/bg_ul'))

                hist2ratios[key2] = _tmp


            ratio = float(hist2ratios['gap'].GetSumOfWeights()/hist2ratios['sb'].GetSumOfWeights())

            print key2, 'ratio=', ratio
            hist2_.Scale(ratio)

        hists2norm.append(copy.deepcopy(hist2_))

        title = key + '(' + year + ')'

        if key == 'sb':
            title = key + '(extrapolated, ' + year + ')'

        titles2norm.append(title)



        ii += 1


    comparisonPlots_alt(hists, titles, False, False, 'bkgcorr/compare_data_obs_' + year + '.pdf', True, True, 'hpe')

    comparisonPlots_alt(hists2norm, titles2norm, False, False, 'bkgcorr/compare_data_obs_' + year + '_yield.pdf', True, True, 'hpe')

        
    ratio = copy.deepcopy(hists[1])
    ratio.Divide(hists[0])
    ratio.SetName('ratio_' + year)
    ratio.SetTitle('ratio_' + year)
    ratio.GetYaxis().SetRangeUser(0.,2.)
        
    envelope_up = TH1F('envelope_up_' + year, 'envelope_up_' + year, ratio.GetXaxis().GetNbins(), ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())
    envelope_down = TH1F('envelope_down_' + year, 'envelope_down_' + year, ratio.GetXaxis().GetNbins(), ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax())
        
    for ibin in range(1, ratio.GetXaxis().GetNbins()+1):
            
        val = ratio.GetBinContent(ibin) - 1.
            
        envelope_down.SetBinContent(ibin, 1+2*val)
        envelope_up.SetBinContent(ibin, 1-2*val)
            
    ratios.append(copy.deepcopy(ratio))
    envelopes.append(copy.deepcopy(envelope_up))
    envelopes.append(copy.deepcopy(envelope_down))



    canvas = TCanvas('can_'+year)

    frame = TH2F('ratio', 'ratio', ratio.GetXaxis().GetNbins(), ratio.GetXaxis().GetXmin(), ratio.GetXaxis().GetXmax(), 100,0,2)
    frame.GetXaxis().SetTitle('Tau rhomass unrolled bin ID')
    frame.GetYaxis().SetTitle('Ratio (after normalized)')
    frame.GetYaxis().SetNdivisions(504)
    frame.Draw()
    ratio.SetMarkerColor(1)
    ratio.SetLineColor(1)
    ratio.SetLineWidth(3)
    ratio.Draw("same")
    envelope_up.SetMarkerColor(2)
    envelope_down.SetMarkerColor(4)
    envelope_up.SetLineColor(2)
    envelope_down.SetLineColor(4)
    envelope_up.SetLineWidth(3)
    envelope_down.SetLineWidth(3)
    envelope_up.Draw('same')
    envelope_down.Draw('same')

    line = TLine(ratio.GetXaxis().GetXmin(), 1, ratio.GetXaxis().GetXmax(), 1)
    line.SetLineStyle(2)
    line.SetLineWidth(3)
    line.Draw('same')

    l2=add_CMS()
    l2.Draw("same")
    l3=add_Preliminary()
    l3.Draw("same")
    l4=add_lumi(year)
    l4.Draw("same")

    gPad.RedrawAxis();

    canvas.SaveAs('bkgcorr/bkgcompare_' + year +'.pdf')
    canvas.SaveAs('bkgcorr/bkgcompare_' + year + '.gif')





datacards = TFile('bkgcorr/correction.root','recreate')
for ratio in ratios:
    ratio.Write()
for envelope in envelopes:
    envelope.Write()

#ratio.Write()

#envelope_up.Write()
#envelope_down.Write()
datacards.Write()
datacards.Close()
