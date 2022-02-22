import os 
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, Double, TCanvas, kRed, TH1D, TRandom3, kBird, THStack, TGraph, gDirectory, kBlue

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def overflow(hist):
    lastp1 = hist.GetBinContent(hist.GetXaxis().GetNbins()+1)
    last = hist.GetBinContent(hist.GetXaxis().GetNbins())
    lastp1e = hist.GetBinError(hist.GetXaxis().GetNbins()+1)
    laste = hist.GetBinError(hist.GetXaxis().GetNbins())
    hist.SetBinContent(hist.GetXaxis().GetNbins(), last+lastp1)
    hist.SetBinError(hist.GetXaxis().GetNbins(), math.sqrt(math.pow(laste,2)+math.pow(lastp1e,2)))
    hist.SetBinContent(hist.GetXaxis().GetNbins()+1, 0)
    hist.SetBinError(hist.GetXaxis().GetNbins()+1, 0)

    firstp1 = hist.GetBinContent(1)
    first = hist.GetBinContent(0)
    firstp1e = hist.GetBinError(1)
    firste = hist.GetBinError(0)
    hist.SetBinContent(1, first+firstp1)
    hist.SetBinError(1, math.sqrt(math.pow(firste,2)+math.pow(firstp1e,2)))


colours = [1, 2, 4, 6, 8, 1, 46, 13, 15, 1,1,1,1,1,1]
styles = [1, 2, 4, 3, 5, 2,7, 1, 1, 1,1,1,1,1,1]

def applyHistStyle(h, i):

    if h.GetName().find('data')!=-1:
        h.SetLineColor(1)
        h.SetMarkerStyle(20)
        h.SetMarkerColor(1)
    else:
        h.SetLineColor(colours[i+1])
        h.SetMarkerColor(colours[i+1])
        h.SetMarkerSize(0)
        h.SetLineStyle(styles[i+1])

#    for ibin in range(1, h.GetXaxis().GetNbins()+1):
#        h.SetBinError(ibin, 0)

    h.GetXaxis().SetLabelSize(0.05)
    h.GetYaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    h.SetMarkerSize(0.2)
    h.SetLineWidth(3)
    h.SetStats(False)


def applyLegendSettings(leg):
    leg.SetBorderSize(0)
    leg.SetFillColor(10)
    leg.SetLineColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
#    leg.SetTextFont(42)
