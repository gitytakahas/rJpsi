import numpy, os, math, copy
from ROOT import Double, gROOT, gStyle, TH1F, TChain, TFile, TGraph, TCanvas, TH2F, TLegend, Double

from common.DisplayManager import DisplayManager, applyLegendSettings
from common.officialStyle import officialStyle

from optparse import OptionParser, OptionValueError

usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

colours = [1, 2, 4, 6, 8, 9, 46, 13, 15, 1,1,1,1,1,1]
styles = [1, 2, 4, 3, 5, 6,7, 1, 1, 1,1,1,1,1,1]


def comparisonPlots(hists, titles, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he', prefix=None):

    display = DisplayManager(pname, isLog, isRatio, 0.6, 0.7, doption)
    display.draw_legend = isLegend
#    display.isEff = isEff
    
    display.Draw(hists, titles, prefix)


def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def overflow(hist):
#    import pdb; pdb.set_trace()
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


def rocCurve(hS, hB, title, option=None):
    ''' Create a ROC TGraph from two input histograms.'''
    
    maxBin = hS.GetNbinsX()
    
    if hS.Integral() == 0.:
        print('ROC curve creator, hist', hS.GetName(), 'has zero entries')
        return
        
    if option==None:
        effsS = [hS.Integral(0, nBin)/hS.Integral(0, maxBin+1) for nBin in range(0, maxBin + 1)]
        rejB = [hB.Integral(0, nBin)/hB.Integral(0, maxBin+1) for nBin in range(0, maxBin + 1)]
    else:
        effsS = [hS.Integral(nBin, maxBin+1)/hS.Integral(0, maxBin+1) for nBin in range(0, maxBin + 1)]
        rejB = [hB.Integral(nBin, maxBin+1)/hB.Integral(0, maxBin+1) for nBin in range(0, maxBin + 1)]
        

    print maxBin, len(effsS), len(rejB)
    rocCurve = TGraph(maxBin, numpy.asarray(effsS), numpy.asarray(rejB))


    rocCurve.SetName(title)

    return rocCurve


rocs = []


for version in ['v1', 'v2']:

    graph = TGraph()
    graph.SetName(version)
    graph.SetTitle(version)

    idx = 0

    for line in open('model_multiple_' + version + '/training_results_roc_csv_None.csv', 'r').readlines():
        
        line = line.rstrip()
        
        if line.find('thresholds')!=-1: continue
        
#        print line 
        fake = line.split(',')[1]
        eff = line.split(',')[2]
    
        graph.SetPoint(idx, Double(fake), Double(eff))

        idx += 1

    rocs.append(copy.deepcopy(graph))



ensureDir('Plots')

can_roc = TCanvas()
can_roc.SetLogx()

frame = TH2F('frame_roc', 'frame_roc', 10000,0.00001,1,10000,0, 1)
frame.GetXaxis().SetTitle('False Alarm Rate')
frame.GetYaxis().SetTitle('Signal efficiency')
frame.Draw()

leg = TLegend(0.2, 0.7, 0.5, 0.9)
applyLegendSettings(leg)

for iroc, roc in enumerate(rocs):
    roc.SetLineColor(colours[iroc])
    roc.SetLineStyle(styles[iroc])
    roc.SetLineWidth(3)
    roc.SetMarkerSize(0)
    
    leg.AddEntry(roc, roc.GetName(), 'l')
    roc.Draw('plsame')

leg.Draw()

can_roc.SaveAs('Plots/roc.gif')

#ofile.Write()
#ofile.Close()
