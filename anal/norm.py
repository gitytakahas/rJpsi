from ROOT import TFile, TH2D, TH1D, gROOT, gStyle, TCanvas, TPaveText, TLatex, TH2F, TH1F, TF1
from array import array
import copy

from common.officialStyle import officialStyle

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)

#gStyle.SetPadRightMargin (0.14)


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

file = TFile('2018_inclusive_None/datacard/xgbs.root')
bg = copy.deepcopy(file.Get('inclusive/bg_ul'))
data = copy.deepcopy(file.Get('inclusive/data_obs'))



frame = TH2F('frame', 'frame', 100, -10,7, 100,0.00001,0.05)
frame.GetXaxis().SetTitle('BDT score')
frame.GetYaxis().SetTitle('a.u.')
frame.GetYaxis().SetNdivisions(504)


canvas = TCanvas('can')
#canvas.SetLogy()
frame.Draw()

hists = [data, bg]

hist2draw = []
for ii, hist in enumerate(hists):
#    hist.Sumw2()
#    hist.Scale(1./hist.GetSumOfWeights())

    hist_ = TH1F(hist.GetName(), hist.GetName(), hist.GetXaxis().GetNbins(), hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
    for ibin in range(1, hist.GetXaxis().GetNbins()+1):
        hist_.SetBinContent(ibin, float(hist.GetBinContent(ibin)/hist.GetSumOfWeights()))
        hist_.SetBinError(ibin, float(hist.GetBinError(ibin)/hist.GetSumOfWeights()))
        

    hist_.SetLineColor(ii+1)
    hist_.SetMarkerColor(ii+1)
    hist2draw.append(copy.deepcopy(hist_))


for hist in hist2draw:
    hist.Draw('lepsame')

l2=add_CMS()
l2.Draw("same")
l3=add_Preliminary()
l3.Draw("same")

canvas.SaveAs('compare.pdf')
canvas.SaveAs('compare.gif')

upper_boundary = 3.

rcanvas = TCanvas('can_ratio')

rframe = TH2F('ratio', 'ratio', 100, -10,7, 100,0.5,2.5)
rframe.GetXaxis().SetTitle('BDT score')
rframe.GetYaxis().SetTitle('Data / Bkg. MC')
rframe.GetYaxis().SetNdivisions(504)

rframe.Draw()

ratio = copy.deepcopy(data)
ratio.Divide(copy.deepcopy(bg))
ratio.Draw('same')

func = TF1('func', '[0]+[1]*x+[2]*x*x', 0, upper_boundary)
ratio.Fit('func', '', '', 0, upper_boundary)

fitline = ratio.GetFunction('func')
fitline.SetLineColor(2)
fitline.SetLineWidth(4)
fitline.Draw('same')


p1 = fitline.GetParameter(0)
p2 = fitline.GetParameter(1)
p3 = fitline.GetParameter(2)

fitted_all = TF1('fitted', str(p1) + '+ ' + str(p2) +' *x + ' + str(p3) + '*x*x', upper_boundary,7)
fitted_all.SetLineStyle(2)
fitted_all.SetLineColor(2)
fitted_all.SetLineWidth(2)

fitted_all.Draw('same')


l2=add_CMS()
l2.Draw("same")
l3=add_Preliminary()
l3.Draw("same")

rcanvas.SaveAs('ratio_compare.pdf')
rcanvas.SaveAs('ratio_compare.gif')


