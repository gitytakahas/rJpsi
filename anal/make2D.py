from ROOT import TFile, TH2D, TH1D, gROOT, gStyle, TCanvas, TPaveText, TLatex
from array import array
import copy
import json

from common.officialStyle import officialStyle

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)

gStyle.SetPadRightMargin (0.14)


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



import ConfigParser
init = ConfigParser.SafeConfigParser()
init.read("./settings.ini")
binning_rhomass1 = json.loads(init.get('common','binning_rhomass1'))
binning_rhomass2 = json.loads(init.get('common','binning_rhomass2'))

print len(binning_rhomass1)
print len(binning_rhomass2)

filenames = [
    ('data', '/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_2018/Data/data.root'),
    ('sig','/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_2018/BcJpsiTau_inclusive/sig.root'),
    ('bkg', '/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_2018/BJpsiX/bkg.root')
]

hists = []

for name, filename in filenames:
    
    file = TFile(filename)
    tree = file.Get('tree')

    hist = TH2D(name, name, len(binning_rhomass2)-1, array('d',binning_rhomass2), len(binning_rhomass1)-1, array('d', binning_rhomass1))
#    hist = TH2D(name, name, 20, 0.2, 1.4, 20, 0.2, 1.4)
    hist.GetXaxis().SetTitle('#rho_{2} mass (GeV)')
    hist.GetYaxis().SetTitle('#rho_{1} mass (GeV)')

#    tree.Draw("tau_rhomass1:tau_rhomass2 >> " + hist.GetName(), "tau_pt > 3. && mu1_isLoose==1 && mu2_isLoose==1 && xgbs < 3.5 && xgbs > 2.5")
    tree.Draw("tau_rhomass1:tau_rhomass2 >> " + hist.GetName(), "tau_pt > 3. && mu1_isLoose==1 && mu2_isLoose==1 && xgbs > 4.")
    
    hists.append(copy.deepcopy(hist))


    canvas = TCanvas(name)
    hist.SetMarkerSize(1.)
    hist.Draw('colztext')

    
#    l1=add_lumi(lumi)
#    l1.Draw("same")
    l2=add_CMS()
    l2.Draw("same")
    l3=add_Preliminary()
    l3.Draw("same")

    canvas.SaveAs(name+'_2d.pdf')


    hname_1d = name + '_1d'

    nbins_var = hist.GetXaxis().GetNbins()*hist.GetYaxis().GetNbins()

    hist_unrolled = TH1D(hname_1d, hname_1d, nbins_var, 0, nbins_var)

    idx_var = 1
    for iy in range(1, hist.GetYaxis().GetNbins()+1):
        for ix in range(1, hist.GetXaxis().GetNbins()+1):
            
            hist_unrolled.SetBinContent(idx_var, hist.GetBinContent(ix, iy))
            hist_unrolled.SetBinError(idx_var, hist.GetBinError(ix, iy))
                       
            idx_var += 1

    hists.append(copy.deepcopy(hist_unrolled))




print hists

ofile=TFile('out.root', 'recreate')

for hist in hists:
    hist.Write()




ofile.Write()
ofile.Close()


