from ROOT import TFile, TH2D, TH1D, gROOT, gStyle, TCanvas, TPaveText, TLatex
from array import array
import copy
import json

from common.officialStyle import officialStyle
from common.helper import *

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

basic = init.get('common', 'basic')

xgbs_sr = 'xgbs > ' + init.get('common', 'sr_low') #+ '&& b_mass > 5.5'
xgbs_sb = 'xgbs > ' + init.get('common', 'sb_low') + ' && xgbs < ' + init.get('common', 'sb_high') #+ '&& b_mass > 5.5'

print len(binning_rhomass1)
print len(binning_rhomass2)

dirname = '/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_2018_MuonPhys'

dirs = {
#    ('data', dirname + '/Data/data.root'),
    'sig_dm10_sr':{'filename':dirname + '/BcJpsiTau_inclusive/sig.root', 'selection':basic + '&&' + xgbs_sr + '&& tau_isRight_3prong == 1 && decayid==10'},
    'sig_dm11_sr':{'filename':dirname + '/BcJpsiTau_inclusive/sig.root', 'selection':basic + '&&' + xgbs_sr + '&& tau_isRight_3prong == 1 && decayid==11'},
    'bkg_sr':{'filename':dirname + '/BJpsiX/bkg.root', 'selection':basic + '&&' + xgbs_sr},
    'bkg_sb':{'filename':dirname + '/BJpsiX/bkg.root', 'selection':basic + '&&' + xgbs_sb},
}

ensureDir('homework')

hists = []

for name, val in dirs.items():

   
    file = TFile(val['filename'])

    print val['filename']

    tree = file.Get('tree')

    hist_coarse = TH2D(name + '_coarse', name + '_coarse', len(binning_rhomass2)-1, array('d',binning_rhomass2), len(binning_rhomass1)-1, array('d', binning_rhomass1))
    hist_coarse.GetXaxis().SetTitle('#rho_{2} mass (GeV)')
    hist_coarse.GetYaxis().SetTitle('#rho_{1} mass (GeV)')


    hist = TH2D(name, name, 30, 0.2, 1.6, 30, 0.2, 1.6)
    hist.GetXaxis().SetTitle('#rho_{2} mass (GeV)')
    hist.GetYaxis().SetTitle('#rho_{1} mass (GeV)')
    
    selstr = None
#
#    if name=='sig' and postfix.find('dm10')!=-1:
#        selstr = basic + "&&"+  xgbs_sr + '&& tau_isRight_3prong == 1 && decayid==10'
#    elif name=='sig' and postfix.find('dm11')!=-1:
#        selstr = basic + "&&"+  xgbs_sr + '&& tau_isRight_3prong == 1 && decayid==11'
#    else:
#        selstr = basic + "&&"+  xgbs_sr + '&& tau_isRight_3prong == 1 && decayid==11'

#    selstr = basic + '&& tau_isRight_3prong == 0'

    print "tau_rhomass1:tau_rhomass2 >> " + hist.GetName(), val['selection']
    tree.Draw("tau_rhomass1:tau_rhomass2 >> " + hist.GetName(), val['selection'])
    tree.Draw("tau_rhomass1:tau_rhomass2 >> " + hist_coarse.GetName(), val['selection'])
    
    hists.append(copy.deepcopy(hist))
    hists.append(copy.deepcopy(hist_coarse))

    canvas = TCanvas(name)
    hist.SetMarkerSize(1.)
    hist.Draw('colztext')

    
#    l1=add_lumi(lumi)
#    l1.Draw("same")
    l2=add_CMS()
    l2.Draw("same")
    l3=add_Preliminary()
    l3.Draw("same")

    canvas.SaveAs('homework/' + name+'_2d.pdf')
    canvas.SaveAs('homework/' + name+'_2d.gif')


    hname_1d = name + '_1d'

    nbins_var = hist_coarse.GetXaxis().GetNbins()*hist_coarse.GetYaxis().GetNbins()

    hist_unrolled = TH1D(hname_1d, hname_1d, nbins_var, 0, nbins_var)
    hist_unrolled.GetXaxis().SetTitle('Unrolled #rho_{1} v.s. #rho_{2}')
    hist_unrolled.GetYaxis().SetTitle('a.u.')


    idx_var = 1
    for iy in range(1, hist_coarse.GetYaxis().GetNbins()+1):
        for ix in range(1, hist_coarse.GetXaxis().GetNbins()+1):
            
            hist_unrolled.SetBinContent(idx_var, hist_coarse.GetBinContent(ix, iy))
            hist_unrolled.SetBinError(idx_var, hist_coarse.GetBinError(ix, iy))
                       
            idx_var += 1

    hist_unrolled.GetYaxis().SetRangeUser(0, hist_unrolled.GetMaximum()*1.3)

    hists.append(copy.deepcopy(hist_unrolled))


    canvas2 = TCanvas(name + '_1d')
    hist_unrolled.Draw('hep')
    
    l22=add_CMS()
    l22.Draw("same")
    l32=add_Preliminary()
    l32.Draw("same")

    canvas2.SaveAs('homework/' + name+'_1d.pdf')
    canvas2.SaveAs('homework/' + name+'_1d.gif')




print hists

ofile=TFile('homework/hist.root', 'recreate')

for hist in hists:
    hist.Write()

ofile.Write()
ofile.Close()


