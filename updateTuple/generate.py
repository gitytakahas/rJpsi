import copy, math, os, collections, sys
from numpy import array
#from CMGTools.H2TauTau.proto.plotter.categories_TauMu import cat_Inc
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, Double, TCanvas, kRed, TH1D, TRandom3
from DisplayManager import DisplayManager
from officialStyle import officialStyle
from array import array
import MultiDraw
import numpy as np

gRandom = TRandom3()

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)


def overflow(hist):
    lastp1 = hist.GetBinContent(hist.GetXaxis().GetNbins()+1)
    last = hist.GetBinContent(hist.GetXaxis().GetNbins())
    lastp1e = hist.GetBinError(hist.GetXaxis().GetNbins()+1)
    laste = hist.GetBinError(hist.GetXaxis().GetNbins())
    hist.SetBinContent(hist.GetXaxis().GetNbins(), last+lastp1)
    hist.SetBinError(hist.GetXaxis().GetNbins(), math.sqrt(math.pow(laste,2)+math.pow(lastp1e,2)))

    firstp1 = hist.GetBinContent(1)
    first = hist.GetBinContent(0)
    firstp1e = hist.GetBinError(1)
    firste = hist.GetBinError(0)
    hist.SetBinContent(1, first+firstp1)
    hist.SetBinError(1, math.sqrt(math.pow(firste,2)+math.pow(firstp1e,2)))

def datarize(hist):
    
    for ibin in range(1, hist.GetXaxis().GetNbins()+1):
        
        # First, fluctuate based on poisson
        entry = gRandom.Poisson(hist.GetBinContent(ibin))
#        print hist.GetBinContent(ibin), entry
        
        hist.SetBinContent(ibin, entry)
        hist.SetBinError(ibin, math.sqrt(entry))


def applyHistStyle(h, i):
#    print h, i
    h.SetLineColor(colours[i])
    h.SetMarkerColor(colours[i])
    h.SetMarkerSize(0.2)
    h.SetLineStyle(styles[i])
    h.SetLineWidth(3)
    h.SetStats(False)

def normalize(hist):

    if hist.GetSumOfWeights()==0:
        print(hist.GetName(), 'does not have any entries ...')
    else:
#        print vkey, hist.GetName(), '->', hist.GetSumOfWeights()
        hist.Scale(1./hist.GetSumOfWeights())
        hist.SetMaximum(hist.GetBinContent(hist.GetMaximumBin())*1.2)

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def comparisonPlots(hists, titles, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='ep', prefix=None):

    display = DisplayManager(pname, isLog, isRatio, 0.2, 0.68, doption)
    display.draw_legend = isLegend
#    display.isEff = isEff
    
    display.Draw(hists, titles, prefix)


def sproducer(vkey, rootfile, ivar, sel):

    hist = TH1F('h_' + vkey, 
                'h_' + vkey, 
                ivar['nbin'], ivar['xmin'], ivar['xmax'])

    hist.Sumw2()
    exp = '(' + sel + ')'
        
    tree = rootfile.Get('tree')

    var2draw = vkey
    if 'var' in ivar:
        var2draw = ivar['var']

#    print ivar['var'] + ' >> ' + hist.GetName(), exp
    
    tree.Draw(var2draw + ' >> ' + hist.GetName(), exp)

#    overflow(hist)

    hist.GetXaxis().SetTitle(ivar['xtitle'])
    hist.GetYaxis().SetTitle('a.u.')


#    num = tree.GetEntries(selection + ' && ' + ivar['var'] + ' < 5.5')
        
    return copy.deepcopy(hist)


ensureDir('toys')

from optparse import OptionParser, OptionValueError
usage = "usage: python compare.py" 
parser = OptionParser(usage) 

parser.add_option("-s", "--sig", dest="sig", default="Myroot_sig_all_0.root")
parser.add_option("-b", "--bkg", dest="bkg", default="Myroot_data.root")

(options, args) = parser.parse_args()

colours = [1, 2, 4, 6, 8, 9, 46, 13, 15, 1,1,1,1,1,1]
styles = [1, 2, 4, 3, 5, 1, 1, 1,1,1,1,1,1]
ptbin = [0,1,2,3,4,5,6,7,8,9,10, 15, 20]


prefix = 'final_root/'

ddir = {}

ddir['data_sr'] = {'file':prefix + '/' + options.bkg, 'sel':'tau_mass > 1 && tau_mass < 1.4'}
ddir['data_sb'] = {'file':prefix + '/' + options.bkg, 'sel':'!(tau_mass > 1 && tau_mass < 1.4)'}
ddir['sig_sr'] = {'file':prefix + '/' + options.sig, 'sel':'tau_mass > 1 && tau_mass < 1.4'}
ddir['sig_sb'] = {'file':prefix + '/' + options.sig, 'sel':'!(tau_mass > 1 && tau_mass < 1.4)'}

vardir = {}
vardir["xgbs"] = {'tree':'tree',  'nbin':12, 'xmin':7., 'xmax':9., 'xtitle':'XGB output score', 'var':'xgbs'}



hists = []
hists_wn = {}
titles = []

for vkey, ivar in vardir.items():

    for type, dvar in ddir.items():

#        print vkey, dvar
        sfile = TFile(dvar['file'])

        hist = sproducer(vkey, sfile, ivar, dvar['sel'])
        hist.SetName(type)
        hist.SetTitle(type)

        hists.append(copy.deepcopy(hist))
        hists_wn[type] = copy.deepcopy(hist)
        titles.append(type)


    for ii, ihist in enumerate(hists):
        applyHistStyle(ihist, ii)
        normalize(ihist)

comparisonPlots(hists, titles, False, 'Plots/' + vkey + '.pdf', True, True, 'HE')




        


#for exp in np.arange(1000, 10001, 1000).tolist():
for exp in np.arange(1000, 8001, 1000).tolist():

    print('generating signal yield with', exp)

    for num in range(0, 100): # number of toys to be generated ... 
#    for num in range(0, 2): # number of toys to be generated ... 

#        print 'toy Nr', num

        ofile = TFile('toys/toy_' + str(exp) + '_Nr' + str(num) + '.root', 'recreate')
    
        ofile.mkdir('mt')
        ofile.cd('mt')

#        for hname, hist in hists_wn.items():
#            hist.Write()
    

        data_sr = copy.deepcopy(hists_wn['data_sr'])
        data_sb = copy.deepcopy(hists_wn['data_sb'])
        data_sb_save = copy.deepcopy(hists_wn['data_sb'])
        sig_sr = copy.deepcopy(hists_wn['sig_sr'])
        sig_sb = copy.deepcopy(hists_wn['sig_sb'])

#        print('signal yield sr/sb', sig_sr.GetSumOfWeights(), sig_sb.GetSumOfWeights())
        sf_sig = exp/sig_sr.GetSumOfWeights()
#        print('sf for sig=', sf_sig)
        sig_sr.Scale(sf_sig)
        sig_sb.Scale(sf_sig)
#        print('signal yield_after sr/sb', sig_sr.GetSumOfWeights(), sig_sb.GetSumOfWeights())
        
#        print('data sb, before', data_sb.GetSumOfWeights(), 'subtract', sig_sb.GetSumOfWeights())
        data_sb.Add(sig_sb, -1.)
#        print('data sr/sb, after', data_sr.GetSumOfWeights(), data_sb.GetSumOfWeights())

        # scale BG dist first ... 
        sf_bg = data_sr.GetSumOfWeights()/data_sb.GetSumOfWeights()
#        print('sf for the sideband data =',sf_bg)
        data_sb.Scale(sf_bg)
#        print('data sb, after norm.', data_sb.GetSumOfWeights())

        
        
#        sig = copy.deepcopy(hists_wn['sig_sr'])
#        print('sig_before =', sig.GetSumOfWeights(),'exp, num=', exp, num)
#        sig.Scale(exp/sig.GetSumOfWeights())
#        print('sig_after =', sig.GetSumOfWeights())
        data_obs = copy.deepcopy(data_sb)
#        print('data_obs yield_before', data_obs.GetSumOfWeights())
        data_obs.Add(sig_sr)
#        print('data_obs yield_after', data_obs.GetSumOfWeights())
        data_obs.SetName('data_obs')
        data_obs.SetTitle('data_obs')
        
#        test = copy.deepcopy(data_obs)
#        test.SetName('test')
#        test.SetTitle('test')
#        test.Write()

        datarize(data_obs)
        data_obs.Write()

        sig_sr.Write()
        data_sb.Write()
        data_sb_save.Scale(sf_bg)
        data_sb_save.SetName('data_sb_orig')
        data_sb_save.SetTitle('data_sb_orig')
        data_sb_save.Write()

        
###        ratio = copy.deepcopy(data_sb)
###        ratio.Divide(data_obs)
###
###        ratio2 = copy.deepcopy(data_sb_save)
###        ratio2.Divide(data_obs)
###
###        ratio.SetName('ratio')
###        ratio2.SetName('ratio_orig')
###        ratio.SetTitle('ratio')
###        ratio2.SetTitle('ratio_orig')
###
###        ratio.Write()
###        ratio2.Write()
#        hists_wn['sig_sr'].Write()
#        hists_wn['data_sb'].Write()


        
        # set uncertainty as if it is data ...

        ofile.Write()
        ofile.Close()
    
