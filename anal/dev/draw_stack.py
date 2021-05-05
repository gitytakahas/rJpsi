import copy, os, math, sys
from numpy import array
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, Double, TCanvas, kRed, TH1D, TRandom3, kBird, THStack
from officialStyle import officialStyle
from array import array
import common.MultiDraw
import numpy as np
from varConfig import vardir
from DataMCPlot import *
from DisplayManager_postfit import DisplayManager

lumi=59.6

gRandom = TRandom3()

#gROOT.SetBatch(False)
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
#gStyle.SetOptStat(0)


gROOT.ProcessLine(".L ~/tool//MultiDraw.cc+");

from optparse import OptionParser, OptionValueError
usage = "usage: python compare.py" 
parser = OptionParser(usage) 

parser.add_option("-c", "--channel", default="sr", type="string", dest="channel")
#parser.add_option('-f', '--ff', action="store_true", default=False, dest='ff')

(options, args) = parser.parse_args() 

colours = [1, 2, 4, 6, 8, 1, 46, 13, 15, 1,1,1,1,1,1]
styles = [1, 2, 4, 3, 5, 2,7, 1, 1, 1,1,1,1,1,1]

isxgbs = True

if not isxgbs:
    for vkey, ivar in vardir.items():
        if vkey.find('xgbs')!=-1:
            vardir.pop(vkey)


print '-'*80

for vkey, ivar in vardir.items():
    print '->', vkey.ljust(20), ivar

print '-'*80


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


def returnTuples(prefix, vardir):

    var_tuples = []
    multihists = {}

    for varname, var in vardir.items():        

        hname = prefix + '_' + varname

        hist_register = TH1D(hname, hname, var['nbin'], var['xmin'], var['xmax'])

        hist_register.GetXaxis().SetTitle(var['xtitle'])
        hist_register.GetYaxis().SetTitle('a.u.')
        hist_register.Sumw2()
        hist_register.SetTitle(varname)
        
        multihists[hname] = hist_register

        var2draw = varname
        if 'var' in var:
            var2draw = var['var']


        var_tuples.append('{var} >> {hist}'.format(var=var2draw, hist=hname))
            
    return var_tuples, multihists



def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def comparisonPlots(hist, pname='sync.pdf', clabel='', isRatio=True):

    display = DisplayManager(pname, isRatio, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)


ensureDir('plots/' + options.channel)
ensureDir('datacard/' + options.channel)

multihists = {}

ddir = {}



phiveto = "!(tau_rhomass1_kk < 1.04) && !(tau_rhomass2_kk < 1.04)"
basic = 'tau_pt > 3. && mu1_isLoose==1 && mu2_isLoose==1 && (pi1_trigMatch==1 || pi2_trigMatch==1 || pi3_trigMatch==1)'
#dnncut = 'perEVT_mc > 0.7'
xgbscut = 'xgbs > 8.5'
#xgbscut = 'xgbs > 8.6'
#rhomass = '((tau_rhomass1 > 0.65 && tau_rhomass1 < 0.89) || (tau_rhomass2 > 0.65 && tau_rhomass2 < 0.89)) '
rhomass = '((tau_rhomass1 > 0.69 && tau_rhomass1 < 0.85) || (tau_rhomass2 > 0.69 && tau_rhomass2 < 0.85)) '

if options.channel in ['sb', 'cr2']:
    rhomass = '!' + rhomass

if options.channel in ['cr1', 'cr2']:
#    xgbscut = 'xgbs > 7.25 && xgbs < 8.5'
    xgbscut = 'xgbs > 7. && xgbs < 8.5'


#cut = '&&'.join([basic, rhomass, dnncut, phiveto, xgbscut])
cut = '&&'.join([basic, rhomass, phiveto, xgbscut])
#cut = '&&'.join([basic, phiveto])


print '-'*80
print 'channel = ', options.channel
print '-'*80
print 'cut = ', cut
print '-'*80

#sys.exit(1)

prefix ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_mass_pt/'

datastr="Data_2018/data.root"
sigstr="BcJpsiTau_inclusive_ul_all_2018/signal.root"
bkgstr="BcJpsiX_ul_2018/bkg.root"


ddir['data_obs'] = {'file':prefix + datastr, 'acut':cut, 'weight':'1', 'scale':1, 'order':2999}
#ddir['bg_ul'] = {'file':prefix + bkgstr, 'acut':cut, 'weight':'puweight', 'scale':7/0.8, 'order':4}
ddir['sig_3p'] = {'file':prefix + sigstr, 'acut':cut + '&& n_occurance==1 && tau_isRight_3prong==1 && tau_isRight_3prong_pi0==0', 'weight':'puweight', 'scale':0.45/(3*0.8), 'order':1}
ddir['sig_3pp'] = {'file':prefix + sigstr, 'acut':cut + '&& n_occurance==1 && tau_isRight_3prong==1 && tau_isRight_3prong_pi0==1', 'weight':'puweight',  'scale':ddir['sig_3p']['scale'], 'order':2}
ddir['sig_others'] = {'file':prefix + sigstr, 'acut':cut + '&& n_occurance==1 && tau_isRight_3prong==0', 'weight':'puweight', 'scale':ddir['sig_3p']['scale'], 'order':3}
ddir['bg_bc'] = {'file':prefix + sigstr, 'acut':cut + '&& n_occurance==0', 'weight':'puweight', 'scale':ddir['sig_3p']['scale'], 'order':4}


# list of variables to produce datacards
finaldiscriminant = ['q2', 'q2_simple', 'b_mass', 'b_mcorr', 'xgbs_zoom', 'estar', 'estar_simple', 'mm2', 'mm2_simple', 'tau_mass_zoom', 'xgbs_fine']
#finaldiscriminant = ['xgbs_fine']

            
for type, ivar in ddir.items():
    
    file = TFile(ivar['file'])
    tree = file.Get('tree')

    print ivar['file'], 'is read ...'

    var_tuples, hists = returnTuples(type, vardir)


    cut = '(' + ivar['acut'] + ')*' + ivar['weight']
    print(type, cut)

    tree.MultiDraw(var_tuples, cut)
    
    multihists.update(copy.deepcopy(hists))


print('write to datacards ...')

writes = {}

for vkey, ivar in vardir.items():

    gStyle.SetPalette(kBird)
    Histo = DataMCPlot(vkey)

    print '-'*80
    print vkey 

    for type, ivar2 in ddir.items():

        multihists[type + '_' + vkey].Scale(Double(ivar2['scale']))

        hist2add = multihists[type + '_' + vkey]

        if type.find('data')!=-1:
            hist2add.SetFillStyle(0)
            hist2add.Sumw2(False)
            hist2add.SetBinErrorOption(1)

        hist2add.SetName(type)
        hist2add.SetTitle(type)
        hist2add.GetXaxis().SetLabelColor(1)
        hist2add.GetXaxis().SetLabelSize(0.0)

        Histo.AddHistogram(hist2add.GetName(), hist2add, ivar2['order'])

        if type.find('data')!=-1:
            Histo.Hist(hist2add.GetName()).stack = False
    

    Histo._ApplyPrefs()
    print Histo

    comparisonPlots(Histo, 'plots/' + options.channel + '/' + vkey  + '.gif')

    Histo.Group('signal', ['sig_3p', 'sig_others', 'sig_3pp'])

#    print 'writing to datacard', vkey
    if vkey in finaldiscriminant:
        Histo.WriteDataCard('datacard/' + options.channel + '/' + vkey + '.root', True, 'RECREATE', options.channel)
