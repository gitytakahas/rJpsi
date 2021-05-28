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
parser.add_option("-s", "--sys", default="None", type="string", dest="sys")
parser.add_option('-m', '--min', action="store_true", default=False, dest='min')

(options, args) = parser.parse_args() 

colours = [1, 2, 4, 6, 8, 1, 46, 13, 15, 1,1,1,1,1,1]
styles = [1, 2, 4, 3, 5, 2,7, 1, 1, 1,1,1,1,1,1]

#isxgbs = True

#if not isxgbs:
#    for vkey, ivar in vardir.items():
#        if vkey.find('xgbs')!=-1:
#            vardir.pop(vkey)




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
    xgbscut = 'xgbs > 6. && xgbs < 8.5'


#cut = '&&'.join([basic, rhomass, dnncut, phiveto, xgbscut])
cut = '&&'.join([basic, rhomass, phiveto, xgbscut])
#cut = '&&'.join([basic, phiveto])


print '-'*80
print 'channel = ', options.channel
print '-'*80
print 'cut = ', cut
print '-'*80

#sys.exit(1)

#prefix ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_mass_pt/'
prefix ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_hammer/'

datastr =  prefix.replace('_hammer','') + "Data_2018/data.root"
sigstr = prefix + "BcJpsiTau_inclusive_ul_all_2018/signal_xgbs6.root"
#bkgstr = "BcJpsiX_ul_2018/bkg_xgbs6.root"

def_sig_3p = cut + '&& n_occurance==1 && tau_isRight_3prong==1 && tau_isRight_3prong_pi0==0'
def_sig_3pp = cut + '&& n_occurance==1 && tau_isRight_3prong==1 && tau_isRight_3prong_pi0==1'
def_sig_others = cut + '&& n_occurance==1 && tau_isRight_3prong==0'
def_bg_bc = cut + '&& n_occurance==0'

bc_sf = 0.45/(3*0.8)
wstr = 'puweight'

ddir['data_obs'] = {'file':datastr, 'cut':cut, 'weight':'1', 'scale':1, 'order':2999}

ddir['sig_3p'] = {'file':sigstr, 'cut':def_sig_3p, 'weight':wstr, 'scale':bc_sf, 'order':1}
ddir['sig_3pp'] = {'file':sigstr, 'cut':def_sig_3pp, 'weight':wstr,  'scale':bc_sf, 'order':2}
ddir['sig_others'] = {'file':sigstr, 'cut':def_sig_others, 'weight':wstr, 'scale':bc_sf, 'order':3}
ddir['bg_bc'] = {'file':sigstr, 'cut':def_bg_bc, 'weight':wstr, 'scale':bc_sf, 'order':4}

wstr_shape = 'hammer_ebe*weight_ctau*puweight'
ddir['sig_3p_shape'] = {'file':sigstr, 'cut':def_sig_3p, 'weight':wstr_shape, 'scale':bc_sf, 'order':1}
ddir['sig_3pp_shape'] = {'file':sigstr, 'cut':def_sig_3pp, 'weight':wstr_shape,  'scale':bc_sf, 'order':2}
ddir['sig_others_shape'] = {'file':sigstr, 'cut':def_sig_others, 'weight':wstr_shape, 'scale':bc_sf, 'order':3}





# list of variables to produce datacards
#finaldiscriminant = ['q2', 'q2_simple', 'b_mass', 'b_mcorr', 'xgbs_zoom', 'estar', 'estar_simple', 'mm2', 'mm2_simple', 'tau_mass_zoom', 'xgbs_fine']
finaldiscriminant = ['q2_simple']
#finaldiscriminant = ['xgbs_fine']

if options.min:
    for vkey, ivar in vardir.items():
        if vkey not in finaldiscriminant:
            vardir.pop(vkey)


print '-'*80

for vkey, ivar in vardir.items():
    print '->', vkey.ljust(20), ivar

print '-'*80


            
for type, ivar in ddir.items():
    
    wstr = ivar['weight']

    file = TFile(ivar['file'])
    tree = file.Get('tree')

    print ivar['file'], 'is read ...'

    var_tuples, hists = returnTuples(type, vardir)    

    if options.sys.find('hammer')!=-1 and type.find('_shape')!=-1:
        wstr = ivar['weight'].replace('hammer_ebe', options.sys)

    elif options.sys.find('ctau')!=-1 and type.find('_shape')!=-1:
        wstr = ivar['weight'].replace('weight_ctau', options.sys)


    cut = '(' + ivar['cut'] + ')*' + wstr
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

    Histo.Group('signal_ref', ['sig_3p', 'sig_others', 'sig_3pp'])
    Histo.Group('signal_shape', ['sig_3p_shape', 'sig_3pp_shape', 'sig_others_shape'])

#    print 'writing to datacard', vkey
    if vkey in finaldiscriminant:
        if options.sys=="None":
            Histo.WriteDataCard('datacard/' + options.channel + '/' + vkey + '.root', True, 'RECREATE', options.channel)
        else:
            Histo.WriteDataCard('datacard/' + options.channel + '/' + vkey + '_' + options.sys + '.root', True, 'RECREATE', options.channel)
