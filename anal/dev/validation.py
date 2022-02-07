import copy, os,  math, sys
from numpy import array
from officialStyle import officialStyle
from array import array
import common.MultiDraw
import numpy as np
from varConfig import vardir
from DataMCPlot import *
from DisplayManager_postfit import DisplayManager
from common.DisplayManager_compare import DisplayManager_compare
from helper import *
from H2TauStyle import *


lumi=59.6

#gROOT.SetBatch(False)
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
#gStyle.SetOptStat(0)


gROOT.ProcessLine(".L ~/tool//MultiDraw.cc+");


file_hammer = TFile('/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi//job_pt_LEGACY/BcJpsiTau_inclusive_ul_all_2018/Myroot_0.root')
hist_hammer = file_hammer.Get('hammer')


from optparse import OptionParser, OptionValueError
usage = "usage: python compare.py" 
parser = OptionParser(usage) 

parser.add_option("-s", "--sys", default="None", type="string", dest="sys")
parser.add_option('-m', '--min', action="store_true", default=False, dest='min')
parser.add_option('-c', '--create', action="store_true", default=False, dest='create')

(options, args) = parser.parse_args() 

print(options)

if not options.create:
    gROOT.Macro('./functionmacro.C+')

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



def comparisonPlots(hist, pname='sync.pdf', isLog=False, isRatio=True, clabel=''):

    display = DisplayManager(pname, isRatio, isLog, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)


def comparisonPlots_alt(hists, titles, norm=False, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he', prefix=None):

    display = DisplayManager_compare(pname, isLog, isRatio, 0.6, 0.7, doption)
    display.draw_legend = isLegend
    
    display.Draw(hists, titles, norm, prefix)

def drawSignificance(hists, channel, vkey, isRight=True):
    
    data_obs = None
    
    for hist in hists:
        if hist.GetName().find('data')!=-1:
            data_obs = copy.deepcopy(hist)

        if hist.GetName().find('sig_3p')!=-1:
            sig = copy.deepcopy(hist)


    maxBin = data_obs.GetNbinsX()

    plot = TGraph()

    sigs = []
    xvals = []

    for nBin in range(0, maxBin+1):
        
        nsig = sig.Integral(nBin, maxBin+1)
        ndata = data_obs.Integral(nBin, maxBin+1)
        xval = sig.GetBinCenter(nBin)

#        print(nBin, xval, 'sig', nsig, 'data', ndata)
        
        if ndata!=0:
            significance = nsig/math.sqrt(ndata)
        else:
            significance = 0.
        
        sigs.append(significance)
        xvals.append(xval)

        plot.SetPoint(nBin, xval, significance)

#    significance = [sig.Integral(nBin, maxBin+1)/math.sqrt(data_obs.Integral(nBin, maxBin+1)) for nBin in range(0, maxBin + 1) if data_obs.Integral(nBin, maxBin+1)!=0]
#    xvals = [sig.GetBinCenter(nBin) for nBin in range(0, maxBin + 1)]

#    plot = TGraph(maxBin, np.asarray(xvals), np.asarray(significance))
    plot.GetXaxis().SetTitle(vkey)
    plot.GetYaxis().SetTitle('Significance')

#    print('significance', len(sigs), sigs)
#    print('xvals', len(xvals), xvals)

    maxsig = max(sigs)
    index = sigs.index(maxsig)
    maxx = xvals[index]

    print('===> max significance', maxsig, 'is reached at the xvalue of ', maxx)

    canvas = TCanvas('can_' + channel + '_' + vkey)
    
    plot.Draw('apl')
    
    canvas.SaveAs('plots/' + channel + '/sig_' + vkey  + '.gif')
    
        


multihists = {}

##################################################

prefix ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_LEGACY/'
prefix_cr ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_LEGACY_cr/'
#prefix_q3 ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_q3/'

datastr = "Data_2018/data.root"
sigstr  = "BcJpsiTau_inclusive_ul_all_2018/sig.root"
#datastr = "Data_2018/Myroot_training_weightAdded.root"
#sigstr  = "BcJpsiTau_inclusive_ul_all_2018/Myroot_training_weightAdded.root"
#bkgstr  = "BcJpsiX_ul_2018/bkg.root"
bkgstr  = "BcJpsiX_ul_2018_new/bkg.root"

bc_sf = 0.45/(3*0.8)
#sig_sf = 'hammer_ebe*puweight/0.55'
#sig_sf = 'hammer_ebe*puweight/0.55'

#sig_sf = 'weight'

ddir = {}


binid = 2

if options.sys.find('hammer')!=-1:
    print('detecting new hammer weight')
    
    binid = -1 

    for ibin in range(1, hist_hammer.GetXaxis().GetNbins()+1):
        if options.sys.replace('hammer_ebe','num') == hist_hammer.GetXaxis().GetBinLabel(ibin):
            binid = ibin


    if binid == -1:
        print('No Hammer weight is found !!! Set to 2')
        binid = 2


hammer_sf = float(hist_hammer.GetBinContent(binid))/float(hist_hammer.GetBinContent(1))

hammer_weight = 'hammer_ebe/' + str(hammer_sf)

print('hammer weight = ', hammer_weight, 'id=', binid)

# data 
ddir['data_obs'] =  {'file':datastr, 'weight':'1', 'scale':1, 'order':2999, 'color':1, 'addcut':'1'}

# signal + Bc 
ddir['sig_3p'] =     {'file':sigstr, 'weight':'puweight*' + hammer_weight, 'scale':bc_sf, 'order':1, 'color':ttcol_v2, 'addcut':'n_occurance==1 && tau_isRight_3prong==1'}

#ddir['sig_lep'] = {'file':sigstr, 'weight':'puweight*(hammer_ebe/0.55)', 'scale':bc_sf, 'order':3, 'color':ttcol, 'addcut':'n_occurance==1 && tau_isRight_3prong==0 && isJpsiTau2Mu==1'}
#ddir['sig_others'] = {'file':sigstr, 'weight':'puweight*(hammer_ebe/0.55)', 'scale':bc_sf, 'order':3, 'color':wcol, 'addcut':'n_occurance==1 && tau_isRight_3prong==0 && isJpsiTau2Mu==0'}
ddir['sig_others'] = {'file':sigstr, 'weight':'puweight*' + hammer_weight, 'scale':bc_sf, 'order':3, 'color':wcol, 'addcut':'n_occurance==1 && tau_isRight_3prong==0'}


#ddir['bg_bc'] =      {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':qcdcol, 'addcut':'n_occurance==0 && isJpsiMu==0'}
#ddir['bg_norm'] =      {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':jtfake, 'addcut':'n_occurance==0 && isJpsiMu==1'}
ddir['bg_bc'] =      {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':qcdcol, 'addcut':'n_occurance==0'}


# other B backgrounds 
#ddir['bg_ul'] =      {'file':bkgstr, 'weight':'puweight', 'scale':7*0.64/0.8, 'order':5, 'color':dycol, 'addcut':'1'}
ddir['bg_ul'] =      {'file':bkgstr, 'weight':'puweight', 'scale':21*7*0.64/0.8, 'order':5, 'color':dycol, 'addcut':'1'}


##################################################


#basic = 'tau_pt > 3. && mu1_isLoose==1 && mu2_isLoose==1'
basic = 'tau_pt > 3. && mu1_isTight==1 && mu2_isTight==1'
phiveto = "!(tau_rhomass1_kk < 1.04) && !(tau_rhomass2_kk < 1.04)"
#xgbs_sr = 'xgbs > 5.35'
#xgbs_sb = 'xgbs > 4.35 && xgbs < 5.35'
#xgbs_lp = 'xgbs > 3.35 && xgbs < 4.35'
xgbs_sr = 'xgbs > 5.1'
xgbs_sb = 'xgbs > 4.1 && xgbs < 5.1'
xgbs_lp = 'xgbs > 3.1 && xgbs < 4.1'

#rhomass = '((tau_rhomass1 > 0.65 && tau_rhomass1 < 0.89) || (tau_rhomass2 > 0.65 && tau_rhomass2 < 0.89)) '
#rhomass = '((tau_rhomass1 > 0.69 && tau_rhomass1 < 0.85) || (tau_rhomass2 > 0.69 && tau_rhomass2 < 0.85)) '

#taumass = '(tau_mass > 1. && tau_mass < 1.3)'
#taumass = '(tau_mass > 0.93 && tau_mass < 1.4)'

channels = {

    'validation':{'cut':'&&'.join([basic])},

#    'extrapolate':{'cut':'&&'.join([basic, 'xgbs > 3.5'])},

#    'sr':{'cut':'&&'.join([basic, xgbs_sr])},

#    'sb':{'cut':'&&'.join([basic, xgbs_sb])},

#    'lp':{'cut':'&&'.join([basic, xgbs_lp])},


#    'cr_sr':{'cut':'&&'.join([basic, xgbs_sr])},
#    'cr_sb':{'cut':'&&'.join([basic, xgbs_sb])},
#    'cr_lp':{'cut':'&&'.join([basic, xgbs_lp])},

#    'q3_sr':{'cut':'&&'.join([basic, xgbs_sr])},
#    'q3_sb':{'cut':'&&'.join([basic, xgbs_sb])},
#    'q3_lp':{'cut':'&&'.join([basic, 'xgbs > 3.'])},

#    'hp_sr':{'cut':'&&'.join([basic, taumass, xgbs_sr])},
#    'hp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_sr])},
#
#    # lp xgbs sideband
#    'lp_sr':{'cut':'&&'.join([basic, taumass, xgbs_sb])},
#    'lp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_sb])},
##
##    # mp xgbs sideband
#    'slp_sr':{'cut':'&&'.join([basic, taumass, xgbs_lp])},
#    'slp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_lp])},
##    
#    'cr_hp_sr':{'cut':'&&'.join([basic, taumass, xgbs_sr])},
#    'cr_hp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_sr])},
##
##    # lp xgbs sideband
#    'cr_lp_sr':{'cut':'&&'.join([basic, taumass, xgbs_sb])},
#    'cr_lp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_sb])},
##
##    # mp xgbs sideband
#    'cr_slp_sr':{'cut':'&&'.join([basic, taumass, xgbs_lp])},
#    'cr_slp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_lp])},

}


#finaldiscriminant = ['xgbs', 'tau_rhomass_unrolled', 'tau_rhomass_unrolled_coarse']
#finaldiscriminant = ['xgbs', 'xgbs_zoom', 'b_mass', 'b_mass_sf', 'tau_rhomass_unrolled', 'tau_rhomass_unrolled_coarse', 'q2_simple', 'jpsi_kpipi']
finaldiscriminant = ['jpsi_kpipi']


if not options.create:
    vardir.pop('b_mass_sf')

if options.min:
    for vkey, ivar in vardir.items():
        if vkey not in finaldiscriminant: 
            vardir.pop(vkey)

print '-'*80

for vkey, ivar in vardir.items():
    print '->', vkey.ljust(20), ivar

print '-'*80



for channel, dict in channels.iteritems():
    print '-'*80
    print 'channel = ', channel
    print '-'*80
    print 'cut = ', dict['cut']
    print '-'*80

    ensureDir('plots/' + channel)

    for type, ivar in ddir.items():

        wstr = ivar['weight']
    
        filename = None

        if channel.find('cr')!=-1:
            filename = prefix_cr + '/' + ivar['file']
        else:        
            filename = prefix + '/' + ivar['file']

        file = TFile(filename)
        tree = file.Get('tree')

        print(filename, 'is read ...')

        var_tuples, hists = returnTuples(type, vardir)    

#        if not options.create and type.find('bg_ul')!=-1:
#        if type.find('bg_ul')!=-1:
#            wstr += '*getWeight(b_mass)'


        if options.sys.find('hammer')!=-1 and type.find('sig')!=-1: 
            wstr = ivar['weight'].replace('hammer_ebe', options.sys)

        elif options.sys.find('ctau')!=-1 and type in ['sig_3p', 'sig_others', 'bg_bc']:
            wstr = ivar['weight'].replace('weight_ctau', options.sys)

        elif options.sys.find('puweight')!=-1 and type.find('data')==-1:
            wstr = ivar['weight'].replace('puweight', options.sys)


        cut = '(' + dict['cut'] + ' &&' + ivar['addcut'] + ')*' + wstr
        print(type, cut)


        

        tree.MultiDraw(var_tuples, cut)
    
        multihists.update(copy.deepcopy(hists))


#    print('write to datacards ...')



    for vkey, ivar in vardir.items():

        gStyle.SetPalette(kBird)
        Histo = DataMCPlot(vkey)

        print '-'*80
        print vkey 

        hists = []
        titles = []


        for type, ivar2 in ddir.items():

            multihists[type + '_' + vkey].Scale(Double(ivar2['scale']))

            ################################
            ### This will be used later ... 
            hist = copy.deepcopy(multihists[type + '_' + vkey])
            hist.SetMarkerColor(ivar2['color'])
            hist.SetLineColor(ivar2['color'])
            hist.SetFillStyle(0)
            hists.append(hist)
            titles.append(type)
            ################################



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
        print(Histo)
        
        commonname = vkey

        if options.sys!="None":
            commonname = vkey + '_' + options.sys


        comparisonPlots(Histo, 'plots/' + channel + '/' + commonname  + '.gif')

#        if channel.find('q3')==-1:

        comparisonPlots(Histo, 'plots/' + channel + '/' + commonname  + '_log.gif', True)
        
#        Histo.Group('signal_ref', ['sig_3p', 'sig_others', 'sig_3pp'])
            
        ensureDir('datacard/' + channel)

        if vkey in finaldiscriminant:

            print(options.sys)

#            if options.sys=="None":
            Histo.WriteDataCard('datacard/' + channel + '/' + commonname + '.root', True, 'RECREATE', channel)

        if options.sys=='None' and vkey=='xgbs':
            drawSignificance(hists, channel, vkey, True)
        
        comparisonPlots_alt(hists, titles, True, False, 'plots/' + channel + '/' + commonname + '_compare.pdf', False, True, 'pE')


    
