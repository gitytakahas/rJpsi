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

prefix ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_Legacy/'
#prefix_cr ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_vprobfsigcr/'
#prefix_q3 ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_q3/'

datastr = "Data_2018/data.root"
sigstr  = "BcJpsiTau_inclusive_ul_all_2018/sig.root"
#datastr = "Data_2018/Myroot_training_weightAdded.root"
#sigstr  = "BcJpsiTau_inclusive_ul_all_2018/Myroot_training_weightAdded.root"
bkgstr  = "BcJpsiX_ul_2018/bkg.root"

bc_sf = 0.45/(3*0.8)
#sig_sf = 'hammer_ebe*puweight/0.55'
#sig_sf = 'hammer_ebe*puweight/0.55'

#sig_sf = 'weight'

ddir = {}

# data 
ddir['data_obs'] =  {'file':datastr, 'weight':'1', 'scale':1, 'order':2999, 'color':1, 'addcut':'1'}

# signal + Bc 
ddir['sig_3p'] =     {'file':sigstr, 'weight':'hammer_ebe*puweight/0.55', 'scale':bc_sf, 'order':1, 'color':ttcol_v2, 'addcut':'n_occurance==1 && tau_isRight_3prong==1'}
ddir['sig_others'] = {'file':sigstr, 'weight':'hammer_ebe*puweight/0.55', 'scale':bc_sf, 'order':3, 'color':wcol, 'addcut':'n_occurance==1 && tau_isRight_3prong==0'}
ddir['bg_bc'] =      {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':qcdcol, 'addcut':'n_occurance==0'}

# other B backgrounds 
#ddir['bg_ul'] =      {'file':bkgstr, 'weight':'puweight', 'scale':7/0.8, 'order':5, 'color':dycol, 'addcut':'1'}
ddir['bg_ul'] =      {'file':bkgstr, 'weight':'puweight', 'scale':7*0.64/0.8, 'order':5, 'color':dycol, 'addcut':'1'}


##################################################


basic = 'tau_pt > 3. && mu1_isLoose==1 && mu2_isLoose==1'
phiveto = "!(tau_rhomass1_kk < 1.04) && !(tau_rhomass2_kk < 1.04)"
xgbs_sr = 'xgbs > 4.84'
xgbs_sb = 'xgbs > 3.5 && xgbs < 4.84'
xgbs_lp = 'xgbs > 3. && xgbs < 3.5'

#rhomass = '((tau_rhomass1 > 0.65 && tau_rhomass1 < 0.89) || (tau_rhomass2 > 0.65 && tau_rhomass2 < 0.89)) '
#rhomass = '((tau_rhomass1 > 0.69 && tau_rhomass1 < 0.85) || (tau_rhomass2 > 0.69 && tau_rhomass2 < 0.85)) '

#taumass = '(tau_mass > 1. && tau_mass < 1.3)'
#taumass = '(tau_mass > 0.93 && tau_mass < 1.4)'


channels = {

#    'inclusive':{'cut':'&&'.join([basic])},

    'extrapolate':{'cut':'&&'.join([basic, 'xgbs > 3.'])},
#
#    'sb':{'cut':'&&'.join([basic, xgbs_lp])},
#
#    'lp':{'cut':'&&'.join([basic, xgbs_slp])},

#    'cr_sr':{'cut':'&&'.join([basic, xgbs_hp])},
#    'cr_sb':{'cut':'&&'.join([basic, xgbs_lp])},
#    'cr_lp':{'cut':'&&'.join([basic, xgbs_slp])},

#    'hp_sr':{'cut':'&&'.join([basic, taumass, xgbs_hp])},
#    'hp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_hp])},
#
#    # lp xgbs sideband
#    'lp_sr':{'cut':'&&'.join([basic, taumass, xgbs_lp])},
#    'lp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_lp])},
##
##    # mp xgbs sideband
#    'slp_sr':{'cut':'&&'.join([basic, taumass, xgbs_slp])},
#    'slp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_slp])},
##    
#    'cr_hp_sr':{'cut':'&&'.join([basic, taumass, xgbs_hp])},
#    'cr_hp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_hp])},
##
##    # lp xgbs sideband
#    'cr_lp_sr':{'cut':'&&'.join([basic, taumass, xgbs_lp])},
#    'cr_lp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_lp])},
##
##    # mp xgbs sideband
#    'cr_slp_sr':{'cut':'&&'.join([basic, taumass, xgbs_slp])},
#    'cr_slp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_slp])},

}


finaldiscriminant = ['xgbs_zoom']

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

        cut = '(' + dict['cut'] + ' &&' + ivar['addcut'] + ')*' + wstr
        print(type, cut)


        

        tree.MultiDraw(var_tuples, cut)
    
        multihists.update(copy.deepcopy(hists))



    for vkey, ivar in vardir.items():

        print '-'*80
        print vkey 

#        hists = []
#        titles = []


#        for type, ivar2 in ddir.items():

            

            ################################
            ### This will be used later ... 
#            hist = copy.deepcopy(multihists[type + '_' + vkey])
#            hist.SetMarkerColor(ivar2['color'])
#            hist.SetLineColor(ivar2['color'])
#            hist.SetFillStyle(0)
#            hists.append(hist)
#            titles.append(type)
            ################################

        

        
        canvas = TCanvas('can_' + channel + '_' + vkey)
        canvas.SetLogy()

        data = copy.deepcopy(multihists['data_obs_' + vkey])
        mc = copy.deepcopy(multihists['bg_ul_' + vkey])


#        for proc in ['bg_bc', 'sig_others','sig_3p']:
##        for proc in ['bg_bc']:
#            _hist = copy.deepcopy(multihists[proc + '_' + vkey])
#            print(proc, _hist.GetSumOfWeights(), 'subtracted')
#            data.Add(_hist, -1)


#        print('data scaling by', data.GetSumOfWeights())
#        data.Scale(1./data.GetSumOfWeights())

        print('data yield=', data.GetSumOfWeights())

        for ibin in range(1, data.GetXaxis().GetNbins()+1):
            if data.GetBinCenter(ibin) > 3.5:
                data.SetBinContent(ibin, -9)


#        data.GetYaxis().SetRangeUser(0, 0.3)
        data.GetYaxis().SetRangeUser(0.1, data.GetBinContent(1)*10)
#        data.Draw('ep')
#        data.Fit('expo', '','',3, 4.5)

#        fitorig = copy.deepcopy(data.GetFunction('expo'))
#        fitorig.SetLineColor(kBlue)
        
#        fitfunc = copy.deepcopy(data.GetFunction('expo'))
#        fitfunc.SetRange(3, 7)
#        fitfunc.SetLineStyle(2)
#        fitfunc.SetLineColor(kBlue)

        ####### 
#        mc.Scale(1./mc.GetSumOfWeights())
        mcsf = data.GetBinContent(1)/mc.GetBinContent(1)
        mc.Scale(mcsf)
        print('mc scaled by', mcsf)
        mc.SetMarkerColor(dycol)
        mc.SetLineColor(dycol)
        mc.SetMarkerStyle(22)
#        mc.Draw()
#        mc.Fit('expo', '','',3, 7)

        algo = 'expo'

        mc.Fit(algo, '','',3, 8)

        fitorig_mc = copy.deepcopy(mc.GetFunction(algo))
        fitorig_mc.SetLineColor(kRed)
        fitorig_mc.Draw('same')
        

        data.Draw('ep')
        mc.Draw('epsame')
#        fitorig.Draw('same')
#        fitfunc.Draw('same')
        fitorig_mc.Draw('E3same')

        leg = TLegend(0.4, 0.7, 0.9, 0.9)
        applyLegendSettings(leg)

        leg.AddEntry(data, 'data (xgbs < 4.8)', 'lep')
        leg.AddEntry(mc, 'bkg. mc', 'lep')
        leg.AddEntry(fitorig_mc, algo, 'l')

        leg.Draw()
        
        canvas.SaveAs('plots/' + channel + '/extrapolate_' + vkey + '.gif')
        

        nsr = float(fitorig_mc.Integral(4.84, 10))
        nsb = float(fitorig_mc.Integral(3.5, 4.5))
        ratio = float(nsr/nsb)
        
        print('mc integral in SR',  nsr)
        print('mc integral in SB', nsb)
        print('ratio=', ratio)
    
