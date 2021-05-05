import copy, os, math, sys
from numpy import array
#from CMGTools.H2TauTau.proto.plotter.categories_TauMu import cat_Inc
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, Double, TCanvas, kRed, TH1D, TRandom3
from common.DisplayManager import DisplayManager
from common.officialStyle import officialStyle
from array import array
import common.MultiDraw
import numpy as np
from varConfig import vardir



gRandom = TRandom3()

gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)


gROOT.ProcessLine(".L ~/tool//MultiDraw.cc+");

from optparse import OptionParser, OptionValueError
usage = "usage: python compare.py" 
parser = OptionParser(usage) 

parser.add_option("-c", "--channel", default="BcJpsiTauNu", type="string", dest="channel")
parser.add_option('-f', '--ff', action="store_true", default=False, dest='ff')


(options, args) = parser.parse_args() 

colours = [1, 2, 4, 6, 8, 1, 46, 13, 15, 1,1,1,1,1,1]
styles = [1, 2, 4, 3, 5, 2,7, 1, 1, 1,1,1,1,1,1]
#varbin = [6, 6.15, 6.3, 6.45, 6.6, 6.75, 6.9, 7.15, 7.3, 7.45, 7.6, 7.75, 7.9, 8.05, 8.2, 8.35, 8.5, 8.65, 8.8] # 10% bias ...
varbin = [7, 7.15, 7.3, 7.45, 7.6, 7.75, 7.9, 8.05, 8.2, 8.35, 8.5, 8.8]



def calc(nsig, ndata):


    sig = 0
    sig_err = 0
    frac = 0
    frac_err = 0


    if ndata!=0:

        sig = nsig/math.sqrt(ndata)
        sig_err = math.sqrt(nsig/ndata + nsig*nsig/(4*ndata*ndata))       

        frac = nsig/ndata
        frac_err = math.sqrt(nsig/(ndata*ndata) + nsig*nsig/(ndata*ndata*ndata))

#        h_significance.SetBinContent(ibin, nsig/math.sqrt(ndata))
#        h_significance.SetBinError(ibin, err_sig)
        
#        h_sig.SetBinContent(ibin, nsig/ndata)
#        h_sig.SetBinError(ibin, err)
        
#    else:
#        h_significance.SetBinContent(ibin, 0)
#        h_significance.SetBinError(ibin, 0)
#        h_sig.SetBinContent(ibin, 0)
#        h_sig.SetBinError(ibin, 0)
        

    return sig, sig_err, frac, frac_err




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

def datarize(hist):
    
    for ibin in range(1, hist.GetXaxis().GetNbins()+1):
        
        # First, fluctuate based on poisson
        entry = gRandom.Poisson(hist.GetBinContent(ibin))
#        print hist.GetBinContent(ibin), entry
        
        hist.SetBinContent(ibin, entry)
        hist.SetBinError(ibin, math.sqrt(entry))


def returnTuples(prefix, vardir):

    var_tuples = []
    multihists = {}

    for varname, var in vardir.items():        

        hname = prefix + '_' + varname

        if 'isVar' in var:
            hist_register = TH1D(hname, hname, len(varbin)-1, array('d', varbin))
        else:
            hist_register = TH1D(hname, hname, var['nbin'], var['xmin'], var['xmax'])

        hist_register.GetXaxis().SetTitle(var['xtitle'])
        hist_register.GetYaxis().SetTitle('a.u.')
#        hist_register.GetXaxis().SetLabelSize(0)
        hist_register.Sumw2()
        hist_register.SetTitle(varname)
        

        multihists[hname] = hist_register

        var2draw = varname
        if 'var' in var:
            var2draw = var['var']


        var_tuples.append('{var} >> {hist}'.format(var=var2draw, hist=hname))
            
    return var_tuples, multihists



def applyHistStyle(h, i):
#    print h, i

    if h.GetName().find('data_sr')!=-1:
        h.SetLineColor(1)
#        h.SetMarkerSize(1)
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


def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def comparisonPlots(hists, titles, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he', prefix=None):

    display = DisplayManager(pname, isLog, isRatio, 0.6, 0.7, doption)
    display.draw_legend = isLegend
#    display.isEff = isEff
    
    display.Draw(hists, titles, prefix)






ensureDir('Plots/' + options.channel + '/sb')
ensureDir('Plots/' + options.channel + '/comp')





multihists = {}

ddir = {}



phiveto = "!(tau_rhomass1_kk < 1.04) && !(tau_rhomass2_kk < 1.04)"
psiveto = "!(3.62 < tau_psimass1 && tau_psimass1 < 3.75) && !(3.62 < tau_psimass2 && tau_psimass2 < 3.75)"
ksveto  = "!(tau_rhomass1_pk > 0.84 && tau_rhomass1_pk < 0.94) && !(tau_rhomass2_pk > 0.84 && tau_rhomass2_pk < 0.94)"
Bveto   = "!(tau_psimass1_kp > 5.15 && tau_psimass1_kp < 5.38) && !(tau_psimass2_kp > 5.15 && tau_psimass2_kp < 5.38) "
Bveto2  = "!(5 < b_mass && b_mass < 5.35)"
#xgbscut = "xgbs > 9."
#xgbscut = "xgbs > 9.6"
#xgbscut = "xgbs > 9.6 && tau_fls3d > 7." 
#xgbscut = "xgbs > 2." 
#xgbscut = "xgbs > 11." 
#xgbscut = "xgbs > 10."


#cut = '&&'.join([phiveto, psiveto, xgbscut])
rhomass1 = "tau_rhomass1 > 0.65 && tau_rhomass1 < 0.85"
rhomass2 = "tau_rhomass2 > 0.65 && tau_rhomass2 < 0.85"

#cut = xgbscut + '&& tau_index == 0'
#cut =  'perEVT_dnn > 0.95' #xgbscut
#cut =  'tau_vprob > 0.1 && tau_fls3d > 3. && xgbs > 10.' #xgbscut
#cut = 'xgbs > 8.'
#cut = 'tau_sumofdnn > 2.5 && tau_q==1'
#cut = 'tau_pt > 3. && mu1_isLoose==1 && mu2_isLoose==1 && xgbs > 9.5 && (pi1_trigMatch==1 || pi2_trigMatch==1 || pi3_trigMatch==1)'
cut = 'tau_pt > 3. && mu1_isLoose==1 && mu2_isLoose==1 && (pi1_trigMatch==1 || pi2_trigMatch==1 || pi3_trigMatch==1)'

cut = '&&'.join([cut, phiveto])

print 'cut = ', cut

#prefix = '../updateTuple/final_root_pt_save/'
#prefix_pnfs ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_multiple/'
#prefix ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_multiple/'
#prefix ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_mass_pt/'
prefix ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt/'
#prefix_pnfs ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_dnn/'

#taumass_sr = '&&'.join([rhomass1, rhomass2])
#taumass_cr = '!(' + taumass_sr + ')'

#print 'taumass_sr', taumass_sr
#print 'taumass_cr', taumass_cr


#ddir['data'] = {'file':prefix_pnfs + 'Data_2018/Myroot.root', 'acut':cut, 'weight':'1', 'lname':'data'}
#ddir['bg_ul'] = {'file':prefix_pnfs + 'BcJpsiX_ul_2018/Myroot.root', 'acut':cut, 'weight':'1', 'lname':'bg_ul'}
#ddir['sig'] = {'file':prefix_pnfs + 'BcJpsiTau_inclusive_ul_all_2018/Myroot.root', 'acut':cut, 'weight':'0.115', 'lname':'sig'}
#ddir['sig_truth'] = {'file':prefix_pnfs + 'BcJpsiTau_inclusive_ul_all_2018/Myroot.root', 'acut':cut + '&& tau_isRight_3prong', 'weight':'0.115', 'lname':'truth'}

ddir['data'] = {'file':prefix + '/Data_2018/data.root', 'acut':cut, 'weight':'1', 'lname':'data', 'scale':1}
ddir['bg_ul'] = {'file':prefix + '/BcJpsiX_ul_2018/bkg.root', 'acut':cut, 'weight':'puweight', 'lname':'bg_ul', 'scale':7/0.8}
ddir['sig'] = {'file':prefix + '/BcJpsiTau_inclusive_ul_all_2018/signal.root', 'acut':cut + '&& tau_isRight_3prong==1', 'weight':'puweight', 'lname':'sig', 'scale':0.45/(3.*0.8)}
#ddir['sig_3p'] = {'file':prefix + '/BcJpsiTau_inclusive_ul_all_2018/signal.root', 'acut':cut + '&& tau_isRight_3prong==1', 'weight':'puweight', 'lname':'truth', 'scale':ddir['sig']['scale']}
#ddir['sig_3pp'] = {'file':prefix + 'signal.root', 'acut':cut + '&& tau_isRight_3prong==1 && tau_isRight_3prong_pi0==1', 'weight':'puweight', 'lname':'truth', 'scale':ddir['sig_3p']['scale']}
#ddir['sig_others'] = {'file':prefix + 'signal.root', 'acut':cut + '&& tau_isRight_3prong==0', 'weight':'puweight', 'lname':'truth', 'scale':ddir['sig_3p']['scale']}

#ddir['data'] = {'file':prefix + 'Myroot_data_2018_0.root', 'acut':cut, 'weight':'1', 'lname':'data', 'scale':1}
#ddir['bg_ul'] = {'file':prefix + 'Myroot_bkg_ul_2018_0.root', 'acut':cut, 'weight':'puweight', 'lname':'bg_ul', 'scale':7/0.8}
#ddir['sig'] = {'file':prefix + 'Myroot_sig_ul_0.root', 'acut':cut, 'weight':'puweight*0.115', 'lname':'sig', 'scale':0.45/(3.*0.8)}
#ddir['sig_truth'] = {'file':prefix + 'Myroot_sig_ul_0.root', 'acut':cut + '&& tau_isRight_3prong', 'weight':'puweight*0.115', 'lname':'truth', 'scale':ddir['sig']['scale']}






finaldiscriminant = 'tau_rhomass1'

isquick = False
if isquick:
    for vkey, ivar in vardir.items():
        if not vkey in ['xgbs', 'perEVT_dnn']:
            vardir.pop(vkey)


isxgbs = True

if not isxgbs:
    for vkey, ivar in vardir.items():
        if vkey.find('xgbs')!=-1:
            vardir.pop(vkey)


print '-'*80

for vkey, ivar in vardir.items():
    print '->', vkey.ljust(20), ivar

print '-'*80

            
for type, ivar in ddir.items():
    
    file = TFile(ivar['file'])
    tree = file.Get('tree')

    print ivar['file'], 'is read ...'

    var_tuples, hists = returnTuples(type, vardir)


    cut = '(' + ivar['acut'] + ')*' + ivar['weight']
    print(type, cut)

    tree.MultiDraw(var_tuples, cut)
    
    multihists.update(copy.deepcopy(hists))



writes = {}

for vkey, ivar in vardir.items():

    hists = []
    titles = []

    for type, ivar2 in ddir.items():

#        print(vkey, type)

        multihists[type + '_' + vkey].Scale(Double(ivar2['scale']))        

        if vkey in [finaldiscriminant]: 

            __hist = copy.deepcopy(multihists[type + '_' + vkey])
            writes[type +'_'+vkey] = __hist

        if type.find('ff')==-1: 
#            print multihists[type + '_' + vkey].GetName(), multihists[type + '_' + vkey].GetEntries(), multihists[type + '_' + vkey].GetSumOfWeights(), Double(ivar2['scale'])

            hists.append(copy.deepcopy(multihists[type + '_' + vkey]))
#            titles.append(ivar2['lname'] + ' (' + str(int(multihists[type + '_' + vkey].GetSumOfWeights())) + ', ' + str(int(multihists[type + '_' + vkey].GetEntries())) + ')')
            titles.append(ivar2['lname'])


    hist_mc = None

    if hists[0].GetSumOfWeights()==0:
        print vkey, ' is skipped!'
        continue

    for ii, ihist in enumerate(hists):
        applyHistStyle(ihist, ii)
#        overflow(ihist)

        if ihist.GetSumOfWeights()==0:
            print('!!!!', vkey, ihist.GetName(), 'does not have any entries ...')
        else:
            ihist.Scale(1./ihist.GetSumOfWeights())
            ihist.SetMaximum(ihist.GetBinContent(ihist.GetMaximumBin())*1.2)

            if (ihist.GetName().find('bg_ul')!=-1 or ihist.GetName().find('sig')!=-1) and ihist.GetName().find('sig_truth')==-1:
                print ihist.GetName(), 'is added ...'
                if hist_mc == None:
                    hist_mc = copy.deepcopy(ihist)
                else:
                    hist_mc.Add(copy.deepcopy(ihist))


    applyHistStyle(hist_mc, len(hists))
#    hists.append(hist_mc)
#    titles.append('total_mc')

    comparisonPlots(hists, titles, ivar.has_key('isLog'), 'Plots/' + options.channel + '/comp/' + vkey + '.pdf', False, True, 'HpE')

    # make significance plots 

    continue

    if vkey!='xgbs_fine': continue

    h_data = copy.deepcopy(multihists['data_' + vkey])

#    h_sig = copy.deepcopy(multihists['sig_' + vkey])
#    h_significance = copy.deepcopy(h_sig)

    h_sig_truth = copy.deepcopy(multihists['truth_' + vkey])
    h_significance_truth = copy.deepcopy(h_sig_truth)

    h_sig_truth.SetMarkerColor(kRed)
    h_sig_truth.SetLineColor(kRed)
    h_significance_truth = copy.deepcopy(h_sig_truth)
    h_significance_truth.SetMarkerColor(kRed)
    h_significance_truth.SetLineColor(kRed)
    

    for ibin in range(1, h_data.GetXaxis().GetNbins()+1):
#        ndata = h_data.GetBinContent(ibin)
#        nsig = h_sig.GetBinContent(ibin)

#        if ivar['isRight']:
#            ndata = h_data.Integral(ibin, h_data.GetXaxis().GetNbins()+1)
#            nsig = h_sig.Integral(ibin, h_data.GetXaxis().GetNbins()+1)
#            ntruth = h_sig_truth.Integral(ibin, h_data.GetXaxis().GetNbins()+1)
#        else:
        ndata = h_data.Integral(ibin, h_data.GetXaxis().GetNbins()+1)
#        nsig = h_sig.Integral(ibin, h_data.GetXaxis().GetNbins()+1)
        ntruth = h_sig_truth.Integral(ibin, h_data.GetXaxis().GetNbins()+1)

#        ndata -= ntruth

#        sig, sig_err, frac, frac_err = calc(nsig, ndata)
#
#
#        h_significance.SetBinContent(ibin, sig)
#        h_significance.SetBinError(ibin, sig_err)
#        
#        h_sig.SetBinContent(ibin, frac)
#        h_sig.SetBinError(ibin, frac_err)

        sig, sig_err, frac, frac_err = calc(ntruth, ndata)

        h_significance_truth.SetBinContent(ibin, sig)
        h_significance_truth.SetBinError(ibin, sig_err)
        
        h_sig_truth.SetBinContent(ibin, frac)
        h_sig_truth.SetBinError(ibin, frac_err)


#        if ndata!=0:

#            err_sig = math.sqrt(nsig/ndata + nsig*nsig/(4*ndata*ndata))
#            err = math.sqrt(nsig/(ndata*ndata) + nsig*nsig/(ndata*ndata*ndata))

#            h_significance.SetBinContent(ibin, nsig/math.sqrt(ndata))
#            h_significance.SetBinError(ibin, err_sig)
#            
#            h_sig.SetBinContent(ibin, nsig/ndata)
#            h_sig.SetBinError(ibin, err)
#
#        else:
#            h_significance.SetBinContent(ibin, 0)
#            h_significance.SetBinError(ibin, 0)
#            h_sig.SetBinContent(ibin, 0)
#            h_sig.SetBinError(ibin, 0)

#    h_sig.Divide(h_data)

    h_sig.GetYaxis().SetTitle('S/(S+B)')
    comparisonPlots([h_sig, h_sig_truth], ['sig', 'truth'], ivar.has_key('isLog'), 'Plots/' + options.channel + '/sb/' + vkey + '.pdf', False, True, 'pzle')

    h_significance.GetYaxis().SetTitle('S/sqrt(S+B)')
    comparisonPlots([h_significance, h_significance_truth], ['sig', 'truth'], ivar.has_key('isLog'), 'Plots/' + options.channel + '/sb/sig_' + vkey + '.pdf', False, True, 'pzle')




#if not isxgbs:
#    sys.exit(1)
#sys.exit(1)

print writes

ofile = TFile('datacard.root', 'recreate')

#_var =  vardir['xgbs']

ofile.cd()
ofile.mkdir('sr')
ofile.cd('sr')

#nsb = writes['data_sb_xgbs'].GetSumOfWeights()
#nsr = writes['data_sr_xgbs'].GetSumOfWeights()


for name, hist in writes.items():

#    overflow(hist)

    print hist.GetName()

#    if hist.GetName().find('data')!=-1:
##        print writes['cr_' + finaldiscriminant], writes['data_sr_' + finaldiscriminant].GetSumOfWeights()
#        h_bg = copy.deepcopy(writes['data_' + finaldiscriminant])
#        h_bg.Add(copy.deepcopy(writes['sig_' + finaldiscriminant]), -1)
#        h_bg.SetTitle('bg')
#        h_bg.SetName('bg')
#        h_bg.Write()
        

    if hist.GetName().find('data')!=-1: 
        hist.SetName('data_obs')
        hist.SetTitle('data_obs')
    else:
        hist.SetName(hist.GetName().split('_')[0])
        hist.SetTitle(hist.GetName().split('_')[0])


    hist.Write()
       

    # Define dummy 1 bin histograms ... 
#    for ibin in range(1, _var['nbin']+1):
#    for ibin in range(1, len(varbin)):
#    print "nbins = ", writes['sig_sr_xgbs'].GetXaxis().GetNbins()+1
#    for ibin in range(1, writes['sig_sr_xgbs'].GetXaxis().GetNbins()+1):

#        dummy = TH1F('bg_bin' + str(ibin), 'bg_bin' + str(ibin), len(varbin)-1, array('d', varbin))
#        dummy = TH1F('bg_bin' + str(ibin), 'bg_bin' + str(ibin), writes['sig_sr_xgbs'].GetXaxis().GetNbins(), writes['sig_sr_xgbs'].GetXaxis().GetXmin(), writes['sig_sr_xgbs'].GetXaxis().GetXmax())


#        dummy.SetBinContent(ibin, 1)
        
#        dummy.SetBinError(ibin, 0)
    
#        dummy.Write()

    
ofile.Write()
ofile.Close()
