import copy, os,  math, sys, shutil
from numpy import array
from common.officialStyle import officialStyle
from array import array
import common.MultiDraw
import numpy as np
from varConfig import vardir
from common.DataMCPlot import *
from common.DisplayManager_postfit import DisplayManager
from common.DisplayManager_compare import DisplayManager_compare
from common.helper import *
from common.H2TauStyle import *
import ConfigParser
from optparse import OptionParser, OptionValueError

#gROOT.SetBatch(False)
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)

gROOT.ProcessLine(".L ~/tool//MultiDraw.cc+");
gROOT.Macro('common/functionmacro.C+')

usage = "usage: python compare.py" 
parser = OptionParser(usage) 

parser.add_option("-s", "--sys", default="None", type="string", dest="sys")
parser.add_option('-m', '--min', action="store_true", default=False, dest='min')
parser.add_option('-b', '--blind', action="store_true", default=False, dest='blind')
parser.add_option("-y", "--year", default="2018", type="string", dest="year")
parser.add_option("-o", "--outdir", default=".", type="string", dest="outdir")
#parser.add_option('-c', '--create', action="store_true", default=False, dest='create')

(options, args) = parser.parse_args() 

#print(options)

print '-'*80
print 'sys=', options.sys
print 'min=', options.min
print 'blind=', options.blind
print 'year=', options.year
print '-'*80


init = ConfigParser.SafeConfigParser()
init.read("./settings.ini")
lumi = init.get('common', 'lumi_' + options.year)
#filedir = init.get('common', 'filedir')

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

def drawSignificance(hists, channel, vkey, pathname, isRight=True):
    
    data_obs = None
    
    for hist in hists:
        if hist.GetName().find('data')!=-1:
            data_obs = copy.deepcopy(hist)

        if hist.GetName().find('bc_jpsi_tau_3p')!=-1:
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
    
    canvas.SaveAs(pathname + '/sig_' + vkey  + '.gif')
    
        


multihists = {}

##################################################

#Just for past reference:

#prefix_yuta ='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi_Legacy_decayBc/job_pt/'
#prefix ='/pnfs/psi.ch/cms/trivcat/store/user/cgalloni/RJpsi_Legacy_decayBc_FromYuta_20220317/job_pt/'
##prefix ='/pnfs/psi.ch/cms/trivcat/store/user/cgalloni/RJpsi_Legacy_vprob_fls3d_invertedOR_v3/job_pt/'

#datastr = "Data_2018/data.root"
#sigstr  = "BcJpsiTau_inclusive_ul_all_2018/sig_20220324.root"
##sigstr  = "BcJpsiTau_inclusive_ul_all_2018/sig.root"    
##datastr = "Data_2018/Myroot_training_weightAdded.root"
##sigstr  = "BcJpsiTau_inclusive_ul_all_2018/Myroot_training_weightAdded.root"
#bkgstr  = "BJpsiX_ul_2018/bkg.root"
#=======
prefix = init.get('common', 'filedir') + '/job_inv_pt_' + options.year + '/'

datastr = init.get('common', 'data_prefix') + '/data.root'
sigstr  = init.get('common', 'sig_prefix') + '/sig.root'
bkgstr  = init.get('common', 'bkg_prefix') + '/bkg.root'


file_hammer = TFile(prefix + '/BcJpsiTau_inclusive/Myroot_0.root')
hist_hammer = file_hammer.Get('hammer')


bc_sf =  init.get('common', 'bc_sf')

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


mu_ID_weight =  "*mu1_SFID*mu2_SFID"
mu_Reco_weight =  "*mu1_SFReco*mu2_SFReco"
mu_weight =mu_ID_weight+mu_Reco_weight


ddir['data_obs'] =  {'file':datastr, 'weight':'1', 'scale':1, 'order':2999, 'color':1, 'addcut':'1'}

ddir['bc_jpsi_tau_3p'] =     {'file':sigstr, 'weight':'puweight*' + hammer_weight + mu_weight +'*weight_ctau', 'scale':bc_sf, 'order':1, 'color':taucol_1, 'addcut':'n_occurance==1 && tau_isRight_3prong==1 && gen_sig_decay==6 && weight_ctau < 100.'}

ddir['bc_jpsi_tau_N3p'] =     {'file':sigstr, 'weight':'puweight*' + hammer_weight + mu_weight+'*weight_ctau', 'scale':bc_sf, 'order':1, 'color':lob_1, 'addcut':'gen_sig_decay==6 && isJpsiTau2Mu!=1 &&  tau_isRight_3prong==0 && weight_ctau < 100.'}

ddir['bc_jpsi_dst'] = {'file':sigstr, 'weight':'puweight'+ mu_weight+'*weight_ctau', 'scale':bc_sf, 'order':4, 'color':jpsi_3, 'addcut':'(gen_sig_decay==10||gen_sig_decay==20)  && weight_ctau < 100.'}

ddir['bc_others'] = {'file':sigstr, 'weight':'puweight'+ mu_weight+'*weight_ctau', 'scale':bc_sf, 'order':4, 'color':others, 'addcut':'(!((n_occurance==1 && tau_isRight_3prong==1 && gen_sig_decay==6) || (gen_sig_decay==6 && isJpsiTau2Mu!=1 &&  tau_isRight_3prong==0) ||gen_sig_decay==10 ||gen_sig_decay==20 ))  && weight_ctau < 100.'}

ddir['bg_ul'] =      {'file':bkgstr, 'weight':'puweight'+ mu_weight+'*genWeightBkgB', 'scale':7*0.64/0.8, 'order':5, 'color':dycol, 'addcut':'1'}


# signal + Bc : All possible channel decays
#ddir['sig_lep'] = {'file':sigstr, 'weight':'puweight*(hammer_ebe/0.55)', 'scale':bc_sf, 'order':3, 'color':ttcol, 'addcut':'n_occurance==1 && tau_isRight_3prong==0 && isJpsiTau2Mu==1'}
#ddir['sig_others'] = {'file':sigstr, 'weight':'puweight*(hammer_ebe/0.55)', 'scale':bc_sf, 'order':3, 'color':wcol, 'addcut':'n_occurance==1 && tau_isRight_3prong==0 && isJpsiTau2Mu==0'}
#ddir['sig_others'] = {'file':sigstr, 'weight':'puweight*' + hammer_weight, 'scale':bc_sf, 'order':3, 'color':wcol, 'addcut':'n_occurance==1 && tau_isRight_3prong==0'}
#ddir['bg_bc'] =      {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':qcdcol, 'addcut':'n_occurance==0 && isJpsiMu==0'}
#ddir['bg_norm'] =      {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':jtfake, 'addcut':'n_occurance==0 && isJpsiMu==1'}
#ddir['bc_jpsi_tau_mu'] =     {'file':sigstr, 'weight':'puweight*' + hammer_weight, 'scale':bc_sf, 'order':1, 'color':taucol_1, 'addcut':'tau_isRight_3prong==0 && gen_sig_decay==6 && isJpsiTau2Mu==1'} # merged in Bc others

#ddir['bc_jpsi_mu'] =     {'file':sigstr, 'weight':'puweight*' + hammer_weight, 'scale':bc_sf, 'order':1, 'color':pink_1, 'addcut':'gen_sig_decay==0'} # merged in Bc others  

#ddir['bc_charmonium_mu'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':bc_mu_1, 'addcut':'gen_sig_decay>1&&gen_sig_decay<6'} #in SR 0
#ddir['bc_psi2s_mu'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':bc_mu_1, 'addcut':'gen_sig_decay==1'} #0
#ddir['bc_chic0_mu'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':bc_mu_2, 'addcut':'gen_sig_decay==2'} #0
#ddir['bc_chic1_mu'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':bc_mu_3, 'addcut':'gen_sig_decay==3'} #0
#ddir['bc_chic2_mu'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':bc_mu_4, 'addcut':'gen_sig_decay==4'} #0
#ddir['bc_hc_mu'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':bc_mu_5, 'addcut':'gen_sig_decay==5'}

#ddir['bc_chic0_tau'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':bc_mu_2, 'addcut':'gen_sig_decay==12'} #0
#ddir['bc_chic1_tau'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':bc_mu_3, 'addcut':'gen_sig_decay==13'} #0
#ddir['bc_chic2_tau'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':bc_mu_4, 'addcut':'gen_sig_decay==14'} #0
#ddir['bc_hc_tau'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':bc_mu_5, 'addcut':'gen_sig_decay==15'}  #0          

#ddir['bc_psi2s_tau'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':tau_1, 'addcut':'gen_sig_decay==7'} #0
#ddir['bc_jpsi_pi'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':jpsi_1, 'addcut':'gen_sig_decay==8'}
#ddir['bc_jpsi_3pi'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':jpsi_2, 'addcut':'gen_sig_decay==9'}
#ddir['bc_jpsi_ds'] = {'file':sigstr, 'weight':'puweight'+ mu_weight+'* weight_ctau', 'scale':bc_sf, 'order':4, 'color':jpsi_3, 'addcut':'gen_sig_decay==10'}
#ddir['bc_jpsi_dsp'] = {'file':sigstr, 'weight':'puweight'+ mu_weight+'* weight_ctau', 'scale':bc_sf, 'order':4, 'color':jpsi_4, 'addcut':'gen_sig_decay==20'}
#ddir['bc_jpsi_5pi'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':jpsi_4, 'addcut':'gen_sig_decay==11'}
#ddir['bc_jpsi_pions'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':lob_1, 'addcut':'(gen_sig_decay==8||gen_sig_decay==9||gen_sig_decay==11)'}
#ddir['bc_others'] = {'file':sigstr, 'weight':'puweight', 'scale':bc_sf, 'order':4, 'color':others, 'addcut':'gen_sig_decay==16'}
#((gen_sig_decay==8||gen_sig_decay==9||gen_sig_decay==11)||(gen_sig_decay>1&&gen_sig_decay<6)||gen_sig_decay==7|| gen_sig_decay==16 ||( tau_isRight_3prong==0 && gen_sig_decay==6 && isJpsiTau2Mu==1)||gen_sig_decay==0)'}

# other B backgrounds 


##################################################


basic = init.get('common', 'basic')

xgbs_sr = 'xgbs > ' + init.get('common', 'sr_low')
xgbs_sr_xl = 'xgbs > ' + init.get('common', 'sr_low_xl')
xgbs_sb = 'xgbs > ' + init.get('common', 'sb_low') + ' && xgbs < ' + init.get('common', 'sb_high')
xgbs_lp = 'xgbs > ' + init.get('common', 'lp_low') + ' && xgbs < ' + init.get('common', 'lp_high')



channels = {

    'inclusive':{'cut':'&&'.join([basic])},
    'sr':{'cut':'&&'.join([basic, xgbs_sr])},
    'sb':{'cut':'&&'.join([basic, xgbs_sb])},
    'lp':{'cut':'&&'.join([basic, xgbs_lp])}, 


#    'inclusive_blind':{'cut':'&&'.join([basic])},
#    'inclusive_blind_all':{'cut':'&&'.join([basic])},
#    'inclusive_blind_all_bptg10':{'cut':'&&'.join([basic,'b_pt>10'])},
#    'inclusive_bptg10':{'cut':'&&'.join([basic,'b_pt>10'])},
#    'inclusive_psi2s':{'cut':'&&'.join([basic, psi2s])},
#    'extrapolate':{'cut':'&&'.join([basic, 'xgbs > 3.5'])},
#    'sr_bptg20':{'cut':'&&'.join([basic, xgbs_sr ,'b_pt>20' ])},
#    'sb_bptg20':{'cut':'&&'.join([basic, xgbs_sb ,'b_pt>20' ])},
#    'sr_highMass':{'cut':'&&'.join([basic, xgbs_sr, 'b_mass>6.5'])},
#    'sb_highMass':{'cut':'&&'.join([basic, xgbs_sb, 'b_mass>6.5'])},
#    'lp_highMass':{'cut':'&&'.join([basic, xgbs_lp, 'b_mass>6.5'])},      
#    'inclusive_highMass':{'cut':'&&'.join([basic, 'b_mass>6.5'])}, 
#    'sr_highMassV2':{'cut':'&&'.join([basic, xgbs_sr, 'b_mass>6'])},
#    'sb_highMassV2':{'cut':'&&'.join([basic, xgbs_sb, 'b_mass>6'])},
#    'lp_highMassV2':{'cut':'&&'.join([basic, xgbs_lp, 'b_mass>6'])},
#    'inclusive_highMassV2':{'cut':'&&'.join([basic, 'b_mass>6'])},
#    'sr_lowMassV2':{'cut':'&&'.join([basic, xgbs_sr,'b_mass>3.5','b_mass<4'])},
#    'sb_lowMassV2':{'cut':'&&'.join([basic, xgbs_sb,'b_mass>3.5', 'b_mass<4'])},
#    'lp_lowMassV2':{'cut':'&&'.join([basic, xgbs_lp,'b_mass>3.5', 'b_mass<4'])},
#    'inclusive_lowMassV2':{'cut':'&&'.join([basic,'b_mass>3.5', 'b_mass<4'])},
#    'sr_lowMass':{'cut':'&&'.join([basic, xgbs_sr, 'b_mass<4'])},
#    'sb_lowMass':{'cut':'&&'.join([basic, xgbs_sb, 'b_mass<4'])},
#    'lp_lowMass':{'cut':'&&'.join([basic, xgbs_lp, 'b_mass<4'])},
#    'inclusive_lowMass':{'cut':'&&'.join([basic, 'b_mass<4'])},
#    'sr_mediumMassV2':{'cut':'&&'.join([basic, xgbs_sr,'b_mass>4','b_mass<6'])},   
#    'sb_mediumMassV2':{'cut':'&&'.join([basic, xgbs_sb,'b_mass>4', 'b_mass<6'])},
#    'lp_mediumMassV2':{'cut':'&&'.join([basic, xgbs_lp,'b_mass>4', 'b_mass<6'])},
#    'inclusive_mediumMassV2':{'cut':'&&'.join([basic,'b_mass>4', 'b_mass<6'])},       
#    'sr_mediumMass':{'cut':'&&'.join([basic, xgbs_sr,'b_mass>4','b_mass<6.5'])},
#    'sb_mediumMass':{'cut':'&&'.join([basic, xgbs_sb,'b_mass>4', 'b_mass<6.5'])},
#    'lp_mediumMass':{'cut':'&&'.join([basic, xgbs_lp,'b_mass>4', 'b_mass<6.5'])},
#    'inclusive_mediumMass':{'cut':'&&'.join([basic,'b_mass>4', 'b_mass<6.5'])},
#    'lp':{'cut':'&&'.join([basic, xgbs_lp])},
    'sr_xl':{'cut':'&&'.join([basic, xgbs_sr_xl])},  
#    'sr_xs':{'cut':'&&'.join([basic, xgbs_sr_xs])},  
#    'cr_sr':{'cut':'&&'.join([basic, xgbs_sr])},
#    'cr_sb':{'cut':'&&'.join([basic, xgbs_sb])},
#    'cr_lp':{'cut':'&&'.join([basic, xgbs_lp])},
#    'q3_sr':{'cut':'&&'.join([basic, xgbs_sr])},
#    'q3_sb':{'cut':'&&'.join([basic, xgbs_sb])},
#    'q3_lp':{'cut':'&&'.join([basic, 'xgbs > 3.'])},
#    'hp_sr':{'cut':'&&'.join([basic, taumass, xgbs_sr])},
#    'hp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_sr])},
#    # lp xgbs sideband
#    'lp_sr':{'cut':'&&'.join([basic, taumass, xgbs_sb])},
#    'lp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_sb])},
##    # mp xgbs sideband
#    'slp_sr':{'cut':'&&'.join([basic, taumass, xgbs_lp])},
#    'slp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_lp])},
#    'cr_hp_sr':{'cut':'&&'.join([basic, taumass, xgbs_sr])},
#    'cr_hp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_sr])},
##    # lp xgbs sideband
#    'cr_lp_sr':{'cut':'&&'.join([basic, taumass, xgbs_sb])},
#    'cr_lp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_sb])},
##    # mp xgbs sideband
#    'cr_slp_sr':{'cut':'&&'.join([basic, taumass, xgbs_lp])},
#    'cr_slp_sb':{'cut':'&&'.join([basic, '!' + taumass, xgbs_lp])},

}


#finaldiscriminant = ['xgbs', 'xgbs_zoom', 'xgbs_sigscan', 'b_mass', 'b_mass_sf', 'tau_rhomass_unrolled', 'tau_rhomass_unrolled_coarse', 'q2_simple', 'jpsi_kpipi']
#finaldiscriminant = ['tau_rhomass_unrolled', 'tau_rhomass_unrolled_coarse', 'tau_rhomass_unrolled_coarse_16', 'q2_simple', 'xgbs', 'xgbs_zoom', 'xgbs_sigscan', 'b_mass' ]
finaldiscriminant = ['tau_rhomass_unrolled_coarse_16']

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


    dirname = options.outdir + '/' + options.year + '_' + channel + '_' + options.sys
    ensureDir(dirname + '/plots/')   
    shutil.copyfile('index.php', dirname + '/index.php')
    shutil.copyfile('index.php', dirname + '/plots/index.php') 
    ensureDir(dirname + '/datacard/')

    for type, ivar in ddir.items():

        wstr = ivar['weight']
        waddcut = ivar['addcut']
        filename = None    

        if channel.find('cr')!=-1:
            filename = prefix_cr + '/' + ivar['file']
        else:        
            filename = prefix + '/' + ivar['file']

        file = TFile(filename)
        tree = file.Get('tree')

        print(filename, 'is read ...')

        var_tuples, hists = returnTuples(type, vardir)    

        if channel=='sr' and options.blind==True and type =='data_obs': 
            waddcut += '&&0'

#        if channel.find('inclusive_blind_all')!=-1 :#   and type =='data_obs':
#            waddcut += '&&xgbs<4.3'        
#
#        if channel.find('inclusive_blind')!=-1  and type =='data_obs':
#            waddcut += '&&xgbs<4.3'

        if options.sys.find('hammer')!=-1 and type.find('bc_jpsi_tau')!=-1: 
            wstr = wstr.replace('hammer_ebe', options.sys)

        elif options.sys.find('ctau')!=-1 and type.find('bc')!=-1:
            wstr = wstr.replace('weight_ctau', options.sys)
            print "CTAU_SYST, wieght is " , wstr

        elif options.sys.find('puweight')!=-1 and type.find('data')==-1:
            wstr = wstr.replace('puweight', options.sys)

        elif options.sys.find('muSFID')!=-1 and type.find('data')==-1:
            wstr = wstr.replace('mu1_SFID', 'mu1_SFID_up' if options.sys.find('up')!=-1 else 'mu1_SFID_down') 
            wstr = wstr.replace('mu2_SFID', 'mu2_SFID_up' if options.sys.find('up')!=-1 else 'mu2_SFID_down')

        elif options.sys.find('muSFReco')!=-1 and type.find('data')==-1:     
            wstr = wstr.replace('mu1_SFReco', 'mu1_SFReco_up' if options.sys.find('up')!=-1   else 'mu1_SFReco_down') 
            wstr = wstr.replace('mu2_SFReco', 'mu2_SFReco_up' if options.sys.find('up')!=-1   else 'mu2_SFReco_down') 

        elif options.sys.find('br_BcJpsiDst')!=-1 and type.find('bc_jpsi_ds')!=-1:
            wstr += '*1.38' if  options.sys.find('up')!=-1 else '*0.62'

        elif options.sys.find('tauBr')!=-1 and type.find('bc_jpsi_tau_3p')!=-1:
            wstr += '*getTauBrWeight_up(gen_dipion_unrolled)' if options.sys.find('up')!=-1 else '*getTauBrWeight_down(gen_dipion_unrolled)'
        
        elif options.sys.find('tauReco')!=-1 and type.find('bc')!=-1:
            wstr += '*1.03' if options.sys.find('up')!=-1 else '*0.97'

        elif options.sys.find('xgbsEff')!=-1 and type.find('bc')!=-1:
            wstr += '*1.05' if options.sys.find('up')!=-1 else '*0.95'

        elif options.sys.find('BcPt')!=-1 and type.find('bc')!=-1:
            wstr += '*getBcWeight(B_pt_gen,1)'if options.sys.find('up')!=-1 else '*getBcWeight(B_pt_gen,-1)'


        cut = '(' + dict['cut'] + ' &&' + waddcut + ')*' + wstr
        print(type, cut)        

        tree.MultiDraw(var_tuples, cut)
    
        multihists.update(copy.deepcopy(hists))




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
       
        comparisonPlots(Histo, dirname + '/plots/' + vkey + '.gif')

        comparisonPlots(Histo, dirname + '/plots/' + vkey + '_log.gif', True)
            

        if vkey in finaldiscriminant:

            Histo.WriteDataCard(dirname + '/datacard/' + vkey + '.root', True, 'RECREATE', channel)

        if options.sys=='None' and vkey.find('xgbs')!=-1:
            drawSignificance(hists, channel, vkey, dirname + '/plots/', True)
        
        comparisonPlots_alt(hists, titles, True, False, dirname + '/plots/' + vkey + '_compare.pdf', False, True, 'pE')


    
