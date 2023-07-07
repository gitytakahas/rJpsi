import copy, os,  math, sys, shutil 
from numpy import array
from common.officialStyle import officialStyle
from array import array
import numpy as np
from varConfig import vardir
from common.DataMCPlot import *
from common.DisplayManager_postfit import DisplayManager
from common.DisplayManager_postfit_wide import DisplayManager_wide
from common.DisplayManager_compare import DisplayManager_compare
from common.helper import *
from common.H2TauStyle import *
import ConfigParser
from optparse import OptionParser, OptionValueError

usage = "usage: python compare.py" 
parser = OptionParser(usage) 
#parser.add_option('-m', '--min', action="store_true", default=False, dest='min')
parser.add_option("-y", "--year", default="all", type="string", dest="year")
parser.add_option('-s', '--scale', action="store_true", default=False, dest='scale')

(options, args) = parser.parse_args() 

process_convention={
    'bc_jpsi_tau':'jpsi_tau',
    'bc_jpsi_dst':'jpsi_hc'
}

nuisance_convention={
    'hammer_ebe':'bglvar',
    'weight_ctau':'ctau',
    'BcPt':'bccorr',
    'puweight':'puWeight',
    'br_BcJpsiDst':'br_jpsi_hc_over_mu',
    'muSFID':'sfIdJpsi',
    'muSFReco':'sfReco'
}


#gROOT.SetBatch(False)
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)
#gStyle.SetOptStat(0)

sf_signal = 0.6389

prefix='tauhad'
systs_hammer = []

for hammer in range(0, 10):
    systs_hammer.append('hammer_ebe_e' + str(hammer) + '_up')
    systs_hammer.append('hammer_ebe_e' + str(hammer) + '_down')

systs_mc = []
#systs_name_mc=['puweight', 'trigger', 'muSFID', 'muSFReco', 'weight_ctau', 'br_BcJpsiDst', 'tauBr',  'tauReco', 'xgbsEff', 'BcPt']
systs_name_mc=['puweight', 'weight_ctau', 'tauBr',  'BcPt', 'BcOthersShape']

for syst in systs_name_mc:
    for ud in ['up', 'down']:
        systs_mc.append(syst + '_' + ud)


#systs_bkg = []
#systs_name_bkg = ["bkgExtra","bkgExtra"]
#
#for syst in systs_name_bkg:
#    for ud in ['up', 'down']:
#        systs_bkg.append(syst + '_' + ud)

datacardpath = '/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/results_simultaneous'
#datacardpath = '/work/cgalloni/Rjpsi_analysis/CMSSW_10_2_10/src/rJpsi/anal/plots_inv_dir_Yuta/'
categories = ['sr','sb']

#finaldiscriminant = ['tau_rhomass_unrolled', 'tau_rhomass_unrolled_coarse']
finaldiscriminant = ['tau_rhomass_unrolled_var']
#finaldiscriminant = ['tau_rhomass_unrolled_coarse']

vardir["tau_rhomass_unrolled_var"] = {'tree':'tree'} # dummy to let the loop below will go through unrolled variable binning distribution


for vkey, ivar in vardir.items():
    if vkey not in finaldiscriminant: 
        vardir.pop(vkey)

print '-'*80

for vkey, ivar in vardir.items():
    print '->', vkey.ljust(20), ivar

print '-'*80



#ratio = 0.092
#fitCat = 'sr'
output='combine_simultaneous'
ensureDir(output)


init = ConfigParser.SafeConfigParser()
init.read("./settings.ini")


eras = ['2016', '2017', '2018']
if options.year!='all':
    eras = [options.year]


def applyHists(hists):
    
    colors = [1, 2, 4]

    for idx, hist in enumerate(hists):
        hist.SetFillStyle(0)
        hist.SetFillColor(0)
        hist.SetLineColor(colors[idx])
        hist.SetMarkerColor(colors[idx])
        hist.SetMarkerSize(0.5)
        hist.SetLineWidth(idx+1)
        hist.SetLineStyle(idx+1)
        hist.GetXaxis().SetTitle('Bin ID of 2-dim. #rho_{1} vs #rho_{2}')
        hist.GetXaxis().SetTitleSize(0.05)
        hist.GetXaxis().SetTitleOffset(0.9)
#        hist.GetXaxis().SetLabelSize(0.05)
#        hist.GetXaxis().SetLabelSize(0.1)



def comparisonPlots(hist, lumi, pname='sync.pdf', isLog=False, isRatio=True, clabel=''):

    display = DisplayManager(pname, isRatio, isLog, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)

def comparisonPlots_wide(hist, lumi, pname='sync.pdf', isLog=False, isRatio=True, clabel=''):

    display = DisplayManager_wide(pname, isRatio, isLog, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)

def comparisonPlots_alt(hists, titles, norm=False, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he'):

    display = DisplayManager_compare(pname, isLog, isRatio, 0.4, 0.7, doption)
    display.draw_legend = isLegend
    
    display.Draw(hists, titles, norm)


def getHist(year, vkey, channel, target, sys='None'):

    filename = datacardpath + '/' + year + '_' + channel + '_' + sys + '/datacard/' + vkey + '.root'

    if not os.path.isfile(filename):
        print 'This file is corrupted!!!'

    file = TFile(filename)

    _hist = copy.deepcopy(file.Get(channel + '/' + target))

    if _hist.GetEntries()==0:
        print vkey, channel, 'is not present!!!'

    file.Close()
    
    return _hist


def setNameTitle(hist, name):
    hist.SetTitle(name)
    hist.SetName(name)

def draw(year, vkey, channels, target, sys='None', subtract=False, saveFig=False, opt='None', sf = 0.265):


    hists = []
    titles = []
    
    for ii, channel in enumerate(channels):

        _hist = getHist(year, vkey, channel, target, sys)
        print( " draw _hist vkey ", vkey, " channel " ,channel, " target " ,target, " sys ", sys , "Nbins ", _hist.GetNbinsX()) 

        if subtract:

            for proc in ['bc_jpsi_tau', 'bc_jpsi_dst','bc_others']: 

                _hist2 = getHist(year, vkey, channel, proc, sys) #copy.deepcopy(file.Get(channel + '/' + proc))

                if proc in ['bc_jpsi_tau']:
                    _hist2.Scale(sf)

                _hist.Add(_hist2, -1)
        
        _hist.SetFillStyle(0)
        _hist.SetMarkerColor(ii+1)
        _hist.SetLineColor(ii+1)
        setNameTitle(_hist, target + '_' + channel)
        hists.append(_hist)
        titles.append(channel)




    hists2return = copy.deepcopy(hists)    


    _dirname = output + '/' + '_'.join(channels) + '_' + target
    ensureDir(_dirname)  
    shutil.copyfile('index.php',output+'/index.php')      
    shutil.copyfile('index.php',_dirname+'/index.php')

    if saveFig:

        comparisonPlots_alt(hists, titles, True, False, _dirname + '/' + vkey +   '_' + year + '.pdf', True, True,   opt+'hpe' )

    file_output = TFile(_dirname + '/' + vkey +  '.root', 'recreate')
    for _hist in hists2return:
        _hist.Write()
    file_output.Write()
    file_output.Close()

        

    return hists2return



for vkey, ivar in vardir.items():

    filename_new = output + '/' + vkey + '.root'

    if options.scale:
        print 'this is scaled!'
        filename_new = filename_new.replace('.root', '_scaled.root')

    file_new = TFile(filename_new, 'recreate')
        
    for year in eras:

        lumi = init.get('common', 'lumi_' + year)

        for fitCat in categories:
   
#            draw(year, vkey, ['lp',  'sb'], 'data_obs', 'None', True, True,'None')

#            draw(year, vkey, ['sb',  'sr'], 'data_obs', 'None', True, True, 'None')



            hists2write = []


#            print("lp-sb comparison")    
#            draw(year, vkey, ['lp',  'sb'], 'data_obs', 'None', True, True, 'ks')
#            print("sb-sr comparison")
#            draw(year, vkey, ['sb',  'sr'], 'data_obs', 'None', True, True, 'ksratio')
#            print("sb-sr_xl comparison")
#            draw(year, vkey, ['sb',  'sr_xl'], 'data_obs', 'None', True, True, 'ksratio')
#            #print("sb-sr comparison hist preparation bkg ")
#            hists4ddbkg_bgmc = draw(year, vkey, ['sb', 'sr'], 'bg_ul', 'None', False, True, 'ks')
#            print("lp-sb comparison hist preparation data")
#            hists4ddbkg_vr = draw(year, vkey, ['lp', 'sb'], 'data_obs', 'None', True, True, 'ks')
        
            file = TFile(datacardpath + '/' + year + '_' + fitCat + '_None/datacard/' + vkey + '.root')
            file.cd(fitCat)
        
            listofprocs = [key.GetName() for key in gDirectory.GetListOfKeys()]

            print listofprocs

    
            for proc in listofprocs:
                _tmp = file.Get(fitCat + '/' + proc)

                if options.scale and _tmp.GetName().find('jpsi_tau')!=-1:
                    _tmp.Scale(sf_signal)
                

                hists2write.append(copy.deepcopy(_tmp))
        
#            hists4ddbkg_sr = draw(year, vkey, ['sr', 'sb'], 'data_obs', 'None', True, False, 'None')
#            bkgHist = copy.deepcopy(hists4ddbkg_sr[1])

#            hist_sb = getHist(year, vkey, 'sb', 'bg_ul') 
#            hist_sr = getHist(year, vkey, 'sr', 'bg_ul') 

        #    ratio from stefanos
        # SR with weight 0.110209
        # SR without weight 0.0967332
        # ratio = 0.110209/0.0967332
        #    ratio = 1.1393089446*

#        ratio_sf = 0.110209/0.0967332
#        
#        if year=='2016':
#            ratio_sf = 0.096054/0.0819467
#        elif year=='2017':
#            ratio_sf = 0.123838/0.0996171
#
#        ratio = ratio_sf*float(hist_sr.GetSumOfWeights())/float(hist_sb.GetSumOfWeights())


##            ratio_file = TFile('correction/'+year + '.root')
##            ratio_graph = ratio_file.Get('Graph')
##
##            xgbs_sr = init.get('common', 'sr_low')
##        
##            if year=='2016':
##                xgbs_sr = init.get('common', 'sr_low_' + year)
##
##            ratio = ratio_graph.Eval(float(xgbs_sr))
##        
##            print '-'*80
##            print year, ': BDT >', xgbs_sr, ', SR/SB extrapolation ratio =', ratio
##            print '-'*80
##
##            bkgHist.Scale(ratio)
##            setNameTitle(bkgHist, 'dd_bkg')
##            hists2write.append(bkgHist)


##        ### draw histograms 
##
##        hists2draw = copy.deepcopy(hists2write)
##        print ( "Calling    Histo = DataMCPlot(vkey)    for  vkey ", vkey )
##        Histo = DataMCPlot(vkey)
##
##        for _hist in hists2draw:
##
##            _name = _hist.GetName()
##
##            if _name=='bg_ul': continue
##
##            if _name.find('data')!=-1:
##                _hist.SetFillStyle(0)
##                _hist.Sumw2(False)
##                _hist.SetBinErrorOption(1)
##            
##
##            _hist.SetName(_name)
##            _hist.SetTitle(_name)
##            _hist.GetXaxis().SetLabelColor(1)
###            _hist.GetXaxis().SetLabelSize(0.0)
##
##            Histo.AddHistogram(_hist.GetName(), _hist)
##            print ( " Histo = DataMCPlot(vkey) with _hist ", _hist.GetName() ," bins ", _hist.GetNbinsX() )
##            if _name.find('data')!=-1:
##                Histo.Hist(_hist.GetName()).stack = False



##        Histo._ApplyPrefs()
##        print(Histo)
##        comparisonPlots_wide(Histo, lumi, output + '/' + vkey + '_' + year +'.gif')
##        comparisonPlots_wide(Histo, lumi, output + '/' + vkey + '_log_' + year + '.gif', True)
###        comparisonPlots(Histo, lumi, output + '/' + vkey + '_' + year +'.pdf')
###        comparisonPlots(Histo, lumi, output + '/' + vkey + '_log_' + year + '.pdf', True)
##
##
##        ### shape variation 
##
##        hists4ddbkg_up = draw(year, vkey, ['sr', 'sb'], 'data_obs', 'None', True, False, 'None', 0.95)
##        bkgHist_up = copy.deepcopy(hists4ddbkg_up[1])
##        bkgHist_up.Scale(ratio)
##        setNameTitle(bkgHist_up, 'dd_bkg_shapeUp')
##        hists2write.append(bkgHist_up)
##    
##        hists4ddbkg_down = draw(year, vkey, ['sr', 'sb'], 'data_obs', 'None', True, False,'None',  0.2)
##        bkgHist_down = copy.deepcopy(hists4ddbkg_down[1])
##        bkgHist_down.Scale(ratio)
##        setNameTitle(bkgHist_down, 'dd_bkg_shapeDown')
##        hists2write.append(bkgHist_down)


            for sys in systs_hammer:

                name_sys = sys.replace('_up','Up').replace('_down', 'Down')
            
  #          hists_sys = draw(year, vkey, ['sr', 'sb'], 'data_obs', sys, True, False, 'None')


#            bkgHist_sys = hists_sys[1]
#            bkgHist_sys.Scale(ratio)
#            setNameTitle(bkgHist_sys, 'dd_bkg_' + name_sys)
#            hists2write.append(bkgHist_sys)
        
                sig_sys = getHist(year, vkey, fitCat, 'bc_jpsi_tau', sys)
            
                setNameTitle(sig_sys, 'bc_jpsi_tau_' + name_sys)

                if options.scale:
                    sig_sys.Scale(sf_signal)
            
                hists2write.append(sig_sys)




            for sys in systs_mc:

                name_sys = sys.replace('_up','Up').replace('_down', 'Down')
            
#            hists_sys = draw(year, vkey, ['sr', 'sb'], 'data_obs', sys, True, False, 'None')
            
#            bkgHist_sys = hists_sys[1]
#            bkgHist_sys.Scale(ratio)

#            setNameTitle(bkgHist_sys, 'dd_bkg_' + name_sys)
#            hists2write.append(bkgHist_sys)

                        
                sig_sys = getHist(year, vkey, fitCat, 'bc_jpsi_tau', sys)
                bc_others = getHist(year, vkey, fitCat, 'bc_others', sys)
                bc_jpsi_dst = getHist(year, vkey, fitCat, 'bc_jpsi_dst', sys)
        
                setNameTitle(sig_sys, 'bc_jpsi_tau_' + name_sys)
                setNameTitle(bc_others, 'bc_others_' + name_sys)
                setNameTitle(bc_jpsi_dst, 'bc_jpsi_dst_' + name_sys)     

                if options.scale:
                    sig_sys.Scale(sf_signal)

            
                hists2write.append(sig_sys)
                hists2write.append(bc_others)
                hists2write.append(bc_jpsi_dst)
        


#        for sys in systs_bkg:
#            continue
#            if vkey == 'tau_rhomass_unrolled' : 
#                continue
#            name_sys = sys.replace('_up','Up').replace('_down', 'Down')
#            hists_sys = draw(vkey, ['sr', 'sb'], 'data_obs', 'None', True, False, 'None')
#            bkgHist_sys = hists_sys[1]
#            bkgHist_sys.Scale(ratio)
#            #To be updated 
#            p0 = 0.871 if bkgHist_sys.GetNbinsX()< 40 else  0.850166  
#            p1 = 0.009 if bkgHist_sys.GetNbinsX()< 40 else  0.0018
#            for bin in range(0, bkgHist_sys.GetNbinsX()):
#                #To be updated          
#                file_ratio =  TFile("syst_bkg/"+vkey+"_"+year+"_ratio.root", 'open')
#                ratio_hist=file_ratio.Get("data_obs_sr")
#                if name_sys.find("Up")!=-1:   
#                    bkgHist_sys.SetBinContent(bin, bkgHist_sys.GetBinContent(bin) *(ratio_hist.GetBinContent(ratio_hist.GetXaxis().FindBin(bkgHist_sys.GetXaxis().GetBinCenter(bin)))))
#                else:
#                    bkgHist_sys.SetBinContent(bin, bkgHist_sys.GetBinContent(bin) *(2-ratio_hist.GetBinContent(ratio_hist.GetXaxis().FindBin(bkgHist_sys.GetXaxis().GetBinCenter(bin)))))
#                file_ratio.Close()
#            setNameTitle(bkgHist_sys, 'dd_bkg_' + name_sys)
#            hists2write.append(bkgHist_sys)
            

  
            file_new.mkdir(prefix + '_' + fitCat + '_' + year)
            file_new.cd(prefix + '_' + fitCat + '_' + year)

            for hist2write in hists2write:
                histname = hist2write.GetName()

                for pkey, pvar in process_convention.items():
                    histname = histname.replace(pkey, pvar)

                for nkey, nvar in nuisance_convention.items():
                    histname = histname.replace(nkey, nvar)

                if hist2write.GetName()!=histname:
                    print 'rename:', hist2write.GetName(),'---->', histname
                    setNameTitle(hist2write, histname)

                hist2write.Write()

    file_new.Write()
    file_new.Close()



    ### shape comparison!
    print 'shape comparison!!'
    file = TFile(filename_new)

    for year in eras:

        for fitCat in categories:

            dirname = prefix + '_' + fitCat + '_' + year

            file.cd(dirname)
        
            listofprocs = [key.GetName() for key in gDirectory.GetListOfKeys()]
            
            print listofprocs

            for proc in listofprocs:
            
                if proc.find('Up')==-1: continue


                #            strs = proc.split('_')
                tmplist =  proc.split('_')
                procname = tmplist[0] + '_' + tmplist[1]
                sysname = ((proc.replace(procname+"_", "")).replace('Up', ''))

                _cent = dirname + '/' + procname
                _down = dirname + '/' + proc.replace('Up','Down')
                _up = dirname + '/' + proc

                print '-'*80
                print '\t', procname, _cent, _down, _up

#                bc_others tauhad_sb_2018/bc_others tauhad_sb_2018/bc_others_sfRecoDown tauhad_sb_2018/bc_others_sfRecoUp
#                jpsi_hc tauhad_sb_2018/jpsi_hc tauhad_sb_2018/jpsi_hc_br_jpsi_hc_over_muDown tauhad_sb_2018/jpsi_hc_br_jpsi_hc_over_muUp
#                jpsi_hc tauhad_sb_2018/jpsi_hc tauhad_sb_2018/jpsi_tau_br_jpsi_hc_over_muDown tauhad_sb_2018/jpsi_tau_br_jpsi_hc_over_muUp


                #            import pdb; pdb.set_trace()
                
                cent = copy.deepcopy(file.Get(_cent))
                down = copy.deepcopy(file.Get(_down))
                up = copy.deepcopy(file.Get(_up))
                

                hists = [cent, down, up]
                titles = [procname, proc, proc.replace('Up','Down')]
                
                print hists, titles
                applyHists(hists)
                
                _dirname = output + '/syscompare/' + vkey + '/' + year + '/' + fitCat
                ensureDir(_dirname)
                shutil.copyfile('index.php',_dirname+'/index.php')
                comparisonPlots_alt(hists, titles, False, False, _dirname + '/' + vkey + '_' + procname + '_' + sysname + '.pdf', True, True, 'hpe')


        




