import copy, os,  math, sys, shutil 
from numpy import array
from common.officialStyle import officialStyle
from array import array
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
#gStyle.SetOptStat(0)


usage = "usage: python compare.py" 
parser = OptionParser(usage) 

parser.add_option("-y", "--year", default="2018", type="string", dest="year")
#parser.add_option('-c', '--create', action="store_true", default=False, dest='create')

(options, args) = parser.parse_args() 



systs_hammer = []

for hammer in range(0, 10):
    systs_hammer.append('hammer_ebe_e' + str(hammer) + '_up')
    systs_hammer.append('hammer_ebe_e' + str(hammer) + '_down')

systs_mc = []

others=['puweight', 'muSFID', 'muSFReco', 'weight_ctau', 'br_BcJpsiDst', 'tauBr', 'tauReco', 'xgbsEff', 'BcPt']

for other in others:
    for ud in ['up', 'down']:
        systs_mc.append(other + '_' + ud)



datacardpath = '/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/results'


finaldiscriminant = ['tau_rhomass_unrolled', 'tau_rhomass_unrolled_coarse']

for vkey, ivar in vardir.items():
    if vkey not in finaldiscriminant: 
        vardir.pop(vkey)

print '-'*80

for vkey, ivar in vardir.items():
    print '->', vkey.ljust(20), ivar

print '-'*80


ratio = 0.07765
fitCat = 'sr'
output='combine/' + options.year 


init = ConfigParser.SafeConfigParser()
init.read("./settings.ini")
lumi = init.get('common', 'lumi_' + options.year)




processes = ["data_obs","dd_bkg","bc_jpsi_dst","bc_others", "bc_jpsi_tau_N3p", "bc_jpsi_tau_3p"]



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


def comparisonPlots(hist, pname='sync.pdf', isLog=False, isRatio=True, clabel=''):

    display = DisplayManager(pname, isRatio, isLog, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)

def comparisonPlots_alt(hists, titles, norm=False, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he'):

    display = DisplayManager_compare(pname, isLog, isRatio, 0.4, 0.7, doption)
    display.draw_legend = isLegend
    
    display.Draw(hists, titles, norm)


def getHist(vkey, channel, target, sys='None'):

    filename = datacardpath + '/' + options.year + '_' + channel + '_' + sys + '/datacard/' + vkey + '.root'

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

def draw(vkey, channels, target, sys='None', subtract=False, saveFig=False, sf = 0.265):

    hists = []
    titles = []
    
    for ii, channel in enumerate(channels):

        _hist = getHist(vkey, channel, target, sys)

        if subtract:

            for proc in ['bc_jpsi_tau_3p', 'bc_jpsi_tau_N3p','bc_jpsi_dst','bc_others']: 

                _hist2 = getHist(vkey, channel, proc, sys) #copy.deepcopy(file.Get(channel + '/' + proc))

                if proc in ['bc_jpsi_tau_3p', 'bc_jpsi_tau_N3p']:
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
    shutil.copyfile('index.php',_dirname+'/index.php')

    if saveFig:
        comparisonPlots_alt(hists, titles, True, False, _dirname + '/' + vkey +  '.pdf', True, True, 'hpe')


    file_output = TFile(_dirname + '/' + vkey +  '.root', 'recreate')
    for _hist in hists2return:
        _hist.Write()
    file_output.Write()
    file_output.Close()

        

    return hists2return



for vkey, ivar in vardir.items():

   
    draw(vkey, ['lp',  'sb'], 'data_obs', 'None', True, True)

    draw(vkey, ['sb',  'sr'], 'data_obs', 'None', True, True)



    hists2write = []

    file = TFile(datacardpath + '/' + options.year + '_' + fitCat + '_None/datacard/' + vkey + '.root')
    file.cd(fitCat)
        
    listofprocs = [key.GetName() for key in gDirectory.GetListOfKeys()]

    print listofprocs
    
    for proc in listofprocs:
       _tmp = file.Get(fitCat + '/' + proc)
       hists2write.append(copy.deepcopy(_tmp))
        
    hists4ddbkg_sr = draw(vkey, ['sr', 'sb'], 'data_obs', 'None', True, False)
    bkgHist = copy.deepcopy(hists4ddbkg_sr[1])
    bkgHist.Scale(ratio)
    setNameTitle(bkgHist, 'dd_bkg')
    hists2write.append(bkgHist)

    ### draw histograms 

    hists2draw = copy.deepcopy(hists2write)

    Histo = DataMCPlot(vkey)

    for _hist in hists2draw:

        _name = _hist.GetName()

        if _name=='bg_ul': continue

        if _name.find('data')!=-1:
            _hist.SetFillStyle(0)
            _hist.Sumw2(False)
            _hist.SetBinErrorOption(1)
            

        _hist.SetName(_name)
        _hist.SetTitle(_name)
        _hist.GetXaxis().SetLabelColor(1)
        _hist.GetXaxis().SetLabelSize(0.0)

        Histo.AddHistogram(_hist.GetName(), _hist)

        if _name.find('data')!=-1:
            Histo.Hist(_hist.GetName()).stack = False


    Histo._ApplyPrefs()
    print(Histo)
    comparisonPlots(Histo, output + '/' + vkey + '.gif')
    comparisonPlots(Histo, output + '/' + vkey + '_log.gif', True)


    ### shape variation 

    hists4ddbkg_up = draw(vkey, ['sr', 'sb'], 'data_obs', 'None', True, False, 0.95)
    bkgHist_up = copy.deepcopy(hists4ddbkg_up[1])
    bkgHist_up.Scale(ratio)
    setNameTitle(bkgHist_up, 'dd_bkg_shapeUp')
    hists2write.append(bkgHist_up)
    
    hists4ddbkg_down = draw(vkey, ['sr', 'sb'], 'data_obs', 'None', True, False, 0.2)
    bkgHist_down = copy.deepcopy(hists4ddbkg_down[1])
    bkgHist_down.Scale(ratio)
    setNameTitle(bkgHist_down, 'dd_bkg_shapeDown')
    hists2write.append(bkgHist_down)


    for sys in systs_hammer:

        name_sys = sys.replace('_up','Up').replace('_down', 'Down')
            
        hists_sys = draw(vkey, ['sr', 'sb'], 'data_obs', sys, True, False)

        bkgHist_sys = hists_sys[1]
        bkgHist_sys.Scale(ratio)
        setNameTitle(bkgHist_sys, 'dd_bkg_' + name_sys)
        hists2write.append(bkgHist_sys)
        
        sig_3p_sys = getHist(vkey, fitCat, 'bc_jpsi_tau_3p', sys)
        sig_others_sys = getHist(vkey, fitCat, 'bc_jpsi_tau_N3p', sys)
        
        setNameTitle(sig_3p_sys, 'bc_jpsi_tau_3p_' + name_sys)
        setNameTitle(sig_others_sys, 'bc_jpsi_tau_N3p_' + name_sys)
        
        hists2write.append(sig_3p_sys)
        hists2write.append(sig_others_sys)


    for sys in systs_mc:

        name_sys = sys.replace('_up','Up').replace('_down', 'Down')

        hists_sys = draw(vkey, ['sr', 'sb'], 'data_obs', sys, True, False)

        bkgHist_sys = hists_sys[1]
        bkgHist_sys.Scale(ratio)

        setNameTitle(bkgHist_sys, 'dd_bkg_' + name_sys)
        hists2write.append(bkgHist_sys)

                        
        sig_3p_sys = getHist(vkey, fitCat, 'bc_jpsi_tau_3p', sys)
        sig_others_sys = getHist(vkey, fitCat, 'bc_jpsi_tau_N3p', sys)
        bc_others = getHist(vkey, fitCat, 'bc_others', sys)
        bc_jpsi_dst = getHist(vkey, fitCat, 'bc_jpsi_dst', sys)
        
        setNameTitle(sig_3p_sys, 'bc_jpsi_tau_3p_' + name_sys)
        setNameTitle(sig_others_sys, 'bc_jpsi_tau_N3p_' + name_sys)
        setNameTitle(bc_others, 'bc_others_' + name_sys)
        setNameTitle(bc_jpsi_dst, 'bc_jpsi_dst_' + name_sys)     

            
        hists2write.append(sig_3p_sys)
        hists2write.append(sig_others_sys)
        hists2write.append(bc_others)
        hists2write.append(bc_jpsi_dst)
        

    filename_new = output + '/' + vkey + '.root'
    file_new = TFile(filename_new, 'recreate')
    file_new.mkdir(fitCat)
    file_new.cd(fitCat)


    for hist2write in hists2write:
        hist2write.Write()

    file_new.Write()
    file_new.Close()



    ### shape comparison!
    file = TFile(filename_new)
    file.cd(fitCat)
        
    listofprocs = [key.GetName() for key in gDirectory.GetListOfKeys()]

    print ("the list of processes is ", listofprocs)

    for proc in listofprocs:

        if proc.find('Up')==-1: continue

        strs = proc.split('_')
        print ("Printing proc", proc) 
        procname = [ processes[x] for x in range(len(processes)) if processes[x] in proc][0]
        print ("procname : ", procname)
        print ("tentative sysname : ", proc.replace(procname+"_", ""))
        sysname = ((proc.replace(procname+"_", "")).replace('Up', ''))

        print(procname, sysname)
        
        _cent = fitCat + '/' + procname
        _down = fitCat + '/' + proc.replace('Up','Down')
        _up = fitCat + '/' + proc

        print('comparing,', _cent, _down, _up)

        cent = copy.deepcopy(file.Get(_cent))
        down = copy.deepcopy(file.Get(_down))
        up = copy.deepcopy(file.Get(_up))
            
        hists = [cent, down, up]
        titles = [procname, proc, proc.replace('Up','Down')]

        applyHists(hists)

        _dirname = output + '/syscompare/' + vkey
        ensureDir(_dirname)
        shutil.copyfile('index.php',_dirname+'/index.php')
        comparisonPlots_alt(hists, titles, False, False, _dirname + '/' + vkey + '_' + procname + '_' + sysname + '.pdf', True, True, 'hpe')


        




