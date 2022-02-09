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


systs_hammer = []

for hammer in range(0, 9):
    systs_hammer.append('hammer_ebe_e' + str(hammer) + '_up')
    systs_hammer.append('hammer_ebe_e' + str(hammer) + '_down')

systs_pu = ['puweight_up', 'puweight_down']


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


def comparisonPlots_alt(hists, titles, norm=False, isLog=False, pname='sync.pdf', isRatio=False, isLegend=False, doption='he'):

    display = DisplayManager_compare(pname, isLog, isRatio, 0.4, 0.7, doption)
    display.draw_legend = isLegend
    
    display.Draw(hists, titles, norm)


def getHist(vkey, channel, target, sys=None):

    filename = 'datacard/' + channel + '/' + vkey + '.root'
        
    if sys!=None:
        filename = 'datacard/' + channel + '/' + vkey + '_' + sys + '.root'

    file = TFile(filename)

    _hist = copy.deepcopy(file.Get(channel + '/' + target))

#    file.Close()
    
    return _hist


def setNameTitle(hist, name):
    hist.SetTitle(name)
    hist.SetName(name)

def draw(vkey, channels, target, sys=None, subtract=False, saveFig=False, sf = 0.265):

    hists = []
    titles = []
    
    for ii, channel in enumerate(channels):

        _hist = getHist(vkey, channel, target, sys)

#        filename = 'datacard/' + channel + '/' + vkey + '.root'
#        
#        if sys!=None:
#            filename = 'datacard/' + channel + '/' + vkey + '_' + sys + '.root'
#
#        file = TFile(filename)
#
#        _hist = copy.deepcopy(file.Get(channel + '/' + target))

        if subtract:

#            ['bg_bc', 'sig_others', 'sig_3p', 'bg_ul', 'data_obs']

            #for proc in ['bg_bc', 'sig_others', 'sig_3p']:
            for proc in ['bc_jpsi_tau_3p', 'bc_jpsi_tau_N3p','bc_jpsi_ds','bc_others']: 
                    
            
#                if channel=='sr' and target=='data_obs' and proc=='sig_3p':
#                    print('THIS IS NOT CONSIDERED !!!!')
#                    continue

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

#        file.Close()



    hists2return = copy.deepcopy(hists)    


    if saveFig:
        _dirname = 'plots/compare/' + '_'.join(channels) + '_' + target
        ensureDir(_dirname)

        comparisonPlots_alt(hists, titles, True, False, _dirname + '/' + vkey +  '.pdf', True, True, 'hpe')


    _dirname = 'datacard/compare/' + '_'.join(channels) + '_' + target
    ensureDir(_dirname)

    file_output = TFile(_dirname + '/' + vkey +  '.root', 'recreate')
    for _hist in hists2return:
        _hist.Write()
    file_output.Write()
    file_output.Close()

        

    return hists2return


finaldiscriminant = ['tau_rhomass_unrolled', 'tau_rhomass_unrolled_coarse']
#finaldiscriminant = ['q2_simple']


for vkey, ivar in vardir.items():
    if vkey not in finaldiscriminant: 
        vardir.pop(vkey)

print '-'*80

for vkey, ivar in vardir.items():
    print '->', vkey.ljust(20), ivar

print '-'*80


ratio = 0.27

fitCat = 'sr'
ensureDir('syscompare/' + fitCat)

##Declare here the list of processes 
processes = ["data_obs","dd_bkg","bc_jpsi_ds","bc_others", "bc_jpsi_tau_N3p", "bc_jpsi_tau_3p"]

for vkey, ivar in vardir.items():

    # first compare the mc in sr and sb ... 
#    ratio1 = draw(vkey, ['hp_sb', 'hp_sr'], 'bg_ul')
#
#    ratio2 = draw(vkey, ['lp_sb', 'lp_sr'], 'data_obs', True)
#
#    ratio3 = draw(vkey, ['slp_sb', 'slp_sr'], 'data_obs', True)
#
#    ratio4 = draw(vkey, ['cr_hp_sb', 'cr_hp_sr'], 'data_obs', True)
#
#    ratio5 = draw(vkey, ['cr_lp_sb', 'cr_lp_sr'], 'data_obs', True)
#
#    ratio6 = draw(vkey, ['cr_slp_sb', 'cr_slp_sr'], 'data_obs', True)
#
#    ratio7 = draw(vkey, ['lp_sb', 'slp_sb', 'hp_sb'], 'data_obs', True)
#
#    ratio8 = draw(vkey, ['cr_lp_sb', 'cr_slp_sb'], 'data_obs', True)
#
#    ratio9 = draw(vkey, ['cr_lp_sr', 'cr_slp_sr'], 'data_obs', True)
#
#    ratio10 = draw(vkey, ['lp_sb', 'cr_lp_sb'], 'data_obs', True)
#
#    ratio11 = draw(vkey, ['lp_sr', 'cr_lp_sr'], 'data_obs', True)
#
#    ratio12 = draw(vkey, ['hp_sb'], 'data_obs', True, ratio2)

#    hists4ratio = draw(vkey, ['sr', 'sb'], 'bg_ul')
    
    draw(vkey, ['cr_lp',  'cr_sb'], 'data_obs', None, True, True)

    draw(vkey, ['cr_sb',  'cr_sr'], 'data_obs', None, True, True)






    hists4ddbkg_bgmc = draw(vkey, ['sb', 'sr'], 'bg_ul', None, False, True)

    hists4ddbkg_vr = draw(vkey, ['lp', 'sb'], 'data_obs', None, True, True)
    

#    if vkey=='tau_rhomass_unrolled':
    if vkey in finaldiscriminant:

        hists2write = []


        file = TFile('datacard/' + fitCat + '/' + vkey + '.root')
        file.cd(fitCat)
        
        listofprocs = [key.GetName() for key in gDirectory.GetListOfKeys()]
    
        for proc in listofprocs:
            _tmp = file.Get(fitCat + '/' + proc)
            hists2write.append(copy.deepcopy(_tmp))

        hists4ddbkg_sr = draw(vkey, ['sr', 'sb'], 'data_obs', None, True, False)
        bkgHist = copy.deepcopy(hists4ddbkg_sr[1])
        bkgHist.Scale(ratio)
        setNameTitle(bkgHist, 'dd_bkg')
        hists2write.append(bkgHist)

        ### shape variation 

        hists4ddbkg_up = draw(vkey, ['sr', 'sb'], 'data_obs', None, True, False, 0.95)
        bkgHist_up = copy.deepcopy(hists4ddbkg_up[1])
        bkgHist_up.Scale(ratio)
        setNameTitle(bkgHist_up, 'dd_bkg_shapeUp')
        hists2write.append(bkgHist_up)

        hists4ddbkg_down = draw(vkey, ['sr', 'sb'], 'data_obs', None, True, False, 0.2)
        bkgHist_down = copy.deepcopy(hists4ddbkg_down[1])
        bkgHist_down.Scale(ratio)
        setNameTitle(bkgHist_down, 'dd_bkg_shapeDown')
        hists2write.append(bkgHist_down)



#        bkgHist_alt = copy.deepcopy(hists4ddbkg_vr[0])
#        bkgHist_alt.Scale(hists4ddbkg[1].GetSumOfWeights()/bkgHist_alt.GetSumOfWeights())
#        bkgHist_alt.Scale(ratio)
#        setNameTitle(bkgHist_alt, 'dd_bkg_shape')


        ### reference
#        sig_3p_cent = getHist(vkey, fitCat, 'sig_3p')
#        sig_others_cent = getHist(vkey, fitCat, 'sig_others')
#        dd_bkg_cent = copy.deepcopy(bkgHist)


        #### shape variation

        for sys in systs_hammer:

            name_sys = sys.replace('_up','Up').replace('_down', 'Down')

            hists_sys = draw(vkey, ['sr', 'sb'], 'data_obs', sys, True, False)

            bkgHist_sys = hists_sys[1]
            bkgHist_sys.Scale(ratio)

            setNameTitle(bkgHist_sys, 'dd_bkg_' + name_sys)
            hists2write.append(bkgHist_sys)

            sig_3p_sys = getHist(vkey, fitCat, 'sig_3p', sys)
            sig_others_sys = getHist(vkey, fitCat, 'sig_others', sys)
            bc_others = getHist(vkey, fitCat, 'bc_others', sys)
            bc_jpsi_dst = getHist(vkey, fitCat, 'bc_jpsi_dst', sys)

            setNameTitle(sig_3p_sys, 'sig_3p_' + name_sys)
            setNameTitle(sig_others_sys, 'sig_others_' + name_sys)
            setNameTitle(bc_others, 'bc_others_' + name_sys)
            setNameTitle(bc_jpsi_dst, 'bc_jpsi_dst_' + name_sys)
            
            hists2write.append(sig_3p_sys)
            hists2write.append(sig_others_sys)
            hists2write.append(bc_others)
            hists2write.append(bc_jpsi_dst)


        for sys in systs_pu:

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
            
        print ("systs are:", systs)            
#            comparisonPlots(hists, titles, cat + ', ' + syst.replace('_up',''), 'syscompare/' + cat + '/' + syst + '.gif')
            

        ### CREATE shape comparisons 

        filename_new = 'datacard/'+ fitCat + '/' + vkey + '_new.root'
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
            print ("tentative sysname : ",proc.lstrip(procname+"_"))
            sysname = ((proc.lstrip(procname+"_")).replace('Up', ''))

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

            _dirname = 'plots/syscompare'
            ensureDir(_dirname)
            comparisonPlots_alt(hists, titles, False, False, _dirname + '/' + vkey + '_' + procname + '_' + sysname + '.pdf', True, True, 'hpe')


        




