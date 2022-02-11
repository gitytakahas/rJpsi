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

systs_mc = ['puweight_up', 'puweight_down', 'muSFID_up', 'muSFID_down', 'muSFReco_up', 'muSFReco_down','weight_ctau_up','weight_ctau_down']

datacardpath = 'datacard_MUSF_blind/'

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

    filename = datacardpath + '/' + channel + '/' + vkey + '.root'
        
    if sys!=None:
        filename = datacardpath + '/' + channel + '/' + vkey + '_' + sys + '.root'

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
            for proc in ['bc_jpsi_tau_3p', 'bc_jpsi_tau_N3p','bc_jpsi_dst','bc_others']: 
                    
            
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
        _dirname = 'plots_MUSF_blind/compare/' + '_'.join(channels) + '_' + target
        ensureDir(_dirname)

        comparisonPlots_alt(hists, titles, True, False, _dirname + '/' + vkey +  '.pdf', True, True, 'hpe')


    _dirname = datacardpath + '/compare/' + '_'.join(channels) + '_' + target
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

#fitCat = 'sr'

##Declare here the list of processes 
processes = ["data_obs","bc_jpsi_dst","bc_others", "bc_jpsi_tau_N3p", "bc_jpsi_tau_3p"]

categories = ['sr','sb']

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
    
#    draw(vkey, ['lp',  'sb'], 'data_obs', None, True, True)

#    draw(vkey, ['sb',  'sr'], 'data_obs', None, True, True)






#    hists4ddbkg_bgmc = draw(vkey, ['sb', 'sr'], 'bg_ul', None, False, True)

#    hists4ddbkg_vr = draw(vkey, ['lp', 'sb'], 'data_obs', None, True, True)
    


#    if vkey=='tau_rhomass_unrolled':
    if vkey in finaldiscriminant:


        datacards = {}

        for fitCat in categories:

            hists2write = []

            file = TFile(datacardpath + '/' + fitCat + '/' + vkey + '.root')
            file.cd(fitCat)
        
            listofprocs = [key.GetName() for key in gDirectory.GetListOfKeys()]
    
            for proc in listofprocs:
                _tmp = file.Get(fitCat + '/' + proc)
                hists2write.append(copy.deepcopy(_tmp))


            #### shape variation

            for sys in systs_hammer:

                name_sys = sys.replace('_up','Up').replace('_down', 'Down')

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
                

            for sys in systs_mc:

                name_sys = sys.replace('_up','Up').replace('_down', 'Down')

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
            

            datacards[fitCat] = copy.deepcopy(hists2write)


        filename_new = datacardpath + '/'+ vkey + '_new.root'
        file_new = TFile(filename_new, 'recreate')

        for fitCat, hists_ in datacards.iteritems():

            print fitCat, hists_

            file_new.cd()
            file_new.mkdir(fitCat)
            file_new.cd(fitCat)
        
            for hist_ in hists_:
                hist_.Write()
            
        file_new.Write()
        file_new.Close()



        ### shape comparison!

        file = TFile(filename_new)

        for fitCat in categories:

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

                _dirname = 'plots_MUSF_blind/syscompare_simultaneous/' + fitCat
                ensureDir(_dirname)
                comparisonPlots_alt(hists, titles, False, False, _dirname + '/' + vkey + '_' + procname + '_' + sysname + '.pdf', True, True, 'hpe')


        




