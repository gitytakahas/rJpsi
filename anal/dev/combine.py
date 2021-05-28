import copy, os, math, sys
from numpy import array
#from CMGTools.H2TauTau.proto.plotter.categories_TauMu import cat_Inc
from ROOT import TFile, TH1F, TH2F, TTree, gROOT, gStyle, Double, TCanvas, kRed, TH1D, TRandom3, kBird, THStack, kBlue
#from common.DisplayManager_stack import DisplayManager
from officialStyle import officialStyle
from array import array
import common.MultiDraw
import numpy as np
from varConfig import vardir
from DataMCPlot import *
#from DisplayManager_postfit import DisplayManager

lumi=59.6
gROOT.SetBatch(True)
officialStyle(gStyle)
gStyle.SetOptTitle(0)

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


#def comparisonPlots(hist, pname='sync.pdf', clabel='', isRatio=True):

#    display = DisplayManager(pname, isRatio, lumi, clabel, 0.42, 0.65)
#    display.Draw(hist)


from common.DisplayManager_compare import DisplayManager

def comparisonPlots(hists, titles, textitle, pname='sync.pdf', ratio=True):
 
    display = DisplayManager(pname, True)
    display.Draw(hists, titles, textitle)



processes = ['signal_ref', 'signal_shape', 'bg_bc', 'data_obs']

finaldiscriminant = ['q2_simple']

#systs = ['None', 'hammer_ebe_a0_up', 'hammer_ebe_a0_down', 'weight_ctau_up', 'weight_ctau_down']
systs = ['None']

for hammer in ['a0', 'a1', 'a2', 'b0', 'b1', 'b2', 'c1', 'c2', 'd0', 'd1', 'd2']:
    systs.append('hammer_ebe_' + hammer + '_up')
    systs.append('hammer_ebe_' + hammer + '_down')


ensureDir('plots/combine')
ensureDir('datacard/combine')

#rjpsi=1


def reject(proc, syst):
    if proc.find('signal_ref')!=-1: return True
    if syst.find('hammer')!=-1 and proc.find('signal')==-1: return True
    if syst.find('ctau')!=-1 and proc.find('signal')==-1: return True

    return False
    
    


for discriminant in finaldiscriminant:

    hists = {}

    for dir in ['sr', 'sb', 'cr1', 'cr2']:
    

        for syst in systs:

            filename = 'datacard/' + dir + '/' + discriminant + '.root'
            
            if syst!='None':
                filename = 'datacard/' + dir + '/' + discriminant + '_' + syst + '.root'

            file = TFile(filename)

            for proc in processes:

                hist = file.Get(dir + '/' + proc)

                hname = dir + '_' + proc + '_' + syst + '_' + discriminant
                hist.SetName(hname)
                hist.SetTitle(hname)
                hists[dir + '_' + proc + '_' + syst] = copy.deepcopy(hist)




    hists4datacard = {}

    for dir in ['sr', 'sb', 'cr1', 'cr2']:

        hists4datacard_ = []

        for syst in systs:

            for proc in processes:
                
                if reject(proc, syst): continue

                hist2add = copy.deepcopy(hists[dir + '_' + proc + '_' + syst])

                hname = proc

                if proc.find('signal_shape')!=-1:
                    hist2add_ref = copy.deepcopy(hists[dir + '_' + proc.replace('_shape', '_ref') + '_' + syst])
                    sf = hist2add_ref.GetSumOfWeights() / hist2add.GetSumOfWeights()
                    hist2add.Scale(sf)
                    
#                    print dir, syst, proc, 'signal_shape has been scaled to ', hist2add_ref.GetSumOfWeights(), hist2add.GetSumOfWeights(), sf

                if syst!='None':
                    hname = proc + '_' + syst.replace('_up','Up').replace('_down', 'Down')

                hist2add.SetName(hname.replace('_shape',''))
                hist2add.SetTitle(hname.replace('_shape',''))
                hist2add.GetXaxis().SetLabelColor(1)
                hist2add.GetXaxis().SetLabelSize(0.0)
                
                hists4datacard_.append(hist2add)
                
        # add BG histograms ... 
        for ibin in range(1, hist2add.GetXaxis().GetNbins()+1):

            dummy = TH1F('bg_bin' + str(ibin), 'bg_bin' + str(ibin),
                         hist2add.GetXaxis().GetNbins(), 
                         hist2add.GetXaxis().GetXmin(),
                         hist2add.GetXaxis().GetXmax())


            dummy.SetBinContent(ibin, 1)
        
            dummy.SetBinError(ibin, 0)
            
            hists4datacard_.append(copy.deepcopy(dummy))
            del dummy


        hists4datacard[dir] = hists4datacard_




    #write to combined datacards

    ofile = TFile('datacard/combine/' + discriminant + '.root', 'recreate')

    for cat, lhist in hists4datacard.items():

        ensureDir('syscompare/' + cat)

        ofile.mkdir(cat)
        ofile.cd(cat)        

        for ihist in lhist:
            ihist.Write()


        ############## making comparison plots ####################

        for syst in systs: # = ['None', 'hammer_ebe_a0_up', 'hammer_ebe_a0_down']


            if syst.find('up')==-1: continue


            hists = []
            titles = []

#            print syst
            
            for ihist in lhist:

#                print '\t ', ihist.GetName()
                
                hname = syst.replace('_up','Up')

                hist2shape = copy.deepcopy(ihist)
                hist2shape.SetFillStyle(0)
                hist2shape.SetMarkerSize(0)
#                hist2shape.SetLineWidth(3)
#                hist2shape.SetLineStyle(ii+1)
                
                
#                print ihist.GetName(), syst
#signal_hammer_ebe_a0Down
                if ihist.GetName() == 'signal':
                    hist2shape.SetLineWidth(3)
                    hists.append(hist2shape)
                    titles.append('default')

                if ihist.GetName().find(hname)!=-1:
                    hist2shape.SetLineColor(kBlue)
                    hist2shape.SetLineStyle(2)
                    hist2shape.SetLineWidth(2)
                    hists.append(hist2shape)
                    titles.append('Up')

                if ihist.GetName().find(hname.replace('Up', 'Down'))!=-1:
                    hist2shape.SetLineColor(kRed)
                    hist2shape.SetLineStyle(3)
                    hist2shape.SetLineWidth(2)
                    hists.append(hist2shape)
                    titles.append('Down')


#            print hists, titles

            comparisonPlots(hists, titles, cat + ', ' + syst.replace('_up',''), 'syscompare/' + cat + '/' + syst + '.gif')


    ofile.Write()
    ofile.Close()
