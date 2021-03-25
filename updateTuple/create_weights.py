#!/usr/bin/env python2

from ROOT import TH2F, TFile, gROOT, gStyle, TChain
#import numpy as np
#import matplotlib.pyplot as plt
#from root_numpy import tree2array, array2root
#from officialStyle import officialStyle
import argparse
import os
import copy

gROOT.SetBatch(True)
#officialStyle(gStyle)
gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)


parser = argparse.ArgumentParser()

parser.add_argument('--bkg_file', help='file to draw the background distribution from', required=True)
parser.add_argument('--sig_file', help='file where to weights are adapted too', required=True)
parser.add_argument('--out_dir', help='output directory', default="weights")

# parser.print_help()



args = parser.parse_args()

if not os.path.exists(args.out_dir):
    os.makedirs(args.out_dir)



def sproducer(name, cut, chain, weight='1'):

    hist = TH2F(name, 
                name, 
                20,0,2.5,20,0,60)

    hist.Sumw2()
    
    selection ='(' + cut + ')*' + weight
#    exp = '(' +  + ')'
        
#    tree = rootfile.Get('tree')

    chain.Draw('b_pt:abs(b_eta) >> ' + hist.GetName(), selection)
    hist.GetXaxis().SetTitle('eta')
    hist.GetYaxis().SetTitle('pt')


#    num = tree.GetEntries(selection + ' && ' + ivar['var'] + ' < 5.5')
        
    return copy.deepcopy(hist)


def returnChain(filelist):

    files = filelist.split(',')

    chain = TChain('tree', 'tree')

    for inputfile in files:
        print('Adding ...', inputfile)
        chain.AddFile(inputfile)
    
    return chain 


#root_file_bkg = TFile(args.bkg_file, 'READ')
#root_tree_bkg = root_file_bkg.Get("tree")

#root_file_sig = TFile(args.sig_file, 'READ')
#root_tree_sig = root_file_sig.Get("tree")

chain_sig = returnChain(args.sig_file)
chain_bg = returnChain(args.bkg_file)

print('# of sig =', chain_sig.GetEntries())
print('# of bkg =', chain_bg.GetEntries())

#sig = sproducer('sig', 'tau_index==0', chain_sig, 'puweight*hammer_ebe')
#bkg = sproducer('bkg', '1.', chain_bg, 'puweight')

sig = sproducer('sig', 'tau_isRight_3prong==1', chain_sig, 'puweight')
bkg = sproducer('bkg', '1.', chain_bg, 'puweight')

#weight = TH2F('weight', 'weight', 20,0, 2.5, 20, 2, 30)
weight = copy.deepcopy(sig)
weight.Divide(bkg)

weight.SetName('ratio')
weight.SetTitle('ratio')

#for ix in range(0, sig.GetXaxis().GetNbins()+2):
#    for iy in range(0, sig.GetYaxis().GetNbins()+2):
#        nsig = sig.GetBinContent(ix, iy)
#        nbkg = bkg.GetBinContent(ix, iy)
#
#        print('(x,y) = ', ix, iy, '(sig, bkg) = ', nsig, nbkg)
#
#        ratio = 1
#
#        if nbkg!=0:
#            ratio = nsig/nbkg
#        else:
#            print(ix, iy, '---> This is dangerous !!!')
#
#        weight.SetBinContent(ix, iy, ratio)


ofile = TFile(args.out_dir + '/weight.root', 'recreate')
weight.Write()
sig.Write()
bkg.Write()
ofile.Write()
ofile.Close()
