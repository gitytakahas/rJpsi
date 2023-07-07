from ROOT import TFile, TH2D, TH1D, gROOT, gStyle, TCanvas, TPaveText, TLatex, TH2F, TH1F, TF1, TGraphErrors, kRed, TVirtualFitter
from array import array
import copy

from common.officialStyle import officialStyle

gROOT.SetBatch(True)
#gROOT.SetBatch(False)
officialStyle(gStyle)
gStyle.SetOptTitle(0)


hist2save = []

for year in ['2016', '2017', '2018', 'inv_2016', 'inv_2017', 'inv_2018']:
#for year in ['2016', '2017', '2018']:

    file = TFile(year + '_inclusive_None/datacard/xgbs.root')
    bg = copy.deepcopy(file.Get('inclusive/bg_ul'))
    data = copy.deepcopy(file.Get('inclusive/data_obs'))
    bc_others = copy.deepcopy(file.Get('inclusive/bc_others'))
    bc_jpsi_dst = copy.deepcopy(file.Get('inclusive/bc_jpsi_dst'))
    bc_sig = copy.deepcopy(file.Get('inclusive/bc_jpsi_tau'))

    norm = float(data.GetSumOfWeights() - bc_others.GetSumOfWeights() - bc_jpsi_dst.GetSumOfWeights() - bc_sig.GetSumOfWeights())/float(bg.GetSumOfWeights())

    print year, norm
    
