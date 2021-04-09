import os, numpy, math, copy, math
from array import array
from ROOT import gStyle, TCanvas, TLegend, Double
from officialStyle import officialStyle
from DisplayManager_postfit import DisplayManager
from DataMCPlot import *

from optparse import OptionParser, OptionValueError

usage = "usage: python draw.py [channel]"
parser = OptionParser(usage)


parser.add_option(
    "-c", "--channel", 
    default="tautau",
    type="string",
    help="channel",
    dest="channel"
    )

parser.add_option(
    "-p", "--pp", 
    default="postfit",
    type="string",
    help="prefit_or_postfit",
    dest="pp"
    )

parser.add_option(
    "-n", "--nodata", 
    default=False,
    action="store_true",
    help="nodata",
    dest="nodata"
    )


parser.add_option(
    "-f", "--freeze", 
    default=None,
    type="string",
    help="freeze",
    dest="freeze"
    )


(options, args) = parser.parse_args()


gROOT.SetBatch(True)
#gROOT.SetBatch(False)
officialStyle(gStyle)
gStyle.SetOptTitle(0)

lumi=35.9

labeldir = {
    'etau':'e#tau_{h}',
    'mutau':'#mu#tau_{h}',
    'tautau':'#tau_{h}#tau_{h}',
    'emu':'e#mu',
    }

def display(proc, numlist, errlist):

    total = 0.
    err = 0.

    for ii in numlist:
        total += ii

    for ii in errlist:
        err += ii*ii

#    print proc, ' = ', '{0:.3f}'.format(total), '$\\pm$', '{0:.3f}'.format(math.sqrt(err))
    print proc, ' = ', '{0:.1f}'.format(total), '$\\pm$', '{0:.1f}'.format(math.sqrt(err))
    return 1



def comparisonPlots(hist, pname='sync.pdf', clabel=labeldir['mutau'], isRatio=True):

    display = DisplayManager(pname, isRatio, lumi, clabel, 0.42, 0.65)
    display.Draw(hist)



isTail = True
#isTail = False
isTail_start = 500.


channel = None

if options.channel=='mutau':
    channel = 'mt'

elif options.channel=='etau':
    channel = 'et'

elif options.channel=='emu':
    channel = 'em'

elif options.channel=='tautau':
    channel = 'tt'

else:
    print 'Invalid channel name', options.channel


print 'channel=', channel


filename = 'output/sm_cards/LIMITS/postfit_1000_shapes.root'

if options.freeze!=None:
    filename = 'output/sm_cards/LIMITS/postfit_1000_' + options.freeze  + '_shapes.root'

#file = TFile('output/sm_cards/LIMITS/postfit_700_shapes.root')
filename = 'output/sm_cards/LIMITS/postfit_1000_tautau_signal_nominal_btag_os_tight_shapes.root'
file = TFile(filename)
#file = TFile('output/sm_cards/LIMITS/postfit_200_shapes.root')
#file = TFile('output/sm_cards/LIMITS/postfit_200_mutau_signal_os_nominal_btag_shapes.root')
#leptoquark_mt_900_LQ.input.root

process = {
    'QCD':{'name':'QCD', 'isSignal':0, 'order':3},
    'TTT':{'name':'TTT', 'isSignal':0, 'order':4},
    'TTJ':{'name':'TTJ', 'isSignal':0, 'order':2},
    'VV':{'name':'VV', 'isSignal':0, 'order':5},
    'W':{'name':'W', 'isSignal':0, 'order':6},
    'STT':{'name':'STT', 'isSignal':0, 'order':9},
    'STJ':{'name':'STJ', 'isSignal':0, 'order':10},
    'ZJ':{'name':'ZJ', 'isSignal':0, 'order':7},
    'ZL':{'name':'ZL', 'isSignal':0, 'order':8},
    'ZTT':{'name':'ZTT', 'isSignal':0, 'order':1},
    'data_obs':{'name':'data_obs', 'isSignal':0, 'order':2999},
    'TotalSig':{'name':'TotalSig', 'isSignal':1, 'order':3001},
#    'Signal_M800':{'name':'Signal_M800', 'isSignal':0, 'order':3001},
}


yields = {}
yields_err = {}


hist = DataMCPlot('h_mass_' + channel)        
hist.legendBorders = 0.5, 0.45, 0.88, 0.85
    
for ii, val in process.iteritems():

    if ii=='TotalSig':
        _h_ = file.Get(channel + '_prefit/' + val['name'])
    else:
        _h_ = file.Get(channel + '_' + options.pp + '/' + val['name'])

    if _h_ == None:
        print ii, 'does not exist!'
        continue

    if channel == 'tt' and ii in ['W', 'ZJ']: continue

    if ii=='data_obs':
#        print 'Poisson error !!!'
        _h_.Sumw2(False)
        _h_.SetBinErrorOption(1)

        if options.nodata and options.channel!="emu":
            print 'No data is selected ... scale 0'
            _h_.Scale(-1.)
        

#        err_low = _h_.GetBinErrorLow(_h_.GetXaxis().GetNbins());
#        err_up = _h_.GetBinErrorUp(_h_.GetXaxis().GetNbins());
        
#        print err_low, err_up


    

#https://root.cern.ch/doc/master/TH1_8h_source.html#l00061
#    // enumeration specifying type of statistics for bin errors
#   61    enum  EBinErrorOpt {
#   62          kNormal = 0,    ///< errors with Normal (Wald) approximation: errorUp=errorLow= sqrt(N)
#   63          kPoisson = 1 ,  ///< errors from Poisson interval at 68.3% (1 sigma)
#   64          kPoisson2 = 2   ///< errors from Poisson interval at 95% CL (~ 2 sigma)
#   65    };


    _h_.SetName(val['name'])
    _h_.GetXaxis().SetLabelColor(1)
    _h_.GetXaxis().SetLabelSize(0.0)
    _h_.GetXaxis().SetTitle('S_{T} (GeV)')
    
    hist.AddHistogram(_h_.GetName(), _h_, val['order'])
    
    if val['name'] == 'data_obs' or val['name'].find('TotalSig')!=-1:
        hist.Hist(_h_.GetName()).stack = False


    start_bin = 0
    
    if isTail:
        start_bin = _h_.FindBin(isTail_start)


    err_ = Double(-1.)
    _h_.IntegralAndError(start_bin, _h_.GetXaxis().GetNbins()+1, err_)
    yields[val['name']] = _h_.Integral(start_bin, _h_.GetXaxis().GetNbins()+1)
    yields_err[val['name']] = err_


        

#hist.Group('Electroweak', ['W', 'ZL', 'ZJ', 'ZTT', 'VV'])
#hist.Group('ZTT_', ['ZL', 'ZJ', 'ZTT'])

if channel == 'tt':
    hist.Group('QCD_', ['QCD'])
#    hist.Group('ZTT_', ['ZL', 'ZTT'])
    hist.Group('EWK', ['VV', 'ZL', 'ZTT'])
else:
#    hist.Group('ZTT_', ['ZL', 'ZJ','ZTT'])
    hist.Group('EWK', ['VV', 'W', 'ZL', 'ZJ', 'ZTT'])

hist.Group('TT', ['TTJ', 'TTT'])
hist.Group('ST', ['STT', 'STJ'])


## for HepData

#import pdb; pdb.set_trace()
hist.DrawStack('HIST', None, None, None, None, 50)  
_total_ = hist.returnTotal()
for ibin in range(1, _total_.weighted.GetXaxis().GetNbins()+1):
    print 'XX ', channel,  '[(', int(hist.Hist('data_obs').obj.GetBinLowEdge(ibin)), ', ', int(hist.Hist('data_obs').obj.GetBinLowEdge(ibin) + hist.Hist('data_obs').obj.GetBinWidth(ibin)), '), ', int(hist.Hist('data_obs').obj.GetBinContent(ibin)), ', ({0:.2f}'.format(_total_.weighted.GetBinContent(ibin)), ', ', '{0:.2f}'.format(_total_.weighted.GetBinError(ibin)), ')]'


## for HepData end.

canvas = TCanvas()
hist.DrawStack('HIST', None, None, None, None, 2)

print hist

if options.freeze==None:
    if options.nodata:
        comparisonPlots(hist, 'plots/' + options.channel + '_' + options.pp  + '_nodata.gif', labeldir[options.channel])
    else:
        comparisonPlots(hist, 'plots/' + options.channel + '_' + options.pp  + '.gif', labeldir[options.channel])


if channel == 'tt':
    print display('TT', [yields['TTT'], yields['TTJ']], [yields_err['TTT'], yields_err['TTJ']])
    print display('single top', [yields['STT'], yields['STJ']], [yields_err['STT'], yields_err['STJ']])
#    print display('diboson', [], [])
    print display('EWK', [yields['VV'], yields['ZL'], yields['ZTT']], [yields_err['VV'], yields_err['ZL'], yields_err['ZTT']])
    print display('QCD', [yields['QCD']], [yields_err['QCD']])
    print display('Total', [yields['TTT'], yields['TTJ'], yields['VV'], yields['STT'], yields['STJ'], yields['ZL'], yields['ZTT'], yields['QCD']], [yields_err['TTT'], yields_err['TTJ'], yields_err['VV'], yields_err['STT'], yields_err['STJ'], yields_err['ZL'], yields_err['ZTT'], yields_err['QCD']])
#print display('Signal(400)', [yields['Signal_M400']], [yields_err['Signal_M400']])
    print display('Signal(700)', [yields['TotalSig']], [yields_err['TotalSig']])

    print display('data_obs', [yields['data_obs']], [yields_err['data_obs']])

elif channel =='em':
    print display('TT', [yields['TTT'], yields['TTJ']], [yields_err['TTT'], yields_err['TTJ']])
    print display('single top', [yields['STT'], yields['STJ']], [yields_err['STT'], yields_err['STJ']])
#    print display('diboson', [yields['VV']], [yields_err['VV']])
    print display('EWK', [yields['VV'], yields['W'], yields['ZL']], [yields_err['VV'], yields_err['W'],  yields_err['ZL']])
    print display('QCD', [yields['QCD']], [yields_err['QCD']])
    print display('Total', [yields['TTT'], yields['TTJ'], yields['VV'], yields['STT'], yields['STJ'], yields['W'], yields['ZL'], yields['QCD']], [yields_err['TTT'], yields_err['TTJ'], yields_err['VV'], yields_err['STT'], yields_err['STJ'], yields_err['W'], yields_err['ZL'], yields_err['QCD']])
#print display('Signal(400)', [yields['Signal_M400']], [yields_err['Signal_M400']])
    print display('Signal(700)', [yields['TotalSig']], [yields_err['TotalSig']])

    print display('data_obs', [yields['data_obs']], [yields_err['data_obs']])


else:
    print display('TT', [yields['TTT'], yields['TTJ']], [yields_err['TTT'], yields_err['TTJ']])
    print display('single top', [yields['STT'], yields['STJ']], [yields_err['STT'], yields_err['STJ']])
#    print display('diboson', [yields['VV']], [yields_err['VV']])
    print display('EWK', [yields['VV'], yields['W'], yields['ZJ'], yields['ZL'], yields['ZTT']], [yields_err['VV'], yields_err['W'], yields_err['ZJ'], yields_err['ZL'], yields_err['ZTT']])
    print display('QCD', [yields['QCD']], [yields_err['QCD']])
    print display('Total', [yields['TTT'], yields['TTJ'], yields['VV'], yields['STT'], yields['STJ'], yields['W'], yields['ZJ'], yields['ZL'], yields['ZTT'], yields['QCD']], [yields_err['TTT'], yields_err['TTJ'], yields_err['VV'], yields_err['STT'], yields_err['STJ'], yields_err['W'], yields_err['ZJ'], yields_err['ZL'], yields_err['ZTT'], yields_err['QCD']])
#print display('Signal(400)', [yields['Signal_M400']], [yields_err['Signal_M400']])
    print display('Signal(700)', [yields['TotalSig']], [yields_err['TotalSig']])

    print display('data_obs', [yields['data_obs']], [yields_err['data_obs']])
