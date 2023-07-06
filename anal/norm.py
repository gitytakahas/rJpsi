from ROOT import TFile, TH2D, TH1D, gROOT, gStyle, TCanvas, TPaveText, TLatex, TH2F, TH1F, TF1, TGraphErrors, kRed, TVirtualFitter
from array import array
import copy

from common.officialStyle import officialStyle

gROOT.SetBatch(True)
#gROOT.SetBatch(False)
officialStyle(gStyle)
gStyle.SetOptTitle(0)

#gStyle.SetPadRightMargin (0.14)


#pnfs_prefix='/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/results_simultaneous/'
 
def add_lumi(year):
    lowX=0.7
    lowY=0.842
    lumi  = TPaveText(lowX, lowY+0.06, lowX+0.65, lowY+0.16, "NDC")
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.SetTextSize(0.055)
    lumi.SetTextFont (   42 )
#    lumi.AddText("2016, " + str(luminumber) + " fb^{-1} (13TeV)")
    lumi.AddText(year)
    return lumi


def add_CMS():
    lowX=0.15
    lowY=0.87
    lumi  = TPaveText(lowX, lowY, lowX+0.15, lowY+0.16, "NDC")
    lumi.SetTextFont(61)
    lumi.SetTextSize(0.05)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
    lumi.AddText("CMS")
    return lumi


def add_channel(cl):
    lowX=0.30
    lowY=0.81
#    lumi  = ROOT.TLatex(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
#    lumi  = ROOT.TLatex(lowX, lowY,  cl)
    lumi  = TLatex(lowX, lowY,  cl)
    lumi.SetNDC(True)
#    lumi.SetTextFont(61)
    lumi.SetTextSize(0.08)
#    lumi.SetBorderSize(   0 )
#    lumi.SetFillStyle(    0 )
#    lumi.SetTextAlign(   12 )
#    lumi.SetTextColor(    1 )
#    lumi.AddText(cl)
    return lumi


def add_Preliminary():
    lowX=0.26
    lowY=0.832
    lumi  = TPaveText(lowX, lowY+0.06, lowX+0.15, lowY+0.16, "NDC")
#    lumi.SetTextFont(52)
    lumi.SetTextSize(0.05)
    lumi.SetBorderSize(   0 )
    lumi.SetFillStyle(    0 )
    lumi.SetTextAlign(   12 )
    lumi.SetTextColor(    1 )
#    lumi.AddText("#it{Simulation} Preliminary")
    lumi.AddText("#it{Preliminary}")
    return lumi


hist2save = []

#for year in ['2016', '2017', '2018', 'inv_2016', 'inv_2017', 'inv_2018']:
for year in ['2016', '2017', '2018']:

    hists = []

    file = TFile(year + '_inclusive_None/datacard/xgbs.root')
    bg = copy.deepcopy(file.Get('inclusive/bg_ul'))
    data = copy.deepcopy(file.Get('inclusive/data_obs'))

    frame = TH2F('frame_'+ year, 'frame_' + year, 100, -10,7, 100,0.00001,0.1)
    frame.GetXaxis().SetTitle('BDT score')
    frame.GetYaxis().SetTitle('a.u.')
    frame.GetYaxis().SetNdivisions(504)


    canvas = TCanvas('can_'+year)
    frame.Draw()

    hists = [data, bg]

    hist2draw = []
    for ii, hist in enumerate(hists):

        print hist.GetName()

        hist_ = TH1F(hist.GetName(), hist.GetName(), hist.GetXaxis().GetNbins(), hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
        for ibin in range(1, hist.GetXaxis().GetNbins()+1):
            hist_.SetBinContent(ibin, float(hist.GetBinContent(ibin)/hist.GetSumOfWeights()))
            hist_.SetBinError(ibin, float(hist.GetBinError(ibin)/hist.GetSumOfWeights()))
        

        hist_.SetLineColor(ii+1)
        hist_.SetMarkerColor(ii+1)
        hist2draw.append(copy.deepcopy(hist_))


    for hist in hist2draw:
        hist.Draw('lepsame')

    l2=add_CMS()
    l2.Draw("same")
    l3=add_Preliminary()
    l3.Draw("same")
    l4=add_lumi(year)
    l4.Draw("same")

    canvas.SaveAs('plots/compare_' + year +'.pdf')
    canvas.SaveAs('plots/compare_' + year + '.gif')

    upper_boundary = 3.5
    lower_boundary = 0.

    if year.find('inv')!=-1:
        lower_boundary = -2.

    if year=='2016':
        upper_boundary = 3.

    rcanvas = TCanvas('can_ratio_' + year)

    miny = 0
    maxy = 3
    rframe = TH2F('ratio', 'ratio', 100, -10,7, 100,miny, maxy)
    rframe.GetXaxis().SetTitle('BDT score')
    rframe.GetYaxis().SetTitle('Data / Bkg. MC')
    rframe.GetYaxis().SetNdivisions(504)
    
    rframe.Draw()
    
#    ratio = copy.deepcopy(data)
    print 'CHECK ====>', hist2draw[0].GetName(), hist2draw[1].GetName()
    ratio = copy.deepcopy(hist2draw[0])
    ratio.SetTitle('ratio_' + year)
    ratio.SetName('ratio_' + year)
#    ratio.Divide(copy.deepcopy(bg))
    ratio.Divide(copy.deepcopy(hist2draw[1]))
    ratio.Draw('same')
    
    func = TF1('func', '[0]+[1]*x+[2]*x*x', lower_boundary, upper_boundary)
    ratio.Fit('func', '', '', lower_boundary, upper_boundary)
    
    fitline = ratio.GetFunction('func')
    fitline.SetLineColor(2)
    fitline.SetLineWidth(4)
    fitline.Draw('same')

    ratio.GetYaxis().SetRangeUser(miny, maxy)
    hist2save.append(copy.deepcopy(ratio))
    
    ### add 95% line 

#    graph = TGraphErrors(ratio.GetNbinsX())



#    result = ratio.Fit('func', 'RNS', '', lower_boundary, upper_boundary)

#    lower_index = -1 
#    lower_index_flag = False

#    for ii in range(1, ratio.GetNbinsX()+1):
#        graph.SetPoint(ii-1, ratio.GetBinCenter(ii), ratio.GetBinContent(ii))
        
#        if ratio.GetBinCenter(ii) > lower_index and lower_index_flag==False:
#            lower_index = ii-1
#            lower_index_flag=True

#    result = graph.Fit('func', 'RN S', '', lower_boundary, upper_boundary)
#    import pdb; pdb.set_trace()

#    result.GetConfidenceIntervals(graph)


#    interval = TGraphErrors(len(values))
#
#    print 'values=', len(values), graph.GetN()
#    
#    for i in range(len(values)):
#        interval.SetPoint(i, graph.GetX()[lower_index + i], func.Eval(graph.GetX()[lower_index+i] ))
#        interval.SetPointError(i, 0, values[i] )


    ###


#    result.GetConfidenceIntervals(0.95, False)
    #interval = TGraphErrors()

    #interval.SetFillColor(r.kRed-7)
    #interval.Draw("3same")

    

    p1 = fitline.GetParameter(0)
    p2 = fitline.GetParameter(1)
    p3 = fitline.GetParameter(2)
    
    fitted_all = TF1('fitted', str(p1) + '+ ' + str(p2) +' *x + ' + str(p3) + '*x*x', upper_boundary,7)
    fitted_all.SetLineStyle(2)
    fitted_all.SetLineColor(2)
    fitted_all.SetLineWidth(2)

    fitted_all.Draw('same')
    #result.Draw("3same")

    l2=add_CMS()
    l2.Draw("same")
    l3=add_Preliminary()
    l3.Draw("same")
    l4=add_lumi(year)
    l4.Draw("same")

    
#    result.SetFillColor(kRed-7)
#    graph.Draw("3same")
    
    rcanvas.SaveAs('plots/ratio_compare_' + year +'.pdf')
    rcanvas.SaveAs('plots/ratio_compare_' + year + '.gif')



ofile = TFile('plots/ratio.root', 'recreate')
for hist in hist2save:
    hist.Write()

ofile.Write()
ofile.Close()
