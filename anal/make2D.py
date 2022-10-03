from ROOT import TFile, TH2D, TH1D
from array import array
import copy

binning = [0.2, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 1.5]
nbins_var = int((len(binning)-1)*(len(binning)-1))

print len(binning)

filenames = [
    ('data', '/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_2018/Data/data.root'),
    ('sig','/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_2018/BcJpsiTau_inclusive/sig.root'),
    ('bkg', '/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_2018/BJpsiX/bkg.root')
]

hists = []

for name, filename in filenames:
    
    file = TFile(filename)
    tree = file.Get('tree')

    hist = TH2D(name, name, len(binning)-1, array('d',binning), len(binning)-1, array('d', binning))

    tree.Draw("tau_rhomass1:tau_rhomass2 >> " + hist.GetName(), "tau_pt > 3. && mu1_isLoose==1 && mu2_isLoose==1 && xgbs > 4.")
    
    hists.append(copy.deepcopy(hist))

    hname_1d = name + '_1d'

    hist_unrolled = TH1D(hname_1d, hname_1d, nbins_var, 0, nbins_var)

    idx_var = 1
    for ix in range(1, hist.GetXaxis().GetNbins()+1):
        for iy in range(1, hist.GetYaxis().GetNbins()+1):
            
            hist_unrolled.SetBinContent(idx_var, hist.GetBinContent(ix, iy))
            hist_unrolled.SetBinError(idx_var, hist.GetBinError(ix, iy))
                       
            idx_var += 1

    hists.append(copy.deepcopy(hist_unrolled))




print hists

ofile=TFile('out.root', 'recreate')

for hist in hists:
    hist.Write()




ofile.Write()
ofile.Close()


