#import os, math, sys
#from ROOT import TFile, TH1F, gROOT, TTree, Double, TChain, TLorentzVector, TVector3
#import numpy as num

from TreeProducerBcJpsiTauNu import *
from correction.PileupTool import *
from correction.ScaleFactorMuonTool import *
from DeltaR import deltaR
import copy
import random
from scipy.constants import c as speed_of_light
import numpy as np

# from https://github.com/scikit-hep/particle/
# installed via 
# >>> pip install particle --user
### 
from particle import Particle

#gROOT.SetBatch(True)

# Bc lifetime
# https://pdglive.lbl.gov/DataBlock.action?node=S091T&home=MXXX049

ctau_pdg    = 0.510e-12 * speed_of_light * 1000. # in mm
ctau_actual = 0.1358 #0.507e-12 * speed_of_light * 1000. # in mm          
ctau_up     = (0.510+0.009)*1e-12 * speed_of_light * 1000. # in mm
ctau_down   = (0.510-0.009)*1e-12 * speed_of_light * 1000. # in mm

def weight_to_new_ctau(old_ctau, new_ctau, ct):
    '''
    Returns an event weight based on the ratio of the normalised lifetime distributions.
    old_ctau: ctau used for the sample production
    new_ctau: target ctau
    ct      : per-event lifetime
    '''
    #print " weight_to_new_ctau:  old_ctau ", old_ctau , ", new_ctau ", new_ctau ," , ct ", ct
    weight = 1.
    weight = old_ctau/new_ctau * np.exp( (1./old_ctau - 1./new_ctau) * ct )
    return weight


# compute the distance between primary and secondary vtx
#sv = jpsi.vertex()
#pv = first_ancestor.vertex()


# compute the distance between primary and secondary vtx
#sv = jpsi.vertex()
#pv = first_ancestor.vertex()







def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def returnMass3(tlvs, m1, m2, m3):
    
    ms = [m1, m2, m3]

    tautlv = ROOT.TLorentzVector()
    
    for itlv, tlv in enumerate(tlvs):
        newtlv = ROOT.TLorentzVector()
        newtlv.SetPtEtaPhiM(tlv.Pt(), tlv.Eta(), tlv.Phi(), ms[itlv])
        tautlv += newtlv


#    print ms, ' ========> original = ', (tlvs[0] + tlvs[1] + tlvs[2]).M(), 'new = ', tautlv.M()

    return tautlv.M()

def returnMass2(tlvs, m1, m2):
    
    ms = [m1, m2]

    tautlv = ROOT.TLorentzVector()
    
    for itlv, tlv in enumerate(tlvs):
        newtlv = ROOT.TLorentzVector()
        newtlv.SetPtEtaPhiM(tlv.Pt(), tlv.Eta(), tlv.Phi(), ms[itlv])
        tautlv += newtlv


#    print ms, ' ========> original = ', (tlvs[0] + tlvs[1] + tlvs[2]).M(), 'new = ', tautlv.M()

    return tautlv.M()

def returnName(pid):
    
    addstr = None
    
    try:
        addstr = Particle.from_pdgid(pid).name    
    except:
        addstr = 'None'

    return addstr
        



from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)

parser.add_option("-o", "--out", default='Myroot.root', type="string", help="output filename", dest="out")
parser.add_option("-p", "--priority", default='pt', type="string", help="priority", dest="priority")
parser.add_option("-t", "--type", default='bkg', type="string", help="type", dest="type")
parser.add_option("-y", "--year", default='UL2017', type="string", help="year", dest="year")
parser.add_option("-f", "--file", default='root://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas//HbToJPsiMuMu_Legacy_2018_20220208/HbToJPsiMuMu_TuneCP5_13TeV-pythia8-evtgen/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/220208_162201/0000/flatTuple_9.root', type="string", help="file", dest="file")
parser.add_option('-c', '--create', action="store_true", default=False, dest='create')



(options, args) = parser.parse_args()

print(options)

if not options.create:
    json_open = open('json_jpsitanu.json', 'r')
    idtable = json.load(json_open)
    print('JSON decay table is read !')

    print idtable 


out = TreeProducerBcJpsiTauNu(options.out, options.type)

files = options.file.split(',')

chain = ROOT.TChain('ntuplizer/tree', 'tree')

hist_hammer = None
hist_hammer_lattice = None

#flag_fill = False

for inputfile in files:
    print('Adding ...', inputfile)
    chain.AddFile(inputfile)

    _file = TFile.Open(inputfile)
    _hist = copy.deepcopy(_file.Get('ntuplizer/cutflow'))

#    if out.hist == None:
#        print('cutflow histo created')
#        out.hist = copy.deepcopy(_hist)
##        out.hist = _hist
#    else:
#    print('cutflow histo added')
#    for ibin in range(1, _hist.GetXaxis().GetNbins()+1):
#        out.hist.GetXaxis().SetBinLabel(ibin, _hist.GetXaxis().GetBinLabel(ibin))
#        print('bin label', _hist.GetXaxis().GetBinLabel(ibin), out.hist.GetXaxis().GetBinLabel(ibin))


    out.hist.Add(_hist)
    
#    if options.type=='signal' and not flag_fill:


    if hist_hammer==None and options.type=='signal':
        hist_hammer = copy.deepcopy(_file.Get('ntuplizer/hammer_width'))

        for ibin in range(1, hist_hammer.GetXaxis().GetNbins()+1):

            out.hist_hammer.SetBinContent(ibin, hist_hammer.GetBinContent(ibin))
            out.hist_hammer.GetXaxis().SetBinLabel(ibin, hist_hammer.GetXaxis().GetBinLabel(ibin))


        
    if hist_hammer_lattice==None and options.type=='signal':
        hist_hammer_lattice = copy.deepcopy(_file.Get('ntuplizer/hammer_width_lattice'))

        for ibin in range(1, hist_hammer_lattice.GetXaxis().GetNbins()+1):

            out.hist_hammer_lattice.SetBinContent(ibin, hist_hammer_lattice.GetBinContent(ibin))
            out.hist_hammer_lattice.GetXaxis().SetBinLabel(ibin, hist_hammer_lattice.GetXaxis().GetBinLabel(ibin))

#        out.hist_hammer_lattice = hist_hammer_lattice


#    _tmp = _file.Get('ntuplizer/hammer_width')
#    out.hist_hammer = _tmp
#            
#    _tmp = _file.Get('ntuplizer/hammer_width_lattice')
#    out.hist_hammer_lattice = _tmp
#
#    print('This is it!', out.hist_hammer, out.hist_hammer_lattice)

#flag_fill = True


chain.SetBranchStatus('*', 0)


chain.SetBranchStatus('JpsiTau_B_pt*', 1)
chain.SetBranchStatus('*simple*', 1)
chain.SetBranchStatus('JpsiTau_B_eta', 1)
chain.SetBranchStatus('JpsiTau_B_phi', 1)
chain.SetBranchStatus('JpsiTau_B_mass', 1)
chain.SetBranchStatus('JpsiTau_B_mcorr', 1)
chain.SetBranchStatus('JpsiTau_B_alpha', 1)
chain.SetBranchStatus('JpsiTau_B_fls3d', 1)
chain.SetBranchStatus('JpsiTau_B_fl3d', 1)
chain.SetBranchStatus('JpsiTau_B_vprob', 1)

chain.SetBranchStatus('JpsiTau_tau_delta*', 1)
chain.SetBranchStatus('JpsiTau_tau_iso*', 1)
chain.SetBranchStatus('JpsiTau_tau_vweight', 1)
chain.SetBranchStatus('JpsiTau_tau_iso_ntracks*', 1)
chain.SetBranchStatus('JpsiTau_tau_iso_mindoca*', 1)
chain.SetBranchStatus('JpsiTau_B_q2', 1)
chain.SetBranchStatus('JpsiTau_B_mm2', 1)
chain.SetBranchStatus('JpsiTau_B_ptmiss', 1)
chain.SetBranchStatus('JpsiTau_B_Es', 1)

chain.SetBranchStatus('JpsiTau_B_pvip', 1)
chain.SetBranchStatus('JpsiTau_B_pvips', 1)
chain.SetBranchStatus('JpsiTau_B_lips', 1)
chain.SetBranchStatus('JpsiTau_B_mindoca', 1)
chain.SetBranchStatus('JpsiTau_B_vx', 1)
chain.SetBranchStatus('JpsiTau_B_vy', 1)
chain.SetBranchStatus('JpsiTau_B_vz', 1)
chain.SetBranchStatus('JpsiTau_nch', 1)
chain.SetBranchStatus('JpsiTau_nch_before', 1)

chain.SetBranchStatus('JpsiTau_PV_vx', 1)
chain.SetBranchStatus('JpsiTau_PV_vy', 1)
chain.SetBranchStatus('JpsiTau_PV_vz', 1)

chain.SetBranchStatus('JpsiTau_bbPV_vx', 1)
chain.SetBranchStatus('JpsiTau_bbPV_vy', 1)
chain.SetBranchStatus('JpsiTau_bbPV_vz', 1)


chain.SetBranchStatus('JpsiTau_tau_pt', 1)
chain.SetBranchStatus('JpsiTau_tau_eta', 1)
chain.SetBranchStatus('JpsiTau_tau_phi', 1)
chain.SetBranchStatus('JpsiTau_tau_mass', 1)
chain.SetBranchStatus('JpsiTau_tau_q', 1)
chain.SetBranchStatus('JpsiTau_tau_rhomass1', 1)
chain.SetBranchStatus('JpsiTau_tau_rhomass2', 1)
chain.SetBranchStatus('JpsiTau_tau_rhomass_ss', 1)
chain.SetBranchStatus('JpsiTau_tau_vprob', 1)
chain.SetBranchStatus('JpsiTau_tau_fls3d*', 1)
chain.SetBranchStatus('JpsiTau_tau_fl3d*', 1)
chain.SetBranchStatus('JpsiTau_tau_alpha*', 1)
chain.SetBranchStatus('JpsiTau_tau_sumofdnn*', 1)

chain.SetBranchStatus('JpsiTau_tau_max_dr_3prong', 1)
chain.SetBranchStatus('JpsiTau_tau_lips', 1)
chain.SetBranchStatus('JpsiTau_tau_pvips', 1)
chain.SetBranchStatus('JpsiTau_tau_refit*', 1)
chain.SetBranchStatus('JpsiTau_tau_pi*', 1)
chain.SetBranchStatus('JpsiTau_tau_v*', 1)
chain.SetBranchStatus('JpsiTau_Jpsi_*', 1)
chain.SetBranchStatus('JpsiTau_mu*', 1)
chain.SetBranchStatus('JpsiTau_ptbal', 1)
chain.SetBranchStatus('JpsiTau_jpsi_tau_alpha', 1)
chain.SetBranchStatus('JpsiTau_perEVT_*', 1)




chain.SetBranchStatus('JpsiTau_nCandidates', 1)
chain.SetBranchStatus('PV_N', 1)
chain.SetBranchStatus('EVENT_*', 1)

#putool = None
#putool_up = None
#putool_down = None

puhist = None
puhist_up = None
puhist_down = None

if options.type!='data':
    chain.SetBranchStatus('*genParticle*', 1)
    chain.SetBranchStatus('JpsiTau_*gen*', 1)
    chain.SetBranchStatus('JpsiTau_genPV_vx', 1)
    chain.SetBranchStatus('JpsiTau_genPV_vy', 1)
    chain.SetBranchStatus('JpsiTau_genPV_vz', 1)
    chain.SetBranchStatus('JpsiTau_genSV_vx', 1)
    chain.SetBranchStatus('JpsiTau_genSV_vy', 1)
    chain.SetBranchStatus('JpsiTau_genSV_vz', 1)
    chain.SetBranchStatus('nPuVtxTrue', 1)
    chain.SetBranchStatus('JpsiTau_st_gentau*',1)
    chain.SetBranchStatus('JpsiTau_st_n_occurance',1)
    chain.SetBranchStatus('JpsiTau_st_decayid',1)
    chain.SetBranchStatus('JpsiTau_isJpsi*', 1)

if options.type=='signal':
    chain.SetBranchStatus('truth_tau_dipion1_mass', 1)
    chain.SetBranchStatus('truth_tau_dipion2_mass', 1)
#    putool = PileupWeightTool(options.year, 'central')
#    putool_up = PileupWeightTool(options.year, 'up')
#    putool_down = PileupWeightTool(options.year, 'down')


if options.type in ['signal','bg']:
#    pufile = TFile('/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/RJPsi_mc_pu_2021Dec08_111.root')

    puname_data = None
    puname_mc = None

    if options.year == '2016':
        puname_data = 'PileupData_GoldenJSON_Full2016.root'
        puname_mc = 'pileup_profile_Summer16.root'
    elif options.year == '2017':
        puname_data = 'PileupHistogram-goldenJSON-13tev-2017-99bins_withVar.root'
        puname_mc = 'mcPileup2017.root'
    elif options.year == '2018':
        puname_data = 'PileupHistogram-goldenJSON-13tev-2018-100bins_withVar.root'
        puname_mc = 'mcPileup2018.root'

    pufile_data = TFile('nanoAOD-tools/python/postprocessing/data/pileup/' + puname_data)
    pufile_mc = TFile('nanoAOD-tools/python/postprocessing/data/pileup/' + puname_mc)

    puhist_data = pufile_data.Get('pileup')
    puhist_data_plus = pufile_data.Get('pileup_plus')
    puhist_data_minus = pufile_data.Get('pileup_minus')

    puhist_mc = pufile_mc.Get('pu_mc')


    if(puhist_data.GetNcells() != puhist_mc.GetNcells()):
        print("ERR: Numerator and Denominator have different number of bins !!!!")

    for _hist in [puhist_data, puhist_data_plus, puhist_data_minus, puhist_mc]:
        _hist.Scale(1./_hist.Integral())

    puhist = copy.deepcopy(puhist_mc)
    puhist_up = copy.deepcopy(puhist_mc)
    puhist_down = copy.deepcopy(puhist_mc)    

    #does the ratio of the histos
    for i in range(puhist_data.GetNcells()):

        val_mc = puhist_mc.GetBinContent(i)
        val_data = puhist_data.GetBinContent(i)
        val_data_plus = puhist_data_plus.GetBinContent(i)
        val_data_minus = puhist_data_minus.GetBinContent(i)

        weight = 1
        weight_up = 1
        weight_down = 1

        if val_mc!=0:
            weight = val_data/val_mc
            weight_up = val_data_plus/val_mc
            weight_down = val_data_minus/val_mc


#        print i, weight, weight_up, weight_down
        puhist.SetBinContent(i, weight)
        puhist_up.SetBinContent(i, weight_up)
        puhist_down.SetBinContent(i, weight_down)


#    puhist = pufile.Get('hweights')
#    puhist_up = pufile.Get('hweights_up')
#    puhist_down = pufile.Get('hweights_down')

    print(puhist, puhist_up, puhist_down, 'is made ...')

    yearstr = options.year
    if options.year == '2016':
        if options.file.find('preVFP')!=-1:
            yearstr = '2016preVFP'
        else:
            yearstr = '2016postVFP'
    print 'muon reco file', yearstr 

    SF_ID = ScaleFactorMuonTool('central', fileName='Efficiency_muon_trackerMuon_Run' + yearstr + '_UL_ID.json', keyName='NUM_LooseID_DEN_TrackerMuons');
    SF_Reco = ScaleFactorMuonTool('central', fileName='Efficiency_muon_generalTracks_Run' + yearstr + '_UL_trackerMuon.json', keyName='NUM_TrackerMuons_DEN_genTracks');
  
if options.type=='bg':
    chain.SetBranchStatus('genWeightBkgB',1)

if options.type=='signal':
    chain.SetBranchStatus('JpsiTau_hammer_*', 1)

Nevt = chain.GetEntries()

print('Total Number of events = ', Nevt)
evtid = 0


mk = 0.493677
mp = 0.139571
mproton = 0.938272

dicts = {}
dict_counter = 0

for evt in xrange(Nevt):
    chain.GetEntry(evt)

    if evt%10000==0: print('{0:.2f}'.format(ROOT.Double(evt)/ROOT.Double(Nevt)*100.), '% processed')

#    print(chain.JpsiTau_isJpsiMu, bool(chain.JpsiTau_isJpsiMu))
#    import pdb; pdb.set_trace()
#    if bool(chain.JpsiTau_isJpsiMu)==True:
#        print(type(chain.JpsiTau_isJpsiMu), bool(chain.JpsiTau_isJpsiMu))


    out.multi.Fill(len(chain.JpsiTau_tau_pt))
    out.filt.Fill(0)

    if len(chain.JpsiTau_tau_pt)==0: 
        print('This is not possible !!')
        continue

    out.filt.Fill(1)

    tlv_jpsi = ROOT.TLorentzVector()
    tlv_jpsi.SetPtEtaPhiM(chain.JpsiTau_Jpsi_pt, chain.JpsiTau_Jpsi_eta, chain.JpsiTau_Jpsi_phi, chain.JpsiTau_Jpsi_mass)
    
    
    tlv_mu1 = ROOT.TLorentzVector()
    tlv_mu1.SetPtEtaPhiM(chain.JpsiTau_mu1_pt, chain.JpsiTau_mu1_eta, chain.JpsiTau_mu1_phi, chain.JpsiTau_mu1_mass)

    tlv_mu2 = ROOT.TLorentzVector()
    tlv_mu2.SetPtEtaPhiM(chain.JpsiTau_mu2_pt, chain.JpsiTau_mu2_eta, chain.JpsiTau_mu2_phi, chain.JpsiTau_mu2_mass)


    tindex_ = -1

    ########################################################
    ## 
    ## Pre-selection here.
    ## 
    ## Please modify selections below, if you want
    ## 
    ## below is an example of filtering vertex prob. > 10% and the flgith sig. > 3 sigma
    ## 
    ########################################################

    ntaucand = 0
    
    if options.priority=='pt':
           
        for itau in range(len(chain.JpsiTau_tau_pt)):
            if chain.JpsiTau_tau_vprob[itau] < 0.1: continue
            if chain.JpsiTau_tau_fls3d[itau] < 3.: continue
#            if  chain.JpsiTau_tau_fls3d[itau] > 3. and chain.JpsiTau_tau_vprob[itau] > 0.1: continue
            if chain.JpsiTau_tau_mass[itau] > 1.7: continue
#            if abs(chain.JpsiTau_tau_q[itau])!=1: continue
#            if bool(chain.JpsiTau_tau_pi1_trigMatch[itau])==False and bool(chain.JpsiTau_tau_pi2_trigMatch[itau])==False and bool(chain.JpsiTau_tau_pi3_trigMatch[itau])==False: 
#                print 'trigger matching was not satisifed ...'
#                continue

            dr_pi1_mu1 = deltaR(chain.JpsiTau_tau_pi1_eta[itau], chain.JpsiTau_tau_pi1_phi[itau], chain.JpsiTau_mu1_eta, chain.JpsiTau_mu1_phi)
            dr_pi1_mu2 = deltaR(chain.JpsiTau_tau_pi1_eta[itau], chain.JpsiTau_tau_pi1_phi[itau], chain.JpsiTau_mu2_eta, chain.JpsiTau_mu2_phi)

            dr_pi2_mu1 = deltaR(chain.JpsiTau_tau_pi2_eta[itau], chain.JpsiTau_tau_pi2_phi[itau], chain.JpsiTau_mu1_eta, chain.JpsiTau_mu1_phi)
            dr_pi2_mu2 = deltaR(chain.JpsiTau_tau_pi2_eta[itau], chain.JpsiTau_tau_pi2_phi[itau], chain.JpsiTau_mu2_eta, chain.JpsiTau_mu2_phi)

            dr_pi3_mu1 = deltaR(chain.JpsiTau_tau_pi3_eta[itau], chain.JpsiTau_tau_pi3_phi[itau], chain.JpsiTau_mu1_eta, chain.JpsiTau_mu1_phi)
            dr_pi3_mu2 = deltaR(chain.JpsiTau_tau_pi3_eta[itau], chain.JpsiTau_tau_pi3_phi[itau], chain.JpsiTau_mu2_eta, chain.JpsiTau_mu2_phi)

#            print(dr_pi1_mu1, dr_pi1_mu2, dr_pi2_mu1, dr_pi2_mu2, dr_pi3_mu1, dr_pi3_mu2)


            if not (    (bool(chain.JpsiTau_tau_pi1_trigMatch[itau])==True and chain.JpsiTau_tau_pi1_pt[itau] > 1.2 and chain.JpsiTau_tau_pi1_d0sig[itau] > 2. and float(chain.JpsiTau_tau_pi1_mumuTrk_chi2[itau])/float(chain.JpsiTau_tau_pi1_mumuTrk_ndof[itau]) < 10. and dr_pi1_mu1 > 0.005 and dr_pi1_mu2 > 0.005) 
                        or
                        (bool(chain.JpsiTau_tau_pi2_trigMatch[itau])==True and chain.JpsiTau_tau_pi2_pt[itau] > 1.2 and chain.JpsiTau_tau_pi2_d0sig[itau] > 2. and float(chain.JpsiTau_tau_pi2_mumuTrk_chi2[itau])/float(chain.JpsiTau_tau_pi2_mumuTrk_ndof[itau]) < 10. and dr_pi2_mu1 > 0.005 and dr_pi2_mu2 > 0.005) 
                        or
                        (bool(chain.JpsiTau_tau_pi3_trigMatch[itau])==True and chain.JpsiTau_tau_pi3_pt[itau] > 1.2 and chain.JpsiTau_tau_pi3_d0sig[itau] > 2. and float(chain.JpsiTau_tau_pi3_mumuTrk_chi2[itau])/float(chain.JpsiTau_tau_pi3_mumuTrk_ndof[itau]) < 10. and dr_pi3_mu1 > 0.005 and dr_pi3_mu2 > 0.005) 
            ): 
                print 'trigger matching was not satisfied ...'
                continue

            ntaucand += 1

            tindex_ = itau
            break

#    if options.priority=='pt':
#        tindex_ = 0
            
#    print evt, 'tindex = ', tindex_ 

    out.multi_tau.Fill(ntaucand)

    if options.priority!='multiple' and tindex_ == -1: continue

    out.filt.Fill(2)

    if len(chain.JpsiTau_B_pt) != len(chain.JpsiTau_tau_pt): continue

    out.filt.Fill(3)

    ### ctau calculation ! 

    if options.type == 'signal':
        #print " len(chain.JpsiTau_genPV_vx) ", len(chain.JpsiTau_genPV_vx)
        lxyz = math.sqrt((chain.JpsiTau_genPV_vx - chain.JpsiTau_genSV_vx)**2 + (chain.JpsiTau_genPV_vy - chain.JpsiTau_genSV_vy)**2 + (chain.JpsiTau_genPV_vz - chain.JpsiTau_genSV_vz)**2)


        tlv_B = ROOT.TLorentzVector()
        ##        tlv_B.SetPtEtaPhiM(chain.JpsiTau_B_pt_gen, chain.JpsiTau_B_eta_gen, chain.JpsiTau_B_phi_gen, chain.JpsiTau_B_mass_gen)
        for igen in range(len(chain.genParticle_pdgId)):
            if abs(chain.genParticle_pdgId[igen])==541 and chain.genParticle_status[igen]==2:
                tlv_B.SetPtEtaPhiM(chain.genParticle_pt[igen], chain.genParticle_eta[igen], chain.genParticle_phi[igen], 6.274)
#        print "Bmass ",  chain.JpsiTau_B_mass_gen 
#        print chain.JpsiTau_B_pt_gen, chain.JpsiTau_B_eta_gen, chain.JpsiTau_B_phi_gen, chain.JpsiTau_B_mass_gen
        #    beta_e = tlv_B.E()
        #    beta_p = tlv_B.P()
        beta = tlv_B.Beta()
        gamma = tlv_B.Gamma()
        
    
        # lorentz boost of the B
        #beta = first_ancestor.p4().Beta()
        #gamma = first_ancestor.p4().Gamma()
        
        # now, lifetime L = beta * gamma * c * t ===> t = (L)/(beta*gamma*c)
    
        #        print lxyz, beta, gamma, weight_central, weight_up, weight_down 
        
        ct = lxyz / (beta * gamma)
        weight_central =1
        weight_up =1
        weight_down =1

        weight_central = weight_to_new_ctau(ctau_actual, ctau_pdg , ct*10.)
        weight_up = weight_to_new_ctau(ctau_actual, ctau_up , ct*10.)
        weight_down = weight_to_new_ctau(ctau_actual, ctau_down , ct*10.)

        out.lxyz[0] = lxyz
        out.beta[0] = beta
        out.gamma[0] = gamma
        out.ctau[0] = ct 
        out.weight_ctau[0] = weight_central
        out.weight_ctau_up[0] = weight_up
        out.weight_ctau_down[0] = weight_down
        if  ct>5:
            #print "B mass ",  chain.JpsiTau_B_mass_gen
            print "new_ctau:  ctau_actual ", ctau_actual, "ctau_pdg " , ctau_pdg ,  " , ct ", ct, ", up ", ctau_up ,", down ", ctau_down 
            print "exponential factor ", (1./ ctau_actual - 1./ctau_pdg)            
            print "ctau weights : ", weight_central, " ", weight_up, " ", weight_down 
            print lxyz, beta, gamma, ct 
            #print chain.JpsiTau_B_pt_gen, chain.JpsiTau_B_eta_gen, chain.JpsiTau_B_phi_gen, chain.JpsiTau_B_mass_gen  
            for igen in range(len(chain.genParticle_pdgId)):
                if abs(chain.genParticle_pdgId[igen])==541:
                    print "pdg status , pt, eta, phi", chain.genParticle_status[igen], chain.genParticle_pt[igen], chain.genParticle_eta[igen], chain.genParticle_phi[igen]
                    for idau in range(len(chain.genParticle_dau[igen])):
                        print '   -> ', (chain.genParticle_dau[igen][idau])           
            print " PV ", chain.JpsiTau_genPV_vx, chain.JpsiTau_genPV_vy, chain.JpsiTau_genPV_vz, " SV ",  chain.JpsiTau_genSV_vx, chain.JpsiTau_genSV_vy, chain.JpsiTau_genSV_vz,
                        #        import pdb; pdb.set_trace()
        out.isJpsiMu[0] = bool(chain.JpsiTau_isJpsiMu)
        out.isJpsiTau2Mu[0] = bool(chain.JpsiTau_isJpsiTau2Mu)
      





    tindex_multiple = 0

    for tindex in range(len(chain.JpsiTau_tau_pt)):


#        if abs(chain.JpsiTau_tau_q[tindex])!=1: continue
#        if abs(chain.JpsiTau_tau_q[tindex])!=1: continue

        if options.priority in ['pt']:
            if tindex != tindex_: continue


        out.filt.Fill(4)

        isRight_3prong = False
        isRight_3prong_pi0 = False

        if options.type != 'data':

            if bool(chain.JpsiTau_tau_pi1_isSignal[tindex]) and chain.JpsiTau_tau_pi1_nprong[tindex]==3 and bool(chain.JpsiTau_tau_pi2_isSignal[tindex]) and chain.JpsiTau_tau_pi2_nprong[tindex]==3 and bool(chain.JpsiTau_tau_pi3_isSignal[tindex]) and chain.JpsiTau_tau_pi3_nprong[tindex]==3:
                isRight_3prong = True

#                print chain.JpsiTau_tau_pi1_nprong_pi0[tindex], chain.JpsiTau_tau_pi2_nprong_pi0[tindex], chain.JpsiTau_tau_pi3_nprong_pi0[tindex], chain.JpsiTau_tau_pi1_nprong[tindex], chain.JpsiTau_tau_pi2_nprong[tindex], chain.JpsiTau_tau_pi3_nprong[tindex]


                if chain.JpsiTau_tau_pi1_nprong_pi0[tindex]>=1 and chain.JpsiTau_tau_pi2_nprong_pi0[tindex]>=1 and chain.JpsiTau_tau_pi3_nprong_pi0[tindex]>=1:
                    isRight_3prong_pi0 = True





        if options.priority=='multiple':

#            if chain.JpsiTau_tau_vprob[itau] < 0.1: continue
#            if chain.JpsiTau_tau_fls3d[itau] < 3.: continue
#            if chain.JpsiTau_tau_mass[itau] > 1.7: continue
            if bool(chain.JpsiTau_tau_pi1_trigMatch[tindex])==False and bool(chain.JpsiTau_tau_pi2_trigMatch[tindex])==False and bool(chain.JpsiTau_tau_pi3_trigMatch[tindex])==False: 
#                print 'trigger matching was not satisifed ...'
                continue

            tindex_multiple+=1

            if tindex_multiple >= 6: break
            # > 10 ---> 0.97316682 efficiency of catching the siganl
            # > 5 ---> 0.95029727 efficiency of catching the siganl


        if options.type in ['bg', 'signal']:
#            print bool(chain.JpsiTau_tau_pi1_isBdecay[tindex]), bool(chain.JpsiTau_tau_pi2_isBdecay[tindex]), bool(chain.JpsiTau_tau_pi3_isBdecay[tindex])
#            print chain.JpsiTau_tau_pi1_isBdecayppdg[tindex], chain.JpsiTau_tau_pi2_isBdecayppdg[tindex], chain.JpsiTau_tau_pi3_isBdecayppdg[tindex]
#            print chain.JpsiTau_tau_pi1_isBdecaypdg[tindex], chain.JpsiTau_tau_pi2_isBdecaypdg[tindex], chain.JpsiTau_tau_pi3_isBdecaypdg[tindex]
            parents = [abs(chain.JpsiTau_tau_pi1_isBdecayppdg[tindex]), abs(chain.JpsiTau_tau_pi2_isBdecayppdg[tindex]), abs(chain.JpsiTau_tau_pi3_isBdecayppdg[tindex])]
            daughters = [abs(chain.JpsiTau_tau_pi1_isBdecaypdg[tindex]), abs(chain.JpsiTau_tau_pi2_isBdecaypdg[tindex]), abs(chain.JpsiTau_tau_pi3_isBdecaypdg[tindex])]
#            print parents 
#            print daughters


            result = all(elem == parents[0] for elem in parents)

            pid = -1

            if result:
                pid = parents[0]

#            print '# of Kaon = ', daughters.count(321)
#            print '# of Pion = ', daughters.count(211)
#            print '# of PU = ', daughters.count(999)

            out.nkaon[0] = daughters.count(321)
            out.npion[0] = daughters.count(211)
            out.npu[0] = daughters.count(999)
            out.pid[0] = pid

            
#        print('test4')
        out.mu1_pt[0] = tlv_mu1.Pt()
        out.mu1_eta[0] = tlv_mu1.Eta()
        out.mu1_phi[0] = tlv_mu1.Phi()
        out.mu1_mass[0] = tlv_mu1.M()
        out.mu1_q[0] = chain.JpsiTau_mu1_q
        out.mu1_isLoose[0] = chain.JpsiTau_mu1_isLoose
        out.mu1_isTight[0] = chain.JpsiTau_mu1_isTight
        out.mu1_isPF[0] = chain.JpsiTau_mu1_isPF
        out.mu1_isGlobal[0] = chain.JpsiTau_mu1_isGlobal
        out.mu1_isTracker[0] = chain.JpsiTau_mu1_isTracker
        out.mu1_isSoft[0] = chain.JpsiTau_mu1_isSoft
        out.mu1_dbiso[0] = chain.JpsiTau_mu1_dbiso

        if options.type!='data':
            central=SF_Reco.getSF(math.fabs(tlv_mu1.Eta()),tlv_mu1.Pt())['value']
            uncertainty =SF_Reco.getSF(math.fabs(tlv_mu1.Eta()),tlv_mu1.Pt())['error']
            out.mu1_SFReco[0]=central
            out.mu1_SFReco_up[0]=central+uncertainty
            out.mu1_SFReco_down[0]=central-uncertainty
            central = SF_ID.getSF(math.fabs(tlv_mu1.Eta()),tlv_mu1.Pt())['value']
            uncertainty = SF_ID.getSF(math.fabs(tlv_mu1.Eta()),tlv_mu1.Pt())['error']
            out.mu1_SFID[0]=central
            out.mu1_SFID_up[0]=central+uncertainty
            out.mu1_SFID_down[0]=central-uncertainty
        else :
            out.mu1_SFReco[0]=1
            out.mu1_SFReco_up[0]=1 
            out.mu1_SFReco_down[0]=1
            out.mu1_SFID[0]=1
            out.mu1_SFID_up[0]=1
            out.mu1_SFID_down[0]=1

        out.mu2_pt[0] = tlv_mu2.Pt()
        out.mu2_eta[0] = tlv_mu2.Eta()
        out.mu2_phi[0] = tlv_mu2.Phi()
        out.mu2_mass[0] = tlv_mu2.M()
        out.mu2_q[0] = chain.JpsiTau_mu2_q
        out.mu2_isLoose[0] = chain.JpsiTau_mu2_isLoose
        out.mu2_isTight[0] = chain.JpsiTau_mu2_isTight
        out.mu2_isPF[0] = chain.JpsiTau_mu2_isPF
        out.mu2_isGlobal[0] = chain.JpsiTau_mu2_isGlobal
        out.mu2_isTracker[0] = chain.JpsiTau_mu2_isTracker
        out.mu2_isSoft[0] = chain.JpsiTau_mu2_isSoft
        out.mu2_dbiso[0] = chain.JpsiTau_mu2_dbiso
        if options.type!='data':
            central=SF_Reco.getSF(math.fabs(tlv_mu2.Eta()),tlv_mu2.Pt())['value']
            uncertainty =SF_Reco.getSF(math.fabs(tlv_mu2.Eta()),tlv_mu2.Pt())['error']
            out.mu2_SFReco[0]=central
            out.mu2_SFReco_up[0]=central+uncertainty
            out.mu2_SFReco_down[0]=central-uncertainty
            central = SF_ID.getSF(math.fabs(tlv_mu2.Eta()),tlv_mu2.Pt())['value']
            uncertainty = SF_ID.getSF(math.fabs(tlv_mu2.Eta()),tlv_mu2.Pt())['error']
            out.mu2_SFID[0]=central
            out.mu2_SFID_up[0]=central+uncertainty
            out.mu2_SFID_down[0]=central-uncertainty
        else :
            out.mu2_SFReco[0]=1
            out.mu2_SFReco_up[0]=1
            out.mu2_SFReco_down[0]=1
            out.mu2_SFID[0]=1
            out.mu2_SFID_up[0]=1
            out.mu2_SFID_down[0]=1
  
      
        out.tau_pt[0] = chain.JpsiTau_tau_pt[tindex]
        out.tau_eta[0] = chain.JpsiTau_tau_eta[tindex]
        out.tau_phi[0] = chain.JpsiTau_tau_phi[tindex]
        out.tau_mass[0] = chain.JpsiTau_tau_mass[tindex]
        out.tau_q[0] = chain.JpsiTau_tau_q[tindex]
        out.tau_sumofdnn[0] = chain.JpsiTau_tau_sumofdnn[tindex]
        out.tau_sumofdnn_1prong[0] = chain.JpsiTau_tau_sumofdnn_1prong[tindex]
        out.tau_sumofdnn_otherB[0] = chain.JpsiTau_tau_sumofdnn_otherB[tindex]
        out.tau_sumofdnn_pu[0] = chain.JpsiTau_tau_sumofdnn_pu[tindex]

        out.tau_max_dr[0] = chain.JpsiTau_tau_max_dr_3prong[tindex]
        out.tau_lips[0] = chain.JpsiTau_tau_lips[tindex]
        out.tau_pvips[0] = chain.JpsiTau_tau_pvips[tindex]
        out.tau_dr_jpsi[0] = deltaR(chain.JpsiTau_tau_eta[tindex], chain.JpsiTau_tau_phi[tindex], 
                                    chain.JpsiTau_Jpsi_eta, chain.JpsiTau_Jpsi_phi)

        out.evtid[0] = evtid
        out.tau_eid[0] = evt
        out.tau_index[0] = tindex
        
        out.tau_rhomass1[0] = chain.JpsiTau_tau_rhomass1[tindex]
        out.tau_rhomass2[0] = chain.JpsiTau_tau_rhomass2[tindex]
        
        out.tau_vprob[0] = chain.JpsiTau_tau_vprob[tindex]
        out.tau_fls3d[0] = chain.JpsiTau_tau_fls3d[tindex]
        out.tau_fl3d[0] = chain.JpsiTau_tau_fl3d[tindex]


        out.tau_isRight_3prong[0] = isRight_3prong
        out.tau_isRight_3prong_pi0[0] = isRight_3prong_pi0

        out.ptbal[0] = chain.JpsiTau_ptbal[tindex]
        out.jpsi_tau_alpha[0] = chain.JpsiTau_jpsi_tau_alpha[tindex]

        out.delta_chi2[0] = chain.JpsiTau_tau_delta_chi2[tindex]
        out.delta_n_ch[0] = chain.JpsiTau_tau_delta_n_ch[tindex]
        out.delta_n_mu[0] = chain.JpsiTau_tau_delta_n_mu[tindex]
        out.vweight[0] = chain.JpsiTau_tau_vweight[tindex]
        out.tau_fl3d_wjpsi[0] = chain.JpsiTau_tau_fl3d_wjpsi[tindex]
        out.tau_fls3d_wjpsi[0] = chain.JpsiTau_tau_fls3d_wjpsi[tindex]
        out.tau_alpha_wjpsi[0] = chain.JpsiTau_tau_alpha_wjpsi[tindex]

#        out.tau_refit_vx[0] = chain.JpsiTau_tau_refit_vx[tindex]
#        out.tau_refit_vy[0] = chain.JpsiTau_tau_refit_vy[tindex]
#        out.tau_refit_vz[0] = chain.JpsiTau_tau_refit_vz[tindex]
        out.tau_refit_chi2[0] = chain.JpsiTau_tau_refit_chi2[tindex]
        out.tau_refit_ndof[0] = chain.JpsiTau_tau_refit_ndof[tindex]
        out.tau_refit_rho[0] = chain.JpsiTau_tau_refit_rho[tindex]


        if len(chain.JpsiTau_B_pt)!=0:
#            print len(chain.JpsiTau_B_pt), len(chain.JpsiTau_tau_pt), tindex
            out.b_pt[0] = chain.JpsiTau_B_pt[tindex]
            out.b_eta[0] = chain.JpsiTau_B_eta[tindex]
            out.b_phi[0] = chain.JpsiTau_B_phi[tindex]
            out.b_mass[0] = chain.JpsiTau_B_mass[tindex]
            out.b_mcorr[0] = chain.JpsiTau_B_mcorr[tindex]
            out.b_alpha[0] = chain.JpsiTau_B_alpha[tindex]
            out.b_fls3d[0] = chain.JpsiTau_B_fls3d[tindex]
            out.b_fl3d[0] = chain.JpsiTau_B_fl3d[tindex]
            out.b_vprob[0] = chain.JpsiTau_B_vprob[tindex]
            out.b_pvip[0] = chain.JpsiTau_B_pvip[tindex]
            out.b_pvips[0] = chain.JpsiTau_B_pvips[tindex]
            out.b_lips[0] = chain.JpsiTau_B_lips[tindex]
            out.b_mindoca[0] = chain.JpsiTau_B_mindoca[tindex]
            
            out.b_vz[0] = chain.JpsiTau_B_vz[tindex]
            out.q2[0] = chain.JpsiTau_B_q2[tindex]
            out.mm2[0] = chain.JpsiTau_B_mm2[tindex]
            out.ptmiss[0] = chain.JpsiTau_B_ptmiss[tindex]
            out.estar[0] = chain.JpsiTau_B_Es[tindex]
            out.B_ptback[0] = chain.JpsiTau_B_ptback[tindex]

            _bdx = chain.JpsiTau_B_vx[0] - chain.JpsiTau_tau_refit_vx[tindex]
            _bdy = chain.JpsiTau_B_vy[0] - chain.JpsiTau_tau_refit_vy[tindex]
            _bdz = chain.JpsiTau_B_vz[0] - chain.JpsiTau_tau_refit_vz[tindex]
            
            out.dx_b_pv[0] = _bdx
            out.dy_b_pv[0] = _bdy
            out.dz_b_pv[0] = _bdz
            out.dr_b_pv[0] = math.sqrt(_bdx*_bdx + _bdy*_bdy + _bdz*_bdz)

            out.filt.Fill(5)

        out.tau_iso_0p5[0] = chain.JpsiTau_tau_iso[tindex][0]
        out.tau_iso_ntracks_0p5[0] = chain.JpsiTau_tau_iso_ntracks[tindex][0]

        out.tau_iso_0p7[0] = chain.JpsiTau_tau_iso[tindex][2]
        out.tau_iso_ntracks_0p7[0] = chain.JpsiTau_tau_iso_ntracks[tindex][2]

        out.tau_iso_0p9[0] = chain.JpsiTau_tau_iso[tindex][4]
        out.tau_iso_ntracks_0p9[0] = chain.JpsiTau_tau_iso_ntracks[tindex][4]

        out.tau_iso_1p1[0] = chain.JpsiTau_tau_iso[tindex][6]
        out.tau_iso_ntracks_1p1[0] = chain.JpsiTau_tau_iso_ntracks[tindex][6]

        out.tau_iso_1p3[0] = chain.JpsiTau_tau_iso[tindex][8]
        out.tau_iso_ntracks_1p3[0] = chain.JpsiTau_tau_iso_ntracks[tindex][8]


        pi1_pt = chain.JpsiTau_tau_pi1_pt[tindex]
        pi1_eta = chain.JpsiTau_tau_pi1_eta[tindex]
        pi1_phi = chain.JpsiTau_tau_pi1_phi[tindex]
        pi1_mass = chain.JpsiTau_tau_pi1_mass[tindex]
        pi1_q = chain.JpsiTau_tau_pi1_q[tindex]
        pi1_dnn = chain.JpsiTau_tau_pi1_dnn[tindex]
        pi1_dnn_1prong = chain.JpsiTau_tau_pi1_dnn_1prong[tindex]
        pi1_dnn_otherB = chain.JpsiTau_tau_pi1_dnn_otherB[tindex]
        pi1_dnn_pu = chain.JpsiTau_tau_pi1_dnn_pu[tindex]
        pi1_trigMatch = bool(chain.JpsiTau_tau_pi1_trigMatch[tindex])
        pi1_trigMatch_dr = chain.JpsiTau_tau_pi1_trigMatch_dr[tindex]
        pi1_d0sig = chain.JpsiTau_tau_pi1_d0sig[tindex]
        pi1_mumuTrk_vprob = chain.JpsiTau_tau_pi1_mumuTrk_vprob[tindex]
        pi1_mumuTrk_chi2 = chain.JpsiTau_tau_pi1_mumuTrk_chi2[tindex]
        pi1_mumuTrk_ndof = chain.JpsiTau_tau_pi1_mumuTrk_ndof[tindex]


        pi2_pt = chain.JpsiTau_tau_pi2_pt[tindex]
        pi2_eta = chain.JpsiTau_tau_pi2_eta[tindex]
        pi2_phi = chain.JpsiTau_tau_pi2_phi[tindex]
        pi2_mass = chain.JpsiTau_tau_pi2_mass[tindex]
        pi2_q = chain.JpsiTau_tau_pi2_q[tindex]
        pi2_dnn = chain.JpsiTau_tau_pi2_dnn[tindex]
        pi2_dnn_1prong = chain.JpsiTau_tau_pi2_dnn_1prong[tindex]
        pi2_dnn_otherB = chain.JpsiTau_tau_pi2_dnn_otherB[tindex]
        pi2_dnn_pu = chain.JpsiTau_tau_pi2_dnn_pu[tindex]
        pi2_trigMatch = bool(chain.JpsiTau_tau_pi2_trigMatch[tindex])
        pi2_trigMatch_dr = chain.JpsiTau_tau_pi2_trigMatch_dr[tindex]
        pi2_d0sig = chain.JpsiTau_tau_pi2_d0sig[tindex]
        pi2_mumuTrk_vprob = chain.JpsiTau_tau_pi2_mumuTrk_vprob[tindex]
        pi2_mumuTrk_chi2 = chain.JpsiTau_tau_pi2_mumuTrk_chi2[tindex]
        pi2_mumuTrk_ndof = chain.JpsiTau_tau_pi2_mumuTrk_ndof[tindex]


        pi3_pt = chain.JpsiTau_tau_pi3_pt[tindex]
        pi3_eta = chain.JpsiTau_tau_pi3_eta[tindex]
        pi3_phi = chain.JpsiTau_tau_pi3_phi[tindex]
        pi3_mass = chain.JpsiTau_tau_pi3_mass[tindex]
        pi3_q = chain.JpsiTau_tau_pi3_q[tindex]
        pi3_dnn = chain.JpsiTau_tau_pi3_dnn[tindex]
        pi3_dnn_1prong = chain.JpsiTau_tau_pi3_dnn_1prong[tindex]
        pi3_dnn_otherB = chain.JpsiTau_tau_pi3_dnn_otherB[tindex]
        pi3_dnn_pu = chain.JpsiTau_tau_pi3_dnn_pu[tindex]
        pi3_trigMatch = bool(chain.JpsiTau_tau_pi3_trigMatch[tindex])
        pi3_trigMatch_dr = chain.JpsiTau_tau_pi3_trigMatch_dr[tindex]
        pi3_d0sig = chain.JpsiTau_tau_pi3_d0sig[tindex]
        pi3_mumuTrk_vprob = chain.JpsiTau_tau_pi3_mumuTrk_vprob[tindex]
        pi3_mumuTrk_chi2 = chain.JpsiTau_tau_pi3_mumuTrk_chi2[tindex]
        pi3_mumuTrk_ndof = chain.JpsiTau_tau_pi3_mumuTrk_ndof[tindex]

        out.pi_forwardest_eta[0] = max(abs(pi1_eta), abs(pi2_eta), abs(pi3_eta))


        out.pi1_pt[0] = pi1_pt
        out.pi1_eta[0] = pi1_eta
        out.pi1_phi[0] = pi1_phi
        out.pi1_dnn[0] = pi1_dnn
        out.pi1_dnn_1prong[0] = pi1_dnn_1prong
        out.pi1_dnn_otherB[0] = pi1_dnn_otherB
        out.pi1_dnn_pu[0] = pi1_dnn_pu
        out.pi1_trigMatch[0] = pi1_trigMatch
        out.pi1_trigMatch_dr[0] = pi1_trigMatch_dr
        out.pi1_d0sig[0] = pi1_d0sig
        out.pi1_mumuTrk_vprob[0] = pi1_mumuTrk_vprob
        out.pi1_mumuTrk_chi2[0] = pi1_mumuTrk_chi2
        out.pi1_mumuTrk_ndof[0] = pi1_mumuTrk_ndof

        out.pi2_pt[0] = pi2_pt
        out.pi2_eta[0] = pi2_eta
        out.pi2_phi[0] = pi2_phi
        out.pi2_dnn[0] = pi2_dnn
        out.pi2_dnn_1prong[0] = pi2_dnn_1prong
        out.pi2_dnn_otherB[0] = pi2_dnn_otherB
        out.pi2_dnn_pu[0] = pi2_dnn_pu
        out.pi2_trigMatch[0] = pi2_trigMatch
        out.pi2_trigMatch_dr[0] = pi2_trigMatch_dr
        out.pi2_d0sig[0] = pi2_d0sig
        out.pi2_mumuTrk_vprob[0] = pi2_mumuTrk_vprob
        out.pi2_mumuTrk_chi2[0] = pi2_mumuTrk_chi2
        out.pi2_mumuTrk_ndof[0] = pi2_mumuTrk_ndof


        out.pi3_pt[0] = pi3_pt
        out.pi3_eta[0] = pi3_eta
        out.pi3_phi[0] = pi3_phi
        out.pi3_dnn[0] = pi3_dnn
        out.pi3_dnn_1prong[0] = pi3_dnn_1prong
        out.pi3_dnn_otherB[0] = pi3_dnn_otherB
        out.pi3_dnn_pu[0] = pi3_dnn_pu
        out.pi3_trigMatch[0] = pi3_trigMatch
        out.pi3_trigMatch_dr[0] = pi3_trigMatch_dr
        out.pi3_d0sig[0] = pi3_d0sig
        out.pi3_mumuTrk_vprob[0] = pi3_mumuTrk_vprob
        out.pi3_mumuTrk_chi2[0] = pi3_mumuTrk_chi2
        out.pi3_mumuTrk_ndof[0] = pi3_mumuTrk_ndof


        tlv1 = ROOT.TLorentzVector()
        tlv1.SetPtEtaPhiM(pi1_pt, pi1_eta, pi1_phi, pi1_mass)
        
        tlv2 = ROOT.TLorentzVector()
        tlv2.SetPtEtaPhiM(pi2_pt, pi2_eta, pi2_phi, pi2_mass)
        
        tlv3 = ROOT.TLorentzVector()
        tlv3.SetPtEtaPhiM(pi3_pt, pi3_eta, pi3_phi, pi3_mass)

        tlv1_k = ROOT.TLorentzVector()
        tlv1_k.SetPtEtaPhiM(pi1_pt, pi1_eta, pi1_phi, mk)
        
        tlv2_k = ROOT.TLorentzVector()
        tlv2_k.SetPtEtaPhiM(pi2_pt, pi2_eta, pi2_phi, mk)
        
        tlv3_k = ROOT.TLorentzVector()
        tlv3_k.SetPtEtaPhiM(pi3_pt, pi3_eta, pi3_phi, mk)
        
        tlvs = [tlv1, tlv2, tlv3]

        mass_kpipi = (tlv_jpsi + tlv1_k + tlv2 + tlv3).M()
#        mass_pikpi = (tlv_jpsi + tlv1 + tlv2_k + tlv3).M()
#        mass_pipik = (tlv_jpsi + tlv1 + tlv2 + tlv3_k).M()

#        mass_b_pm = 5.27934
#        
#        masses_b_pm = [mass_kpipi, mass_pikpi, mass_pipik]
#
#        diffes_b_pm = [abs(m - mass_b_pm) for m in masses_b_pm]
#        
#        min_mass_b_pm = min(diffes_b_pm)
#        index_b_pm = diffes_b_pm.index(min_mass_b_pm)

#        print("-"*50)
#        print("masses_b_pm", masses_b_pm)
#        print("diffes_b_pm", diffes_b_pm)
#        print("index_b_pm", index_b_pm)
#        print("closest_b_pm", masses_b_pm[index_b_pm])


        out.jpsi_kpipi[0] = mass_kpipi
#        out.jpsi_pikpi[0] = mass_pikpi
#        out.jpsi_pipik[0] = mass_pipik
        out.jpsi_isrestos[0] = (pi2_q*pi3_q==-1)
        out.jpsi_kpipi_psi2smass[0] = (tlv_jpsi + tlv2 + tlv3).M()
        out.jpsi_kpipi_rho[0] = (tlv2 + tlv3).M()
        out.jpsi_kpipi_ks1[0] = (tlv1_k + tlv2).M()
        out.jpsi_kpipi_ks2[0] = (tlv1_k + tlv3).M()
#        out.jpsi_kpipi_closest[0] = masses_b_pm[index_b_pm]

        

        out.tau_m12s[0] = math.pow((tlv1+tlv2).M(), 2);
        out.tau_m23s[0] = math.pow((tlv2+tlv3).M(), 2);
        out.tau_m13s[0] = math.pow((tlv1+tlv3).M(), 2);

        out.tau_kpp[0] = returnMass3(tlvs, mk, mp, mp)
        out.tau_pkp[0] = returnMass3(tlvs, mp, mk, mp)
        out.tau_ppk[0] = returnMass3(tlvs, mp, mp, mk)
        out.tau_kkp[0] = returnMass3(tlvs, mk, mk, mp)
        out.tau_kpk[0] = returnMass3(tlvs, mk, mp, mk)
        out.tau_pkk[0] = returnMass3(tlvs, mp, mk, mk)
        out.tau_kkk[0] = returnMass3(tlvs, mk, mk, mk)
        
        m_kp = []
        m_pk = []
        m_kk = []
        m_lambda = []

       
        tlvs_psi = []
        tlvs_psi_pk = []
        tlvs_psi_kp = []
        tlvs_psi_kk = []


        if pi1_q*pi2_q == -1:
            tlvs_12 = [tlv1, tlv2]
        
            m_kp.append(returnMass2(tlvs_12, mk, mp))
            m_pk.append(returnMass2(tlvs_12, mp, mk))
            m_kk.append(returnMass2(tlvs_12, mk, mk))
            m_lambda.append(returnMass2(tlvs_12, mproton, mp))
            tlvs_psi.append((tlv_jpsi + tlv1 + tlv2).M())
            tlvs_psi_pk.append((tlv_jpsi + tlv1 + tlv2_k).M())
            tlvs_psi_kp.append((tlv_jpsi + tlv1_k + tlv2).M())
            tlvs_psi_kk.append((tlv_jpsi + tlv1_k + tlv2_k).M())

        if pi1_q*pi3_q == -1:
            tlvs_13 = [tlv1, tlv3]

            m_kp.append(returnMass2(tlvs_13, mk, mp))
            m_pk.append(returnMass2(tlvs_13, mp, mk))
            m_kk.append(returnMass2(tlvs_13, mk, mk))
            m_lambda.append(returnMass2(tlvs_13, mproton, mp))
            tlvs_psi.append((tlv_jpsi + tlv1 + tlv3).M())
            tlvs_psi_pk.append((tlv_jpsi + tlv1 + tlv3_k).M())
            tlvs_psi_kp.append((tlv_jpsi + tlv1_k + tlv3).M())
            tlvs_psi_kk.append((tlv_jpsi + tlv1_k + tlv3_k).M())
        
        if pi2_q*pi3_q == -1:
            tlvs_23 = [tlv2, tlv3]

            m_kp.append(returnMass2(tlvs_23, mk, mp))
            m_pk.append(returnMass2(tlvs_23, mp, mk))
            m_kk.append(returnMass2(tlvs_23, mk, mk))
            m_lambda.append(returnMass2(tlvs_23, mproton, mp))
            tlvs_psi.append((tlv_jpsi + tlv2 + tlv3).M())        
            tlvs_psi_pk.append((tlv_jpsi + tlv2 + tlv3_k).M())
            tlvs_psi_kp.append((tlv_jpsi + tlv2_k + tlv3).M())
            tlvs_psi_kk.append((tlv_jpsi + tlv2_k + tlv3_k).M())


#        rhomass_ss = -1.
#
#        if pi1_q*pi2_q == 1:
#            rhomass_ss = (tlv1 + tlv2).M()
#
#        if pi1_q*pi3_q == 1:
#            rhomass_ss = (tlv1 + tlv3).M()
#
#        if pi2_q*pi3_q == 1:
#            rhomass_ss = (tlv2 + tlv3).M()

        out.tau_rhomass_ss[0] = chain.JpsiTau_tau_rhomass_ss[tindex]

#        print 'check', rhomass_ss, 
        out.tau_rhomass[0] = random.choice([out.tau_rhomass1[0], out.tau_rhomass2[0]])


        if not (len(m_kp)!=2 or len(m_pk)!=2 or len(m_kk)!=2 or len(m_lambda)!=2):

            out.tau_rhomass1_kp[0] = m_kp[0]
            out.tau_rhomass2_kp[0] = m_kp[1]
            
            out.tau_rhomass1_pk[0] = m_pk[0]
            out.tau_rhomass2_pk[0] = m_pk[1]
            
            out.tau_rhomass1_kk[0] = m_kk[0]
            out.tau_rhomass2_kk[0] = m_kk[1]

            out.tau_lambda1[0] = m_lambda[0]
            out.tau_lambda2[0] = m_lambda[1]

            out.tau_psimass1[0] = tlvs_psi[0]
            out.tau_psimass2[0] = tlvs_psi[1]

            out.tau_psimass1_pk[0] = tlvs_psi_pk[0]
            out.tau_psimass2_pk[0] = tlvs_psi_pk[1]
            
            out.tau_psimass1_kp[0] = tlvs_psi_kp[0]
            out.tau_psimass2_kp[0] = tlvs_psi_kp[1]
            
            out.tau_psimass1_kk[0] = tlvs_psi_kk[0]
            out.tau_psimass2_kk[0] = tlvs_psi_kk[1]

            out.tau_jpsi_pi1[0] = (tlv_jpsi + tlv1).M()
            out.tau_jpsi_k1[0] = (tlv_jpsi + tlv1_k).M()
            
            out.tau_jpsi_pi2[0] = (tlv_jpsi + tlv2).M()
            out.tau_jpsi_k2[0] = (tlv_jpsi + tlv2_k).M()
            
            out.tau_jpsi_pi3[0] = (tlv_jpsi + tlv3).M()
            out.tau_jpsi_k3[0] = (tlv_jpsi + tlv3_k).M()

        # per-event quantity
        out.nch[0] = chain.JpsiTau_nch
        out.nch_before_dnn[0] = chain.JpsiTau_nch_before
        out.bbpv_vz[0] = chain.JpsiTau_bbPV_vz
        out.pv_vz[0] = chain.JpsiTau_PV_vz
        out.ncand[0] = chain.JpsiTau_nCandidates

        _dx = chain.JpsiTau_Jpsi_vx - chain.JpsiTau_tau_vx[tindex]
        _dy = chain.JpsiTau_Jpsi_vy - chain.JpsiTau_tau_vy[tindex]
        _dz = chain.JpsiTau_Jpsi_vz - chain.JpsiTau_tau_vz[tindex]
        
        out.dx_jpsi_tau[0] = _dx
        out.dy_jpsi_tau[0] = _dy
        out.dz_jpsi_tau[0] = _dz
        out.dr_jpsi_tau[0] = math.sqrt(_dx*_dx + _dy*_dy + _dz*_dz)
    
        

        out.jpsi_pt[0] = chain.JpsiTau_Jpsi_pt
        out.jpsi_eta[0] = chain.JpsiTau_Jpsi_eta
        out.jpsi_phi[0] = chain.JpsiTau_Jpsi_phi
        out.jpsi_mass[0] = chain.JpsiTau_Jpsi_mass
        out.jpsi_vprob[0] = chain.JpsiTau_Jpsi_vprob
        out.jpsi_fls3d[0] = chain.JpsiTau_Jpsi_fls3d
        out.jpsi_fl3d[0] = chain.JpsiTau_Jpsi_fl3d
        out.jpsi_alpha2d[0] = chain.JpsiTau_Jpsi_alpha2d
        out.jpsi_fls2d[0] = chain.JpsiTau_Jpsi_fls2d
        out.jpsi_dr[0] = chain.JpsiTau_Jpsi_dr
        out.jpsi_maxdoca[0] = chain.JpsiTau_Jpsi_maxdoca
        out.jpsi_mindoca[0] = chain.JpsiTau_Jpsi_mindoca
        
        out.npv[0] = chain.PV_N

        out.b_pt_simple[0] = chain.JpsiTau_B_pt_simple[tindex]
        out.b_eta_simple[0] = chain.JpsiTau_B_eta_simple[tindex]
        out.b_phi_simple[0] = chain.JpsiTau_B_phi_simple[tindex]
        out.b_mass_simple[0] = chain.JpsiTau_B_mass_simple[tindex]

        out.q2_simple[0] = chain.JpsiTau_B_q2_simple[tindex]
        out.mm2_simple[0] = chain.JpsiTau_B_mm2_simple[tindex]
        out.ptmiss_simple[0] = chain.JpsiTau_B_ptmiss_simple[tindex]
        out.estar_simple[0] = chain.JpsiTau_B_Es_simple[tindex]
        out.B_ptback_simple[0] = chain.JpsiTau_B_ptback_simple[tindex]


        out.perEVT_data[0] = chain.JpsiTau_perEVT_data
        out.perEVT_mc[0] = chain.JpsiTau_perEVT_mc
#        out.perEVT_otherB[0] = chain.JpsiTau_perEVT_otherB
#        out.perEVT_sig[0] = chain.JpsiTau_perEVT_sig
#        out.perEVT_leptonic[0] = chain.JpsiTau_perEVT_leptonic
#        out.perEVT_1prong[0] = chain.JpsiTau_perEVT_1prong

        out.evt[0] = chain.EVENT_event
        out.lumi[0] = chain.EVENT_lumiBlock
        out.run[0] = chain.EVENT_run
        
#        if options.type == 'signal':
#
#            for igen in range(len(chain.genParticle_pdgs)):
#
#                print '-'*80
#                print 'gen = ', igen
#                print '-'*80
#
#                for ipdg in range(len(chain.genParticle_pdgs[igen])):
#            
#                    print '  '*2*int(chain.genParticle_layers[igen][ipdg]), 'pdg  = ', chain.genParticle_pdgs[igen][ipdg], '(',  returnName(chain.genParticle_pdgs[igen][ipdg]) , '), (pt, eta, phi) = ', '({0:.2f}'.format(chain.genParticle_ppt[igen][ipdg]), '{0:.2f}'.format(chain.genParticle_peta[igen][ipdg]), '{0:.2f}'.format(chain.genParticle_pphi[igen][ipdg]), '), isfinal=',  chain.genParticle_isfinal[igen][ipdg]








#        if options.type=='bg' and mass_kpipi > 5.25 and mass_kpipi < 5.3: 
        if options.type=='bg':

#            for idau in range(len(chain.genParticle_dau[igen])):
#
#                print '\t -->', chain.genParticle_dau[igen][idau]
#
#                flag_pi = []
#
#                for ipi, itlv in enumerate([tlv1, tlv2, tlv3]):
#
#                    for ipdg in range(len(chain.genParticle_pdgs[igen])):            
#                                
#                        _dr = deltaR(itlv.Eta(), itlv.Phi(), chain.genParticle_peta[igen][ipdg], chain.genParticle_pphi[igen][ipdg])
#                                
#                        if _dr < 0.1:
#                            flag_pi.append(ipi)
#                                    
#                print '\t', flag_pi




#            for igen in range(len(chain.genParticle_pdgId)):
#                if abs(chain.genParticle_pdgId[igen])>=511 and abs(chain.genParticle_pdgId[igen])<=545:
#
#                    print 'TRACE', chain.genParticle_pdgId[igen]
#                
#                    for idau in range(len(chain.genParticle_dau[igen])):
#                        print '\t -->', chain.genParticle_dau[igen][idau]
                    


            decayTable = []

#            flag_jpsi = False
            
#####            print '-'*80
#####            print 'evtid = ', evtid 
#####            print '-'*80
#####            print 'there are ', len(chain.genParticle_pdgs), 'Bs ...'
#####
#####
#####            for igen in range(len(chain.genParticle_pdgs)):
#####
#####                print ''
#####                print 'gen = ', igen
#####                print ''
#####
#####                for ipdg in range(len(chain.genParticle_pdgs[igen])):
#####            
#####
#####                    print '  '*2*int(chain.genParticle_layers[igen][ipdg]), 'pdg  = ', chain.genParticle_pdgs[igen][ipdg], '(',  returnName(chain.genParticle_pdgs[igen][ipdg]) , '), (pt, eta, phi) = ', '({0:.2f}'.format(chain.genParticle_ppt[igen][ipdg]), '{0:.2f}'.format(chain.genParticle_peta[igen][ipdg]), '{0:.2f}'.format(chain.genParticle_pphi[igen][ipdg]), ')'
#####
#####                    if abs(chain.genParticle_pdgs[igen][ipdg]) == 13:
#####                        dr1 = deltaR(tlv_mu1.Eta(), tlv_mu1.Phi(), chain.genParticle_peta[igen][ipdg], chain.genParticle_pphi[igen][ipdg])
#####                        dr2 = deltaR(tlv_mu2.Eta(), tlv_mu2.Phi(), chain.genParticle_peta[igen][ipdg], chain.genParticle_pphi[igen][ipdg])
#####                        
#####                        if dr1 < 0.05:
######                            flag_mu1 = True
#####                            print '*'
#####                            
#####                        if dr2 < 0.05:
######                            flag_mu2 = True
#####                            print '*'
#####
#####                    else:
#####                        
#####                        for itlv in [tlv1, tlv2, tlv3]:
#####
#####                            _dr = deltaR(itlv.Eta(), itlv.Phi(), chain.genParticle_peta[igen][ipdg], chain.genParticle_pphi[igen][ipdg])
#####                        
#####                            if _dr < 0.1:
#####                                print '*'

                            



            for igen in range(len(chain.genParticle_pdgs)):

                flag_mu1 = False
                flag_mu2 = False


                for ipdg in range(len(chain.genParticle_pdgs[igen])):
            

                    if abs(chain.genParticle_pdgs[igen][ipdg]) == 13:
                        dr1 = deltaR(tlv_mu1.Eta(), tlv_mu1.Phi(), chain.genParticle_peta[igen][ipdg], chain.genParticle_pphi[igen][ipdg])
                        dr2 = deltaR(tlv_mu2.Eta(), tlv_mu2.Phi(), chain.genParticle_peta[igen][ipdg], chain.genParticle_pphi[igen][ipdg])
                        
                        if dr1 < 0.05:
                            flag_mu1 = True
#                            print '*'
                            
                        if dr2 < 0.05:
                            flag_mu2 = True
#                            print '*'

                if not (flag_mu1 and flag_mu2): continue


#                print '-'*80
#                print 'gen = ', igen
#                print '-'*80


#                for ipdg in range(len(chain.genParticle_pdgs[igen])):
            
#                    print '  '*2*int(chain.genParticle_layers[igen][ipdg]), 'pdg  = ', chain.genParticle_pdgs[igen][ipdg], '(',  returnName(chain.genParticle_pdgs[igen][ipdg]) , '), (pt, eta, phi) = ', '({0:.2f}'.format(chain.genParticle_ppt[igen][ipdg]), '{0:.2f}'.format(chain.genParticle_peta[igen][ipdg]), '{0:.2f}'.format(chain.genParticle_pphi[igen][ipdg]), '), isfinal=',  chain.genParticle_isfinal[igen][ipdg]



                # first identify the layer at which B+/- is created

                _layerid = -1 
                addstr = 'None'

                for ipdg in range(len(chain.genParticle_pdgs[igen])):

                    if abs(chain.genParticle_pdgs[igen][ipdg]) in [511, 521, 531, 541]: 
                        addstr = returnName(chain.genParticle_pdgs[igen][ipdg])
                        _layerid = chain.genParticle_layers[igen][ipdg]
                        break

                    
                for ipdg in range(len(chain.genParticle_pdgs[igen])):
                    if chain.genParticle_layers[igen][ipdg]==_layerid+1:
                        addstr += '_' + returnName(chain.genParticle_pdgs[igen][ipdg])


#                print '\t EVENT', chain.EVENT_event, ' ==========>', addstr
                decayTable.append(addstr)


#                for itlv in [tlv1, tlv2, tlv3]:
#
#                    flag_pi = False
#                    addstr = 'PU'
#
#
#                    for ipdg in range(len(chain.genParticle_pdgs[igen])):            
#
#                        _dr = deltaR(itlv.Eta(), itlv.Phi(), chain.genParticle_peta[igen][ipdg], chain.genParticle_pphi[igen][ipdg])
#                        
#                        if _dr < 0.1:
#                            flag_pi = True
#                            addstr = returnName(chain.genParticle_pdgs[igen][ipdg])
#                            break
#                            
#                                
#                    decayTable.append(addstr)


#                    else:
#
#                        flag_pi_other = False
#
#                        for igen2 in range(len(chain.genParticle_pdgs)):
#                            if igen2==igen: continue
#
#                            addstr = returnName(chain.genParticle_pdgs[igen2][0])
#                        
#                            for ipdg in range(len(chain.genParticle_pdgs[igen2])):
#                            
#                                if chain.genParticle_layers[igen2][ipdg]==1:
#                                    addstr += '_' + returnName(chain.genParticle_pdgs[igen2][ipdg])
#
#
#                            _dr = deltaR(itlv.Eta(), itlv.Phi(), chain.genParticle_peta[igen2][ipdg], chain.genParticle_pphi[igen2][ipdg])
#                            
#                            if _dr < 0.1:
#                                flag_pi_other = True
#                                decayTable.append(addstr)
                                

 

            if len(decayTable)!=1:
                out.procid[0] = -9
            else:

                for idt in decayTable:
                    dict_counter += 1
            
                    if idt in dicts: 
#                    print 'This is already there!!!'
                        dicts[idt] += 1
                    else:
                        #                    print 'This is new!'
                        dicts[idt] = 1 


                    if not options.create:

                        if idt in idtable:
                            out.procid[0] = idtable[idt]['id']
                        else:
                            out.procid[0] = -1
                    else:
                        out.procid[0] = -99



#####
#####
######        if options.type=='bg':
######            if (out.tau_rhomass1_kk[0] > 1.04) and (out.tau_rhomass2_kk[0] > 1.04) and (3.62 > out.tau_psimass1[0] or out.tau_psimass1[0] > 3.75) and (3.62 > out.tau_psimass2[0] or out.tau_psimass2[0] > 3.75) and (out.tau_rhomass1_pk[0] < 0.84 or out.tau_rhomass1_pk[0] > 0.94) and (out.tau_rhomass2_pk[0] < 0.84 or out.tau_rhomass2_pk[0] < 0.94) and (out.tau_psimass1_kp[0] < 5.15 or out.tau_psimass1_kp[0] > 5.38) and (out.tau_psimass2_kp[0] < 5.15 or out.tau_psimass2_kp[0] > 5.38) and out.tau_rhomass1[0] > 0.75 and out.tau_rhomass1[0] < 0.79:
######
######                print '-'*80
######
######                for igen in range(len(chain.genParticle_pdgId)):
######                    if abs(chain.genParticle_pdgId[igen])>=511 and abs(chain.genParticle_pdgId[igen])<=545:
######
######                        print 'TRACE2', chain.genParticle_pdgId[igen]
######                
######                        lp = []
######
######                        for idau in range(len(chain.genParticle_dau[igen])):
######                            print 'TRACE2 -->', chain.genParticle_dau[igen][idau]
######                            lp.append(abs(chain.genParticle_dau[igen][idau]))
######                    
######                        if 15 in lp and 443 in lp and abs(chain.genParticle_pdgId[igen])==541:
######                            print 'THIS IS THE SIGNAL2 !!!'
######                            flag_veto = True
#####
#####
#####            print decayTable
#####    
#####            if flag_jpsi: counter_jpsiassociation += 1
#####            else: counter_nojpsiassociation += 1
#####
#####            out.v.clear()
#####            for s in decayTable:
#####                out.v.push_back( s )
######                t.Fill()
#####
######            out.dchain[0] = '='.join(decayTable)

        if options.type != 'data':

            out.gen_vz[0] = chain.JpsiTau_genPV_vz

            out.npv_true[0] = chain.nPuVtxTrue[0] # current BX
            
            out.tau_genpt[0] = chain.JpsiTau_st_gentau_pt
            out.tau_geneta[0] = chain.JpsiTau_st_gentau_eta
            out.tau_genphi[0] = chain.JpsiTau_st_gentau_phi
            out.isjpsimatched[0] = chain.JpsiTau_isgenmatched
        
            out.q2_gen[0] = chain.JpsiTau_q2_gen
            out.B_pt_gen[0] = chain.JpsiTau_B_pt_gen

            out.n_occurance[0] = chain.JpsiTau_st_n_occurance
            out.decayid[0] = chain.JpsiTau_st_decayid

            if options.type == 'bg':
                out.genWeightBkgB[0] = chain.genWeightBkgB

            if options.type in ['signal']:

#                if isRight_3prong:
#                    if len(chain.truth_tau_dipion1_mass)==1:
#                        out.gen_dipion1_mass[0] = chain.truth_tau_dipion1_mass[0]
#                        out.gen_dipion2_mass[0] = chain.truth_tau_dipion2_mass[0]
#                        out.gen_dipion_unrolled[0] = int((min([chain.truth_tau_dipion2_mass[0], 1.3]) - 0.2)/0.11) + 11*int((min([chain.truth_tau_dipion1_mass[0], 1.3]) - 0.2)/0.11)
#                        out.gen_dipion_unrolled_coarse[0] = int((min([chain.truth_tau_dipion2_mass[0], 1.3]) - 0.2)/0.22) + 6*int((min([chain.truth_tau_dipion1_mass[0], 1.3]) - 0.2)/0.22)
#                    else:
#                        out.gen_dipion1_mass[0] = -1
#                        out.gen_dipion2_mass[0] = -1
#                        out.gen_dipion_unrolled[0] = -1 
#                        out.gen_dipion_unrolled_coarse[0] = -1
#                else:
#                        out.gen_dipion1_mass[0] = -9
#                        out.gen_dipion2_mass[0] = -9
#                        out.gen_dipion_unrolled[0] = -9 
#                        out.gen_dipion_unrolled_coarse[0] = -9
 

                if len(chain.truth_tau_dipion1_mass)==1 and chain.JpsiTau_st_n_occurance == 1:
                    out.gen_dipion1_mass[0] = chain.truth_tau_dipion1_mass[0]
                    out.gen_dipion2_mass[0] = chain.truth_tau_dipion2_mass[0]
#                    out.gen_dipion_unrolled[0] = int((min([chain.truth_tau_dipion2_mass[0], 1.3]) - 0.2)/0.11) + 11*int((min([chain.truth_tau_dipion1_mass[0], 1.3]) - 0.2)/0.11)
#                    out.gen_dipion_unrolled_coarse[0] = int((min([chain.truth_tau_dipion2_mass[0], 1.3]) - 0.2)/0.22) + 6*int((min([chain.truth_tau_dipion1_mass[0], 1.3]) - 0.2)/0.22)
                else:
                    out.gen_dipion1_mass[0] = -1
                    out.gen_dipion2_mass[0] = -1
#                    out.gen_dipion_unrolled[0] = -1 
#                    out.gen_dipion_unrolled_coarse[0] = -1


                   

#                import pdb; pdb.set_trace()
            if options.type in ['signal'] and chain.JpsiTau_st_n_occurance == 1:

                out.filt.Fill(6)
                
                if len(chain.JpsiTau_hammer_ebe)==1:


#                    print evtid, chain.JpsiTau_hammer_ebe[0], len(chain.JpsiTau_hammer_ebe)
                    
                    if chain.JpsiTau_hammer_ebe[0]==-1:
                        continue

                    out.hammer_ebe[0] = chain.JpsiTau_hammer_ebe[0]
                    out.hammer_wratio[0] = ROOT.Double(hist_hammer.GetBinContent(2))/ROOT.Double(hist_hammer.GetBinContent(1))
                    
                    for ii in range(0, 15):
                        for ud in ['up', 'down']:
                            
                            hammer_weight = getattr(chain, 'JpsiTau_hammer_ebe_e' + str(ii) + '_' + ud)[0]
                            if hammer_weight==-1: continue
                            getattr(out, 'hammer_ebe_e' + str(ii) + '_' + ud)[0] = hammer_weight

                        
#                    out.hammer_ebe_a0_up[0] = chain.JpsiTau_hammer_ebe_a0_up[0]
#                    out.hammer_ebe_a0_down[0] = chain.JpsiTau_hammer_ebe_a0_down[0]
#                    out.hammer_ebe_a1_up[0] = chain.JpsiTau_hammer_ebe_a1_up[0]
#                    out.hammer_ebe_a1_down[0] = chain.JpsiTau_hammer_ebe_a1_down[0]
#                    out.hammer_ebe_a2_up[0] = chain.JpsiTau_hammer_ebe_a2_up[0]
#                    out.hammer_ebe_a2_down[0] = chain.JpsiTau_hammer_ebe_a2_down[0]
#
#                    out.hammer_ebe_b0_up[0] = chain.JpsiTau_hammer_ebe_b0_up[0]
#                    out.hammer_ebe_b0_down[0] = chain.JpsiTau_hammer_ebe_b0_down[0]
#                    out.hammer_ebe_b1_up[0] = chain.JpsiTau_hammer_ebe_b1_up[0]
#                    out.hammer_ebe_b1_down[0] = chain.JpsiTau_hammer_ebe_b1_down[0]
#                    out.hammer_ebe_b2_up[0] = chain.JpsiTau_hammer_ebe_b2_up[0]
#                    out.hammer_ebe_b2_down[0] = chain.JpsiTau_hammer_ebe_b2_down[0]
#
#                    out.hammer_ebe_c1_up[0] = chain.JpsiTau_hammer_ebe_c1_up[0]
#                    out.hammer_ebe_c1_down[0] = chain.JpsiTau_hammer_ebe_c1_down[0]
#                    out.hammer_ebe_c2_up[0] = chain.JpsiTau_hammer_ebe_c2_up[0]
#                    out.hammer_ebe_c2_down[0] = chain.JpsiTau_hammer_ebe_c2_down[0]
#
#                    out.hammer_ebe_d0_up[0] = chain.JpsiTau_hammer_ebe_d0_up[0]
#                    out.hammer_ebe_d0_down[0] = chain.JpsiTau_hammer_ebe_d0_down[0]
#                    out.hammer_ebe_d1_up[0] = chain.JpsiTau_hammer_ebe_d1_up[0]
#                    out.hammer_ebe_d1_down[0] = chain.JpsiTau_hammer_ebe_d1_down[0]
#                    out.hammer_ebe_d2_up[0] = chain.JpsiTau_hammer_ebe_d2_up[0]
#                    out.hammer_ebe_d2_down[0] = chain.JpsiTau_hammer_ebe_d2_down[0]




#                    for iham in range(len(chain.JpsiTau_hammer_ebe_toy[0])):
            
#                        out.hammer_ebe_toy[iham] = ROOT.Double(chain.JpsiTau_hammer_ebe_toy[0][iham])

                else:
                    continue


#                if chain.JpsiTau_hammer_ebe[0]==-1: continue
            


                if len(chain.JpsiTau_hammer_ebe_lattice)==1:
#                    print evtid, chain.JpsiTau_hammer_ebe[0], len(chain.JpsiTau_hammer_ebe)

                    if chain.JpsiTau_hammer_ebe_lattice[0]==-1:
                        continue

                    out.hammer_ebe_lattice[0] = chain.JpsiTau_hammer_ebe_lattice[0]
                    out.hammer_wratio_lattice[0] = ROOT.Double(hist_hammer_lattice.GetBinContent(2))/ROOT.Double(hist_hammer_lattice.GetBinContent(1))

                    for ii in range(0, 15):
                        for ud in ['up', 'down']:
                            hammer_weight = getattr(chain, 'JpsiTau_hammer_ebe_e' + str(ii) + '_' + ud + '_lattice')[0]
                            if hammer_weight==-1: continue
                            getattr(out, 'hammer_ebe_e' + str(ii) + '_' + ud + '_lattice')[0] = hammer_weight
 

#                    out.hammer_ebe_a0_up_lattice[0] = chain.JpsiTau_hammer_ebe_a0_up_lattice[0]
#                    out.hammer_ebe_a0_down_lattice[0] = chain.JpsiTau_hammer_ebe_a0_down_lattice[0]
#                    out.hammer_ebe_a1_up_lattice[0] = chain.JpsiTau_hammer_ebe_a1_up_lattice[0]
#                    out.hammer_ebe_a1_down_lattice[0] = chain.JpsiTau_hammer_ebe_a1_down_lattice[0]
#                    out.hammer_ebe_a2_up_lattice[0] = chain.JpsiTau_hammer_ebe_a2_up_lattice[0]
#                    out.hammer_ebe_a2_down_lattice[0] = chain.JpsiTau_hammer_ebe_a2_down_lattice[0]
#
#                    out.hammer_ebe_b0_up_lattice[0] = chain.JpsiTau_hammer_ebe_b0_up_lattice[0]
#                    out.hammer_ebe_b0_down_lattice[0] = chain.JpsiTau_hammer_ebe_b0_down_lattice[0]
#                    out.hammer_ebe_b1_up_lattice[0] = chain.JpsiTau_hammer_ebe_b1_up_lattice[0]
#                    out.hammer_ebe_b1_down_lattice[0] = chain.JpsiTau_hammer_ebe_b1_down_lattice[0]
#                    out.hammer_ebe_b2_up_lattice[0] = chain.JpsiTau_hammer_ebe_b2_up_lattice[0]
#                    out.hammer_ebe_b2_down_lattice[0] = chain.JpsiTau_hammer_ebe_b2_down_lattice[0]
#
#                    out.hammer_ebe_c1_up_lattice[0] = chain.JpsiTau_hammer_ebe_c1_up_lattice[0]
#                    out.hammer_ebe_c1_down_lattice[0] = chain.JpsiTau_hammer_ebe_c1_down_lattice[0]
#                    out.hammer_ebe_c2_up_lattice[0] = chain.JpsiTau_hammer_ebe_c2_up_lattice[0]
#                    out.hammer_ebe_c2_down_lattice[0] = chain.JpsiTau_hammer_ebe_c2_down_lattice[0]
#
#                    out.hammer_ebe_d0_up_lattice[0] = chain.JpsiTau_hammer_ebe_d0_up_lattice[0]
#                    out.hammer_ebe_d0_down_lattice[0] = chain.JpsiTau_hammer_ebe_d0_down_lattice[0]
#                    out.hammer_ebe_d1_up_lattice[0] = chain.JpsiTau_hammer_ebe_d1_up_lattice[0]
#                    out.hammer_ebe_d1_down_lattice[0] = chain.JpsiTau_hammer_ebe_d1_down_lattice[0]
#                    out.hammer_ebe_d2_up_lattice[0] = chain.JpsiTau_hammer_ebe_d2_up_lattice[0]
#                    out.hammer_ebe_d2_down_lattice[0] = chain.JpsiTau_hammer_ebe_d2_down_lattice[0]




#                    for iham in range(len(chain.JpsiTau_hammer_ebe_toy[0])):
            
#                        out.hammer_ebe_toy[iham] = ROOT.Double(chain.JpsiTau_hammer_ebe_toy[0][iham])

                else:
                    continue
#                    out.hammer_ebe_lattice[0] = -1

#                if chain.JpsiTau_hammer_ebe[0]==-1: continue
            
                out.filt.Fill(7)



#            out.puweight[0] = ROOT.Double(putool.getWeight(chain.nPuVtxTrue[0]))
#            out.puweight_up[0] = ROOT.Double(putool_up.getWeight(chain.nPuVtxTrue[0]))
#            out.puweight_down[0] = ROOT.Double(putool_down.getWeight(chain.nPuVtxTrue[0]))

            pubin = max(1, min(puhist.GetNbinsX(), puhist.FindBin(chain.nPuVtxTrue[0])))

            if chain.nPuVtxTrue[0] < puhist.GetNbinsX():
                out.puweight[0] = puhist.GetBinContent(pubin)
                out.puweight_up[0] = puhist_up.GetBinContent(pubin)
                out.puweight_down[0] = puhist_down.GetBinContent(pubin)
            else:
                print('This cannot happen !!!!')
                out.puweight[0] = 1
                out.puweight_up[0] = 1
                out.puweight_down[0] = 1
                


#            print(evt, 'pubin=', pubin, out.puweight[0],out.puweight_up[0],out.puweight_down[0])



        if options.type =='signal':
            #GEN weight for the signal 
#            print ' NEW Event '
            final_daus = [] 

            for igen in range(len(chain.genParticle_pdgId)):
                if abs(chain.genParticle_pdgId[igen])==541 or abs(chain.genParticle_pdgId[igen])==543  and chain.genParticle_status[igen]==2:
#                    print 'Event found a Mother ---- of status ', chain.genParticle_status[igen] , ' pdgId ' ,  chain.genParticle_pdgId[igen]
#                    print ' Bc daughers are : '
                    final_daus = []
                    for idau in range(len(chain.genParticle_dau[igen])):
#                        print '   -> ', (chain.genParticle_dau[igen][idau])
                        if abs(chain.genParticle_dau[igen][idau]) not in [22, 541]: #avoid propagation and photon emission
                            final_daus.append(abs(chain.genParticle_dau[igen][idau]))
                        dau_number = [x for x in range(len(chain.genParticle_pdgId)) if  chain.genParticle_pdgId[x] == chain.genParticle_dau[igen][idau]][0]
                        #for iidau in range(len(chain.genParticle_dau[dau_number])):
                            #print ' nephew  tentative --> ', (chain.genParticle_dau[dau_number][iidau])
#                    print "Bc final_daus: ", final_daus

            final_daus.sort()
#            print " Event final_daus: ", final_daus
            if final_daus==[13, 14, 443]:
                out.gen_sig_decay[0] = 0 ##Bc -> Jpsi mu nu  
            elif final_daus==[13, 14, 100443]:
                out.gen_sig_decay[0] = 1  ##Bc -> psi(2s) mu nu 
            elif final_daus==[13, 14, 10441]:
                out.gen_sig_decay[0] = 2##Bc -> chic0 mu nu
            elif final_daus==[13, 14, 20443]:
                out.gen_sig_decay[0] = 3 ##Bc -> chic1 mu nu              
            elif final_daus==[13, 14, 445]:
                out.gen_sig_decay[0] = 4 ##Bc -> chic2 mu nu             
            elif final_daus==[13, 14, 10443]:
                out.gen_sig_decay[0] = 5 ##Bc -> hc mu nu               
            elif final_daus==[15, 16, 443]:
                out.gen_sig_decay[0] = 6 ##Bc -> Jpsi tau nu            
            elif final_daus==[15, 16, 100443]:
                out.gen_sig_decay[0] = 7 ##Bc -> Psi(2s) tau nu            
            elif final_daus==[211, 443]:
                out.gen_sig_decay[0] = 8 ##Bc -> Jpsi pi
            elif final_daus==[211, 211, 211, 443]:
                out.gen_sig_decay[0] = 9 ##Bc -> Jpsi  3pi  
            elif final_daus==[431, 443]:
                out.gen_sig_decay[0] = 10 ##Bc -> Jpsi Ds+
            elif final_daus==[211, 211, 211, 211, 211, 443]:
                out.gen_sig_decay[0] = 11 ##Bc -> Jpsi  5pi  
            elif final_daus==[15, 16, 10441]:
                out.gen_sig_decay[0] = 12 ##Bc -> chic0 tau nu                                                                                                                              
            elif final_daus==[15, 16, 20443]:
                out.gen_sig_decay[0] = 13 ##Bc -> chic1 tau nu                                                                                                 
            elif final_daus==[15, 16, 445]:
                out.gen_sig_decay[0] = 14 ##Bc -> chic2 tau nu                                                                                                             
            elif final_daus==[15, 16, 10443]:
                out.gen_sig_decay[0] = 15 ##Bc -> hc tau nu 
            elif final_daus==[311, 413, 443]:
                out.gen_sig_decay[0] = 16 ##[B_c+ -> K0 (D*+ -> (D+ -> anti-K0 mu+ nu_mu) pi0) (J/psi -> mu+ mu-)]cc
            elif final_daus==[443, 10433] or final_daus==[443, 20433]:
                out.gen_sig_decay[0] = 17 ### [B_c+ -> (D_s1+ -> (D_s*+ -> (D_s+ -> eta mu+ nu_mu) gamma) pi0) (J/psi -> mu+ mu-)]cc 
            elif final_daus==[321, 423, 443]:
                out.gen_sig_decay[0] = 18 ##Bc -> K+ D*(2007)0 Jpsi
            elif final_daus==[321, 421, 443]:
                out.gen_sig_decay[0] = 19 ##Bc -> K+ D0 Jpsi
            elif final_daus==[433, 443]:
                out.gen_sig_decay[0] = 20 ##Bc -> Jpsi Ds*+     
            elif final_daus==[413, 443]:
                out.gen_sig_decay[0] = 21 ##Bc -> Jpsi D*(2010)+
            elif final_daus==[211, 321, 321, 443]:
                out.gen_sig_decay[0] = 22 ##Bc -> Jpsi pi K K       
            else :
                #print "final_daus not recognized: ", final_daus            
                out.gen_sig_decay[0] = 23 ##Bc no more decays 
    
        out.filt.Fill(8)
        evtid += 1
        out.tree.Fill()
    
    out.filt.Fill(9)

print(Nevt, 'evt processed.', evtid, 'evt has matching')

from collections import Counter
x = Counter(dicts)


if options.create:
    dict2save = {}

id_counter = 0

for k,l in sorted([(j,i) for i,j in x.items()], reverse=True):
    print str(id_counter).ljust(10), l.ljust(50), str(k).ljust(10), 'counts out of', str(dict_counter).ljust(10), ' frac = ({0:2f}'.format(float(k)/float(dict_counter)), ')'


    if options.create:
        dict2save[l] = {'count':k, 'frac':float(k)/float(dict_counter), 'id':id_counter}
    else:
        if l in idtable:
            print 'id =', idtable[l]


    id_counter += 1


#for key, idt in dicts.iteritems():
#    print key.ljust(40), dicts[key], 'counts'

if options.create:

    with open('json_jpsitanu.json', 'w') as outfile:
        json.dump(dict2save, outfile)

out.endJob()

