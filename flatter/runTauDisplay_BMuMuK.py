
#import os, math, sys
#from ROOT import TFile, TH1F, gROOT, TTree, Double, TChain, TLorentzVector, TVector3
#import numpy as num

from TreeProducerBMuMuK import *
from correction.PileupTool import *
from correction.ScaleFactorMuonTool import *
from DeltaR import deltaR
import copy
import random
from scipy.constants import c as speed_of_light
import numpy as np
from particle import Particle
import json

# from https://github.com/scikit-hep/particle/
# installed via 
# >>> pip install particle --user
### 
#from particle import Particle

#gROOT.SetBatch(True)

# Bc lifetime
# https://pdglive.lbl.gov/DataBlock.action?node=S091T&home=MXXX049
ctau_pdg    = 0.510e-12 * speed_of_light * 1000. # in mm
ctau_actual = 0.507e-12 * speed_of_light * 1000. # in mm
ctau_up     = (0.510+0.009)*1e-12 * speed_of_light * 1000. # in mm
ctau_down   = (0.510-0.009)*1e-12 * speed_of_light * 1000. # in mm


def weight_to_new_ctau(old_ctau, new_ctau, ct):
    '''
    Returns an event weight based on the ratio of the normalised lifetime distributions.
    old_ctau: ctau used for the sample production
    new_ctau: target ctau
    ct      : per-event lifetime
    '''
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

#def returnName(pid):
#    
#    addstr = None
#    
#    try:
#        addstr = Particle.from_pdgid(pid).name    
#    except:
#        addstr = 'None'
#
#    return addstr
    

def returnName(pid):
    
    addstr = None
    
    try:
        addstr = Particle.from_pdgid(pid).name    
    except:
        addstr = 'None'

    return addstr

def build_uncertainties(sf):
    keys = ["nominal"]
    keys += ["up_syst", "down_syst"] if "syst" in sf else []
    keys += ["up_stat", "down_stat"] if "stat" in sf else []
    
    content = [sf["value"]]
    content += [sf["value"] + sf["syst"], sf["value"] - sf["syst"]] if "syst" in sf else []
    content += [sf["value"] + sf["error"], sf["value"] - sf["error"]] if "error" in sf else []
    
    return ({
        "keys": keys,
        "content": content
    })


def parse_str(key, prefix="abseta:"):
    if not key.startswith(prefix + "["):
        raise ValueError
    lo, hi = map(float, key[len(prefix + "["):-1].split(","))
    return lo, hi
    

def edges_pts(sf):
    # Could happen that higher pt bin edge comes
    # lexicographically before lower one, so sort bins first
    sf_sorted_data = {}
    sf_sorted_hi = {}
    for binstr, data in sf.items():
        if not binstr.startswith("pt:["):
            raise ValueError
        lo, hi = map(float, binstr[len("pt:["):-1].split(","))
        sf_sorted_data[lo] = data
        sf_sorted_hi[lo] = hi

    edges = []
    content = []
    for i in sorted(sf_sorted_data):
        lo = i
        data = sf_sorted_data[i]
        if len(edges) == 0:
            edges.append(lo)
        if edges[-1] != lo:
            raise ValueError
        edges.append(sf_sorted_hi[lo])
        content.append(build_uncertainties(data))
    
    # print " Pt edges", edges, "content", content 
    return ({
        #"nodetype": "binning",
        "edges": edges,
        "content": content,
    })


def build_etas(sf):
    bins = [parse_str(s, "abseta:") for s in sf]
    edges = sorted(set(edge for bin in bins for edge in bin))
      
    content = [None] * (len(edges) - 1)
    for s, data in sf.items():
        lo, hi = parse_str(s, "abseta:")
        found = False
        for i, bin in enumerate(bins):
            if bin[0] >= lo and bin[1] <= hi:
                content[i] = build_pts(data)
                found = True

    print " Eta edges", edges, "content", content
    return ({
    #    "nodetype": "binning",
        "edges": edges,
        "content": content,
    })
    #return 1

def find_etas_pt(sf, eta,pt):
    bins = [parse_str(s, "abseta:") for s in sf]
    print "bins ", bins
    edges_eta = sorted(set(edge for bin in bins for edge in bin))
    print " Eta edges ", edges_eta
    eta_limit={} 
    for i in range(len(edges_eta)):
        print "edges_eta [i] ",i," ",edges_eta [i] 
        if edges_eta[i] > eta:  
            eta_limit=sorted({edges_eta[i-1],edges_eta[i]}) 
            break
    bins_pt = [parse_str(s, "pt:") for s in sf["abseta:[{:.2f},{:.2f}]".format(eta_limit[0], eta_limit[1])]] 
    edges_pt = sorted(set(edge_pt for bin_pt in bins_pt for edge_pt in bin_pt))
    print " Pt edges", edges_pt
    pt_limit={}
    for i in range(len(edges_pt)):
        #print "edges_pt [i] ",i," ",edges_pt [i]                                                                                                                                                                                                                                           
        if edges_pt[i] > pt:
            pt_limit=sorted({edges_pt[i-1],edges_pt[i]})
            break
    if pt_limit=={0,0}: pt_limit=sorted({edges_pt[-2],edges_pt[-1]})
    print "eta_limit ", eta_limit, " ;pt_limit ", pt_limit 
    return (eta_limit,pt_limit)

def find_eta_pt_bins(sf):
    bins = [parse_str(s, "abseta:") for s in sf]
    print "bins ", bins
    edges_eta = sorted(set(edge for bin in bins for edge in bin))
    print " Eta edges ", edges_eta
    bins_pt = [parse_str(s, "pt:") for s in sf["abseta:[{:.2f},{:.2f}]".format(edges_eta[0], edges_eta[1])]]
    edges_pt = sorted(set(edge_pt for bin_pt in bins_pt for edge_pt in bin_pt))
    print " Pt edges", edges_pt
    return ({'eta': edges_eta,'pt': edges_pt})

def find_eta_boundaries(edges, eta):
    eta_string=""
    for i in range(len(edges['eta'])):
        print "edges_eta [i] ",i," ",edges['eta'] [i]
        if edges['eta'][i] > eta:
            eta_string="abseta:[{:.2f},{:.2f}]".format(edges['eta'][i-1],edges['eta'][i])
            break;
    if eta_string=="": eta_string="abseta:[{:.2f},{:.2f}]".format(edges['eta'][-2],edges['eta'][-1])
    return eta_string
def find_pt_boundaries(edges, pt):
    string=""
    for i in range(len(edges['pt'])):
        print "edges_pt [i] ",i," ",edges['pt'] [i]
        if edges['pt'][0]> pt: 
            ValueError("pt value too small,no sf available")
        if edges['pt'][i] > pt:
            string="pt:[{:.2f},{:.2f}]".format(edges['pt'][i-1],edges['pt'][i])
            break;
    if string=="": string="pt:[{:.2f},{:.2f}]".format(edges['pt'][-2],edges['pt'][-1])
    return string


from optparse import OptionParser, OptionValueError
usage = "usage: python runTauDisplay_BsTauTau.py"
parser = OptionParser(usage)

parser.add_option("-o", "--out", default='Myroot.root', type="string", help="output filename", dest="out")
parser.add_option("-p", "--priority", default='pt', type="string", help="priority", dest="priority")
parser.add_option("-t", "--type", default='bg', type="string", help="type", dest="type")

parser.add_option("-y", "--year", default='2018', type="string", help="year", dest="year")
#parser.add_option("-f", "--file", default='root://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas/BcToJPsiMuMu_Legacy_2018_20210430/BcToJPsiMuMu_inclusive_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/210430_140808/0000/flatTuple_100.root', type="string", help="file", dest="file")
#parser.add_option("-f", "--file", default='root://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas/JpsiX_Legacy_2018_20210826/JpsiX_MuMu_J_211119/cgalloni-Autumn18_10_2_9_miniAOD-39a089a8e7301f392b8b059e430f83ef/210826_070112/0000/flatTuple_5.root', type="string", help="file", dest="file")   
parser.add_option("-f", "--file", default="root://storage01.lcg.cscs.ch/pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas//DIGI_JpsiK_multiple_2018_20210907/DIGI_JpsiX_MuMu_20210905/cgalloni-UL18_MINIAOD_v1_noDuplCheck-cb0e7829b1b4ca4eee675686c1769096/210907_055844/0000/flatTuple_1.root", type="string", help="file", dest="file") 
parser.add_option('-c', '--create', action="store_true", default=False, dest='create')

(options, args) = parser.parse_args()

print(options)

if not options.create:
    json_open = open('json_jpsik.json', 'r')
    idtable = json.load(json_open)
    print('JSON decay table is read !')

    print idtable 


#with open('correction/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.json') as f:
#    sf_file = json.load(f)
#print "loaded json SF:  ", sf_file
#build_etas(sf_file["NUM_TrackerMuons_DEN_genTracks"]["abseta_pt"])
#edges = find_eta_pt_bins(sf_file["NUM_TrackerMuons_DEN_genTracks"]["abseta_pt"]);#,1.5,5)
#print "edges" , edges
#print "edges[eta]" , edges['eta']
#print "edges[pt]" , edges['pt']
#print "eta label ", find_eta_boundaries(edges, 1.5)
#print "pt label ", find_pt_boundaries(edges,5)
#print "SF ", sf_file["NUM_TrackerMuons_DEN_genTracks"]["abseta_pt"][find_eta_boundaries(edges,1.5)] [find_pt_boundaries(edges,55)]

out = TreeProducerBcJpsiKNu(options.out, options.type)

files = options.file.split(',')

chain = ROOT.TChain('ntuplizer/tree', 'tree')

hist_hammer = None

for inputfile in files:
    print('Adding ...', inputfile)
    chain.AddFile(inputfile)

    _file = TFile.Open(inputfile)
    _hist = _file.Get('ntuplizer/cutflow')

    if out.hist == None:
        print('cutflow histo created')
        out.hist = copy.deepcopy(_hist)
    else:
        print('cutflow histo added')
        out.hist.Add(copy.deepcopy(_hist))

    if hist_hammer==None:
        hist_hammer = _file.Get('ntuplizer/hammer_width')

chain.SetBranchStatus('*', 0)


chain.SetBranchStatus('JpsiK_B_pt*', 1)
chain.SetBranchStatus('*simple*', 1)
chain.SetBranchStatus('JpsiK_B_eta', 1)
chain.SetBranchStatus('JpsiK_B_phi', 1)
chain.SetBranchStatus('JpsiK_B_mass', 1)
chain.SetBranchStatus('JpsiK_B_alpha', 1)
chain.SetBranchStatus('JpsiK_B_fls3d', 1)
chain.SetBranchStatus('JpsiK_B_fl3d', 1)
chain.SetBranchStatus('JpsiK_B_vprob', 1)


chain.SetBranchStatus('JpsiK_B_q2', 1)
chain.SetBranchStatus('JpsiK_B_mm2', 1)
chain.SetBranchStatus('JpsiK_B_ptmiss', 1)
chain.SetBranchStatus('JpsiK_B_Es', 1)

chain.SetBranchStatus('JpsiK_B_pvip', 1)
chain.SetBranchStatus('JpsiK_B_pvips', 1)
chain.SetBranchStatus('JpsiK_B_lips', 1)
chain.SetBranchStatus('JpsiK_B_mindoca', 1)
chain.SetBranchStatus('JpsiK_B_vx', 1)
chain.SetBranchStatus('JpsiK_B_vy', 1)
chain.SetBranchStatus('JpsiK_B_vz', 1)
chain.SetBranchStatus('JpsiK_nch', 1)
chain.SetBranchStatus('JpsiK_nch_before', 1)

chain.SetBranchStatus('JpsiK_PV_vx', 1)
chain.SetBranchStatus('JpsiK_PV_vy', 1)
chain.SetBranchStatus('JpsiK_PV_vz', 1)

chain.SetBranchStatus('JpsiK_bbPV_vx', 1)
chain.SetBranchStatus('JpsiK_bbPV_vy', 1)
chain.SetBranchStatus('JpsiK_bbPV_vz', 1)


chain.SetBranchStatus('JpsiK_pi_pt', 1)
chain.SetBranchStatus('JpsiK_pi_eta', 1)
chain.SetBranchStatus('JpsiK_pi_phi', 1)
chain.SetBranchStatus('JpsiK_pi_mass', 1)
chain.SetBranchStatus('JpsiK_pi_q', 1)
chain.SetBranchStatus('JpsiK_pi_vprob', 1)
chain.SetBranchStatus('JpsiK_pi_fls3d*', 1)
chain.SetBranchStatus('JpsiK_pi_fl3d*', 1)
chain.SetBranchStatus('JpsiK_pi_trigMatch',1)

#chain.SetBranchStatus('JpsiK_pi_max_dr_3prong', 1)
chain.SetBranchStatus('JpsiK_pi_lips', 1)
chain.SetBranchStatus('JpsiK_pi_pvips', 1)
chain.SetBranchStatus('JpsiK_pi_refit*', 1)
chain.SetBranchStatus('JpsiK_pi_pi*', 1)
chain.SetBranchStatus('JpsiK_pi_v*', 1)
chain.SetBranchStatus('JpsiK_Jpsi_*', 1)
chain.SetBranchStatus('JpsiK_mu*', 1)
chain.SetBranchStatus('JpsiK_jpsi_pi_alpha', 1)
chain.SetBranchStatus('JpsiK_perEVT_*', 1)


chain.SetBranchStatus('JpsiK_nCandidates', 1)
chain.SetBranchStatus('PV_N', 1)
chain.SetBranchStatus('EVENT_*', 1)


putool = None

if options.type!='data':
    chain.SetBranchStatus('*genParticle*', 1)
    chain.SetBranchStatus('JpsiK_*gen*', 1)
    chain.SetBranchStatus('JpsiK_genPV_vx', 1)
    chain.SetBranchStatus('JpsiK_genPV_vy', 1)
    chain.SetBranchStatus('JpsiK_genPV_vz', 1)
    chain.SetBranchStatus('JpsiK_genSV_vx', 1)
    chain.SetBranchStatus('JpsiK_genSV_vy', 1)
    chain.SetBranchStatus('JpsiK_genSV_vz', 1)
    chain.SetBranchStatus('JpsiK_hammer_*', 1)
    chain.SetBranchStatus('nPuVtxTrue', 1)
    chain.SetBranchStatus('JpsiK_st_gentau*',1)
    chain.SetBranchStatus('JpsiK_st_n_occurance',1)
    chain.SetBranchStatus('JpsiK_st_decayid',1)

    putool = PileupWeightTool(options.year, 'central')
    SF_ID = ScaleFactorMuonTool('central', fileName='Efficiency_muon_trackerMuon_Run2018_UL_ID.json', keyName='NUM_LooseID_DEN_TrackerMuons');
    SF_Reco = ScaleFactorMuonTool('central', fileName='Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.json', keyName='NUM_TrackerMuons_DEN_genTracks');
    #print SF_ID.getSF(1.5, 55)['value']
    #print SF_ID.getSF(1.1, 25)['value']
    #print SF_Reco.getSF(1.5, 55)['value']

if options.type=='bg':
    chain.SetBranchStatus('genWeightBkgB',1)


Nevt = chain.GetEntries()

print('Total Number of events = ', Nevt)
evtid = 0


mk = 0.493677
mp = 0.139571


dicts = {}
dict_counter = 0


for evt in xrange(Nevt):
    chain.GetEntry(evt)
#    if evt>100 : break

    if evt%10000==0: print('{0:.2f}'.format(ROOT.Double(evt)/ROOT.Double(Nevt)*100.), '% processed')

    out.multi.Fill(len(chain.JpsiK_pi_pt))

    if len(chain.JpsiK_pi_pt)==0: 
#        print('This is not possible !!')
        continue


    tlv_jpsi = ROOT.TLorentzVector()
    tlv_jpsi.SetPtEtaPhiM(chain.JpsiK_Jpsi_pt, chain.JpsiK_Jpsi_eta, chain.JpsiK_Jpsi_phi, chain.JpsiK_Jpsi_mass)
    
    
    tlv_mu1 = ROOT.TLorentzVector()
    tlv_mu1.SetPtEtaPhiM(chain.JpsiK_mu1_pt, chain.JpsiK_mu1_eta, chain.JpsiK_mu1_phi, chain.JpsiK_mu1_mass)

    tlv_mu2 = ROOT.TLorentzVector()
    tlv_mu2.SetPtEtaPhiM(chain.JpsiK_mu2_pt, chain.JpsiK_mu2_eta, chain.JpsiK_mu2_phi, chain.JpsiK_mu2_mass)


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
    
    if options.priority=='pt':
            
        for itau in range(len(chain.JpsiK_pi_pt)):
            if chain.JpsiK_B_vprob[itau] < 0.1: continue
#            if chain.JpsiK_pi_fls3d[itau] < 3.: continue
#            if chain.JpsiK_pi_mass[itau] > 1.7: continue
            if bool(chain.JpsiK_pi_trigMatch[itau])==False: 
#                print 'trigger matching was not satisifed ...'
                continue
            
            # you can add tau mass cuts here

            tindex_ = itau
            break

#    if options.priority=='pt':
#        tindex_ = 0
            
    if options.priority!='multiple' and tindex_ == -1: continue


    if len(chain.JpsiK_B_pt) != len(chain.JpsiK_pi_pt): continue



    ### ctau calculation ! 

    if options.type == 'signal':
        lxyz = math.sqrt((chain.JpsiK_genPV_vx - chain.JpsiK_genSV_vx)**2 + (chain.JpsiK_genPV_vy - chain.JpsiK_genSV_vy)**2 + (chain.JpsiK_genPV_vz - chain.JpsiK_genSV_vz)**2)


        tlv_B = ROOT.TLorentzVector()
        tlv_B.SetPtEtaPhiM(chain.JpsiK_B_pt_gen, chain.JpsiK_B_eta_gen, chain.JpsiK_B_phi_gen, chain.JpsiK_B_mass_gen)
        
#        print chain.JpsiK_B_pt_gen, chain.JpsiK_B_eta_gen, chain.JpsiK_B_phi_gen, chain.JpsiK_B_mass_gen
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
        
        weight_central = weight_to_new_ctau(ctau_actual, ctau_pdg , ct*10.)
        weight_up = weight_to_new_ctau(ctau_actual, ctau_up , ct*10.)
        weight_down = weight_to_new_ctau(ctau_actual, ctau_down , ct*10.)
        
        out.weight_ctau[0] = weight_central
        out.weight_ctau_up[0] = weight_up
        out.weight_ctau_down[0] = weight_down

      








    for tindex in range(len(chain.JpsiK_pi_pt)):


#        if abs(chain.JpsiK_pi_q[tindex])!=1: continue ## Charge 1 

        if options.priority in ['pt']:
            if tindex != tindex_: continue


        isRight_3prong = False
        isRight_3prong_pi0 = False

        
          




        if options.priority=='multiple':

            if tindex > 5: break
            # > 10 ---> 0.97316682 efficiency of catching the siganl
            # > 5 ---> 0.95029727 efficiency of catching the siganl


     

            

        out.mu1_pt[0] = tlv_mu1.Pt()
        out.mu1_eta[0] = tlv_mu1.Eta()
        out.mu1_phi[0] = tlv_mu1.Phi()
        out.mu1_mass[0] = tlv_mu1.M()
        out.mu1_q[0] = chain.JpsiK_mu1_q
        out.mu1_isLoose[0] = chain.JpsiK_mu1_isLoose
        out.mu1_isTight[0] = chain.JpsiK_mu1_isTight
        out.mu1_isPF[0] = chain.JpsiK_mu1_isPF
        out.mu1_isGlobal[0] = chain.JpsiK_mu1_isGlobal
        out.mu1_isTracker[0] = chain.JpsiK_mu1_isTracker
        out.mu1_isSoft[0] = chain.JpsiK_mu1_isSoft
        out.mu1_dbiso[0] = chain.JpsiK_mu1_dbiso
        if options.type!='data':
            out.mu1_SFReco[0]=SF_Reco.getSF(math.fabs(tlv_mu1.Eta()),tlv_mu1.Pt())['value']
            out.mu1_SFID[0]=SF_ID.getSF(math.fabs(tlv_mu1.Eta()),tlv_mu1.Pt())['value']
        else :
            out.mu1_SFReco[0]=1                                                                                                                              
            out.mu1_SFID[0]=1       

        out.mu2_pt[0] = tlv_mu2.Pt()
        out.mu2_eta[0] = tlv_mu2.Eta()
        out.mu2_phi[0] = tlv_mu2.Phi()
        out.mu2_mass[0] = tlv_mu2.M()
        out.mu2_q[0] = chain.JpsiK_mu2_q
        out.mu2_isLoose[0] = chain.JpsiK_mu2_isLoose
        out.mu2_isTight[0] = chain.JpsiK_mu2_isTight
        out.mu2_isPF[0] = chain.JpsiK_mu2_isPF
        out.mu2_isGlobal[0] = chain.JpsiK_mu2_isGlobal
        out.mu2_isTracker[0] = chain.JpsiK_mu2_isTracker
        out.mu2_isSoft[0] = chain.JpsiK_mu2_isSoft
        out.mu2_dbiso[0] = chain.JpsiK_mu2_dbiso
        if options.type!='data':
            out.mu2_SFReco[0]=SF_Reco.getSF(math.fabs(tlv_mu2.Eta()),tlv_mu2.Pt())['value']
            out.mu2_SFID[0]=SF_ID.getSF(math.fabs(tlv_mu2.Eta()),tlv_mu2.Pt())['value']
        else :                                                                                                                                                                                                                                                                             
            out.mu2_SFReco[0]=1
            out.mu2_SFID[0]=1


        out.mu1_mu2_dr[0] = deltaR(tlv_mu1.Eta(),tlv_mu1.Phi(),tlv_mu2.Eta(),tlv_mu2.Phi())

        out.pi_pt[0] = chain.JpsiK_pi_pt[tindex]
        out.pi_eta[0] = chain.JpsiK_pi_eta[tindex]
        out.pi_phi[0] = chain.JpsiK_pi_phi[tindex]
        out.pi_mass[0] = chain.JpsiK_pi_mass[tindex]
        out.pi_q[0] = chain.JpsiK_pi_q[tindex]
       
        out.pi_dr_jpsi[0] = deltaR(chain.JpsiK_pi_eta[tindex], chain.JpsiK_pi_phi[tindex], 
                                    chain.JpsiK_Jpsi_eta, chain.JpsiK_Jpsi_phi)

        out.evtid[0] = evtid


      

#        out.pi_refit_vx[0] = chain.JpsiK_pi_refit_vx[tindex]
#        out.pi_refit_vy[0] = chain.JpsiK_pi_refit_vy[tindex]
#        out.pi_refit_vz[0] = chain.JpsiK_pi_refit_vz[tindex]


        if len(chain.JpsiK_B_pt)!=0:
#            print len(chain.JpsiK_B_pt), len(chain.JpsiK_pi_pt), tindex
            out.b_pt[0] = chain.JpsiK_B_pt[tindex]
            out.b_eta[0] = chain.JpsiK_B_eta[tindex]
            out.b_phi[0] = chain.JpsiK_B_phi[tindex]
            out.b_mass[0] = chain.JpsiK_B_mass[tindex]
            out.b_alpha[0] = chain.JpsiK_B_alpha[tindex]
            out.b_fls3d[0] = chain.JpsiK_B_fls3d[tindex]
            out.b_fl3d[0] = chain.JpsiK_B_fl3d[tindex]
            out.b_vprob[0] = chain.JpsiK_B_vprob[tindex]
            out.b_pvip[0] = chain.JpsiK_B_pvip[tindex]
            out.b_pvips[0] = chain.JpsiK_B_pvips[tindex]
            out.b_lips[0] = chain.JpsiK_B_lips[tindex]
            
            out.b_vz[0] = chain.JpsiK_B_vz[tindex]



        
       

        tlv_k = ROOT.TLorentzVector()
        k_mass= 0.493677
        tlv_k.SetPtEtaPhiM(chain.JpsiK_pi_pt[tindex], chain.JpsiK_pi_eta[tindex], chain.JpsiK_pi_phi[tindex], k_mass)
       
        tlv_B_k = ROOT.TLorentzVector()
       
        tlv_B_k = tlv_k + tlv_mu1 + tlv_mu2
       
        out.k_pt[0]=tlv_k.Pt()
        out.B_k_mass[0]=tlv_B_k.M()
        out.B_k_pt[0]=tlv_B_k.Pt()
     
        # per-event quantity
        out.bbpv_vz[0] = chain.JpsiK_bbPV_vz
        out.pv_vz[0] = chain.JpsiK_PV_vz
        out.ncand[0] = chain.JpsiK_nCandidates
        #out.nch[0] = chain.JpsiK_nch
    
        

        out.jpsi_pt[0] = chain.JpsiK_Jpsi_pt
        out.jpsi_eta[0] = chain.JpsiK_Jpsi_eta
        out.jpsi_phi[0] = chain.JpsiK_Jpsi_phi
        out.jpsi_mass[0] = chain.JpsiK_Jpsi_mass
        out.jpsi_vprob[0] = chain.JpsiK_Jpsi_vprob
        out.jpsi_fls3d[0] = chain.JpsiK_Jpsi_fls3d
        out.jpsi_fl3d[0] = chain.JpsiK_Jpsi_fl3d
        
        out.npv[0] = chain.PV_N



        out.evt[0] = chain.EVENT_event
        out.lumi[0] = chain.EVENT_lumiBlock
        out.run[0] = chain.EVENT_run

        init_float=[-99 for i in range(12)]
        init_bool=[0 for i in range(12)]
        #out.JpsiK_st_doca3d.clear()
        #print " type of chain.JpsiK_st_doca3d) ", type(chain.JpsiK_st_doca3d)
        #out.JpsiK_st_doca3d=init_float
        #out.JpsiK_st_doca3d=chain.JpsiK_st_doca3d
  

        if options.type=='bg':

            decayTable = []


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


#                print '\t ==========>', addstr
                decayTable.append(addstr)




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




#####        flag_veto = False
#####
#####        if options.type == 'bg' and bool(chain.JpsiK_tau_isRight[tindex]):
#####            print 'TRACE', '-'*80
#####
#####            for igen in range(len(chain.genParticle_pdgId)):
#####                if abs(chain.genParticle_pdgId[igen])>=511 and abs(chain.genParticle_pdgId[igen])<=545:
#####
#####                    print 'TRACE', chain.genParticle_pdgId[igen]
#####                
#####                    lp = []
#####
#####                    for idau in range(len(chain.genParticle_dau[igen])):
#####                        print 'TRACE -->', chain.genParticle_dau[igen][idau]
#####                        lp.append(abs(chain.genParticle_dau[igen][idau]))
#####                    
#####                    if 15 in lp and 443 in lp and abs(chain.genParticle_pdgId[igen])==541:
#####                        print 'THIS IS THE SIGNAL !!!'
#####                        flag_veto = True
#####
#####        if options.type=='bg':
#####            out.isveto[0] = flag_veto
#####
#####            
#####
#####        if options.type=='bg':
#####
#####            decayTable = []
#####
#####            flag_jpsi = False
#####            
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
#####                    print '  '*int(chain.genParticle_layers[igen][ipdg]), 'pdg  = ', chain.genParticle_pdgs[igen][ipdg], '(',  returnName(chain.genParticle_pdgs[igen][ipdg]) , '), (pt, eta, phi) = ', '({0:.2f}'.format(chain.genParticle_ppt[igen][ipdg]), '{0:.2f}'.format(chain.genParticle_peta[igen][ipdg]), '{0:.2f}'.format(chain.genParticle_pphi[igen][ipdg]), ')'
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
#####
#####                            
#####
#####
#####
#####            for igen in range(len(chain.genParticle_pdgs)):
#####
######                print ''
######                print 'gen = ', igen
######                print ''
#####
#####                flag_mu1 = False
#####                flag_mu2 = False
#####
#####                addstr = returnName(chain.genParticle_pdgs[igen][0])
#####
#####                for ipdg in range(len(chain.genParticle_pdgs[igen])):
#####            
#####                    if chain.genParticle_layers[igen][ipdg]==1:
#####                        addstr += '_' + returnName(chain.genParticle_pdgs[igen][ipdg])
#####
#####
######                    print '  '*int(chain.genParticle_layers[igen][ipdg]), 'pdg  = ', chain.genParticle_pdgs[igen][ipdg], '(',  returnName(chain.genParticle_pdgs[igen][ipdg]) , '), (pt, eta, phi) = ', '({0:.2f}'.format(chain.genParticle_ppt[igen][ipdg]), '{0:.2f}'.format(chain.genParticle_peta[igen][ipdg]), '{0:.2f}'.format(chain.genParticle_pphi[igen][ipdg]), ')'
#####
#####
#####                    if abs(chain.genParticle_pdgs[igen][ipdg]) == 13:
#####                        dr1 = deltaR(tlv_mu1.Eta(), tlv_mu1.Phi(), chain.genParticle_peta[igen][ipdg], chain.genParticle_pphi[igen][ipdg])
#####                        dr2 = deltaR(tlv_mu2.Eta(), tlv_mu2.Phi(), chain.genParticle_peta[igen][ipdg], chain.genParticle_pphi[igen][ipdg])
#####                        
#####                        if dr1 < 0.05:
#####                            flag_mu1 = True
######                            print '*'
#####                            
#####                        if dr2 < 0.05:
#####                            flag_mu2 = True
######                            print '*'
#####
#####                if not (flag_mu1 and flag_mu2): continue
#####
#####                flag_jpsi = True 
#####                decayTable.append(addstr)
#####
#####
#####                for itlv in [tlv1, tlv2, tlv3]:
#####
#####                    flag_pi = False
#####
#####                    for ipdg in range(len(chain.genParticle_pdgs[igen])):            
#####
#####                        _dr = deltaR(itlv.Eta(), itlv.Phi(), chain.genParticle_peta[igen][ipdg], chain.genParticle_pphi[igen][ipdg])
#####                        
#####                        if _dr < 0.1:
#####                            flag_pi = True
######                            print '*'
#####                                
#####                    if flag_pi:
#####                        decayTable.append(addstr)
#####
#####                    else:
#####
#####                        flag_pi_other = False
#####
#####                        for igen2 in range(len(chain.genParticle_pdgs)):
#####                            if igen2==igen: continue
#####
#####                            addstr = returnName(chain.genParticle_pdgs[igen2][0])
#####                        
#####                            for ipdg in range(len(chain.genParticle_pdgs[igen2])):
#####                            
#####                                if chain.genParticle_layers[igen2][ipdg]==1:
#####                                    addstr += '_' + returnName(chain.genParticle_pdgs[igen2][ipdg])
#####
#####
#####                            _dr = deltaR(itlv.Eta(), itlv.Phi(), chain.genParticle_peta[igen2][ipdg], chain.genParticle_pphi[igen2][ipdg])
#####                            
#####                            if _dr < 0.1:
#####                                flag_pi_other = True
#####                                decayTable.append(addstr)
#####                                
#####
#####                        if not flag_pi_other:
#####                            decayTable.append('PU')
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

            out.gen_vz[0] = chain.JpsiK_genPV_vz

            out.npv_true[0] = chain.nPuVtxTrue[0] # current BX
            
            out.B_pt_gen[0] = chain.JpsiK_B_pt_gen


            if options.type == 'bg':
                out.genWeightBkgB[0] = chain.genWeightBkgB

#            if options.type in ['signal'] and chain.JpsiK_st_n_occurance == 1:

#                print evtid, chain.JpsiK_hammer_ebe[0]
                if options.type in ['signal']:
                    out.hammer_ebe[0] = chain.JpsiK_hammer_ebe[0]
                    out.hammer_wratio[0] = ROOT.Double(hist_hammer.GetBinContent(2))/ROOT.Double(hist_hammer.GetBinContent(1))
                    out.hammer_ebe_a0_up[0] = chain.JpsiK_hammer_ebe_a0_up[0]
                    out.hammer_ebe_a0_down[0] = chain.JpsiK_hammer_ebe_a0_down[0]
                    out.hammer_ebe_a1_up[0] = chain.JpsiK_hammer_ebe_a1_up[0]
                    out.hammer_ebe_a1_down[0] = chain.JpsiK_hammer_ebe_a1_down[0]
                    out.hammer_ebe_a2_up[0] = chain.JpsiK_hammer_ebe_a2_up[0]
                    out.hammer_ebe_a2_down[0] = chain.JpsiK_hammer_ebe_a2_down[0]

                    out.hammer_ebe_b0_up[0] = chain.JpsiK_hammer_ebe_b0_up[0]
                    out.hammer_ebe_b0_down[0] = chain.JpsiK_hammer_ebe_b0_down[0]
                    out.hammer_ebe_b1_up[0] = chain.JpsiK_hammer_ebe_b1_up[0]
                    out.hammer_ebe_b1_down[0] = chain.JpsiK_hammer_ebe_b1_down[0]
                    out.hammer_ebe_b2_up[0] = chain.JpsiK_hammer_ebe_b2_up[0]
                    out.hammer_ebe_b2_down[0] = chain.JpsiK_hammer_ebe_b2_down[0]

                    out.hammer_ebe_c1_up[0] = chain.JpsiK_hammer_ebe_c1_up[0]
                    out.hammer_ebe_c1_down[0] = chain.JpsiK_hammer_ebe_c1_down[0]
                    out.hammer_ebe_c2_up[0] = chain.JpsiK_hammer_ebe_c2_up[0]
                    out.hammer_ebe_c2_down[0] = chain.JpsiK_hammer_ebe_c2_down[0]

                    out.hammer_ebe_d0_up[0] = chain.JpsiK_hammer_ebe_d0_up[0]
                    out.hammer_ebe_d0_down[0] = chain.JpsiK_hammer_ebe_d0_down[0]
                    out.hammer_ebe_d1_up[0] = chain.JpsiK_hammer_ebe_d1_up[0]
                    out.hammer_ebe_d1_down[0] = chain.JpsiK_hammer_ebe_d1_down[0]
                    out.hammer_ebe_d2_up[0] = chain.JpsiK_hammer_ebe_d2_up[0]
                    out.hammer_ebe_d2_down[0] = chain.JpsiK_hammer_ebe_d2_down[0]




                #    for iham in range(len(chain.JpsiK_hammer_ebe_toy[0])):
            
                #        out.hammer_ebe_toy[iham] = ROOT.Double(chain.JpsiK_hammer_ebe_toy[0][iham])

                #else:
                #    out.hammer_ebe[0] = -1

            




            out.puweight[0] = ROOT.Double(putool.getWeight(chain.nPuVtxTrue[0]))


        evtid += 1
        out.tree.Fill()
    

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

    with open('json_jpsik.json', 'w') as outfile:
        json.dump(dict2save, outfile)

out.endJob()


