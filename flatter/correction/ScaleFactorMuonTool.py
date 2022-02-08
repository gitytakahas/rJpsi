# Author: Izaak Neutelings (November 2018)
import os, re, json, numpy 
datadir = os.environ['PWD'] + '/correction/'
from ROOT import TFile

class ScaleFactorMuonTool:
  
  def __init__(self, era, sigma='central', fileName='Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.json', keyName='NUM_TrackerMuons_DEN_genTracks'):
    """Load JSON with SFs"""
    with open(datadir+'/'+fileName) as f:    
      self.sf_file=json.load(f) 
    #print self.sf_file
    self.fileName = fileName
    self.keyName = keyName    
    self.edges = self.find_eta_pt_bins(self.sf_file[self.keyName]["abseta_pt"])


  def build_uncertainties(self, sf):
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
    

  def parse_str(self, key, prefix='abseta:'):
    if not key.startswith(prefix + "["):
      raise ValueError
    lo, hi = map(float, key[len(prefix + "["):-1].split(","))
    return lo, hi
          

  def find_eta_pt_bins(self, sf):
    #print sf
    bins = [self.parse_str(s, 'abseta:') for s in sf]
    #print "bins ", bins
    edges_eta = sorted(set(edge for bin in bins for edge in bin))
    #print " Eta edges ", edges_eta
    bins_pt = [self.parse_str(s, "pt:") for s in sf["abseta:[{:.2f},{:.2f}]".format(edges_eta[0], edges_eta[1])]]
    edges_pt = sorted(set(edge_pt for bin_pt in bins_pt for edge_pt in bin_pt))
    #print " Pt edges", edges_pt
    return ({'eta': edges_eta,'pt': edges_pt})
      
  def find_eta_boundaries(self, eta):
    eta_string=""
    for i in range(len(self.edges['eta'])):
    #print "edges_eta [i] ",i," ",edges['eta'] [i]
      if self.edges['eta'][i] > eta:
        eta_string="abseta:[{:.2f},{:.2f}]".format(self.edges['eta'][i-1],self.edges['eta'][i])
        break;
    if eta_string=="": eta_string="abseta:[{:.2f},{:.2f}]".format(self.edges['eta'][-2],self.edges['eta'][-1])
    return eta_string
    
  def find_pt_boundaries(self, pt):
    string=""
    for i in range(len(self.edges['pt'])):
        #print "edges_pt [i] ",i," ",edges['pt'] [i]
        if self.edges['pt'][0]> pt:
            ValueError("pt value too small,no sf available")
        if self.edges['pt'][i] > pt:
            string="pt:[{:.2f},{:.2f}]".format(self.edges['pt'][i-1],self.edges['pt'][i])
            break;
    if string=="": string="pt:[{:.2f},{:.2f}]".format(self.edges['pt'][-2],self.edges['pt'][-1])
    return string
  
  
  def getSF(self,eta,pt):
    """Get SF for a given eta and pt."""
  
    sf = self.sf_file[self.keyName]["abseta_pt"][self.find_eta_boundaries(eta)][self.find_pt_boundaries(pt)]
    return sf
  
  
