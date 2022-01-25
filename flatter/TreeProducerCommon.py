import numpy as num
import os, math, sys
#from ROOT import TFile, TTree, TH1F, gROOT, TTree, Double, TChain, TLorentzVector, TVector3
#import numpy as num
from array import array
import ROOT

root_dtype = {
  float: 'D',  int: 'I',  bool: 'O',
  'f':   'D',  'i': 'I',  '?':  'O',  'b': 'b', 'v':'D', 's':'s'
}

num_dtype = {
  'D':   'f',  'I': 'i',  'O':  '?',  'b': 'b'
}


class TreeProducerCommon(object):
    """Class to create a custom output file & tree; as well as create and contain branches."""
    
    def __init__(self, name, dataType, **kwargs):
        print('TreeProducerCommon is called for', name)
        
        self.name       = name
        self._isData    = dataType=='data'
        
        # TREE
        self.outputfile = ROOT.TFile(name, 'RECREATE')
        self.tree       = ROOT.TTree('tree','tree')
        self.v = ROOT.std.vector( ROOT.std.string )()
        self.tree._v = self.v


        ###################################
        self.hist = ROOT.TH1F('cutflow', 'cutflow', 15,0,15)
#        self.hist = None
        self.multi = ROOT.TH1F('multi', 'multi', 100,0,100)

        self.filt = ROOT.TH1F('filter', 'filter', 10,0,10)

#        if dataType=='signal':
  #           print('This is for signal !!! ')
        self.hist_hammer = ROOT.TH1F('hammer', 'hammer', 32,0,32)
        self.hist_hammer_lattice = ROOT.TH1F('hammer_lattice', 'hammer_lattice', 32,0,32)

        
        self.addBranch('evt',                  'i')
        self.addBranch('run',                  'i')
        self.addBranch('lumi',                  'i')

        self.addBranch('evtid',                  'i')

        self.addBranch('b_pt',                  'f')
        self.addBranch('b_eta',                  'f')
        self.addBranch('b_phi',                  'f')
        self.addBranch('b_mass',                  'f')
        self.addBranch('b_mcorr',                  'f')
        self.addBranch('b_alpha',                  'f')
        self.addBranch('b_fls3d',                  'f')
        self.addBranch('b_fl3d',                  'f')
        self.addBranch('b_vprob',                  'f')

        self.addBranch('b_pt_simple',                  'f')
        self.addBranch('b_eta_simple',                  'f')
        self.addBranch('b_phi_simple',                  'f')
        self.addBranch('b_mass_simple',                  'f')


        for iso in range(5):
            iso_ = str(0.5 + iso*0.2).replace('.','p')
            self.addBranch('tau_iso_' + str(iso_),                  'f')
            self.addBranch('tau_iso_ntracks_' + str(iso_),                  'i')
#            self.addBranch('b_iso_mindoca_' + str(iso_),                  'f')

#        self.addBranch('b_iso_nocut',                  'f')
#        self.addBranch('b_iso_ntracks_nocut',                  'i')
#        self.addBranch('b_iso_mindoca_nocut',                  'f')

#        self.addBranch('pv_vx',                  'f')
#        self.addBranch('pv_vy',                  'f')
        self.addBranch('pv_vz',                  'f')

#        self.addBranch('bbpv_vx',                  'f')
#        self.addBranch('bbpv_vy',                  'f')
        self.addBranch('bbpv_vz',                  'f')

        self.addBranch('b_pvip',                  'f')
        self.addBranch('b_pvips',                  'f')
        self.addBranch('b_lips',                  'f')
        self.addBranch('b_mindoca',                  'f')
#        self.addBranch('b_vx',                  'f')
#        self.addBranch('b_vy',                  'f')
        self.addBranch('b_vz',                  'f')


        if dataType!='data':
#            self.addBranch('gen_vx',                  'f')
#            self.addBranch('gen_vy',                  'f')
            self.addBranch('gen_vz',                  'f')
            self.addBranch('puweight',                  'f')
            self.addBranch('puweight_up',                  'f')
            self.addBranch('puweight_down',                  'f')
            self.addBranch('npv_true',                  'i')


#        self.addBranch('b_mcol',                  'f')

#        if not self._isData:
#          self.addBranch('pt_genboson',           'f', -1)
    
        #################################

    def addBranch(self, name, dtype='f', default=None):
        """Add branch with a given name, and create an array of the same name as address."""
        if hasattr(self,name):
          print("ERROR! TreeProducerCommon.addBranch: Branch of name '%s' already exists!"%(name))
          exit(1)
        if isinstance(dtype,str):
          if dtype.lower()=='f': # 'f' is only a 'float32', and 'F' is a 'complex64', which do not work for filling float branches
            dtype = float        # float is a 'float64' ('f8')
          elif dtype.lower()=='i': # 'i' is only a 'int32'
            dtype = int            # int is a 'int64' ('i8')


        if dtype=='v':
            setattr(self,name,array('d', 1000*[0]))
#            setattr(self,name,array('d', 1000,dtype='f'))
            self.tree.Branch(name, getattr(self,name), '%s[1000]/%s'%(name,root_dtype[dtype]))
        elif dtype=='s':

            self.tree.Branch( name, self.v )

#            setattr(self,name,None)

#            self.tree.Branch(name, getattr(self,name), '%s'%(name))
        else:
            setattr(self,name,num.zeros(1,dtype=dtype))
            self.tree.Branch(name, getattr(self,name), '%s/%s'%(name,root_dtype[dtype]))

        if default!=None:
          getattr(self,name)[0] = default
        
    
    def endJob(self):
        """Write and close files after the job ends."""
#        self.hist.Write()
#        ROOT.gDirectory.cd()
#        print(self.hist_hammer, self.hist_hammer_lattice)
#        self.hist_hammer.Write()
#        self.hist_hammer_lattice.Write()

        self.outputfile.Write()
        self.outputfile.Close()
        

