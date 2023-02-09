import ROOT
import math 
import numpy as num 
from TreeProducerCommon import *

class TreeProducerJpsiPi(TreeProducerCommon):
    """Class to create a custom output file & tree; as well as create and contain branches."""

    def __init__(self, name, dataType, **kwargs):
        super(TreeProducerJpsiPi, self).__init__(name,dataType,**kwargs)

        self.addBranch('mu1_pt',                  'f')
        self.addBranch('mu1_eta',                  'f')
        self.addBranch('mu1_phi',                  'f')
        self.addBranch('mu1_mass',                  'f')
        self.addBranch('mu1_q', 'i')
        self.addBranch('mu1_isLoose', 'i')
        self.addBranch('mu1_isTight', 'i')
        self.addBranch('mu1_isPF', 'i')
        self.addBranch('mu1_isGlobal', 'i')
        self.addBranch('mu1_isTracker', 'i')
        self.addBranch('mu1_isSoft', 'i')
        self.addBranch('mu1_dbiso', 'f')
        self.addBranch('mu1_SFReco', 'f')
        self.addBranch('mu1_SFID', 'f')
  
        self.addBranch('mu2_pt',                  'f')
        self.addBranch('mu2_eta',                  'f')
        self.addBranch('mu2_phi',                  'f')
        self.addBranch('mu2_mass',                  'f')
        self.addBranch('mu2_q', 'i')
        self.addBranch('mu2_isLoose', 'i')
        self.addBranch('mu2_isTight', 'i')
        self.addBranch('mu2_isPF', 'i')
        self.addBranch('mu2_isGlobal', 'i')
        self.addBranch('mu2_isTracker', 'i')
        self.addBranch('mu2_isSoft', 'i')
        self.addBranch('mu2_dbiso', 'f')
        self.addBranch('mu2_SFReco', 'f')
        self.addBranch('mu2_SFID', 'f')

        self.addBranch('mu1_mu2_dr', 'f')
  
        self.addBranch('pi_pt',                  'f')
        self.addBranch('pi_eta',                  'f')
        self.addBranch('pi_phi',                  'f')
        self.addBranch('pi_mass',                  'f')
        self.addBranch('pi_q',                  'i')
        self.addBranch('pi_dr_jpsi',                  'f')
        
        
        self.addBranch('ptbal',                  'f')
        self.addBranch('jpsi_tau_alpha',                  'f')
        
        self.addBranch('k_pt',                  'f') 
        self.addBranch('B_k_mass',                'f')      
        self.addBranch('B_k_pt',                  'f')      

#        self.addBranch('tau_ppp',                  'f')
        
        self.addBranch('npv',                  'i')

        self.addBranch('ncand',                  'i')

        self.addBranch('dx_jpsi_tau',                  'f')
        self.addBranch('dy_jpsi_tau',                  'f')
        self.addBranch('dz_jpsi_tau',                  'f')
        self.addBranch('dr_jpsi_tau',                  'f')
        
        self.addBranch('dx_b_pv',                  'f')
        self.addBranch('dy_b_pv',                  'f')
        self.addBranch('dz_b_pv',                  'f')
        self.addBranch('dr_b_pv',                  'f')


        self.addBranch('jpsi_pt',                  'f')
        self.addBranch('jpsi_eta',                  'f')
        self.addBranch('jpsi_phi',                  'f')
        self.addBranch('jpsi_mass',                  'f')
        self.addBranch('jpsi_vprob',                  'f')
        self.addBranch('jpsi_fls3d',                  'f')
        self.addBranch('jpsi_fl3d',                  'f')

        self.addBranch('estar',                  'f')
        self.addBranch('q2',                  'f')
        self.addBranch('ptmiss',                  'f')
        self.addBranch('mm2',                  'f')
        self.addBranch('B_ptback',                  'f')

        self.addBranch('estar_simple',                  'f')
        self.addBranch('q2_simple',                  'f')
        self.addBranch('ptmiss_simple',                  'f')
        self.addBranch('mm2_simple',                  'f')
        self.addBranch('B_ptback_simple',                  'f')

        self.addBranch('delta_chi2',                  'f')
        self.addBranch('delta_n_ch',                  'f')
        self.addBranch('delta_n_mu',                  'f')
        self.addBranch('vweight',                  'f')
        self.addBranch('tau_fl3d_wjpsi',                  'f')
        self.addBranch('tau_fls3d_wjpsi',                  'f')

        self.addBranch('tau_refit_chi2',                  'f')
        self.addBranch('tau_refit_ndof',                  'f')
        self.addBranch('tau_refit_rho',                  'f')

#        self.addBranch('estar_vis',                  'f')
#        self.addBranch('q2_vis',                  'f')
#        self.addBranch('ptmiss_vis',                  'f')
#        self.addBranch('mm2_vis',                  'f')

        self.addBranch('perEVT_mc',                  'f')
        self.addBranch('perEVT_data',                  'f')

#        self.addBranch('perEVT_old',                  'f')
#        self.addBranch('perEVT_otherB',                  'f')
#        self.addBranch('perEVT_sig',                  'f')
#        self.addBranch('perEVT_leptonic',                  'f')
#        self.addBranch('perEVT_1prong',                  'f')

#        if dataType=='bg':
        self.addBranch('nkaon',                  'i')
        self.addBranch('npion',                  'i')
        self.addBranch('npu',                  'i')
        self.addBranch('pid',                  'i')

        if dataType!='data':
 #           self.addBranch('tau_isgen3matched',                  '?')
 #           self.addBranch('tau_isgen3',                  '?')
 #           self.addBranch('tau_gendm',                  'i')
 #           self.addBranch('tau_genpt',                  'f')
 #           self.addBranch('tau_geneta',                  'f')
 #           self.addBranch('tau_genphi',                  'f')
#            self.addBranch('tau_genmass',                  'f')
 #           self.addBranch('ngentau',                  'i')

            self.addBranch('genWeightBkgB',                  'f')
            self.addBranch('q2_gen',                  'f')
            self.addBranch('B_pt_gen',                  'f')
            self.addBranch('isjpsimatched',                  '?')
            self.addBranch('tau_genpt',                  'f')
            self.addBranch('tau_geneta',                  'f')
            self.addBranch('tau_genphi',                  'f')
            self.addBranch('n_occurance',                  'i')
            self.addBranch('decayid',                  'i')
            self.addBranch('procid',                  'i')

            if dataType in ['signal']:
                self.addBranch('hammer_ebe',                  'f')
                self.addBranch('hammer_wratio',                  'f')
                self.addBranch('hammer_ebe_toy',                  'v')
                self.addBranch('weight_ctau',                  'f')
                self.addBranch('weight_ctau_up',                  'f')
                self.addBranch('weight_ctau_down',                  'f')

                self.addBranch('hammer_ebe_a0_up',                  'f')
                self.addBranch('hammer_ebe_a0_down',                  'f')
                self.addBranch('hammer_ebe_a1_up',                  'f')
                self.addBranch('hammer_ebe_a1_down',                  'f')
                self.addBranch('hammer_ebe_a2_up',                  'f')
                self.addBranch('hammer_ebe_a2_down',                  'f')
                self.addBranch('hammer_ebe_b0_up',                  'f')
                self.addBranch('hammer_ebe_b0_down',                  'f')
                self.addBranch('hammer_ebe_b1_up',                  'f')
                self.addBranch('hammer_ebe_b1_down',                  'f')
                self.addBranch('hammer_ebe_b2_up',                  'f')
                self.addBranch('hammer_ebe_b2_down',                  'f')
                self.addBranch('hammer_ebe_c1_up',                  'f')
                self.addBranch('hammer_ebe_c1_down',                  'f')
                self.addBranch('hammer_ebe_c2_up',                  'f')
                self.addBranch('hammer_ebe_c2_down',                  'f')
                self.addBranch('hammer_ebe_d0_up',                  'f')
                self.addBranch('hammer_ebe_d0_down',                  'f')
                self.addBranch('hammer_ebe_d1_up',                  'f')
                self.addBranch('hammer_ebe_d1_down',                  'f')
                self.addBranch('hammer_ebe_d2_up',                  'f')
                self.addBranch('hammer_ebe_d2_down',                  'f')
                self.addBranch('gen_sig_decay',                  'i')  

            if dataType=='bg':
                self.addBranch('isveto',                  '?')

#                self.addBranch('sameB',                  '?')
#                self.addBranch('diffB',                  '?')
#                self.addBranch('pu',                  '?')


        self.addBranch('nch',                  'i')
#        self.addBranch('nch_qr',                  'i')
        self.addBranch('nch_before_dnn',                  'i')
 #       self.addBranch('nch_before_dnn',                  'i')


#        self.addBranch('dchain',                  's')



#        self.addBranch('phi_mass',     'f')
#        self.addBranch('phi_pt',     'f')
        
#        if not self._isData:
#          self.addBranch('genPartFlav_1',            'i', -1)

##        self.addBranch('JpsiK_st_size','i') 
##        self.addBranch('JpsiK_st_doca3d','vec12')
###        self.JpsiK_st_doca3d=ROOT.vector('float')()
###        self.tree.Branch('JpsiK_st_doca3d', ROOT.AddressOf(self.JpsiK_st_doca3d),'')
##        self.addBranch('JpsiK_st_doca2d','vec12')
##        self.addBranch('JpsiK_st_doca3de','vec12')
##        self.addBranch('JpsiK_st_doca2de','vec12')
##        self.addBranch('JpsiK_st_doca3ds','vec12')
##        self.addBranch('JpsiK_st_doca2ds','vec12')
##        self.addBranch('JpsiK_st_dz','vec12')
##        
##        self.addBranch('JpsiK_st_isAssociate','vec12')
##        self.addBranch('JpsiK_st_pvAssociationQuality','vec12')
##        self.addBranch('JpsiK_st_pt','vec12')
##        self.addBranch('JpsiK_st_eta','vec12')
##        self.addBranch('JpsiK_st_phi','vec12')
##        self.addBranch('JpsiK_st_mass','vec12')
##        self.addBranch('JpsiK_st_charge','vec12')
##
##        self.addBranch('JpsiK_st_isBdecay','vec12')
##        self.addBranch('JpsiK_st_isBdecaypdg','vec12')
##        self.addBranch('JpsiK_st_isBdecayppdg','vec12')
##        self.addBranch('JpsiK_st_isSignal','vec12')
##        self.addBranch('JpsiK_st_nprong','vec12')
##        self.addBranch('JpsiK_st_nprong_pi0','vec12')
##        self.addBranch('JpsiK_st_near_dz','vec12')
##        self.addBranch('JpsiK_st_dr_jpsi','vec12')
##        self.addBranch('JpsiK_st_trigMatch_dr','vec12')
##        self.addBranch('JpsiK_st_trigMatch','vec12')
        
