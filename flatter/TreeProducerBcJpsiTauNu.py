import ROOT
import math 
import numpy as num 
from TreeProducerCommon import *

class TreeProducerBcJpsiTauNu(TreeProducerCommon):
    """Class to create a custom output file & tree; as well as create and contain branches."""

    def __init__(self, name, dataType, **kwargs):
        print('TreeProducerBsTauTau is called for', name)
        super(TreeProducerBcJpsiTauNu, self).__init__(name,dataType,**kwargs)

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

  
        self.addBranch('tau_pt',                  'f')
        self.addBranch('tau_eta',                  'f')
        self.addBranch('tau_phi',                  'f')
        self.addBranch('tau_mass',                  'f')
        self.addBranch('tau_q',                  'i')
        self.addBranch('tau_dr_jpsi',                  'f')
        self.addBranch('tau_rhomass1',                  'f')
        self.addBranch('tau_rhomass2',                  'f')
        self.addBranch('tau_rhomass_ss',                  'f')
        self.addBranch('tau_rhomass',                  'f')
        self.addBranch('tau_psimass1',                  'f')
        self.addBranch('tau_psimass2',                  'f')

        self.addBranch('tau_psimass1_pk',                  'f')
        self.addBranch('tau_psimass2_pk',                  'f')
        self.addBranch('tau_psimass1_kp',                  'f')
        self.addBranch('tau_psimass2_kp',                  'f')
        self.addBranch('tau_psimass1_kk',                  'f')
        self.addBranch('tau_psimass2_kk',                  'f')

        self.addBranch('tau_jpsi_pi1',                  'f')
        self.addBranch('tau_jpsi_k1',                  'f')

        self.addBranch('tau_jpsi_pi2',                  'f')
        self.addBranch('tau_jpsi_k2',                  'f')

        self.addBranch('tau_jpsi_pi3',                  'f')
        self.addBranch('tau_jpsi_k3',                  'f')

        self.addBranch('tau_m12s',                  'f')
        self.addBranch('tau_m23s',                  'f')
        self.addBranch('tau_m13s',                  'f')
        self.addBranch('tau_vprob',                  'f')
        self.addBranch('tau_fls3d',                  'f')
        self.addBranch('tau_fl3d',                  'f')
        self.addBranch('tau_sumofdnn',                  'f')
        self.addBranch('tau_sumofdnn_1prong',                  'f')
        self.addBranch('tau_sumofdnn_otherB',                  'f')
        self.addBranch('tau_sumofdnn_pu',                  'f')


        self.addBranch('tau_max_dr',                  'f')
        self.addBranch('tau_lips',                  'f')
        self.addBranch('tau_pvips',                  'f')
        self.addBranch('tau_index',                  'i')
        self.addBranch('tau_eid',                  'i')

        self.addBranch('ptbal',                  'f')
        self.addBranch('jpsi_tau_alpha',                  'f')

        self.addBranch('pi1_pt',                  'f')
        self.addBranch('pi1_eta',                  'f')
        self.addBranch('pi1_phi',                  'f')
        self.addBranch('pi1_dnn',                  'f')
        self.addBranch('pi1_dnn_1prong',                  'f')
        self.addBranch('pi1_dnn_otherB',                  'f')
        self.addBranch('pi1_dnn_pu',                  'f')
        self.addBranch('pi1_trigMatch',                  '?')
        self.addBranch('pi1_trigMatch_dr',                  'f')

        self.addBranch('pi2_pt',                  'f')
        self.addBranch('pi2_eta',                  'f')
        self.addBranch('pi2_phi',                  'f')
        self.addBranch('pi2_dnn',                  'f')
        self.addBranch('pi2_dnn_1prong',                  'f')
        self.addBranch('pi2_dnn_otherB',                  'f')
        self.addBranch('pi2_dnn_pu',                  'f')
        self.addBranch('pi2_trigMatch',                  '?')
        self.addBranch('pi2_trigMatch_dr',                  'f')

        self.addBranch('pi3_pt',                  'f')
        self.addBranch('pi3_eta',                  'f')
        self.addBranch('pi3_phi',                  'f')
        self.addBranch('pi3_dnn',                  'f')
        self.addBranch('pi3_dnn_1prong',                  'f')
        self.addBranch('pi3_dnn_otherB',                  'f')
        self.addBranch('pi3_dnn_pu',                  'f')
        self.addBranch('pi3_trigMatch',                  '?')
        self.addBranch('pi3_trigMatch_dr',                  'f')



#        self.addBranch('tau_ppp',                  'f')
        self.addBranch('tau_kpp',                  'f')
        self.addBranch('tau_pkp',                  'f')
        self.addBranch('tau_ppk',                  'f')
        self.addBranch('tau_kkp',                  'f')
        self.addBranch('tau_kpk',                  'f')
        self.addBranch('tau_pkk',                  'f')
        self.addBranch('tau_kkk',                  'f')

        self.addBranch('tau_rhomass1_kp',                  'f')
        self.addBranch('tau_rhomass2_kp',                  'f')

        self.addBranch('tau_rhomass1_pk',                  'f')
        self.addBranch('tau_rhomass2_pk',                  'f')

        self.addBranch('tau_rhomass1_kk',                  'f')
        self.addBranch('tau_rhomass2_kk',                  'f')

        self.addBranch('tau_isRight_3prong',                  '?')
        self.addBranch('tau_isRight_3prong_pi0',                  '?')

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

            if dataType in ['signal']:
                self.addBranch('hammer_ebe',                  'f')
                self.addBranch('hammer_wratio',                  'f')
                self.addBranch('hammer_ebe_lattice',                  'f')
                self.addBranch('hammer_wratio_lattice',                  'f')
                self.addBranch('isJpsiMu',                  '?')
                self.addBranch('isJpsiTau2Mu',                  '?')


                for ii in range(0, 15):
                    for ud in ['up', 'down']:
                        self.addBranch('hammer_ebe_e' + str(ii) + '_' + ud,                  'f')
                        self.addBranch('hammer_ebe_e' + str(ii) + '_' + ud + '_lattice',                  'f')

#                self.addBranch('hammer_ebe_a0_up',                  'f')
#                self.addBranch('hammer_ebe_a0_down',                  'f')
#                self.addBranch('hammer_ebe_a1_up',                  'f')
#                self.addBranch('hammer_ebe_a1_down',                  'f')
#                self.addBranch('hammer_ebe_a2_up',                  'f')
#                self.addBranch('hammer_ebe_a2_down',                  'f')
#                self.addBranch('hammer_ebe_b0_up',                  'f')
#                self.addBranch('hammer_ebe_b0_down',                  'f')
#                self.addBranch('hammer_ebe_b1_up',                  'f')
#                self.addBranch('hammer_ebe_b1_down',                  'f')
#                self.addBranch('hammer_ebe_b2_up',                  'f')
#                self.addBranch('hammer_ebe_b2_down',                  'f')
#                self.addBranch('hammer_ebe_c1_up',                  'f')
#                self.addBranch('hammer_ebe_c1_down',                  'f')
#                self.addBranch('hammer_ebe_c2_up',                  'f')
#                self.addBranch('hammer_ebe_c2_down',                  'f')
#                self.addBranch('hammer_ebe_d0_up',                  'f')
#                self.addBranch('hammer_ebe_d0_down',                  'f')
#                self.addBranch('hammer_ebe_d1_up',                  'f')
#                self.addBranch('hammer_ebe_d1_down',                  'f')
#                self.addBranch('hammer_ebe_d2_up',                  'f')
#                self.addBranch('hammer_ebe_d2_down',                  'f')
#
#                self.addBranch('hammer_ebe_a0_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_a0_down_lattice',                  'f')
#                self.addBranch('hammer_ebe_a1_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_a1_down_lattice',                  'f')
#                self.addBranch('hammer_ebe_a2_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_a2_down_lattice',                  'f')
#                self.addBranch('hammer_ebe_b0_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_b0_down_lattice',                  'f')
#                self.addBranch('hammer_ebe_b1_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_b1_down_lattice',                  'f')
#                self.addBranch('hammer_ebe_b2_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_b2_down_lattice',                  'f')
#                self.addBranch('hammer_ebe_c1_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_c1_down_lattice',                  'f')
#                self.addBranch('hammer_ebe_c2_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_c2_down_lattice',                  'f')
#                self.addBranch('hammer_ebe_d0_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_d0_down_lattice',                  'f')
#                self.addBranch('hammer_ebe_d1_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_d1_down_lattice',                  'f')
#                self.addBranch('hammer_ebe_d2_up_lattice',                  'f')
#                self.addBranch('hammer_ebe_d2_down_lattice',                  'f')



#                self.addBranch('hammer_ebe_toy',                  'v')
                self.addBranch('weight_ctau',                  'f')
                self.addBranch('weight_ctau_up',                  'f')
                self.addBranch('weight_ctau_down',                  'f')


#            if dataType=='bg':
#                self.addBranch('isveto',                  '?')

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
