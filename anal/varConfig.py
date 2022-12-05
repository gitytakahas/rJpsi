import math

vardir = {}

vardir["tau_pt"] = {'tree':'tree',  'nbin':30, 'xmin':0, 'xmax':30, 'xtitle':'tau pT (GeV)'}
vardir["tau_eta"] = {'tree':'tree',  'nbin':30, 'xmin':-2.5, 'xmax':2.5, 'xtitle':'tau eta'}
vardir["tau_phi"] = {'tree':'tree',  'nbin':30, 'xmin':-math.pi, 'xmax':math.pi, 'xtitle':'tau phi'}
vardir["tau_q"] = {'tree':'tree',  'nbin':7, 'xmin':-3, 'xmax':4, 'xtitle':'tau q'}

vardir["tau_vprob"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1, 'xtitle':'Tau vertex prob.'}
vardir["tau_fls3d"] = {'tree':'tree',  'nbin':30, 'xmin':2., 'xmax':30., 'xtitle':'Tau flight length sig.'}
vardir["tau_fl3d"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'Tau flight length'}
vardir["tau_sumofdnn"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Sum of DNN'}
#vardir["tau_sumofdnn_old"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Sum of DNN old'}
vardir["tau_sumofdnn_1prong"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Sum of DNN 1prong'}
vardir["tau_sumofdnn_otherB"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Sum of DNN otherB'}
vardir["tau_sumofdnn_pu"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Sum of DNN pu'}

vardir["tau_sumofdnn_log"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Sum of DNN', 'var':'tau_sumofdnn', 'isLog':True}
#vardir["tau_sumofdnn_old_log"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Sum of DNN old', 'var':'tau_sumofdnn_old', 'isLog':True}
vardir["tau_sumofdnn_1prong_log"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Sum of DNN 1prong', 'var':'tau_sumofdnn_1prong', 'isLog':True}
vardir["tau_sumofdnn_otherB_log"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Sum of DNN otherB', 'var':'tau_sumofdnn_otherB', 'isLog':True}
vardir["tau_sumofdnn_pu_log"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Sum of DNN pu', 'var':'tau_sumofdnn_pu', 'isLog':True}


vardir["tau_index"] = {'tree':'tree',  'nbin':10, 'xmin':0., 'xmax':10., 'xtitle':'tau index'}
#vardir["tau_sumofdnn_zoom"] = {'tree':'tree',  'nbin':60, 'xmin':0, 'xmax':3., 'xtitle':'Sum of DNN zoom', 'var':'tau_sumofdnn'}
#vardir["tau_sumofdnn_ratio"] = {'tree':'tree',  'nbin':20, 'xmin':0., 'xmax':1., 'xtitle':'Sum of DNN ratio'}

#vardir["tau_lips"] = {'tree':'tree',  'nbin':20, 'xmin':-20., 'xmax':20., 'xtitle':'tau lips'}
#vardir["tau_pvips"] = {'tree':'tree',  'nbin':20, 'xmin':0., 'xmax':30., 'xtitle':'tau pvips'}

#vardir["reg_taupt"] = {'tree':'tree',  'nbin':20, 'xmin':0, 'xmax':20, 'xtitle':'reg. tau pT (GeV)'}
#vardir["reg_taueta"] = {'tree':'tree',  'nbin':20, 'xmin':-2.5, 'xmax':2.5, 'xtitle':'reg. tau eta'}
vardir["tau_mass"] = {'tree':'tree',  'nbin':30, 'xmin':0.4, 'xmax':1.8, 'xtitle':'Tau mass (GeV)'}
#vardir["tau_max_dr"] = {'tree':'tree',  'nbin':20, 'xmin':0., 'xmax':2*math.pi, 'xtitle':'Tau max dr'}
vardir["tau_mass_zoom"] = {'tree':'tree',  'nbin':40, 'xmin':0.7, 'xmax':1.7, 'xtitle':'Tau mass wide (GeV)', 'var':'tau_mass'}
vardir["tau_rhomass1"] = {'tree':'tree',  'nbin':25, 'xmin':0.2, 'xmax':1.4, 'xtitle':'Tau rhomass1 (GeV)'}
vardir["tau_rhomass2"] = {'tree':'tree',  'nbin':25, 'xmin':0.2, 'xmax':1.4, 'xtitle':'Tau rhomass2 (GeV)'}
#vardir["tau_rhomass_unrolled"] = {'tree':'tree',  'nbin':625, 'xmin':0, 'xmax':625, 'xtitle':'Tau rhomasses Unrolled bin ID', 'var':'(min(tau_rhomass2, 1.4) - 0.2)/0.048 + 25*(min(tau_rhomass1, 1.4) - 0.2)/0.048'}
vardir["tau_rhomass_unrolled"] = {'tree':'tree',  'nbin':121, 'xmin':0, 'xmax':121, 'xtitle':'Tau rhomasses Unrolled bin ID', 'var':'int((min(tau_rhomass2, 1.3) - 0.2)/0.11) + 11*int((min(tau_rhomass1, 1.3) - 0.2)/0.11)'}

vardir["tau_rhomass_unrolled_coarse"] = {'tree':'tree',  'nbin':36, 'xmin':0, 'xmax':36, 'xtitle':'Tau rhomasses Unrolled bin ID', 'var':'int((min(tau_rhomass2, 1.3) - 0.2)/0.22) + 6*int((min(tau_rhomass1, 1.3) - 0.2)/0.22)'}
#vardir["tau_rhomass_unrolled_coarse"] = {'tree':'tree',  'nbin':25, 'xmin':0, 'xmax':25, 'xtitle':'Tau rhomasses Unrolled bin ID', 'var':'int((min(tau_rhomass2, 1.1) - 0.2)/0.225) + 5*int((min(tau_rhomass1, 1.1) - 0.2)/0.225)'}

#vardir["tau_rhomass_unrolled_var"] = {'tree':'tree',  'nbin':25, 'xmin':0, 'xmax':25, 'xtitle':'Tau rhomasses1', 'ytitle':'Tau rhomasses2', 'var':'int((min(tau_rhomass2, 1.1) - 0.2)/0.225) + 5*int((min(tau_rhomass1, 1.1) - 0.2)/0.225)','isVar':True}


vardir["tau_rhomass_ss"] = {'tree':'tree',  'nbin':30, 'xmin':0.2, 'xmax':1.5, 'xtitle':'Tau rhomass ss (GeV)'}
vardir["tau_rhomass_random"] = {'tree':'tree',  'nbin':30, 'xmin':0.2, 'xmax':1.5, 'xtitle':'Tau rhomass random (GeV)', 'var':'tau_rhomass'}
vardir["tau_rhomass_ave"] = {'tree':'tree',  'nbin':30, 'xmin':0.2, 'xmax':1.5, 'xtitle':'Tau rhomass average (GeV)', 'var':'(tau_rhomass1 + tau_rhomass2)/2.'}

vardir["tau_rhomass1_kp"] = {'tree':'tree',  'nbin':40, 'xmin':0.4, 'xmax':2.4, 'xtitle':'Tau rhomass1 kp (GeV)'}
vardir["tau_rhomass2_kp"] = {'tree':'tree',  'nbin':40, 'xmin':0.4, 'xmax':2.4, 'xtitle':'Tau rhomass2 kp (GeV)'}

vardir["tau_rhomass1_pk"] = {'tree':'tree',  'nbin':40, 'xmin':0.4, 'xmax':2.4, 'xtitle':'Tau rhomass1 pk (GeV)'}
vardir["tau_rhomass2_pk"] = {'tree':'tree',  'nbin':40, 'xmin':0.4, 'xmax':2.4, 'xtitle':'Tau rhomass2 pk (GeV)'}

vardir["tau_rhomass1_kk"] = {'tree':'tree',  'nbin':40, 'xmin':0.7, 'xmax':3., 'xtitle':'Tau rhomass1 kk (GeV)'}
vardir["tau_rhomass2_kk"] = {'tree':'tree',  'nbin':40, 'xmin':0.7, 'xmax':3., 'xtitle':'Tau rhomass2 kk (GeV)'}

vardir["tau_m12s"] = {'tree':'tree',  'nbin':30, 'xmin':0.2, 'xmax':5., 'xtitle':'Tau mass 12^2 (GeV)'}
vardir["tau_m23s"] = {'tree':'tree',  'nbin':30, 'xmin':0.2, 'xmax':5., 'xtitle':'Tau mass 23^2 (GeV)'}
vardir["tau_m13s"] = {'tree':'tree',  'nbin':30, 'xmin':0.2, 'xmax':5., 'xtitle':'Tau mass 13^2 (GeV)'}
#vardir["tau_ppp"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':3., 'xtitle':'Tau ppp (GeV)'}
vardir["tau_ppk"] = {'tree':'tree',  'nbin':40, 'xmin':0., 'xmax':4., 'xtitle':'Tau ppk (GeV)'}
vardir["tau_pkp"] = {'tree':'tree',  'nbin':40, 'xmin':0., 'xmax':4., 'xtitle':'Tau pkp (GeV)'}
vardir["tau_kpp"] = {'tree':'tree',  'nbin':40, 'xmin':0., 'xmax':4., 'xtitle':'Tau kpp (GeV)'}
vardir["tau_pkk"] = {'tree':'tree',  'nbin':40, 'xmin':1., 'xmax':5., 'xtitle':'Tau pkk (GeV)'}
vardir["tau_kpk"] = {'tree':'tree',  'nbin':40, 'xmin':1., 'xmax':5., 'xtitle':'Tau kpk (GeV)'}
vardir["tau_kkp"] = {'tree':'tree',  'nbin':40, 'xmin':1., 'xmax':5., 'xtitle':'Tau kkp (GeV)'}
vardir["tau_kkk"] = {'tree':'tree',  'nbin':40, 'xmin':1., 'xmax':5., 'xtitle':'Tau kkk (GeV)'}
vardir["b_mindoca"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':0.01, 'xtitle':'B mindoca', 'isLog':True}
#vardir["b_iso_mindoca"] = {'tree':'tree',  'nbin':20, 'xmin':0., 'xmax':0.1, 'xtitle':'B iso mindoca'}
    
vardir["b_pt"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':50, 'xtitle':'B pT (GeV)'}
vardir["b_eta"] = {'tree':'tree',  'nbin':30, 'xmin':-3., 'xmax':3., 'xtitle':'B eta'}
vardir["b_phi"] = {'tree':'tree',  'nbin':30, 'xmin':-math.pi, 'xmax':math.pi, 'xtitle':'B phi'}
vardir["b_pt_simple"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':50, 'xtitle':'B pT simple (GeV)'}
vardir["b_eta_simple"] = {'tree':'tree',  'nbin':30, 'xmin':-3., 'xmax':3., 'xtitle':'B eta simple'}
vardir["b_phi_simple"] = {'tree':'tree',  'nbin':30, 'xmin':-math.pi, 'xmax':math.pi, 'xtitle':'B phi simple'}
vardir["tau_iso_0p7"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':200, 'xtitle':'B iso (pT > 0.7 GeV)', 'isLog':True}
#vardir["nch"] = {'tree':'tree',  'nbin':50, 'xmin':0, 'xmax':50, 'xtitle':'Num of charged hadrons (50)'}
#vardir["nch_qr"] = {'tree':'tree',  'nbin':120, 'xmin':0, 'xmax':120, 'xtitle':'Nch qr'}
#vardir["nch_after_dnn"] = {'tree':'tree',  'nbin':25, 'xmin':0, 'xmax':25, 'xtitle':'Nch after DNN'}
#vardir["nch_before_dnn"] = {'tree':'tree',  'nbin':120, 'xmin':0, 'xmax':120, 'xtitle':'Nch before DNN'}
vardir["npv"] = {'tree':'tree',  'nbin':100, 'xmin':0, 'xmax':100, 'xtitle':'Num of PV'}
vardir["tau_iso_ntracks_0p7"] = {'tree':'tree',  'nbin':8, 'xmin':0, 'xmax':8, 'xtitle':'B iso ntracks (pT > 0.7 GeV)'}
vardir["b_vprob"] = {'tree':'tree',  'nbin':30, 'xmin':0, 'xmax':1, 'xtitle':'B vertex prob.', 'isLog':True}
vardir["b_alpha"] = {'tree':'tree',  'nbin':30, 'xmin':0.95, 'xmax':1, 'xtitle':'B cos(alpha)'}
vardir["b_alpha_zoom"] = {'tree':'tree',  'nbin':30, 'xmin':0.99, 'xmax':1, 'xtitle':'B cos(alpha)', 'var':'b_alpha'}
vardir["b_alpha_log"] = {'tree':'tree',  'nbin':30, 'xmin':0.95, 'xmax':1, 'xtitle':'B cos(alpha)', 'var':'b_alpha', 'isLog':True}
vardir["b_alpha_zoom_log"] = {'tree':'tree',  'nbin':30, 'xmin':0.99, 'xmax':1, 'xtitle':'B cos(alpha)', 'var':'b_alpha', 'isLog':True}
vardir["b_pvips"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':20., 'xtitle':'B pvips'}
vardir["b_fls3d"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':30., 'xtitle':'B flight length sig.'}
vardir["b_fl3d"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'B flight length'}
vardir["b_lips"] = {'tree':'tree',  'nbin':30, 'xmin':-10., 'xmax':10., 'xtitle':'B lips'}
vardir["b_mcorr"] = {'tree':'tree',  'nbin':30, 'xmin':4, 'xmax':12, 'xtitle':'B corrected mass (GeV)'}
#vardir["m_corr_bd"] = {'tree':'tree',  'nbin':30, 'xmin':3, 'xmax':15, 'xtitle':'B corrected mass bd (GeV)'}
vardir["b_mass"] = {'tree':'tree',  'nbin':35, 'xmin':4., 'xmax':6, 'xtitle':'B mass (GeV)'}
vardir["b_mass_high"] = {'tree':'tree',  'nbin':150, 'xmin':6.5, 'xmax':14, 'xtitle':'B mass high (GeV)', 'var':'b_mass'}
vardir["b_mass_low"] = {'tree':'tree',  'nbin':20, 'xmin':2., 'xmax':4, 'xtitle':'B mass low (GeV)', 'var':'b_mass'}
vardir["b_mass_simple"] = {'tree':'tree',  'nbin':35, 'xmin':4., 'xmax':6, 'xtitle':'B mass simple (GeV)'}
vardir["b_mass_sf"] = {'tree':'tree', 'nbin':50, 'xmin':3.5, 'xmax':8, 'xtitle':'B mass (GeV)', 'var':'b_mass'}
vardir["jpsi_kpipi"] = {'tree':'tree',  'nbin':50, 'xmin':5., 'xmax':5.5, 'xtitle':'J/#psi + K#pi#pi mass (GeV)'}
vardir["jpsi_k"] = {'tree':'tree', 'nbin':50, 'xmin':5., 'xmax':5.5, 'xtitle':'J/#psi + K mass (GeV)', 'var':'b_mass'}

vardir["dx_b_pv"] = {'tree':'tree',  'nbin':30, 'xmin':-0.4, 'xmax':0.4, 'xtitle':'#Deltax (B, bbPV) (cm)'}
vardir["dy_b_pv"] = {'tree':'tree',  'nbin':30, 'xmin':-0.4, 'xmax':0.4, 'xtitle':'#Deltay (B, bbPV) (cm)'}
vardir["dz_b_pv"] = {'tree':'tree',  'nbin':30, 'xmin':-0.4, 'xmax':0.4, 'xtitle':'#Deltaz (B, bbPV) (cm)'}
vardir["dr_b_pv"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':0.6, 'xtitle':'#Deltar (B, bbPV) (cm)'}

vardir["dx_jpsi_tau"] = {'tree':'tree',  'nbin':30, 'xmin':-0.4, 'xmax':0.4, 'xtitle':'#Deltax (J/psi, tau) (cm)'}
vardir["dy_jpsi_tau"] = {'tree':'tree',  'nbin':30, 'xmin':-0.4, 'xmax':0.4, 'xtitle':'#Deltay (J/psi, tau) (cm)'}
vardir["dz_jpsi_tau"] = {'tree':'tree',  'nbin':30, 'xmin':-0.4, 'xmax':0.4, 'xtitle':'#Deltaz (J/psi, tau) (cm)'}
vardir["dr_jpsi_tau"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':0.4, 'xtitle':'#Deltar (J/psi, tau) (cm)'}

vardir["ncand"] = {'tree':'tree',  'nbin':10, 'xmin':0., 'xmax':10., 'xtitle':'Nr. of tau cand.', 'var':'ncand'}

#vardir["met"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':70., 'xtitle':'MET', 'var':'met'}
#vardir["met_puppi"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':70., 'xtitle':'puppi MET', 'var':'met_puppi'}

vardir["jpsi_pt"] = {'tree':'tree',  'nbin':30, 'xmin':0, 'xmax':40, 'xtitle':'jpsi pT (GeV)'}
vardir["jpsi_eta"] = {'tree':'tree',  'nbin':30, 'xmin':-2.5, 'xmax':2.5, 'xtitle':'jpsi eta'}
vardir["jpsi_phi"] = {'tree':'tree',  'nbin':30, 'xmin':-math.pi, 'xmax':math.pi, 'xtitle':'jpsi phi'}
vardir["jpsi_fls3d"] = {'tree':'tree',  'nbin':30, 'xmin':0, 'xmax':30, 'xtitle':'jpsi flight sig.'}
vardir["jpsi_fl3d"] = {'tree':'tree',  'nbin':30, 'xmin':0, 'xmax':1., 'xtitle':'jpsi flight length'}
#vardir["jpsi_unfit_mass"] = {'tree':'tree',  'nbin':30, 'xmin':2.9, 'xmax':3.3, 'xtitle':'jpsi unfit mass'}
vardir["jpsi_mass"] = {'tree':'tree',  'nbin':30, 'xmin':2.9, 'xmax':3.3, 'xtitle':'jpsi mass'}


#if options.mva:
#vardir["bdt"] = {'tree':'tree',  'nbin':30, 'xmin':-1., 'xmax':1., 'xtitle':'BDT output score', 'var':'bdt'}
#vardir["dnn"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'DNN output score', 'var':'dnn'}

vardir["q2"] = {'tree':'tree',  'nbin':20, 'xmin':4., 'xmax':11., 'xtitle':'q2'}
vardir["mm2"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':6., 'xtitle':'missing mass squared'}
vardir["ptmiss"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':30., 'xtitle':'ptmiss'}
vardir["estar"] = {'tree':'tree',  'nbin':30, 'xmin':0.8, 'xmax':2.4, 'xtitle':'E*'}
vardir["B_ptback"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':70., 'xtitle':'B pT back (GeV)'}

vardir["q2_simple"] = {'tree':'tree',  'nbin':20, 'xmin':4., 'xmax':11., 'xtitle':'q2 simple'}
vardir["mm2_simple"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':6., 'xtitle':'missing mass squared simple'}
vardir["ptmiss_simple"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':30., 'xtitle':'ptmiss simple'}
vardir["estar_simple"] = {'tree':'tree',  'nbin':30, 'xmin':0.8, 'xmax':2.4, 'xtitle':'E* simple'}
vardir["B_ptback_simple"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':70., 'xtitle':'B pT back simple (GeV)'}

vardir["delta_chi2"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':25., 'xtitle':'delta chi2'}
vardir["vweight"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1.5, 'xtitle':'vweight'}
vardir["tau_fls3d_wjpsi"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':15, 'xtitle':'tau fls3d w.r.t j/psi'}
vardir["tau_fl3d_wjpsi"] = {'tree':'tree',  'nbin':30, 'xmin':-0.3, 'xmax':0.7, 'xtitle':'tau fl3d w.r.t j/psi'}

vardir["tau_refit_chi2"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':200., 'xtitle':'refitted PV chi2'}
vardir["tau_refit_ndof"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':200., 'xtitle':'refitted PV ndof'}
vardir["tau_refit_rho"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':2, 'xtitle':'refitted PV rho'}

vardir["pi1_pt"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':15, 'xtitle':'pi1 pt'}
vardir["pi1_eta"] = {'tree':'tree',  'nbin':30, 'xmin':-2.5, 'xmax':2.5, 'xtitle':'pi1 eta'}
vardir["pi1_dnn"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1, 'xtitle':'pi1 DNN'}
vardir["pi1_dnn_1prong"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'pi1 DNN 1prong'}
vardir["pi1_dnn_otherB"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'pi1 DNN otherB'}
vardir["pi1_dnn_pu"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'pi1 DNN pu'}
vardir["pi1_trigMatch"] = {'tree':'tree',  'nbin':2, 'xmin':0., 'xmax':2., 'xtitle':'pi1 trigMatch'}


vardir["pi2_pt"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':15, 'xtitle':'pi2 pt'}
vardir["pi2_eta"] = {'tree':'tree',  'nbin':30, 'xmin':-2.5, 'xmax':2.5, 'xtitle':'pi2 eta'}
vardir["pi2_dnn"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1, 'xtitle':'pi2 DNN'}
vardir["pi2_dnn_1prong"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'pi2 DNN 1prong'}
vardir["pi2_dnn_otherB"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'pi2 DNN otherB'}
vardir["pi2_dnn_pu"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'pi2 DNN pu'}
vardir["pi2_trigMatch"] = {'tree':'tree',  'nbin':2, 'xmin':0., 'xmax':2., 'xtitle':'pi2 trigMatch'}

vardir["pi3_pt"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':15, 'xtitle':'pi3 pt'}
vardir["pi3_eta"] = {'tree':'tree',  'nbin':30, 'xmin':-2.5, 'xmax':2.5, 'xtitle':'pi3 eta'}
vardir["pi3_dnn"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1, 'xtitle':'pi3 DNN'}
vardir["pi3_dnn_1prong"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'pi3 DNN 1prong'}
vardir["pi3_dnn_otherB"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'pi3 DNN otherB'}
vardir["pi3_dnn_pu"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':1., 'xtitle':'pi3 DNN pu'}
vardir["pi3_trigMatch"] = {'tree':'tree',  'nbin':2, 'xmin':0., 'xmax':2., 'xtitle':'pi3 trigMatch'}

vardir["anytrigMatch"] = {'tree':'tree',  'nbin':2, 'xmin':0., 'xmax':2., 'xtitle':'anytrigMatch', 'var':'(pi1_trigMatch==1 || pi2_trigMatch==1 || pi3_trigMatch==1)'}

vardir["ptbal"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':20., 'xtitle':'pT balance'}
vardir["jpsi_tau_alpha"] = {'tree':'tree',  'nbin':30, 'xmin':-1., 'xmax':1., 'xtitle':'alpha(j/psi, tau)'}
vardir["tau_dr_jpsi"] = {'tree':'tree',  'nbin':30, 'xmin':0., 'xmax':2*math.pi, 'xtitle':'#DeltaR(tau, J/psi)'}

#vardir["bdt_zoom"] = {'tree':'tree',  'nbin':30, 'xmin':0.5, 'xmax':1., 'xtitle':'zoom BDT output score', 'var':'bdt'}
#vardir["dnn_zoom"] = {'tree':'tree',  'nbin':30, 'xmin':0.7, 'xmax':1., 'xtitle':'zoom DNN output score', 'var':'dnn'}

#vardir["bdt_zoom2"] = {'tree':'tree',  'nbin':30, 'xmin':0.9, 'xmax':1., 'xtitle':'zoom2 BDT output score', 'var':'bdt'}
#vardir["dnn_zoom2"] = {'tree':'tree',  'nbin':30, 'xmin':0.9, 'xmax':1., 'xtitle':'zoom2 DNN output score', 'var':'dnn'}

#vardir["bdt_zoom3"] = {'tree':'tree',  'nbin':30, 'xmin':0.9, 'xmax':1., 'xtitle':'zoom2 BDT output score', 'var':'bdt'}
#vardir["dnn_zoom3"] = {'tree':'tree',  'nbin':30, 'xmin':0.9, 'xmax':1., 'xtitle':'zoom2 DNN output score', 'var':'dnn'}

#if options.datacard:
#vardir["xgbs"] = {'tree':'tree',  'nbin':30, 'xmin':8.5, 'xmax':13, 'xtitle':'zoom XGB output score', 'var':'xgbs'}
#vardir["xgbs_log"] = {'tree':'tree',  'nbin':30, 'xmin':8.5, 'xmax':13, 'xtitle':'zoom XGB output score', 'var':'xgbs', 'isLog':True}
#vardir["xgbs_wide"] = {'tree':'tree',  'nbin':30, 'xmin':8.5, 'xmax':13, 'xtitle':'XGB output score', 'var':'xgbs'}

#vardir["xgbs_fine"] = {'tree':'tree',  'nbin':30, 'xmin':6, 'xmax':12, 'xtitle':'XGB output score fine', 'var':'xgbs', 'isRight':True}
vardir["xgbs"] = {'tree':'tree',  'nbin':100, 'xmin':-15, 'xmax':9, 'xtitle':'XGB output score', 'var':'xgbs', 'isRight':True}

vardir["xgbs_zoom"] = {'tree':'tree',  'nbin':60, 'xmin':1.1, 'xmax':9.1, 'xtitle':'XGB output score', 'var':'xgbs', 'isRight':True}
vardir["xgbs_zoom_extra"] = {'tree':'tree',  'nbin':60, 'xmin':1.1, 'xmax':9.1, 'xtitle':'XGB output score', 'var':'xgbs', 'isRight':True}
vardir["xgbs_zoom_thight"] = {'tree':'tree',  'nbin':60, 'xmin':3.1, 'xmax':9.1, 'xtitle':'XGB output score', 'var':'xgbs', 'isRight':True}
vardir["xgbs_sigscan"] = {'tree':'tree',  'nbin':600, 'xmin':3.1, 'xmax':7., 'xtitle':'XGB output score', 'var':'xgbs', 'isRight':True}
vardir["xgbs_fit"] = {'tree':'tree',  'nbin':50, 'xmin':2., 'xmax':7., 'xtitle':'XGB output score', 'var':'xgbs', 'isRight':True}

#vardir["xgbs_log"] = {'tree':'tree',  'nbin':30, 'xmin':-15, 'xmax':9, 'xtitle':'XGB output score', 'var':'xgbs', 'isLog':True}
#vardir["xgbs_zoom"] = {'tree':'tree',  'nbin':30, 'xmin':9., 'xmax':11, 'xtitle':'zoom XGB output score', 'var':'xgbs', 'isLog':True}
#vardir["xgbs_zoom"] = {'tree':'tree',  'nbin':30, 'xmin':7., 'xmax':12, 'xtitle':'zoom XGB output score', 'var':'xgbs', 'isLog':True}
#vardir["xgbs_wide"] = {'tree':'tree',  'nbin':30, 'xmin':-15, 'xmax':13, 'xtitle':'XGB output score', 'var':'xgbs'}
#vardir["perEVT_old"] = {'tree':'tree',  'nbin':20, 'xmin':0., 'xmax':1, 'xtitle':'perEVT DNN', 'isRight':True}
#vardir["perEVT_old_log"] = {'tree':'tree',  'nbin':30, 'xmin':0, 'xmax':1, 'xtitle':'perEVT DNN log', 'var':'perEVT_old', 'isLog':True}
#vardir["perEVT_mc"] = {'tree':'tree',  'nbin':20, 'xmin':0., 'xmax':1, 'xtitle':'perEVT mc'}
#vardir["perEVT_data"] = {'tree':'tree',  'nbin':20, 'xmin':0., 'xmax':1, 'xtitle':'perEVT data'}
#vardir["perEVT_otherB"] = {'tree':'tree',  'nbin':20, 'xmin':0., 'xmax':1, 'xtitle':'perEVT DNN otherB'}
#vardir["perEVT_leptonic"] = {'tree':'tree',  'nbin':20, 'xmin':0., 'xmax':1, 'xtitle':'perEVT DNN leptonic'}
#vardir["perEVT_1prong"] = {'tree':'tree',  'nbin':20, 'xmin':0., 'xmax':1, 'xtitle':'perEVT DNN 1prong'}

vardir["mu1_pt"] = {'tree':'tree',  'nbin':30, 'xmin':0, 'xmax':40, 'xtitle':'mu1 pT (GeV)'}
vardir["mu1_eta"] = {'tree':'tree',  'nbin':30, 'xmin':-2.5, 'xmax':2.5, 'xtitle':'mu1 eta'}
vardir["mu1_phi"] = {'tree':'tree',  'nbin':30, 'xmin':-math.pi, 'xmax':math.pi, 'xtitle':'mu1 phi'}
vardir["mu1_isLoose"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu1 isLoose', 'isLog':True}
vardir["mu1_isTight"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu1 isTight', 'isLog':True}
vardir["mu1_isPF"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu1 isPF', 'isLog':True}
vardir["mu1_isGlobal"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu1 isGlobal', 'isLog':True}
vardir["mu1_isTracker"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu1 isTracker', 'isLog':True}
vardir["mu1_isSoft"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu1 isSoft', 'isLog':True}
vardir["mu1_dbiso"] = {'tree':'tree',  'nbin':30, 'xmin':0, 'xmax':30, 'xtitle':'mu1 dbiso'}

vardir["mu2_pt"] = {'tree':'tree',  'nbin':30, 'xmin':0, 'xmax':40, 'xtitle':'mu2 pT (GeV)'}
vardir["mu2_eta"] = {'tree':'tree',  'nbin':30, 'xmin':-2.5, 'xmax':2.5, 'xtitle':'mu2 eta'}
vardir["mu2_phi"] = {'tree':'tree',  'nbin':30, 'xmin':-math.pi, 'xmax':math.pi, 'xtitle':'mu2 phi'}
vardir["mu2_isLoose"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu2 isLoose', 'isLog':True}
vardir["mu2_isTight"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu2 isTight', 'isLog':True}
vardir["mu2_isPF"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu2 isPF', 'isLog':True}
vardir["mu2_isGlobal"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu2 isGlobal', 'isLog':True}
vardir["mu2_isTracker"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu2 isTracker', 'isLog':True}
vardir["mu2_isSoft"] = {'tree':'tree',  'nbin':2, 'xmin':0, 'xmax':2, 'xtitle':'mu2 isSoft', 'isLog':True}
vardir["mu2_dbiso"] = {'tree':'tree',  'nbin':30, 'xmin':0, 'xmax':30, 'xtitle':'mu2 dbiso'}


#vardir["xgbs"] = {'tree':'tree',  'nbin':15, 'xmin':6.5, 'xmax':10., 'xtitle':'zoom XGB output score', 'var':'xgbs'}
#vardir["xgbs"] = {'tree':'tree',  'nbin':12, 'xmin':6.8, 'xmax':8.8, 'xtitle':'zoom XGB output score', 'var':'xgbs'}


