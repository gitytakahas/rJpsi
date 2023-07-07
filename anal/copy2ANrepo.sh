prefix="/work/ytakahas/work/Note/61e532a4361d0d70e5b2cf63/figures/"
pnfs_prefix="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/results_simultaneous/"


####### for inclusive distribution !!! 

for year in 2016 2017 2018
do
    for var in tau_pt tau_eta tau_fl3d tau_fls3d_wjpsi tau_vprob tau_mass tau_rhomass1 tau_rhomass2 tau_sumofdnn tau_sumofdnn_1prong tau_sumofdnn_otherB tau_sumofdnn_pu jpsi_pt jpsi_eta jpsi_mass jpsi_fl3d jpsi_fls3d dr_jpsi_tau b_pt b_eta b_mass b_alpha b_vprob b_lips b_pvips b_mindoca dr_b_pv b_fls3d xgbs xgbs_compare #sig_xgbs_zoom
    do
	echo ${var}
	cp ./${year}_inclusive_None/plots/${var}.pdf ${prefix}/${year}_inclusive_None/plots
    done
done




####### for MVA input variable 
#for year in 2016 2017 2018
#do
#    for var in b_pt b_eta b_alpha b_vprob tau_iso_0p7 b_lips b_pvips b_mindoca dr_b_pv tau_fls3d tau_vprob tau_fls3d_wjpsi tau_sumofdnn tau_sumofdnn_1prong tau_sumofdnn_otherB tau_sumofdnn_pu ncand estar
#    do
#	echo ${var}
#	#    cp 2018_inclusive_None/plots/${var}_compare.pdf ${prefix}/compare_mva/
#    done
#done



for reg in sr sb
do
    for var in bc_others bg_ul data_obs jpsi_hc jpsi_tau
    do
	echo ${var}, ${reg}
	cp display/${var}_${reg}.pdf ${prefix}/shape_comparison

    done
done


for year in 2016 2017 2018
do
    echo $year

    cp combine_inv/sb_sr_data_obs/tau_rhomass_unrolled_var_${year}.pdf ${prefix}/combine_sb3p5_sr4/inv_sb_sr_data_obs    
    cp combine/sb_sr_bg_ul/tau_rhomass_unrolled_var_${year}.pdf ${prefix}/combine_sb3p5_sr4/sb_sr_bg_ul/
    cp combine/lp_sb_data_obs/tau_rhomass_unrolled_var_${year}.pdf ${prefix}/combine_sb3p5_sr4/lp_sb_data_obs/
    cp ${pnfs_prefix}/${year}_gap_None/plots/tau_rhomass_unrolled_var.pdf ${prefix}/${year}_gap_None/plots/ 
    cp ${pnfs_prefix}/${year}_sr_None/plots/tau_rhomass_unrolled_var.pdf ${prefix}/${year}_sr_None/plots/ 
    cp ${pnfs_prefix}/${year}_sb_None/plots/tau_rhomass_unrolled_var.pdf ${prefix}/${year}_sb_None/plots/ 
done



for year in 2016 2017 2018
do
    for var in tau_pt tau_eta tau_fl3d tau_fls3d_wjpsi tau_vprob tau_mass tau_rhomass1 tau_rhomass2 tau_sumofdnn tau_sumofdnn_1prong tau_sumofdnn_otherB tau_sumofdnn_pu jpsi_pt jpsi_eta jpsi_mass jpsi_fl3d jpsi_fls3d dr_jpsi_tau b_pt b_eta b_mass b_alpha b_vprob b_lips b_pvips b_mindoca dr_b_pv b_fls3d 
    do
	echo ${year}, ${var}
	cp ./${year}_sr_None/plots/${var}_coarse.pdf ${prefix}/${year}_sr_None/plots/
	cp ./${year}_sb_None/plots/${var}_coarse.pdf ${prefix}/${year}_sb_None/plots/
    done
done


#comparision between gap region and the 
cp -r bkgcorr/*.pdf ${prefix}/bkg_shape_unc/

cp plots/compare_20*.pdf ${prefix}
cp plots/ratio_compare_20*.pdf ${prefix}
cp plots/compare_inv_20*.pdf ${prefix}
cp plots/ratio_compare_inv_20*.pdf ${prefix}

cp -r combine_simultaneous/syscompare $prefix/


#for year in 2016 2017 2018
#do
#    echo ${year}
#    cp inv_${year}_inclusive_None/plots/xgbs.pdf ${prefix}/${year}_xgbs_inv.pdf
#
#done


