today=`date "+%Y%m%d%H%M"`

#name="nomass_pt"
#name="pt_Legacy"
name="pt"
year="2018"
#name="multiple"

# location of the input files
prefix="/pnfs/psi.ch/cms/trivcat/store/user/${USER}/RJpsi/job_${name}_${year}"


sig_inclusive_all="${prefix}/BcJpsiTau_inclusive/Myroot_training.root" 
#sig_inclusive_ul_all_2018="${prefix}/BcJpsiTau_inclusive_ul_all_2018/Myroot.root" 
bkg_new="${prefix}/BJpsiX_inclusive/Myroot_weightAdded.root"
#data_ul_2018="${prefix}/Data_2018/Myroot_training_weightAdded.root"

##############################################
#
# Build analysis BDT based on xgboos
#
#############################################

source $PWD/setup.sh 

ls -lart ${sig_inclusive_all}
ls -lart ${bkg_new}

python TrainModel_XGB.py -s ${sig_inclusive_all} -b ${bkg_new} --dir ${name}_${year}_ref


#python TrainModel_XGB.py -s ${sig_inclusive_all} -b ${bkg_new} --dir ${name}_${year}_val -o
#python TrainModel_XGB_Dbkg.py -s ${sig_inclusive_all} -b ${bkg_new} --dir ${name}_${year}_test1










#python TrainModel_XGB.py -s ${sig_inclusive_ul_all_2018} -b ${bkg_ul_2018} --dir ${name}
#python TrainModel_XGB.py -s ${sig_inclusive_ul_all_2018} -b ${bkg_ul_2018_new} --dir ${name}_orig_novalid_morestat -o
#python TrainModel_XGB.py -s ${sig_inclusive_ul_all_2018} -b ${bkg_ul_2018_new} --dir ${name}_test
#python TrainModel_XGB.py -s ${sig_inclusive_ul_all_2018} -b ${bkg_ul_2018_new} --dir ${name}_test_sumofdnn2
#python TrainModel_XGB.py -s ${sig_inclusive_ul_all_2018} -b ${bkg_ul_2018_new} --dir ${name}_test_removed_ncand_iso

