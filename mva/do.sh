today=`date "+%Y%m%d%H%M"`

#name="nomass_pt"
name="pt_Legacy"
#name="multiple"

# location of the input files
prefix="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${name}"


sig_inclusive_ul_all_2018="${prefix}/BcJpsiTau_inclusive_ul_all_2018/Myroot_training_weightAdded.root" 
#bkg_ul_2018="${prefix}/BcJpsiX_ul_2018/Myroot_training_default.root"
data_ul_2018="${prefix}/Data_2018/Myroot_training_weightAdded.root"

##############################################
#
# Build analysis BDT based on xgboos
#
#############################################

source $PWD/setup.sh 

ls -lart ${sig_inclusive_ul_all_2018}
#ls -lart ${bkg_ul_2018}
ls -lart ${data_ul_2018}

#python TrainModel_XGB.py -s ${sig_inclusive_ul_all_2018} -b ${bkg_ul_2018} --dir ${name}
python TrainModel_XGB.py -s ${sig_inclusive_ul_all_2018} -b ${data_ul_2018} --dir ${name}_val -o
