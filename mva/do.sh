today=`date "+%Y%m%d%H%M"`

#name="multiple"
name="pt"

prefix="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${name}"


sig_inclusive_ul_all_2018="${prefix}/BcJpsiTau_inclusive_ul_all_2018" 
bkg_ul_2018="${prefix}/BcJpsiX_ul_2018"
data_ul_2018="${prefix}/Data_2018"


##############################################
#
# Build analysis BDT based on xgboos
#
#############################################

source $PWD/setup.sh 

ls -lart ${sig_inclusive_ul_all_2018}/Myroot_training.root
ls -lart ${bkg_ul_2018}/Myroot_training.root

python TrainModel_XGB.py -s ${sig_inclusive_ul_all_2018}/Myroot_training.root -b ${bkg_ul_2018}/Myroot_training.root --dir ${name}
