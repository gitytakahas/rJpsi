today=`date "+%Y%m%d%H%M"`

#name="pt"
#name="pt_vprobfsigcr"
name="pt_Legacy"
#name="nomass_pt"
#name="multiple"

# input files
prefix="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${name}"

sig_inclusive_ul_all_2018="${prefix}/BcJpsiTau_inclusive_ul_all_2018"
bkg_ul_2018="${prefix}/BcJpsiX_ul_2018"
data_ul_2018="${prefix}/Data_2018"

##############################################
#
# 1. update analysis files with signal_all to split into 1000 chunks (depending on FF)
# At this stage, you should finish MVA training !!!!
#
##############################################

model_name="pt_Legacy"

#python application.py --file ${bkg_ul_2018}/Myroot.root --prefix bkg_xgbs --model ${model_name} --outdir ${bkg_ul_2018}
#python application.py --file ${sig_inclusive_ul_all_2018}/Myroot_analysis.root --prefix sig_xgbs --model ${model_name} --outdir ${sig_inclusive_ul_all_2018}
#python getDataset_data.py --path ${data_ul_2018} --odir ${data_ul_2018} --jdir data_application_${today} --name data_xgbs_2018 --chunk 5 --model $model_name


#python application.py --file ${data_ul_2018}/Myroot_analysis.root --prefix data_xgbs --model ${model_name} --outdir ${data_ul_2018}

##############################################
#
# 2. combine
#
##############################################


#for prio in ${name}
#do
#    echo ${prio}
hadd -f ${data_ul_2018}/data.root ${data_ul_2018}/Myroot_*xgbs*.root
hadd -f ${sig_inclusive_ul_all_2018}/sig.root  ${sig_inclusive_ul_all_2018}/Myroot_*xgbs*.root 
hadd -f ${bkg_ul_2018}/bkg.root ${bkg_ul_2018}/Myroot_*xgbs*.root

#    cd -;
#done




#python getDataset_data.py --path ${sig_inclusive_ul_all_2018} --odir ${sig_inclusive_ul_all_2018} --jdir signal_hammer_application_${today} --name sig_xgbs6_hammer --chunk 5 --model $model_name --signal


