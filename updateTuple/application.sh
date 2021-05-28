today=`date "+%Y%m%d%H%M"`

name="pt"
#name="nomass_pt"
#name="multiple"

# input files
prefix="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${name}"

sig_inclusive_ul_all_2018="${prefix}_hammer/BcJpsiTau_inclusive_ul_all_2018"
bkg_ul_2018="${prefix}/BcJpsiX_ul_2018"
data_ul_2018="${prefix}/Data_2018"

##############################################
#
# 1. update analysis files with signal_all to split into 1000 chunks (depending on FF)
# At this stage, you should finish MVA training !!!!
#
##############################################

oname="${name}"

#python application.py --file ${bkg_ul_2018}/Myroot_analysis_default.root --prefix bkg_xgbs6_default --model ${oname} --outdir ${bkg_ul_2018}
#python application.py --file ${sig_inclusive_ul_all_2018}/Myroot_analysis_default.root --prefix sig_xgbs6_default --model ${oname} --outdir ${sig_inclusive_ul_all_2018}



##############################################
#
# 2. combine
#
##############################################


for prio in ${name}
do
    echo ${prio}
#    hadd -f /pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${prio}/Data_2018/data.root /pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${prio}/Data_2018/Myroot_data_xgbs6_2018_*.root
    hadd -f ${sig_inclusive_ul_all_2018}/signal_xgbs6.root  ${sig_inclusive_ul_all_2018}/Myroot_sig_xgbs6_default*.root 
    hadd -f ${bkg_ul_2018}/bkg_xgbs6.root ${bkg_ul_2018}/Myroot_bkg_xgbs6_default*.root

#    cd -;
done


#python getDataset_data.py --path ${data_ul_2018} --odir ${data_ul_2018} --jdir data_application_${today} --name data_xgbs6_2018 --chunk 5 --model $oname
#python getDataset_data.py --path ${sig_inclusive_ul_all_2018} --odir ${sig_inclusive_ul_all_2018} --jdir signal_hammer_application_${today} --name sig_xgbs6_hammer --chunk 5 --model $oname --signal


