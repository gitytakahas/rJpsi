today=`date "+%Y%m%d%H%M"`

#name="mass_pt"
#name="nomass_pt"
name="multiple"

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

oname="${name}"

#python application.py --file ${bkg_ul_2018}/Myroot_analysis.root --prefix bkg_ul_2018 --model ${oname} --outdir $PWD/final_root_${oname}
#python application.py --file ${sig_inclusive_ul_all_2018}/Myroot_analysis.root --prefix sig_ul --model ${oname} --outdir $PWD/final_root_${oname}

#python getDataset_data.py --path ${data_ul_2018} --odir ${data_ul_2018} --jdir data_application_${today} --name data_xgbs6_2018 --chunk 5 --model $oname


##############################################
#
# 2. combine
#
##############################################


for prio in ${name}
do
    echo ${prio}
    cd $PWD/final_root_${prio}
#    hadd -f /pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${prio}/Data_2018/data.root Myroot_data_2018_*
    hadd -f /pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${prio}/BcJpsiTau_inclusive_ul_all_2018/signal.root  Myroot_sig_ul_*.root 
    hadd -f /pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${prio}/BcJpsiX_ul_2018/bkg.root Myroot_bkg_ul_2018_*.root

    cd -;
done




#python getDataset_data.py --path ${sig_inclusive_ul_all_2018} --odir ${sig_inclusive_ul_all_2018} --jdir signal_application_${today} --name signal_xgbs6_2018 --chunk 5 --model $oname
#python getDataset_data.py --path ${bkg_ul_2018} --odir ${bkg_ul_2018} --jdir bg_application_${today} --name bg_xgbs6_2018 --chunk 5 --model $oname
#python application.py --file ${data_ul_2018}/Myroot.root --prefix data_2018 --model ${oname} --outdir $PWD/final_root_${oname}
