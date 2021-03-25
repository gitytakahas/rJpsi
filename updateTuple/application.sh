today=`date "+%Y%m%d%H%M"`

name="pt"
#name="multiple"

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

python application.py --file ${bkg_ul_2018}/Myroot_anal.root --prefix bkg_ul_2018 --model ${oname} --outdir $PWD/final_root_${oname}
python application.py --file ${sig_inclusive_ul_all_2018}/Myroot_anal.root --prefix sig_ul --model ${oname} --outdir $PWD/final_root_${oname}
python application.py --file ${data_ul_2018}/Myroot.root --prefix data_2018 --model ${oname} --outdir $PWD/final_root_${oname}


##############################################
#
# 2. combine
#
##############################################


for prio in ${name}
do
    echo ${prio}
    cd $PWD/final_root_${prio}
#    hadd -f /pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${prio}/Data_2018/Myroot_data_2018.root Myroot_data_2018_*
#    hadd -f /pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${prio}/BcJpsiTau_ul_all_2018/Myroot_sig_2018.root  Myroot_sig_new_-1_*.root 
#    hadd -f /pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_${prio}/BcJpsiX_rereco_2018/Myroot_bkg_rereco_2018.root Myroot_bkg_rereco_2018_*.root

    hadd -f data.root Myroot_data_2018_*.root
    hadd -f signal.root  Myroot_sig_ul_*.root
    hadd -f bkg.root Myroot_bkg_ul_2018_*.root
    cd -;
done


