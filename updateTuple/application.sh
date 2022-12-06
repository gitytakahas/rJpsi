today=`date "+%Y%m%d%H%M"`

name="pt"
year="2016"


# input files
prefix="/pnfs/psi.ch/cms/trivcat/store/user/${USER}/RJpsi/job_${name}_${year}"

sig_inclusive="${prefix}/BcJpsiTau_inclusive"
bkg="${prefix}/BJpsiX"
#bkg_ul_${year}_new="${prefix}/BJpsiX_ul_${year}_new"
data="${prefix}/Data"

##############################################
#
# 1. update analysis files with signal_all to split into 1000 chunks (depending on FF)
# At this stage, you should finish MVA training !!!!
#
##############################################

model_name="${name}_2018_val"

python application.py --file ${bkg}/Myroot.root --prefix bkg_xgbs --model ${model_name} --outdir ${bkg}

if [ $year = "2018" ]
then
    python application.py --file ${sig_inclusive}/Myroot_analysis.root --prefix sig_xgbs --model ${model_name} --outdir ${sig_inclusive}
else
    python application.py --file ${sig_inclusive}/Myroot.root --prefix sig_xgbs --model ${model_name} --outdir ${sig_inclusive}
fi

#python getDataset_data.py --path ${data} --odir ${data} --jdir data_application_${today}_${year} --name data_xgbs_${year} --chunk 1 --model $model_name


##############################################
#
# 2. combine
#
##############################################


#hadd -f ${data}/data.root ${data}/Myroot_*xgbs*.root
hadd -f ${bkg}/bkg.root ${bkg}/Myroot_*xgbs*.root
hadd -f ${sig_inclusive}/sig.root  ${sig_inclusive}/Myroot_*xgbs*.root 

