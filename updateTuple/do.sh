today=`date "+%Y%m%d%H%M"`


# Specify here the output file types 
#name="mass_pt"

name="pt_LEGACY_v2"
#name="pt"
#name="multiple"

prefix="/pnfs/psi.ch/cms/trivcat/store/user/${USER}/RJpsi_Legacy_decayBc/job_${name}"

sig_inclusive_ul_all_2018="${prefix}/BcJpsiTau_inclusive_ul_all_2018"
bkg_ul_2018="${prefix}/BJpsiX_ul_2018"
data_ul_2018="${prefix}/Data_2018"


function hadding (){

    dir=$1
    target=$2
    src=$3
    
    hadd -f ${dir}/$target ${dir}/${src}
}


##############################################
#
# 0. hadd files
#
##############################################

#for dir in $sig_inclusive_ul_all_2018 $bkg_ul_2018 # $data_ul_2018         
for dir in $bkg_ul_2018 # $data_ul_2018
do

    echo "hadding ...", $dir
    hadding $dir Myroot.root Myroot_*.root

done


##############################################
#
# 1. Split BG and signal MC with a fraction of 0.2 and 0.8. 
# 
# split is done based on random number so that the events will be 
# equally distributed
# 
# "0.2 sample" will be used for training analysis BDT (called "training files")
# "0.8 sample" will be used for the analysis (i.e. building the template) (called "analysis files")
# In this way, we can keep orthogonality. 
#
##############################################

#for dir in $sig_inclusive_ul_all_2018 $bkg_ul_2018
for dir in $sig_inclusive_ul_all_2018
do
    echo "splitting into two ...", $dir
    #sh updateTuple_pick.sh ${dir}/Myroot.root ${dir}/Myroot_analysis.root ${dir}/Myroot_training.root 0.2
done

for dir in $data_ul_2018
do
    echo "splitting into two ...", $dir
    #sh updateTuple_pick.sh ${dir}/Myroot.root ${dir}/Myroot_analysis.root ${dir}/Myroot_training.root 0.006
done


##############################################
#
# 3. Create pt/eta map for the TMVA training
# 
# to avoid biasing, we can reweigh based on B pT and eta
#
##############################################
#python create_weights.py --sig_file ${sig_inclusive_ul_all_2018}/Myroot_training.root --bkg_file ${data_ul_2018}/Myroot_training.root --out_dir weight_${name}


##############################################
#
# 4. add pt/eta weight calculated above
#    -> add "weight" variable that contains all weightings (including pu weight, BkgB weights)
# We do it for both training files and analysis files.
# 
##############################################

#sh updateTuple.sh ${sig_inclusive_ul_all_2018}/Myroot_training.root $PWD/weight_${name}/weight.root
#sh updateTuple.sh ${data_ul_2018}/Myroot_training.root $PWD/weight_${name}/weight.root























#sh updateTuple.sh ${bkg_ul_2018}/Myroot_training.root $PWD/weight_${name}/weight.root ${bkg_ul_2018} -1
#sh updateTuple.sh ${bkg_ul_2018}/Myroot_analysis.root $PWD/weight_${name}/weight.root ${bkg_ul_2018} -1
#sh updateTuple.sh ${bkg_ul_2018}/Myroot.root $PWD/weight_${name}/weight.root ${bkg_ul_2018} -1
#sh updateTuple.sh ${sig_inclusive_ul_all_2018}/Myroot_analysis.root $PWD/weight_pt/weight.root ${sig_inclusive_ul_all_2018} -1
#python getDataset.py --path ${sig_inclusive_ul_all_2018}/Myroot_analysis.root --wfile $PWD/weight_${name}/weight.root --odir signal_${today} --chunk 5

