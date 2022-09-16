today=`date "+%Y%m%d%H%M"`


# Specify here the output file types 

name="pt"
year="2016"

prefix="/pnfs/psi.ch/cms/trivcat/store/user/${USER}/RJpsi/job_${name}_${year}"

sig_inclusive_all="${prefix}/BcJpsiTau_inclusive"
bkg="${prefix}/BJpsiX"
bkg_new="${prefix}/BJpsiX_inclusive"
data="${prefix}/Data"


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

for dir in $sig_inclusive_all $bkg $data
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

if [ $year = "2018" ]; then

    hadding $bkg_new Myroot.root Myroot_*.root

    for dir in $sig_inclusive_all
    do
	echo "splitting into two ...", $dir
	sh updateTuple_pick.sh ${dir}/Myroot.root ${dir}/Myroot_analysis.root ${dir}/Myroot_training.root 0.2
    done


    ##############################################
    #
    # 3. Create pt/eta map for the TMVA training
    # 
    # to avoid biasing, we can reweigh based on B pT and eta
    #
    ##############################################
    python create_weights.py --sig_file ${sig_inclusive_all}/Myroot_training.root --bkg_file ${bkg_new}/Myroot.root --out_dir weight_${name}_${year}


    ##############################################
    #
    # 4. add pt/eta weight calculated above
    #    -> add "weight" variable that contains all weightings (including pu weight, BkgB weights)
    # We do it for both training files and analysis files.
    # 
    ##############################################
    sh updateTuple.sh ${bkg_new}/Myroot.root $PWD/weight_${name}_${year}/weight.root

fi




# obsolete

#sh updateTuple.sh ${data}/Myroot_training.root $PWD/weight_${name}/weight.root
#sh updateTuple.sh ${sig_inclusive_all}/Myroot_training.root $PWD/weight_${name}/weight.root
#sh updateTuple.sh ${bkg_2018}/Myroot_training.root $PWD/weight_${name}/weight.root ${bkg_2018} -1
#sh updateTuple.sh ${bkg_2018}/Myroot_analysis.root $PWD/weight_${name}/weight.root ${bkg_2018} -1
#sh updateTuple.sh ${bkg_2018}/Myroot.root $PWD/weight_${name}/weight.root ${bkg_2018} -1
#sh updateTuple.sh ${sig_inclusive_all_2018}/Myroot_analysis.root $PWD/weight_pt/weight.root ${sig_inclusive_all_2018} -1
#python getDataset.py --path ${sig_inclusive_all_2018}/Myroot_analysis.root --wfile $PWD/weight_${name}/weight.root --odir signal_${today} --chunk 5

