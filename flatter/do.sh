####################
# This is the flatter script
####################

# 
# Location where output files will be stored
#
pnfs="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/"


#
# analysis to be run runTauDisplay_XXX.py
# XXX is the name here. 
#
analysis="BcJpsiTauNu"


#
# Write here the "top directory" name of the Ntuplizer files
# The rest of files under this directory will be processed
# 
# You can get it by doing 
# uberftp -ls gsiftp://storage01.lcg.cscs.ch//pnfs/lcg.cscs.ch/cms/trivcat/store/user/ytakahas/
# 
sigmc="BcToJPsiMuMu_Legacy_2018_20210509"
bgmc="HbToJPsiMuMu_legacy_2018_20210331"
#dataset="Charmonium_legacy_2018_20210331"
dataset="Charmonium_Legacy_q3_2018_20210518"


priority="pt"
nchunk_sig=10
nchunk_bg=10
nchunk_data=5



# uncomment here in case of "multpile" options
#nchunk_sig=5
#nchunk_bg=3
#nchunk_data=3
#priority="multiple"





outdir="job_${priority}_q3"

#########################################
# for signal MC
#########################################
for year in 2018
do

    pustr=$year

    if [ $year = "2017" ]; then
	pustr="UL$year"
    fi

    echo "signal", $year

#    python getDataset.py --file ${sigmc} --chunk ${nchunk_sig} --analysis ${analysis} --type signal --name BcJpsiTau_inclusive_ul_all_${year} --select UL --year $pustr --priority ${priority} --odir ${pnfs} --jdir ${outdir}


done


#########################################
# for J/psi + X BG
#########################################
#python getDataset.py --file ${bgmc} --chunk ${nchunk_bg} --analysis ${analysis} --type bg --name BcJpsiX_ul_2018 --year 2018 --priority ${priority} --odir ${pnfs} --jdir ${outdir}

#########################################
# Data (2016)
#########################################

for year in 2018
do
    echo "data", $year
    python getDataset.py --file ${dataset} --chunk ${nchunk_data} --analysis ${analysis} --type data --name Data_${year} --priority ${priority} --odir ${pnfs} --jdir ${outdir}
done




