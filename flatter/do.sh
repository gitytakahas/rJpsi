####################
# This is the flatter script
####################

pnfs="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/"

analysis="BcJpsiTauNu"


sigmc="BcToJPsiMuMu_legacy_mc_2018_20210322"
bgmc="HbToJPsiMuMu_legacy_2018_20210322"
dataset="Charmonium_Data_legacy_2018_20210321"


nchunk_sig=10
nchunk_bg=10
nchunk_data=5
priority="pt"


#nchunk_sig=5
#nchunk_bg=3
#nchunk_data=3
#priority="multiple"

outdir="job_${priority}"

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

    python getDataset.py --file ${sigmc} --chunk ${nchunk_sig} --analysis ${analysis} --type signal --name BcJpsiTau_inclusive_ul_all_${year} --select UL --year $pustr --priority ${priority} --odir ${pnfs} --jdir ${outdir}


done


#########################################
# for J/psi + X BG
#########################################
python getDataset.py --file ${bgmc} --chunk ${nchunk_bg} --analysis ${analysis} --type bg --name BcJpsiX_ul_2018 --year 2018 --priority ${priority} --odir ${pnfs} --jdir ${outdir}

#########################################
# Data (2016)
#########################################

for year in 2018
do
    echo "data", $year
    python getDataset.py --file ${dataset} --chunk ${nchunk_data} --analysis ${analysis} --type data --name Data_${year} --priority ${priority} --odir ${pnfs} --jdir ${outdir}
done




