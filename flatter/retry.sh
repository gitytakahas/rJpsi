for id in 443 447
#for id in 177 200 139 205 318 120 307 118 54 62 45 92 112 157 86 80 121 140 152 135 72 154 109 76 134 79 303 49 156 151 73 123 56 57
do 

    sbatch -p short --account=t3 --error=/work/ytakahas/work/analysis/CMSSW_10_2_10/src/rJpsi/flatter/job_pt/Data_2018/err.${id} --output=/work/ytakahas/work/analysis/CMSSW_10_2_10/src/rJpsi/flatter/job_pt/Data_2018/out.${id} /work/ytakahas/work/analysis/CMSSW_10_2_10/src/rJpsi/flatter/job_pt/Data_2018/job_${id}.sh

done

