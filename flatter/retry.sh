for id in 337 336 356 122 326 33 7 92 147 316 331 174 334 2 118 137 228 189 220 
#for id in 177 200 139 205 318 120 307 118 54 62 45 92 112 157 86 80 121 140 152 135 72 154 109 76 134 79 303 49 156 151 73 123 56 57
do 

    sbatch -p short --account=t3 --error=/work/ytakahas/work/analysis/CMSSW_10_2_10/src/rJpsi/flatter/job_pt_BMuMuK/Data_2018/err.${id} --output=/work/ytakahas/work/analysis/CMSSW_10_2_10/src/rJpsi/flatter/job_pt_BMuMuK/Data_2018/out.${id} /work/ytakahas/work/analysis/CMSSW_10_2_10/src/rJpsi/flatter/job_pt_BMuMuK/Data_2018/job_${id}.sh

done

