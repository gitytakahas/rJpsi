#!/bin/bash
# update each ntuple by adding QCD fake weight and the MVA for ttbar rejection


if [ $# -ne 4 ]; then
  echo "# of specified params: #" 1>&2
  echo "You need to give three param" 1>&2
  exit 1
fi


#basedir="/work/ytakahas/work/BsTauTau/CMSSW_10_2_10/src/job/"

#for dir in BcJpsiTau_large_BcJpsiTauNu_truth_2020-04-24-150429 BJpsiX_BcJpsiTauNu_mc_2020-04-24-150431
#do

#ifile="${basedir}/${dir}/Myroot.root";
#ofile="${basedir}/${dir}/Myroot4mva.root";

ifile="$1"
ofile="$2"
ofile_training="$3"
frac="$4"

#    ifile="/work/ytakahas/work/BsTauTau/CMSSW_10_2_10/src/job/BcJpsiTau_large_BcJpsiTauNu_truth_2020-04-24-150429/Myroot.root";
#    ofile="/work/ytakahas/work/BsTauTau/CMSSW_10_2_10/src/job/BcJpsiTau_large_BcJpsiTauNu_truth_2020-04-24-150429/Myroot.root";
    
root -l -q -b 'pick.C("'${ifile}'", "'${ofile}'", "'${ofile_training}'", "'${frac}'")'


