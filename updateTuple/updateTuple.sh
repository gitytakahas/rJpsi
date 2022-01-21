#!/bin/bash
# update each ntuple by adding QCD fake weight and the MVA for ttbar rejection


if [ $# -ne 2 ]; then
  echo "# of specified params: #" 1>&2
  echo "You need to give 2 param" 1>&2
  exit 1
fi

#basedir="/work/ytakahas/work/BsTauTau/CMSSW_10_2_10/src/job/"

#for dir in BcJpsiTau_large_BcJpsiTauNu_truth_2020-04-24-150429 BJpsiX_BcJpsiTauNu_mc_2020-04-24-150431
#do

#ifile="${basedir}/${dir}/Myroot.root";
#ofile="${basedir}/${dir}/Myroot4mva.root";

ifile="$1"
#ofile="$2"
wfile="$2"
#tdir="$3"
#id="$4"
#    ifile="/work/ytakahas/work/BsTauTau/CMSSW_10_2_10/src/job/BcJpsiTau_large_BcJpsiTauNu_truth_2020-04-24-150429/Myroot.root";
#    ofile="/work/ytakahas/work/BsTauTau/CMSSW_10_2_10/src/job/BcJpsiTau_large_BcJpsiTauNu_truth_2020-04-24-150429/Myroot.root";
    
#root -l -q -b 'addvariable.C("'${ifile}'", "'${wfile}'", "'${tdir}'", "'${id}'")'
root -l -q -b 'addvariable.C("'${ifile}'", "'${wfile}'")'

#if [ $id == "-1" ]; then
#    echo "${id} is -1!"
#    mv -f ${ofile} ${ifile}
#else
#    echo "${id} is not -1!"
#    mv -f ${ofile}_${id}.root ${ifile}_${id}.root
#fi



