path="/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/"

for year in 2016 2017 2018
do
    for typename in "Data/data.root" "BcJpsiTau_inclusive/sig.root" "BJpsiX/bkg.root"
    do
	for postfix in approval MuonPhys
	do
	    filename="${path}/job_pt_${year}_${postfix}/${typename}"

	    echo "==========================="
	    echo $filename

	    root -l $filename <<EOF
tree->GetEntries()
EOF
	done
    done
done
