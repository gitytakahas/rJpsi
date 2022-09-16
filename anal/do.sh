
python draw.py --min -o plots_dir

#python draw.py --min --blind True # to veto data in SR

for h in $(seq 0 14)
do
    echo "test" $h
    python draw.py --sys hammer_ebe_e${h}_up --min -o plots_dir
    python draw.py --sys hammer_ebe_e${h}_down --min -o plots_dir
done


python draw.py --sys puweight_up --min -o plots_dir
python draw.py --sys puweight_down --min -o plots_dir

python draw.py --sys muSFID_up --min -o plots_dir
python draw.py --sys muSFID_down --min -o plots_dir

python draw.py --sys muSFReco_up --min -o plots_dir
python draw.py --sys muSFReco_down --min -o plots_dir

python draw.py --sys weight_ctau_up --min -o plots_dir
python draw.py --sys weight_ctau_down --min -o plots_dir

python draw.py --sys br_BcJpsiDst_up --min -o plots_dir
python draw.py --sys br_BcJpsiDst_down --min -o plots_dir

python draw.py --sys tauBr_up --min -o plots_dir
python draw.py --sys tauBr_down --min  -o plots_dir

python draw.py --sys tauReco_up --min -o plots_dir
python draw.py --sys tauReco_down --min -o plots_dir

python draw.py --sys xgbsEff_up --min -o plots_dir
python draw.py --sys xgbsEff_down --min -o plots_dir

python draw.py --sys BcPt_up --min -o plots_dir
python draw.py --sys BcPt_down --min -o plots_dir
