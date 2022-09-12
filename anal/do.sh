
python draw.py --min

#python draw.py --min --blind True # to veto data in SR

for h in $(seq 0 9)
do
    echo "test" $h
    python draw.py --sys hammer_ebe_e${h}_up --min
    python draw.py --sys hammer_ebe_e${h}_down --min
done


python draw.py --sys puweight_up --min
python draw.py --sys puweight_down --min

python draw.py --sys muSFID_up --min
python draw.py --sys muSFID_down --min

python draw.py --sys muSFReco_up --min
python draw.py --sys muSFReco_down --min

python draw.py --sys weight_ctau_up --min
python draw.py --sys weight_ctau_down --min

python draw.py --sys br_BcJpsiDst_up --min 
python draw.py --sys br_BcJpsiDst_down --min 

python draw.py --sys tauBr_up --min 
python draw.py --sys tauBr_down --min 

python draw.py --sys tauReco_up --min
python draw.py --sys tauReco_down --min

python draw.py --sys xgbsEff_up --min     
python draw.py --sys xgbsEff_down --min

python draw.py --sys BcPt_up --min
python draw.py --sys BcPt_down --min
