
python draw.py --min

python draw.py --min --blind  # to veto data in SR

for h in $(seq 0 9)
do
    echo "test" $h
    python draw.py --sys hammer_ebe_e${h}_up --min
    python draw.py --sys hammer_ebe_e${h}_down --min
done


python draw.py --sys puweight_up --min
python draw.py --sys puweight_down --min
