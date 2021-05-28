run(){
    
    echo "region" $1

    python draw_stack.py --channel ${1} --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_a0_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_a0_down --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_a1_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_a1_down --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_a2_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_a2_down --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_b0_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_b0_down --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_b1_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_b1_down --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_b2_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_b2_down --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_c1_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_c1_down --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_c2_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_c2_down --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_d0_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_d0_down --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_d1_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_d1_down --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_d2_up --min
    python draw_stack.py --channel ${1} --sys hammer_ebe_d2_down --min

#    python draw_stack.py --channel ${1} --sys weight_ctau_up --min
#    python draw_stack.py --channel ${1} --sys weight_ctau_down --min
    
}


for ch in sr sb cr1 cr2
do
    run $ch &
done
