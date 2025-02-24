#!/bin/bash

qec_kind=("BS17" "BS49" "BS97" "BS161")
shor_kind=("weak"
 "strong")
states=("X" "Z")
n_samples=1

for ((i=0; i <= 285; i +=1)); do
    for q in "${qec_kind[@]}"; do
        for k in "${shor_kind[@]}"; do
            for s in "${states[@]}"; do
                PYTHONPATH="src:chp" python src/sim_shor_ec.py \
                --chp_path "chp/chp_extended" \
                --out_dir "out" \
                --n_samples $n_samples \
                --init_state "$s" \
                --qec_kind "$q" \
                --subset_idx "$i" \
                --shor_kind "$k"
            done
        done
    done
done