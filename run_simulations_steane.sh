#!/bin/bash

distances=(3 5 7 9)
states=("X" "Z")
n_samples=1

for ((i=0; i <= 285; i +=1)); do
    for d in "${distances[@]}"; do
        for s in "${states[@]}"; do
            PYTHONPATH="src:chp" python src/sim_steane_ec.py \
                --chp_path "chp/chp_extended" \
                --out_dir "out" \
                --n_samples "$n_samples" \
                --init_state "$s" \
                --distance "$d" \
                --subset_idx "$i"
        done
    done
done