#!/bin/bash

qec_kind=("BS17" "BS49" "BS97" "BS161")
shor_kind=("weak" "strong")
states=("X" "Z")

for q in "${qec_kind[@]}"; do
    for s in "${states[@]}"; do
        for k in "${shor_kind[@]}"; do
            PYTHONPATH="src:" python src/summarize.py \
                "out/shor/$s/$q/$k/probabilities"
            n_1q_gates=$(jq --arg key "$q" '.[$key][0]' gatesNumShorEC.json)
            n_2q_gates=$(jq --arg key "$q" '.[$key][1]' gatesNumShorEC.json)
            n_measurements=$(jq --arg key "$q" '.[$key][2]' gatesNumShorEC.json)
            PYTHONPATH="src:chp" python src/plogical.py \
                --input_dir "out/shor/$s/$q/$k/probabilities" \
                --n_1q_gates $n_1q_gates \
                --n_2q_gates $n_2q_gates \
                --n_measurements $n_measurements
        done
    done
done

