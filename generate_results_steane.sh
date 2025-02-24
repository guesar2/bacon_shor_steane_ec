#!/bin/bash

distances=("3" "5" "7" "9")
states=("X" "Z")

for d in "${distances[@]}"; do
    for s in "${states[@]}"; do
        PYTHONPATH="src:" python src/summarize.py \
            "out/steane/$s/$d/probabilities"
        n_1q_gates=$(jq --arg key "$d" '.[$key][0]' gatesNumSteaneEC.json)
        n_2q_gates=$(jq --arg key "$d" '.[$key][1]' gatesNumSteaneEC.json)
        n_measurements=$(jq --arg key "$d" '.[$key][2]' gatesNumSteaneEC.json)
        PYTHONPATH="src:chp" python src/plogical.py \
            --input_dir "out/steane/$s/$d/probabilities" \
            --n_1q_gates $n_1q_gates \
            --n_2q_gates $n_2q_gates \
            --n_measurements $n_measurements
    done
done

