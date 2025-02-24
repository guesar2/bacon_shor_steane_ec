"""
plogical_steane_ec.py

Calculate logical rates for BS49 Steane error correction
"""
import MC_functions as mc 
import chper_wrapper as wrapper 
import copy
import json
import numpy as np
import sys
import os
import argparse

parser = argparse.ArgumentsParser()
parser.add_argument('--n_1q_gates', type=int, required=True)
parser.add_argument('--n_2q_gates', type=int, required=True)
parser.add_argument('--n_measurements', type=int, required=True)
parser.add_argument('--input_dir', type=str, required=True)
args = parser.parse_args()

n_gates = [args.n_1q_gates, args.n_2q_gates, args.n_measurements]
input_dir = args.input_file

print(n_gates)

def ocurrence_p(subset, p_list):
    """Calculates the probability of occurrence of a given error subset

    Parameters
    ----------
    subset : list
        Number of faulty: 1-qubit gates, 2-qubit gates and measurement gates, respectively
    p_list : list
        Error probability in each kind of gate

    Returns
    -------
    float
        Probaility of occurrence of the given error subset
    """

    subset_list = subset.split("_")
    subset_list = [int(item) for item in subset_list]
    probability = wrapper.prob_for_subset_general(n_gates, subset_list, p_list)
    
    return probability

summary_dict = {}
summary_path = os.path.join(input_dir, 'summary.json')
with open (summary_path, 'r') as json_file:
    summary_dict = json.load(json_file)



p_kinds = summary_dict["0_0_0"].keys()

p_upper_kinds = [i+"_upper" for i in p_kinds]
file_content = 'p ' + 'q ' + ' '.join(p_kinds) + ' ' + ' '.join(p_upper_kinds) + '\n'

probabilities = np.logspace(-6, -1, num=50)

for p in probabilities:
    for q in probabilities:

        p_set = [p,p,q]
        p_dict = dict.fromkeys(p_kinds, 0)
        p_dict_upper = dict.fromkeys(p_upper_kinds, 1)

        for subset, subset_dict in summary_dict.items():

            for i in range(len(p_kinds)):

                p_dict[p_kinds[i]] += ocurrence_p(subset, p_set)*subset_dict[p_kinds[i]]
                p_dict_upper[p_upper_kinds[i]] += ocurrence_p(subset, p_set)*(subset_dict[p_kinds[i]]-1)
                
            
        file_content += "%.15f" % p + " " + "%.15f" % q + " " +  " ".join("%.15f" % p_dict[p_kind] for p_kind in p_kinds) + " " + " ".join("%.15f" % p_dict_upper[p_kind] for p_kind in p_upper_kinds) + "\n"

data_filename = os.path.join(input_dir, 'results_2d.dat') 
with open(data_filename, 'w') as data_file:
    data_file.write(file_content)


    

