from multiprocessing import pool
import circuit 
import chper_wrapper as wrapper 
import MC_functions as mc 
import qcircuit_wrapper as qwrap
from circuits import prep_plus, prep_y, prep_zero, meas_xi, meas_zi
import copy
import json
import correction as cor
import BS17 
import BS49
import BS97
import BS161
import os
import multiprocessing
import argparse
from fullSampler import full_sampler


p = 0.2

parser = argparse.ArgumentParser(
    description="Simulates QEC on Bacon-Shor codes using the Steane syndrome extraction method"
)

parser.add_argument('--chp_path', type=str, required=True)
parser.add_argument('--out_dir', type=str, required=True)
parser.add_argument('--n_samples', type=int, default=1000)
parser.add_argument('--init_state', type=str, default='Z')
parser.add_argument('--distance', type=int, required=True)
parser.add_argument('--subset_idx', type=int, required=True)
parser.add_argument('--verification_kind', type=str, default="")
parser.add_argument('--use_exhaustive_sampling', action="store_true")
args = parser.parse_args()

chp_loc = args.chp_path
n_runs = args.n_samples
init_state_Pauli = args.init_state
index = args.subset_idx
d = args.distance
variant = args.verification_kind
folder_name = args.out_dir
use_exhaustive_samp = args.use_exhaustive_sampling

if d == 3:
    lookupTableX = BS17.Code.lookuptable['Xstabs']
    lookupTableZ = BS17.Code.lookuptable['Zstabs']
    code_stabs = BS17.Code.stabilizers[:]
    log_oper = BS17.Code.logical_opers[init_state_Pauli]
elif d == 5:
    lookupTableX = BS49.Code.lookuptable['Xstabs']
    lookupTableZ = BS49.Code.lookuptable['Zstabs']
    code_stabs = BS49.Code.stabilizers[:]
    log_oper = BS49.Code.logical_opers[init_state_Pauli]
elif d == 7:
    lookupTableX = BS97.Code.lookuptable['Xstabs']
    lookupTableZ = BS97.Code.lookuptable['Zstabs']
    code_stabs = BS97.Code.stabilizers[:]
    log_oper = BS97.Code.logical_opers[init_state_Pauli]
    circ_log_stabs, circ_log_destabs = BS161.Code.init_states[init_state_Pauli]
elif d == 9:
    lookupTableX = BS161.Code.lookuptable['Xstabs']
    lookupTableZ = BS161.Code.lookuptable['Zstabs']
    code_stabs = BS161.Code.stabilizers[:]
    log_oper = BS161.Code.logical_opers[init_state_Pauli]
    circ_log_stabs, circ_log_destabs = BS161.Code.init_states[init_state_Pauli]
    
error_model = 'standard'
p_1q, p_2q, p_meas = p, p, p 
error_dict, Is_after2q, Is_after_1q = wrapper.dict_for_error_model(error_model, p_1q, p_2q, p_meas)
error_info = mc.read_error_info(error_dict) 
faulty_gates_grouped = [['PrepareZPlus', 'H', 'X', 'Z', 'Y'], ['CX', 'CZ'], ['ImZ', 'ImX']] 

#Verifications
if d == 3 or variant == 'nv':
    cnots = []
elif d == 5:
    cnots = [(0,4)]
elif d == 7: 
    if variant == 'a':
        cnots = [(1,4), (2,5)]
    else:
        cnots = [(1,4), (1,5), (3,5)]
elif d == 9: 
    if variant == 'a':
        cnots = [(1,6), (2,7), (4,7)]
    else: 
        cnots = [(1,4), (2,7), (3,5), (4,7)]

num_ancilla = 2*d*len(cnots)

if num_ancilla > 0:
    ver = True
else:
    ver = False

#Measurements circuit
####################################################################################################
circ_meas = circuit.Circuit()
prep_plus(circ_meas, d, cnots, d*d, anc=3*d*d)
prep_zero(circ_meas, d, cnots, 2*d*d, anc=3*d*d+num_ancilla/2)
for i in range(0,d*d):
    circ_meas.add_gate_at([i,d*d+i], 'CX')
    circ_meas.add_gate_at([2*d*d+i,i], 'CX')
for k in range(0,d*d):
    circ_meas.add_gate_at([d*d+k], 'ImZ')
    circ_meas.add_gate_at([d*d+k], 'MeasureZ')
    circ_meas.add_gate_at([2*d*d+k], 'ImX')
    circ_meas.add_gate_at([2*d*d+k], 'MeasureX')
#if ver else 0 # for d=7 4d, for d=9 6d, (d-3)*d for d < 11
for l in range(d*d, 3*d*d+num_ancilla): #Changes if verified
   circ_meas.to_ancilla([l])

circ_meas_copy = [copy.deepcopy(circ_meas)]
####################################################################################################

summary_dict = {}  

#Run logical state preparation circuit
if d < 9:
#Logical state preparation circuit
    circ_log = circuit.Circuit()

    if init_state_Pauli == 'Z':
        prep_zero(circ_log, d)
    elif init_state_Pauli == 'X':
        prep_plus(circ_log, d)
    elif init_state_Pauli == 'Y':
        prep_y(circ_log, d)

    circ_log_copy = [copy.deepcopy(circ_log)]
    qoper_log = qwrap.Quantum_Operation([[],[]], circ_log_copy, chp_loc)
    qoper_log.run_one_circ(0)
    circ_log_stabs, circ_log_destabs = qoper_log.stabs[:], qoper_log.destabs[:]

list_faulty_gates = list(wrapper.gates_list_general([circ_meas], faulty_gates_grouped))

def subset_sim(run):

    circ_meas_copy = [copy.deepcopy(circ_meas)]

    if use_exhaustive_samp:
        faulty_circ = [faulty_circs[run]]
    else:
        errors_dict, carry_run, faulty_circ = wrapper.add_errors_fast_sampler_ion(list_faulty_gates, n_errors, circ_meas_copy, error_info)

    qoper_meas = qwrap.Quantum_Operation([circ_log_stabs, circ_log_destabs], faulty_circ, chp_loc)
    dict = qoper_meas.run_one_circ(0)
    #brow.from_circuit(faulty_circs[0], True)
    circ_meas_stabs, circ_meas_destabs = qoper_meas.stabs[:], qoper_meas.destabs[:]
    #print("Measurments",circ_meas_stabs)
    #print(ver)

    n_ver_fail = 0
    n_ver_total = 0
    #number of ancilla might be needed
    if ver:

        ver_pass_plus = True
        ver_pass_zero = True 

        for anc in range(0,int(num_ancilla/2)):
            if dict[3*d*d+anc][0] == 1:
                ver_pass_plus = False
            if dict[3*d*d+num_ancilla/2 + anc][0] == 1:
                ver_pass_zero = False
                n_ver_fail += 1 
            n_ver_total += 1

        ver_pass_both = ver_pass_zero and ver_pass_plus
        ver_pass_none = not ver_pass_zero and not ver_pass_plus

        if not ver_pass_both: 
            return ver_pass_both, ver_pass_none, ver_pass_plus, ver_pass_zero, 0, 0, -1, -1, n_ver_fail, n_ver_total

    syndrome_x = ""
    syndrome_z = ""

    
    
    for i in range(d-1):
        syndrome_x += str(meas_xi(dict, i, 2*d*d, d))
        syndrome_z += str(meas_zi(dict, i, d*d, d))

    error_z = lookupTableX[syndrome_x]
    error_x = lookupTableZ[syndrome_z]
    xErrorStr = '_'.join([str(_) for _ in error_z])

    #[[1,2], "Z"]

    circ_correct = circuit.Circuit()
    
    for j in range(d*d):
        circ_correct.add_gate_at([j], "I")
    for loc in error_z:
        circ_correct.add_gate_at([loc], "Z")
    for loc in error_x:
        circ_correct.add_gate_at([loc], "X")    

    circ_correct_copy = [copy.deepcopy(circ_correct)]
    ##brow.from_circuit(circ_correct_copy[0], True)

    qoper_correct = qwrap.Quantum_Operation([circ_meas_stabs, circ_meas_destabs], circ_correct_copy, chp_loc)
    qoper_correct.run_one_circ(0)
    circ_correct_stabs, circ_correct_destabs = qoper_correct.stabs[:], qoper_correct.destabs[:]
    #print("Correction: ", circ_correct_stabs)

    #Projection

    corr_circ= cor.Bare_Correct.generate_rep_bare_meas(d*d, code_stabs, 1, False, False,
                                                            False, False, False, True)

    
    corr_circ_list = []
    for supra_gate in corr_circ.gates:
        corr_circ_list += [supra_gate.circuit_list[0]]

    corr_object = qwrap.QEC_d3([circ_correct_stabs[:],circ_correct_destabs[:]], corr_circ_list[:], chp_loc)
    dm1, dm2, X_stab, Z_stab, Z_corr, X_corr, blah_stab, blah_destab = corr_object.run_fullQEC_CSS_perfect("BS"+str(2*d**2-1))
    z_weight = Z_corr.count("Z")
    x_weight = X_corr.count("X")
    corr_stabs = corr_object.stabs[:]
    corr_destabs = corr_object.destabs[:]
    

    #print("Projection", corr_stabs)
    # Determine if a failure has occured (for the lookup table decoder).
    
    #fail = False
    #for stab in corr_stabs[:25]:
    #    if stab[0] != '+':
    
    #        fail = True
    #        break
    #if fail:
    #    n_fail += 1

    # If the code is the surface code, we don't determine whether or not the logical
    # qubit was entangled with a gauge qubit.
    #log_random = False
    
    #print(log_oper)

    log_meas_circ = cor.Bare_Correct.generate_bare_meas(d*d, [log_oper], False, False)
    #brow.from_circuit(log_meas_circ, True)
    
    log_meas_object = qwrap.Quantum_Operation([corr_stabs,corr_destabs],
                                    [log_meas_circ], chp_loc)
    final_dict = log_meas_object.run_one_circ(0)
    log_outcome, log_random = list(final_dict.values())[0][0], list(final_dict.values())[0][1]

    n_fail = 0 
    if log_outcome == 1:
        n_fail = 1
    
        if init_state_Pauli == "Z":
            x_weight = d-x_weight
        elif init_state_Pauli == "X":
            z_weight = d-z_weight

    n_random = 0
    if log_random:
        n_random = 1
    if ver:
        return ver_pass_both, ver_pass_none, ver_pass_plus, ver_pass_zero, n_fail, n_random, x_weight, z_weight, n_ver_fail, n_ver_total
    else:
        return n_fail, n_random, x_weight, z_weight
    
max_weight = int((d-1)/2)
subsets = [[e,f,h] for e in range(max_weight) for f in range(max_weight) for h in range(max_weight) if e+f+h <= max_weight]
n_fail = 0. 
x_weights = [0. for i in range(d+1)]
z_weights = [0. for i in range(d+1)]
n_random = 0.
n_errors = subsets[index % max_weight]
if use_exhaustive_samp:
    faulty_circs = full_sampler(list_faulty_gates, n_errors, circ_meas_copy[0], error_info)
    print(len(faulty_circs))
#print(ver)
if ver: 
    n_ver_pass_zero = 0.  # Changes if verified
    n_ver_pass_plus = 0.
    n_ver_pass_both = 0.
    n_ver_pass_none = 0.

subset_str = str(n_errors[0]) + "_" + str(n_errors[1]) + "_" + str(n_errors[2]) 
runs = [i for i in range(n_runs)]
pool = multiprocessing.Pool()            
results = pool.map(subset_sim, runs)
pool.close()
pool.join()
if ver:
    n_ver_pass_both, n_ver_pass_none, n_ver_pass_plus, n_ver_pass_zero, n_fail, n_random, n_x_weights, n_z_weights, n_ver_fail, n_ver_total = zip(*results)
    n_ver_pass_both = float(n_ver_pass_both.count(True))
    n_ver_pass_none = float(n_ver_pass_none.count(True))
    n_ver_pass_plus = float(n_ver_pass_plus.count(True))
    n_ver_pass_zero = float(n_ver_pass_zero.count(True))
    n_ver_fail = float(sum(n_ver_fail))
    n_ver_total = float(sum(n_ver_total))
else:
    n_fail, n_random, n_x_weights, n_z_weights = zip(*results)
n_fail = float(sum(n_fail))
n_random = float(sum(n_random))
x_weights = [float(n_x_weights.count(i)) for i in range(d+1)]
z_weights = [float(n_z_weights.count(i)) for i in range(d+1)]

if ver:

    run_dict = {
        "n_runs" : n_runs,
        "n_fail": n_fail,
        "n_random": n_random,
        "n_ver_pass_zero": n_ver_pass_zero,
        "n_ver_pass_plus": n_ver_pass_plus,
        "n_ver_pass_both": n_ver_pass_both,
        "n_ver_pass_none": n_ver_pass_none,
        "n_ver_fail": n_ver_fail,
        "n_ver_total": n_ver_total
    }

    if n_ver_pass_both == 0:
        n_ver_pass_both = -1

    p_dict = {
        "p_fail" : n_fail/n_ver_pass_both,
        "p_random" : n_random/n_ver_pass_both,
        "p_ver_pass_zero": n_ver_pass_zero/n_runs,
        "p_ver_pass_plus": n_ver_pass_plus/n_runs,
        "p_ver_pass_both": n_ver_pass_both/n_runs,
        "p_ver_pass_none": n_ver_pass_none/n_runs, 
        "p_ver_fail": n_ver_fail/n_ver_total
    }

else:

    run_dict = {
        "n_runs" : n_runs,
        "n_fail": n_fail,
        "n_random": n_random,
    }

    p_dict = {
        "p_fail" : n_fail/n_runs,
        "p_random" : n_random/n_runs,
    }

for i in range(d+1):
    run_dict["n_x_weight_"+str(i)+"_errors"] = x_weights[i]
    p_dict["p_x_weight_"+str(i)+"_errors"] = x_weights[i]/n_ver_pass_both if ver else x_weights[i]/n_runs
    run_dict["n_z_weight_"+str(i)+"_errors"] = z_weights[i]
    p_dict["p_z_weight_"+str(i)+"_errors"] = z_weights[i]/n_ver_pass_both if ver else z_weights[i]/n_runs


summary_dict[subset_str] = p_dict

json_folder = os.path.join(folder_name, "steane", init_state_Pauli, str(d), variant)
file_name = subset_str + ".json"
run_path = os.path.join(json_folder, "occurrences", file_name)
p_path = os.path.join(json_folder, "probabilities", file_name)
run_json = json.dumps(run_dict, indent=4)
p_json = json.dumps(p_dict, indent=4)


os.makedirs(os.path.join(json_folder, "occurrences"), exist_ok=True)
os.makedirs(os.path.join(json_folder, "probabilities"), exist_ok=True)

with open(run_path, 'w') as json_file:
    json_file.write(run_json)

with open(p_path, 'w') as json_file:
    json_file.write(p_json)          




