from multiprocessing import pool
import sys
import os
import json
import copy
import multiprocessing as mp
import circuit as cir
from circuits import prep_plus, prep_y, prep_zero
import BS17
import BS49
import BS97
import BS161
import surface17 as surf17
import surface49 as surf49
import correction as cor
import chper_wrapper as wrapper
import MC_functions as mc
import qcircuit_wrapper as qwrap
from decoder import find_usable_substring, numNonOverlapping11Subs
import math
import argparse


def perfCorr(lookupTableX, lookupTableZ, syndrome, d, measStabs, measDestabs):
    """
    Corrects X and Z errors, projects onto logical space, determines if the correction was correct and returns the error weight

    Args:
    lookupTableX (list): Look-up table for Z errros
    lookupTableZ (list): Look-up table for X errors
    syndrome (list): List with both X and Z stabilizers syndromes
    d (int): Distance
    measStabs: Measurment circuit final state stabilizers
    measDestabs: Measurment circuit final state destabilizers
    """
    syndromeX = "".join(map(str, syndrome[: int(len(syndrome) / 2)]))
    syndromeZ = "".join(map(str, syndrome[int(len(syndrome) / 2) :]))
    errorZ = lookupTableX[syndromeX]  # Check intervals
    errorX = lookupTableZ[syndromeZ]

    circ_correct = cir.Circuit()

    for i in range(d * d):
        circ_correct.add_gate_at([i], "I")
    for loc in errorZ:
        circ_correct.add_gate_at([loc], "Z")
    for loc in errorX:
        circ_correct.add_gate_at([loc], "X")

    circ_correct_copy = [copy.deepcopy(circ_correct)]

    qoper_correct = qwrap.Quantum_Operation(
        [measStabs, measDestabs], circ_correct_copy, chp_loc
    )
    qoper_correct.run_one_circ(0)
    circ_correct_stabs, circ_correct_destabs = (
        qoper_correct.stabs[:],
        qoper_correct.destabs[:],
    )

    corr_circ = cor.Bare_Correct.generate_rep_bare_meas(
        d * d, code_stabs, 1, False, False, False, False, False, True
    )  # set second-last to one

    corr_circ_list = []
    for supra_gate in corr_circ.gates:
        corr_circ_list += [supra_gate.circuit_list[0]]

    corr_object = qwrap.QEC_d3(
        [circ_correct_stabs[:], circ_correct_destabs[:]], corr_circ_list[:], chp_loc
    )
    dm1, dm2, X_stab, Z_stab, Z_corr, X_corr, blah_stab, blah_destab = (
        corr_object.run_fullQEC_CSS_perfect(QEC_kind)
    )
    z_weight = Z_corr.count("Z")
    x_weight = X_corr.count("X")
    corr_stabs = corr_object.stabs[:]
    corr_destabs = corr_object.destabs[:]

    log_meas_circ = cor.Bare_Correct.generate_bare_meas(d * d, [log_oper], False, False)
    log_meas_object = qwrap.Quantum_Operation(
        [corr_stabs, corr_destabs], [log_meas_circ], chp_loc
    )
    final_dict = log_meas_object.run_one_circ(0)
    log_outcome, log_random = list(final_dict.values())[0][0], list(final_dict.values())[0][1]

    if log_outcome == 1:
        if state == "Z":
            x_weight = d - x_weight
        elif state == "X":
            z_weight = d - z_weight

    return x_weight, z_weight, log_random, log_outcome


def generate_bare_meas(n_data, stabilizer_list):
    """
    Generates stabilizer measurements with bare ancilla, no cat states.

    stabilizer_list is assumed to have the following format:
        [
          [('X',0), ('X',4)], ...

        ]
    """

    n_ancilla = len(stabilizer_list)

    bare_meas_circ = cir.Circuit()

    for i in range(n_data, n_data + n_ancilla):
        bare_meas_circ.add_gate_at([i], "PrepareXPlus")

    for i in range(len(stabilizer_list)):
        gateprefix = "C"
        for gate in stabilizer_list[i]:
            bare_meas_circ.add_gate_at([n_data + i, gate[1]], gateprefix + gate[0])

    for i in range(len(stabilizer_list)):
        bare_meas_circ.add_gate_at([n_data + i], "ImX")
        bare_meas_circ.add_gate_at([n_data + i], "MeasureX")

    bare_meas_circ.to_ancilla(range(n_data, n_data + n_ancilla))

    return bare_meas_circ


def generate_rep_bare_meas(n_data, stabilizer_list, num_rounds):
    """ """

    rep_meas_circ = cir.Circuit()

    for i in range(num_rounds):
        gate_nameX = "SX_%i" % i
        stab_circ = generate_bare_meas(n_data, stabilizer_list[:])
        stab_circ = cir.Encoded_Gate("S%i" % i, [stab_circ]).circuit_wrap()
        rep_meas_circ.join_circuit(stab_circ, False)

    return rep_meas_circ


def run_QEC_adaptive(init_state, QEC_circ_list, kind="strong"):
    """
    kind: 'strong' or 'weak'
    """
    t = math.floor(float(d - 1) / 2)
    if kind == "strong":
        QEC_object = qwrap.QEC_d3(init_state, QEC_circ_list[:], chp_loc)
        s1 = QEC_object.run_one_bare_anc_new(0)
        # print(s1)
        s2 = QEC_object.run_one_bare_anc_new(1)
        # print(s2)
        diff_vec = "0" if s1 == s2 else "1"
        syndromes = []

        if d == 3:
            s1 = QEC_object.run_one_bare_anc_new(0)
            s2 = QEC_object.run_one_bare_anc_new(1)
            diff_vec = "0" if s1 == s2 else "1"
            if s1 == s2:
                syndrome = s2
                n_rounds = 2
            else:
                s3 = QEC_object.run_one_bare_anc_new(2)
                syndrome = s3
                n_rounds = 3
        else:
            end = False
            s1 = QEC_object.run_one_bare_anc_new(0)
            syndromes = [s1]
            diff_vec = ""
            i = 1
            while not end:
                s = QEC_object.run_one_bare_anc_new(i)
                i += 1
                syndromes.append(s)
                diff_vec += "0" if syndromes[-1] == syndromes[-2] else "1"
                usblSubs = find_usable_substring(diff_vec, t)
                num11subs, index = numNonOverlapping11Subs(diff_vec)
                if len(usblSubs) > 0:
                    end = True
                    index = usblSubs[0]
                    syndrome = syndromes[index]

                elif num11subs == t:
                    end = True
                    syndrome = syndromes[index]

            n_rounds = len(syndromes)

    elif kind == "weak":
        QEC_object = qwrap.QEC_d3(init_state, QEC_circ_list[:], chp_loc)
        s1 = QEC_object.run_one_bare_anc_new(0)  # or syndromes
        syndrome = s1
        if d == 3:
            n_rounds = 1
            if any(i != 0 for i in s1):
                s2 = QEC_object.run_one_bare_anc_new(1)
                syndrome = s2
                n_rounds = 2
        else:
            end = False
            s2 = QEC_object.run_one_bare_anc_new(1)
            syndromes = [s1, s2]
            diff_vec = "0" if s1 == s2 else "1"
            if any(i != 0 for i in s1):
                diff_vec = ""
                syndromes = syndromes[1:]
                t = t - 1
            else:
                diff_vec = "0" + diff_vec
                syndromes = [s1] + syndromes
            i = 2
            while not end:
                s = QEC_object.run_one_bare_anc_new(i)
                i += 1
                syndromes.append(s)
                diff_vec += "0" if syndromes[-1] == syndromes[-2] else "1"
                usblSubs = find_usable_substring(diff_vec, t)
                num11subs, index = numNonOverlapping11Subs(diff_vec)
                if len(usblSubs) > 0:
                    end = True
                    index = usblSubs[0]
                    syndrome = syndromes[index]
                elif num11subs == t:
                    end = True
                    syndrome = syndromes[-1]

            n_rounds = len(syndromes)
    x_weight, z_weight, log_random, log_outcome = perfCorr(
        lookupTableX, lookupTableZ, syndrome, d, QEC_object.stabs, QEC_object.destabs
    )

    return x_weight, z_weight, log_outcome, log_random, n_rounds


def run_subset(run):
    errors_dict, carry_run, faulty_circs = wrapper.add_errors_fast_sampler_ion(
        [oneq_gates, twoq_gates, meas_gates], n_errors, QEC_circ_list_copy, error_info
    )

    x_weight, z_weight, log_outcome, log_random, n_rounds = run_QEC_adaptive(
        init_state, faulty_circs, kind
    )
    failed, random = 0, 0
    if log_outcome == 1:
        failed = 1
    if log_random:
        random = 1
    return failed, random, x_weight, z_weight, n_rounds


####################################################################################################
############################################### Main ###############################################
####################################################################################################

parser = argparse.ArgumentParser(
    description="Simulates QEC on Bacon-Shor and Surface codes using the Shor syndrome extraction method"
)
parser.add_argument('--chp_path', type=str, required=True)
parser.add_argument('--out_dir', type=str, required=True)
parser.add_argument('--n_samples', type=int, default=1000)
parser.add_argument('--init_state', type=str, default='Z')
parser.add_argument('--qec_kind', type=str, required=True)
parser.add_argument('--subset_idx', type=int, required=True)
parser.add_argument('--shor_kind', type=str, required=True)
args = parser.parse_args()

n_runs = args.n_samples
state = args.init_state  # 'X' or 'Z' or 'Y'
i = args.subset_idx
QEC_kind = args.qec_kind
folder_name = args.out_dir
kind = args.shor_kind
chp_loc = args.chp_path


error_model = "standard"
p = 0.2
p_1q, p_2q, p_meas = p, p, p
error_dict, Is_after2q, Is_after_1q = wrapper.dict_for_error_model(
    error_model, p_1q, p_2q, p_meas
)
error_info = mc.read_error_info(error_dict)

# Define initial state
if QEC_kind == "BS17":
    code_stabs = BS17.Code.stabilizers[:]
    lookupTableX = BS17.Code.lookuptable["Xstabs"]
    lookupTableZ = BS17.Code.lookuptable["Zstabs"]
    log_oper = BS17.Code.logical_opers[state]
    n_data = 9
    num_rounds = 3
    d = 3
elif QEC_kind == "BS49":
    code_stabs = BS49.Code.stabilizers[:]
    lookupTableX = BS49.Code.lookuptable["Xstabs"]
    lookupTableZ = BS49.Code.lookuptable["Zstabs"]
    log_oper = BS49.Code.logical_opers[state]
    n_data = 25
    num_rounds = 5  # 8 d=7, 12 d=9
    d = 5
elif QEC_kind == "BS97":
    code_stabs = BS97.Code.stabilizers[:]
    lookupTableX = BS97.Code.lookuptable["Xstabs"]
    lookupTableZ = BS97.Code.lookuptable["Zstabs"]
    log_oper = BS97.Code.logical_opers[state]
    n_data = 49
    num_rounds = 8  # 8 d=7, 12 d=9
    d = 7
elif QEC_kind == "BS161":
    code_stabs = BS161.Code.stabilizers[:]
    lookupTableX = BS161.Code.lookuptable["Xstabs"]
    lookupTableZ = BS161.Code.lookuptable["Zstabs"]
    log_oper = BS161.Code.logical_opers[state]
    n_data = 81
    num_rounds = 12  # 8 d=7, 12 d=9
    d = 9
    circ_log_stabs, circ_log_destabs = BS161.Code.init_states[state]


elif QEC_kind == "surface17":
    init_stabs = surf17.Code.stabilizers_CHP[state][:]
    init_destabs = surf17.Code.destabilizers_CHP[state][:]
    init_state = [init_stabs, init_destabs]
    code_stabs = surf17.Code.stabilizers[:]
    lookupTableX = surf17.Code.lookuptable["Xstabs"]
    lookupTableZ = surf17.Code.lookuptable["Zstabs"]
    log_oper = surf17.Code.logical_opers[state]
    n_data = 9
    num_rounds = 3
    d = 3

elif QEC_kind == "surface49":
    init_stabs = surf49.Code.stabilizers_CHP[state][:]
    init_destabs = surf49.Code.destabilizers_CHP[state][:]
    init_state = [init_stabs, init_destabs]
    code_stabs = surf49.Code.stabilizers[:]
    lookupTableX = surf49.Code.lookuptable["Xstabs"]
    lookupTableZ = surf49.Code.lookuptable["Zstabs"]
    log_oper = surf49.Code.logical_opers[state]
    n_data = 25
    num_rounds = 5
    d = 5
else:
    raise NotImplementedError("QEC code not implemented.")
max_weight = int((d-1)/2)
# Define circuit and circuit list
QEC_circ = generate_rep_bare_meas(n_data, code_stabs, num_rounds)
QEC_circ_list = []
for supra_gate in QEC_circ.gates:
    QEC_circ_list += [supra_gate.circuit_list[0]]
# Prepare logical state
circ_log = cir.Circuit()
if state == "X":
    prep_plus(circ_log, d)
elif state == "Z":
    prep_zero(circ_log, d)
elif state == "Y":
    prep_y(circ_log, d)
circ_log_copy = [copy.deepcopy(circ_log)]
qoper_log = qwrap.Quantum_Operation([[], []], circ_log_copy, chp_loc)
qoper_log.run_one_circ(0)
circ_log_stabs, circ_log_destabs = qoper_log.stabs[:], qoper_log.destabs[:]
init_state = [circ_log_stabs, circ_log_destabs]


# Define the list of error-prone 1-q and 2-q gates
faulty_gates_grouped = [
    ["PrepareXPlus", "PrepareZPlus", "H"],
    ["CX", "CZ"],
    ["ImX", "ImZ"],
]
oneq_gates, twoq_gates, meas_gates = wrapper.gates_list_general(
    QEC_circ_list, faulty_gates_grouped
)

subsets = [
    [l, q, r]
    for l in range(max_weight + 1)
    for q in range(max_weight + 1)
    for r in range(max_weight + 1)
    if l + q + r <= max_weight
]

QEC_circ_list_copy = []
for subcirc in QEC_circ_list:
    QEC_circ_list_copy += [copy.deepcopy(subcirc)]

n_errors = subsets[i % max_weight]

runs = [i for i in range(n_runs)]


pool = mp.Pool()
results = pool.map(run_subset, runs)
pool.close()
pool.join()

n_fails, n_random, n_x_weights, n_z_weights, n_rounds = zip(*results)
n_rounds_avg = float(sum(n_rounds)) / n_runs
n_fails = float(sum(n_fails))
n_random = float(sum(n_random))
x_weights = [float(n_x_weights.count(i)) for i in range(d + 1)]
z_weights = [float(n_z_weights.count(i)) for i in range(d + 1)]

run_dict = {
    "n_runs": n_runs,
    "n_fail": n_fails,
    "n_random": n_random,
    "n_rounds_avg": n_rounds_avg,
}

p_dict = {"p_fail": n_fails / n_runs}

for i in range(d + 1):
    run_dict["n_x_weight_" + str(i) + "_errors"] = x_weights[i]
    run_dict["n_z_weight_" + str(i) + "_errors"] = z_weights[i]
    p_dict["p_x_weight_" + str(i) + "_errors"] = x_weights[i] / n_runs
    p_dict["p_z_weight_" + str(i) + "_errors"] = z_weights[i] / n_runs

json_folder = os.path.join(folder_name, "shor", state, QEC_kind, kind)
file_name = "_".join(map(str, n_errors)) + ".json"
run_path = os.path.join(json_folder, "occurrences", file_name)
p_path = os.path.join(json_folder, "probabilities", file_name)
run_json = json.dumps(run_dict, indent=4)
p_json = json.dumps(p_dict, indent=4)

os.makedirs(os.path.join(json_folder, "occurrences"), exist_ok=True)
os.makedirs(os.path.join(json_folder, "probabilities"), exist_ok=True)
with open(run_path, "w") as json_file:
    json_file.write(run_json)

with open(p_path, "w") as json_file:
    json_file.write(p_json)


##################################################################################################
############################################## End ###############################################
##################################################################################################
