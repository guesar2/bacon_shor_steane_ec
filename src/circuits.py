"""
circuits.py
"""
import circuit as chp_circuit

def prep_plus(circ, d, cnots=[], start=0, **kwargs):
    num_ancilla = len(cnots)
    anc = kwargs.get("anc", d * d)
    for k in range(start, start + d):
        for i in range(0, d):
            circ.add_gate_at([d * i + k], "PrepareZPlus")

        for i in range(num_ancilla):
            circ.add_gate_at([anc + i + num_ancilla * (k % d)], "PrepareZPlus")

        circ.add_gate_at([k], "H")

        for i in range(0, d - 1):
            circ.add_gate_at([d * i + k, d * (i + 1) + k], "CX")

        for i in range(num_ancilla):
            l, m = cnots[i]
            circ.add_gate_at([k + d * l, anc + num_ancilla * (k % d) + i], "CX")
            circ.add_gate_at([k + d * m, anc + num_ancilla * (k % d) + i], "CX")
            circ.add_gate_at([anc + num_ancilla * (k % d) + i], "ImZ")
            circ.add_gate_at([anc + num_ancilla * (k % d) + i], "MeasureZ")


def prep_zero(circ, d, cnots=[], start=0, **kwargs):
    num_ancilla = len(cnots)
    anc = kwargs.get("anc", d * d)
    for k in range(0, d):
        for i in range(start, d + start):
            circ.add_gate_at([d * k + i], "PrepareZPlus")

        for i in range(num_ancilla):
            circ.add_gate_at([anc + num_ancilla * k + i], "PrepareZPlus")

        circ.add_gate_at([d * k + start], "H")

        for i in range(start, d - 1 + start):
            circ.add_gate_at([d * k + i, d * k + i + 1], "CX")

        for i in range(num_ancilla):
            l, m = cnots[i]
            circ.add_gate_at([start + d * k + l, anc + num_ancilla * k + i], "CX")
            circ.add_gate_at([start + d * k + m, anc + num_ancilla * k + i], "CX")
            circ.add_gate_at([anc + num_ancilla * k + i], "ImZ")

        for i in range(start, d + start):
            circ.add_gate_at([d * k + i], "H")

        for i in range(num_ancilla):
            circ.add_gate_at([anc + num_ancilla * k + i], "MeasureZ")

def prep_y(circ, d, cnots=[], start=0, **kwargs):
    """
    Prepares logical Y state
    """
    anc = kwargs.get("anc", d * d)
    if cnots:
        raise NotImplementedError()

    prep_plus(circ, d, start=start)
    for i in range(d - 1, 0, -1):
        circ.add_gate_at([start + i, start + i - 1], "CX")
    circ.add_gate_at([0], "P")
    for i in range(1, d):
        circ.add_gate_at([start + i, start + i - 1], "CX")

def log_state_prep_circuit(state, d, **kwargs):
    circuit = chp_circuit.Circuit()
    prep_func_dict = {"X": prep_plus, "Y": prep_y, "Z": prep_zero}
    return prep_func_dict[state](circuit, d, kwargs)

def steane_measurement_circuit(d, verification):
    circuit = chp_circuit.Circuit()
    n_qubits = d * d
    plus_idx = n_qubits
    zero_idx = plus_idx + n_qubits
    plus_anc_idx = zero_idx + n_qubits
    zero_anc_idx = plus_anc_idx + len(verification)
    
    prep_plus(circuit, d, verification, start=plus_idx, anc=plus_anc_idx)
    prep_zero(circuit, d, verification, start=zero_idx, anc=zero_anc_idx)

    for i in range(n_qubits):
        circuit.add_gate_at([i, plus_idx + i], 'CX')
        circuit.add_gate_at([zero_idx + i, i], 'CX')

    for j in range(n_qubits):
        circuit.add_gate_at([plus_idx + j], 'ImZ')
        circuit.add_gate_at([plus_idx + j], 'MeasureZ')
        circuit.add_gate_at([zero_idx + j], 'ImX')
        circuit.add_gate_at([zero_idx + j], 'MeasureX')
    
    for k in range(plus_idx, zero_anc_idx + len(verification)):
        circuit.to_ancilla([k])

    return [circuit]

def generate_bare_meas(n_data, stabilizer_list):
    '''
    Generates stabilizer measurements with bare ancilla, no cat states.

    stabilizer_list is assumed to have the following format:
        [
          [('X',0), ('X',4)], ... 

        ]
    '''

    n_ancilla = len(stabilizer_list)
    
    bare_meas_circ = chp_circuit.Circuit()

    for i in range(n_data, n_data+n_ancilla):
        bare_meas_circ.add_gate_at([i], 'PrepareXPlus')
        
    for i in range(len(stabilizer_list)):
        gateprefix = 'C'
        for gate in stabilizer_list[i]:
            bare_meas_circ.add_gate_at([n_data+i,gate[1]], gateprefix+gate[0])

    for i in range(len(stabilizer_list)):
        bare_meas_circ.add_gate_at([n_data+i], 'ImX')
        bare_meas_circ.add_gate_at([n_data+i], 'MeasureX')
       
    bare_meas_circ.to_ancilla(range(n_data, n_data+n_ancilla))


def generate_rep_bare_meas(n_data, stabilizer_list, num_rounds):
    '''
    '''

    rep_meas_circ = chp_circuit.Circuit()

    for i in range(num_rounds):
        stab_circ = generate_bare_meas(n_data, stabilizer_list[:])
        stab_circ = chp_circuit.Encoded_Gate('S%i'%i, [stab_circ]).circuit_wrap()
        rep_meas_circ.join_circuit(stab_circ, False)

    return rep_meas_circ

def shor_measurement_circuit(d, stabilizers, rounds):
    circuit = generate_rep_bare_meas(d * d, stabilizers, rounds)
    circuit_list = [supra_gate.circuit_list[0] for supra_gate in circuit.gates]
    return circuit_list

def prep_GHZ(circ, d, vers=[]):
    num_ancilla = len(vers)
    for i in range(0, d + num_ancilla):
        circ.add_gate_at([i], "PrepareZPlus")

    circ.add_gate_at([0], "H")

    for i in range(0, d - 1):
        circ.add_gate_at([i, i + 1], "CX")
    
    for i, anc in enumerate(vers):
        for ver in anc:
            circ.add_gate_at([ver, d + i], "CX")

    for i in range(d, d + num_ancilla):
        circ.add_gate_at([i], "ImZ")
    for i in range(0, d + num_ancilla):
        circ.add_gate_at([i], "MeasureZ")
    for i in range(d, d + num_ancilla):
        circ.to_ancilla([i])




def meas_xi(dict, i, start, d):  # i in 0,d-1
    stabm = 0
    for j in range(0, d):
        stabm += dict[start + d * j + i][0] + dict[start + d * j + i + 1][0]
    stabm = stabm % 2
    return stabm


def meas_zi(dict, i, start, d):  # i in 0,4
    stabm = 0
    for j in range(0, d):
        stabm += dict[start + d * i + j][0] + dict[start + d * (i + 1) + j][0]
    stabm = stabm % 2
    return stabm


def mZZ(dict, i):
    return (dict[i][0] + dict[i + 1][0]) % 2


"""
def perfCorr(lookupTableX, lookupTableZ, syndrome, d, measStabs, measDestabs):
    Corrects X and Z errors, projects onto logical space, determines if the correction was correct and returns the error weight
    
    Args:
    lookupTableX (list): Look-up table for Z errros
    lookupTableZ (list): Look-up table for X errors
    syndrome (list): List with both X and Z stabilizers syndromes
    d (int): Distance
    measStabs: Measurment circuit final state stabilizers
    measDestabs: Measurment circuit final state destabilizers
    syndromeX = ''.join(map(str, syndrome[:d-1]))  
    syndromeZ = ''.join(map(str, syndrome[d-1:]))  
    errorZ = lookupTableX[syndromeX] #Check intervals
    errorX = lookupTableZ[syndromeZ]

    circ_correct = circuit.Circuit()

    for i in range(d*d):
        circ_correct.add_gate_at([i], "I")
    for loc in errorZ:
        circ_correct.add_gate_at([loc], "Z")
    for loc in errorX:
        circ_correct.add_gate_at([loc], "X")    

    circ_correct_copy = [copy.deepcopy(circ_correct)]
    ##brow.from_circuit(circ_correct_copy[0], True)

    qoper_correct = qwrap.Quantum_Operation([measStabs, measDestabs], circ_correct_copy, chp_loc)
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
    dm1, dm2, X_stab, Z_stab, Z_corr, X_corr, blah_stab, blah_destab = corr_object.run_fullQEC_CSS_perfect(QEC_kind)
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
    log_outcome, log_random = final_dict.values()[0][0], final_dict.values()[0][1]


    if log_outcome == 1:

        if state == "Z":
            x_weight = d-x_weight
        elif state == "X":
            z_weight = d-z_weight 

    
    return x_weight, z_weight, log_random, log_outcome
"""
