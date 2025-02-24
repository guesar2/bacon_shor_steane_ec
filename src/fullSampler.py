'''
fullSampler.py
''' 
import copy
import circuit 
from visualizer import browser_vis as brow
import chper_wrapper as wrapper
import MC_functions as mc
import operator as op
from functools import reduce
from itertools import combinations

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom  # or / in Python 2

def getCardinality(nErrors, nGates):
    total = 1
    for i in range(len(nErrors)):
        total *= ncr(nGates[i], nErrors[i])
    return total


def getAllErrors(nErrors, nGates):
    allErrors = []
    for i in range(nErrors):
        if i == 0:
            newErrors = [[j] for j in range(nGates)]
        else:
            newErrors = []
            for error in allErrors:
                for k in range(error[-1]):
                    newError = copy.deepcopy(error)
                    newError.append(k)
                    newErrors.append(newError)
        allErrors = newErrors
    return allErrors

def getErrors(nErrors, nGates):
    allErrors =  list(combinations([i for i in range(nGates)], nErrors))
    return allErrors

def full_sampler(gate_indices, n_errors, circ, error_info):

    total_selected_gates = [[] for _ in range(len(n_errors))]

    for i in range(len(n_errors)):
         
        possibleErrors = getAllErrors(n_errors[i], len(gate_indices[i]))
        for error in possibleErrors:
            selected_gates = []
            for index in error:
                selected_gates.append(gate_indices[i][index])
            sorted_selected_gates = sorted(selected_gates, key=lambda gate: gate[0]) 
            total_selected_gates[i] += [sorted_selected_gates]
        if len(possibleErrors) == 0:
            total_selected_gates[i] += [[]]

    possible_gate_selections = []

    for i in total_selected_gates[0]:
        for j in total_selected_gates[1]:
            for k in total_selected_gates[2]:
                possible_gate_selections.append([i, j, k])

    faulty_circs = []
        
    for total_selected_gates in possible_gate_selections: 

        local_gates = []

        for i_gate in range(len(n_errors)):
            for pair in total_selected_gates[i_gate]:
                if 0 == pair[0]:
                    local_gates += [pair[1]]
        
        if len(local_gates) != 0:
            faulty_circs += add_error(circ, error_info, local_gates)
    
    return faulty_circs

        

def add_error(circ, error_info, gate_indices=[]):
    circuit = copy.deepcopy(circ)
    circuits = [circuit]
    gate_indices.sort()  # this step is critical
    for i in gate_indices[::-1]:
        g = circ.gates[i]
        error_ratio_dic = error_info.dic[g.gate_name]['error_ratio']
        errors = [str(x[0]) for x in error_ratio_dic.items()]
        new_circuits = []
        if len(g.qubits) == 1:
            for c in circuits:
                for key in errors: 
                    circuit = copy.deepcopy(c)
                
                    g = circuit.gates[i]
                    new_g = circuit.insert_gate(g, g.qubits, '', key, False) #same below
                    new_g.is_error = True
                    new_circuits.append(circuit)
                    
            # MGA 12/23/19.  New addition for the NN decoder
        else:
            for c in circuits:
                for error in errors:
                    circuit = copy.deepcopy(c)
                    g = circuit.gates[i]
                    for q_index in range(len(g.qubits))[::-1]:
                        new_g = circuit.insert_gate(g, [g.qubits[q_index]], '', error[q_index], False)
                        new_g.is_error = True
                    new_circuits.append(circuit)
        circuits = new_circuits
    return circuits

if __name__ == '__main__':

    import sys


    p = 0.2
    chp_loc = '../chp/chp_extended'
    error_model = 'standard'
    p_1q, p_2q, p_meas = p, p, p 
    error_dict, Is_after2q, Is_after_1q = wrapper.dict_for_error_model(error_model, p_1q, p_2q, p_meas)
    error_info = mc.read_error_info(error_dict) 
    faulty_gates_grouped = [['PrepareZPlus', 'H', 'X', 'Z', 'Y'], ['CX', 'CZ'], ['ImZ', 'ImX']] 


    testCircuit = circuit.Circuit()
    testCircuit.add_gate_at([0], 'PrepareZPlus')
    testCircuit.add_gate_at([1], 'PrepareZPlus')
    testCircuit.add_gate_at([0, 1], 'CX')


    testCircuitCopy = copy.deepcopy(testCircuit)
    list_faulty_gates = list(wrapper.gates_list_general([testCircuit], faulty_gates_grouped))
    #faulty_circs, t, a = wrapper.add_errors_fast_sampler_ion(list_faulty_gates, [1,1,1], [testCircuitCopy], error_info)
    faulty_circs = full_sampler(list_faulty_gates, [2,1,2], testCircuitCopy, error_info)

    #brow.from_circuit(faulty_circs[0], True)

    print(list_faulty_gates)
    print(len(faulty_circs))
    print(getCardinality([4,0,0], [315,414,2166]))
    print("9b: ", getCardinality([4,0,0], [333,450,234]))
    print("9a: ", getCardinality([4,0,0], [315,414,216]))
    


