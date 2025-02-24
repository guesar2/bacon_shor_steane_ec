'''
surface49.py
'''

import sys
import os
import json
from circuit import *
import correction

if __name__ == "__main__":

    from . import chper_wrapper as wrapper
    from . import MC_functions as mc 
    from . import circuit as cir
    from . import qcircuit_wrapper as qwrap
    from .visualizer import browser_vis as brow
    from . import correction as cor

    chp_loc = './chp_extended'

    error_model = 'standard'
    p = 0.2
    p_1q, p_2q, p_meas = p, p, p 
    error_dict, Is_after2q, Is_after_1q = wrapper.dict_for_error_model(error_model, p_1q, p_2q, p_meas)
    error_info = mc.read_error_info(error_dict)

    circ_prep = cir.Circuit()
    prep_CNOT_list = [[1,0,5,6], [2,1], [3,2,7,8], [4,3], [12,7,6,11], [14,9,8,13], [16,11,10,15], [18,13,12,17], [21,20], [22,17,16,21], [23,22], [24,19,18,23]]
    for index_list in prep_CNOT_list:
        circ_prep.add_gate_at([index_list[0]], 'H')
        for i in range(1, len(index_list)):
            circ_prep.add_gate_at([index_list[0], index_list[i]], 'CX')
    #for i in range(25):
    #   circ_prep.add_gate_at([i], 'H')
    circ_prep_copy = [copy.deepcopy(circ_prep)]
    qoper = qwrap.Quantum_Operation([[],[]], circ_prep_copy, chp_loc)
    qoper.run_one_circ(0)
    stabs, destabs = qoper.stabs[:], qoper.destabs[:]
    print(stabs) 
    sys.exit()
    log_meas_circ = cor.Bare_Correct.generate_bare_meas(25, stabilizers, False, False) #added the parallel thing

    log_meas_object = qwrap.Quantum_Operation([stabs,destabs],
                                    [log_meas_circ], chp_loc)
    final_dict = log_meas_object.run_one_circ(0)
    log_outcome, log_random = list(final_dict.values()), list(final_dict.values())
    print(log_outcome)
   #brow.from_circuit(circ_prep_copy[0], True)

class Code:
    '''
    Defines constants and methods useful for dealing with the Surface-49 code.
    '''

    code_name = 'Surface49'

    stabilizers = [
                    [('X',1), ('X',2)],
                    [('X',3), ('X',4)],
                    [('X',0), ('X',1), ('X',5), ('X',6)],
                    [('X',2), ('X',3), ('X',7), ('X',8)],
                    [('X',6), ('X',7), ('X',11), ('X',12)],
                    [('X',8), ('X',9), ('X',13), ('X',14)],
                    [('X',10), ('X',11), ('X',15), ('X',16)],
                    [('X',12), ('X',13), ('X',17), ('X',18)],
                    [('X',16), ('X',17), ('X',21), ('X',22)],
                    [('X',18), ('X',19), ('X',23), ('X',24)],
                    [('X',20), ('X',21)],
                    [('X',22), ('X',23)],
                    [('Z',0), ('Z',5)],
                    [('Z',1), ('Z',6), ('Z',2), ('Z',7)],
                    [('Z',3), ('Z',8), ('Z',4), ('Z',9)],
                    [('Z',5), ('Z',10), ('Z',6), ('Z',11)],
                    [('Z',7), ('Z',12), ('Z',8), ('Z',13)],
                    [('Z',9), ('Z',14)],
                    [('Z',10), ('Z',15)],
                    [('Z',11), ('Z',16), ('Z',12), ('Z',17)],
                    [('Z',13), ('Z',18), ('Z',14), ('Z',19)],
                    [('Z',15), ('Z',20), ('Z',16), ('Z',21)],
                    [('Z',17), ('Z',22), ('Z',18), ('Z',23)],
                    [('Z',19), ('Z',24)]]



    stabilizers_CHP = {
                        'Z': 
                            ['+XXIIIXXIIIIIIIIIIIIIIIIII', '+IXXIIIIIIIIIIIIIIIIIIIIII', '+IIXXIIIXXIIIIIIIIIIIIIIII', '+IIIXXIIIIIIIIIIIIIIIIIIII', '+IIIIIIXXIIIXXIIIIIIIIIIII', '+IIIIIIIIXXIIIXXIIIIIIIIII', '+IIIIIIIIIIXXIIIXXIIIIIIII', '+IIIIIIIIIIIIXXIIIXXIIIIII', '+IIIIIIIIIIIIIIIIXXIIIXXII', '+IIIIIIIIIIIIIIIIIIXXIIIXX', '+IIIIIIIIIIIIIIIIIIIIXXIII', '+IIIIIIIIIIIIIIIIIIIIIIXXI', '+ZZZZZIIIIIIIIIIIIIIIIIIII', '+IZZZZZIIIIIIIIIIIIIIIIIII', '+IIIZZIIZIIIIZIIIIIZIIIIIZ', '+IIIIIZZIIIIIZIIIIIZIIIIIZ', '+IIIIIIIZZIIIZIZIIIZIIIIIZ', '+IIIIIIIIIZIIIIZIIIIIIIIII', '+IIIIIIIIIIZIIIIIZIIIIIZZZ', '+IIIIIIIIIIIZZIIIZIZIIIZZI', '+IIIIIIIIIIIIIZZIIIZIIIIIZ', '+IIIIIIIIIIIIIIIZZIIIIIZZZ', '+IIIIIIIIIIIIIIIIIZZIIIZZI', '+IIIIIIIIIIIIIIIIIIIZIIIIZ', '+IIIIIIIIIIIIIIIIIIIIZZZZZ']
                      }
    
    destabilizers_CHP = {
                        'Z':
                        ['+IZZZZIIIIIIIIIIIIIIIIIIII', '+IIZZZIIIIIIIIIIIIIIIIIIII', '+IIIZZIIIIIIIIIIIIIIIIIIII', '+IIIIZIIIIIIIIIIIIIIIIIIII', '+IIIIIIIIIIIIZIIIIIZIIIIIZ', '+IIIIIIIIIIIIIIZIIIIIIIIII', '+IIIIIIIIIIIIIIIIZIIIIIZZZ', '+IIIIIIIIIIIIIIIIIIZIIIIIZ', '+IIIIIIIIIIIIIIIIIIIIIIZZZ', '+IIIIIIIIIIIIIIIIIIIIIIIIZ', '+IIIIIIIIIIIIIIIIIIIIIZZZZ', '+IIIIIIIIIIIIIIIIIIIIIIIZZ', '+XIIIIIIIIIIIIIIIIIIIIIIII', '+IIIIIXXIIIIIIIIIIIIIIIIII', '+IIIIIIIXXIIIIIIIIIIIIIIII', '+IIIIIIXIIIIIIIIIIIIIIIIII', '+IIIIIIIIXIIIIIIIIIIIIIIII', '+IIIIIIIIIXIIIIIIIIIIIIIII', '+IIIIIIIIIIXIIIIIIIIIIIIII', '+IIIIIIIIIIIXIIIIIIIIIIIII', '+IIIIIIIIIIIIIXIIIIIIIIIII', '+IIIIIIIIIIIIIIIXIIIIIIIII', '+IIIIIIIIIIIIIIIIIXIIIIIII', '+IIIIIIIIIIIIIIIIIIIXIIIII', '+IIIIIIIIIIIIIIIIIIIIXIIII']
                      }


    logical_opers = {'X': [('X',0), ('X',5), ('X',10), ('X',15), ('X',20)],
                     'Z': [('Z',0), ('Z',1), ('Z',2), ('Z',3), ('Z',4)],
                     'Y': [('Y',0), ('Z',1), ('Z',2), ('Z',3), ('Z',4),
                           ('X',5), ('X',10), ('X',15), ('X',20)]
                    }


    # The complete lookup table is imported from the json file.
    # The json file was generated with the script 'complete_lookups_surface49.py'.
    # The idea is to add all the errors in the basic lookup table
    # until we have obtained the 2**12 = 4096 possible syndromes. 
    # Just like for surface17, we have 2 lookup tables, one for the X stabilizers
    # and another one for the Z stabilizers

    lookup_folder = 'chp/lookup_tables_surface49/'
    lookuptable = {}
    stab_kinds = ['X','Z']
    for stab_kind in stab_kinds:
        json_filename = 'lookup_stabs%s.json'%stab_kind
        abs_filename = lookup_folder + json_filename
        json_file = open(abs_filename, 'r')
        local_table = json.load(json_file)
        json_file.close()
        lookuptable['%sstabs'%stab_kind] = local_table


